#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# Step 2: Reference indexing (idempotent)
if [ ! -f "$REF_DIR/chrM.fa.fai" ] || [ ! -f "$REF_DIR/chrM.fa.bwt" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Step 3-5: Per-sample alignment, sorting, indexing
for sample in "${SAMPLES[@]}"; do
    bam_file="$RESULTS_DIR/${sample}.bam"
    bai_file="${bam_file}.bai"
    vcf_file="$RESULTS_DIR/${sample}.vcf"
    vcf_gz_file="${vcf_file}.gz"
    vcf_tbi_file="${vcf_gz_file}.tbi"

    # Skip if BAM and BAI already exist
    if [ -f "$bam_file" ] && [ -f "$bai_file" ]; then
        :
    else
        # Step 3 & 4: Align and sort
        bwa mem -t $THREADS \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "$RAW_DIR/${sample}_1.fq.gz" \
            "$RAW_DIR/${sample}_2.fq.gz" \
            | samtools sort -@ $THREADS -o "$bam_file" -
        
        # Step 5: Index BAM
        samtools index -@ $THREADS "$bam_file"
    fi

    # Step 6: Variant calling with lofreq
    if [ ! -f "$vcf_file" ]; then
        lofreq call-parallel --pp-threads $THREADS --verbose \
            --ref "$REF" \
            --out "$vcf_file" \
            --sig \
            --bonf \
            "$bam_file"
    fi

    # Step 7: VCF compression and indexing
    if [ ! -f "$vcf_gz_file" ] || [ ! -f "$vcf_tbi_file" ]; then
        bgzip -f "$vcf_file"
        tabix -p vcf "$vcf_gz_file"
    fi
done

# Step 8: Collapse step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
needs_collapse=false

for sample in "${SAMPLES[@]}"; do
    vcf_gz_file="$RESULTS_DIR/${sample}.vcf.gz"
    if [ ! -f "$vcf_gz_file" ]; then
        needs_collapse=true
        break
    fi
done

if [ ! -f "$COLLAPSED" ]; then
    needs_collapse=true
fi

if [ "$needs_collapse" = true ]; then
    # Create header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    
    # Append data for each sample
    for sample in "${SAMPLES[@]}"; do
        vcf_gz_file="$RESULTS_DIR/${sample}.vcf.gz"
        if [ -f "$vcf_gz_file" ]; then
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_gz_file" >> "$COLLAPSED"
        fi
    done
fi