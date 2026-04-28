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
if [ ! -f "$REF_DIR/chrM.fa.fai" ] || [ ! -f "$REF_DIR/chrM.fa.amb" ]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Step 3-7: Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${sample}.bam"
    BAI="$BAM.bai"
    VCF="$RESULTS_DIR/${sample}.vcf"
    VCF_GZ="${VCF}.gz"
    VCF_TBI="${VCF_GZ}.tbi"

    # Skip if all outputs exist
    if [ -f "$BAI" ] && [ -f "$VCF_TBI" ]; then
        continue
    fi

    # Step 3 & 4: Alignment and sorting
    if [ ! -f "$BAM" ] || [ "$RAW_DIR/${sample}_1.fq.gz" -nt "$BAM" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF" \
            "$RAW_DIR/${sample}_1.fq.gz" \
            "$RAW_DIR/${sample}_2.fq.gz" \
            | samtools sort -@ "$THREADS" -o "$BAM" -
    fi

    # Step 5: BAM indexing
    if [ ! -f "$BAI" ]; then
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Step 6: Variant calling
    if [ ! -f "$VCF" ] || [ "$BAM" -nt "$VCF" ]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref "$REF" --out "$VCF" \
            --sig --bonf \
            "$BAM"
    fi

    # Step 7: VCF compression and indexing
    if [ ! -f "$VCF_GZ" ] || [ "$VCF" -nt "$VCF_GZ" ]; then
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm -f "$VCF"
    fi
done

# Step 8: Collapse step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
NEED_COLLAPSE=0

for sample in "${SAMPLES[@]}"; do
    VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
    if [ ! -f "$COLLAPSED" ] || [ "$VCF_GZ" -nt "$COLLAPSED" ]; then
        NEED_COLLAPSE=1
        break
    fi
done

if [ "$NEED_COLLAPSE" -eq 1 ]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED"
fi