#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_FA="data/ref/chrM.fa"

mkdir -p "$RESULTS_DIR"

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx "$REF_FA"
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
    bwa index "$REF_FA"
fi

for sample in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${sample}.bam"
    BAI="$RESULTS_DIR/${sample}.bam.bai"
    VCF_RAW="$RESULTS_DIR/${sample}.vcf"
    VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
    TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"

    # Alignment and Sorting
    if [ ! -f "$BAI" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "$REF_FA" \
            data/raw/${sample}_1.fq.gz \
            data/raw/${sample}_2.fq.gz | \
            samtools sort -@ "$THREADS" -o "$BAM"
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Variant Calling
    if [ ! -f "$TBI" ]; then
        lofreq call-parallel --pp-threads "$THREADS" \
            -f "$REF_FA" \
            -o "$VCF_RAW" \
            "$BAM"
        
        # Note: bgzip is part of htslib, usually available if samtools/tabix are. 
        # If not explicitly in inventory but required for .vcf.gz, we use it here.
        bgzip -c "$VCF_RAW" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF_RAW"
    fi
done

# Collapse results
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
# Check if any VCF is newer than the TSV to determine if rebuild is needed
NEEDS_REBUILD=false
for sample in "${SAMPLES[@]}"; do
    if [ "$RESULTS_DIR/${sample}.vcf.gz" -nt "$COLLAPSED" ] && [ -f "$COLLAPSED" ]; then
        NEEDS_REBUILD=true
        break
    fi
done

if [ "$NEEDS_REBUILD" = true ] || [ ! -f "$COLLAPSED" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RESULTS_DIR/${sample}.vcf.gz" >> "$COLLAPSED"
    done
fi