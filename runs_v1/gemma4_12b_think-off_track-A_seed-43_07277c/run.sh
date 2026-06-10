#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_FA="data/ref/chrM.fa"

mkdir -p "$RESULTS_DIR"

# Reference indexing
if [ ! -f "data/ref/chrM.fa.fai" ]; then
    samtools faidx "$REF_FA"
fi
if [ ! -f "data/ref/chrM.fa.amb" ]; then
    bwa index "$REF_FA"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${SAMPLE}.bam"
    BAI="$RESULTS_DIR/${SAMPLE}.bam.bai"
    VCF_RAW="$RESULTS_DIR/${SAMPLE}.vcf"
    VCF_GZ="$RESULTS_DIR/${SAMPLE}.vcf.gz"
    TBI="$RESULTS_DIR/${SAMPLE}.vcf.gz.tbi"

    # Alignment and Sorting
    if [ ! -f "$BAI" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "$REF_FA" \
            "data/raw/${SAMPLE}_1.fq.gz" \
            "data/raw/${SAMPLE}_2.fq.gz" | \
            samtools sort -@ "$THREADS" -o "$BAM"
        
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Variant Calling
    if [ ! -f "$TBI" ]; then
        lofreq call-parallel --pp-threads "$THREADS" \
            --ref "$REF_FA" \
            "$BAM" > "$VCF_RAW"
        
        bgzip -c "$VCF_RAW" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF_RAW"
    fi
done

# Collapsed Table
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
# Check if any input VCF is newer than the TSV
NEEDS_REBUILD=false
for SAMPLE in "${SAMPLES[@]}"; do
    if [ "$COLLAPSED" != "" ] && [ "$(find "$RESULTS_DIR/${SAMPLE}.vcf.gz" -newer "$COLLAPSED" | wc -l)" -gt 0 ]; then
        NEEDS_REBUILD=true
    fi
done

if [ "$NEEDS_REBUILD" = true ] || [ ! -f "$COLLAPSED" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RESULTS_DIR/${SAMPLE}.vcf.gz" >> "$COLLAPSED"
    done
fi