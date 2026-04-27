#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
REF_BAI="data/ref/chrM.fa.bai"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory
mkdir -p "$OUT_DIR"

# Index reference if not present
if [[ ! -f "$REF_BAI" ]]; then
    samtools faidx "$REF"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Alignment with BWA-MEM
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools view -bS - > "$BAM"
    fi

    # Sort and index BAM
    if [[ ! -f "$BAI" ]]; then
        samtools sort -o "$BAM" "$BAM"
        samtools index "$BAM"
    fi

    # Variant calling with LoFreq
    if [[ ! -f "$VCF_GZ" ]]; then
        lofreq call -f "$REF" -o "${OUT_DIR}/${SAMPLE}.vcf" "$BAM"
        bgzip "${OUT_DIR}/${SAMPLE}.vcf"
        rm "${OUT_DIR}/${SAMPLE}.vcf"
    fi

    # Index VCF
    if [[ ! -f "$TBI" ]]; then
        tabix "$VCF_GZ"
    fi
done

# Generate collapsed table
COLLAPSED="${OUT_DIR}/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        # Extract sample, chrom, pos, ref, alt, AF (INFO/AF)
        bcftools query -f "%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" >> "$COLLAPSED"
    done
fi