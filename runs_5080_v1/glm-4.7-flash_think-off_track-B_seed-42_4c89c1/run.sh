#!/usr/bin/env bash
set -euo pipefail

# Directory setup
RAW_DIR="data/raw"
REF_DIR="data/ref"
OUT_DIR="results"

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Sample list
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Index reference genome
if [ ! -f "$REF_DIR/chrM.fa.bai" ]; then
    samtools index "$REF_DIR/chrM.fa"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # Define input and output paths
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    VCF="${OUT_DIR}/${SAMPLE}.vcf.gz"

    # Check if outputs already exist
    if [ -f "$BAM" ] && [ -f "$BAM.bai" ] && [ -f "$VCF" ] && [ -f "$VCF.tbi" ]; then
        continue
    fi

    # Align reads
    bwa mem -t 4 "$REF_DIR/chrM.fa" "$R1" "$R2" | \
        samtools view -b - | \
        samtools sort -@ 2 -o "$BAM"
    samtools index "$BAM"

    # Call variants
    lofreq call-parallel -f "$REF_DIR/chrM.fa" -o "$VCF" -m indel -t 2 "$BAM"
    tabix -p vcf "$VCF"
done

# Collapse VCFs into a single TSV
# Using bcftools to concatenate and then SnpSift to format
bcftools concat -Oz -o "${OUT_DIR}/merged.vcf.gz" "${OUT_DIR}"/*.vcf.gz
tabix -p vcf "${OUT_DIR}/merged.vcf.gz"

# Filter and format
# Using SnpSift to extract relevant fields and calculate allele frequency
java -jar $(which SnpSift.jar) Filter -e "AF > 0" "${OUT_DIR}/merged.vcf.gz" | \
    java -jar $(which SnpSift.jar) Tsv -e "sample,chrom,pos,ref,alt,AF" -o "${OUT_DIR}/collapsed.tsv"

# Cleanup intermediate files
rm "${OUT_DIR}/merged.vcf.gz" "${OUT_DIR}/merged.vcf.gz.tbi"