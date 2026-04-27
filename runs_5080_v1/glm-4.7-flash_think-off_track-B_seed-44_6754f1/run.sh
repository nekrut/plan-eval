#!/usr/bin/env bash
set -euo pipefail

# Directories
RAW_DIR="data/raw"
REF_DIR="data/ref"
OUT_DIR="results"

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# Define samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Index reference if not already indexed
if [ ! -f "$REF_DIR/chrM.fa.bai" ]; then
    samtools index "$REF_DIR/chrM.fa"
fi

# Loop over samples
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    
    # Check if input files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "Error: Input files not found for $SAMPLE"
        exit 1
    fi

    # BAM file path
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    
    # Skip if BAM already exists
    if [ -f "$BAM" ]; then
        continue
    fi

    # 1. Align reads using BWA-MEM
    bwa mem -t 4 "$REF_DIR/chrM.fa" "$R1" "$R2" | \
        samtools view -b - > "$BAM"

    # 2. Sort BAM
    samtools sort -@ 2 -o "${BAM}.sorted" "$BAM"
    mv "${BAM}.sorted" "$BAM"

    # 3. Index BAM
    samtools index "$BAM"

    # 4. Call variants using LoFreq
    lofreq call-parallel -f "$REF_DIR/chrM.fa" -o "${OUT_DIR}/${SAMPLE}.vcf.gz" -m indel -t 2 "$BAM"

    # 5. Index VCF
    tabix -p vcf "${OUT_DIR}/${SAMPLE}.vcf.gz"
done

# 6. Collapse VCFs into a single TSV
#    Using bcftools to concatenate and extract required fields
bcftools concat -Oz -o "${OUT_DIR}/collapsed.vcf.gz" "${OUT_DIR}"/*.vcf.gz
tabix -p vcf "${OUT_DIR}/collapsed.vcf.gz"

# 7. Convert to TSV and filter
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "${OUT_DIR}/collapsed.vcf.gz" > "${OUT_DIR}/collapsed.tsv"

echo "Workflow completed successfully."