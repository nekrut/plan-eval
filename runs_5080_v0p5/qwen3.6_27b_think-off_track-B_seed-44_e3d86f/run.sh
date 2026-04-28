#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4

# Sample list
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Index reference if not already indexed
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${BAM}.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${VCF_GZ}.tbi"

    # Check if BAM exists and is valid (idempotency check)
    if [[ -f "$BAM" ]] && [[ -f "$BAI" ]] && [[ -f "$VCF_GZ" ]] && [[ -f "$TBI" ]]; then
        continue
    fi

    # Align reads
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | \
            samtools view -bS -o "$BAM"
    fi

    # Sort and index BAM if not already done
    if [[ ! -f "$BAI" ]]; then
        samtools sort -o "$BAM" "$BAM"
        samtools index "$BAM"
    fi

    # Call variants using lofreq
    if [[ ! -f "$VCF_GZ" ]]; then
        lofreq call -f "$REF" -o - "$BAM" | \
            bcftools view -Oz -o "$VCF_GZ"
    fi

    # Index VCF if not already done
    if [[ ! -f "$TBI" ]]; then
        tabix -p vcf "$VCF_GZ"
    fi
done

# Create collapsed TSV
COLLAPSED="${OUT_DIR}/collapsed.tsv"

# Check if collapsed.tsv already exists and has content (idempotency)
if [[ -f "$COLLAPSED" ]] && [[ -s "$COLLAPSED" ]]; then
    exit 0
fi

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"

# Process each VCF and append to collapsed TSV
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    
    # Extract variants and format as TSV
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n' "$VCF_GZ" | \
        awk -v sample="$SAMPLE" '{print sample "\t" $0}' >> "$COLLAPSED"
done