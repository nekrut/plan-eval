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
mkdir -p "${OUT_DIR}"

# Index reference if not already indexed
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
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
    if [[ -f "${BAM}" ]] && [[ -f "${BAI}" ]] && [[ -f "${VCF_GZ}" ]] && [[ -f "${TBI}" ]]; then
        continue
    fi

    # Align reads
    bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" | \
        samtools view -bS - | \
        samtools sort -@ "${THREADS}" -o "${BAM}"

    # Index BAM
    samtools index "${BAM}"

    # Call variants using lofreq
    # lofreq call -> filter -> sort -> index
    lofreq call -f "${REF}" -o - "${BAM}" | \
        lofreq filter -f "${REF}" -o - | \
        bcftools sort -Oz -o "${VCF_GZ}"

    # Index VCF
    tabix -p vcf "${VCF_GZ}"
done

# Generate collapsed TSV
COLLAPSED="${OUT_DIR}/collapsed.tsv"

# Check if collapsed.tsv already exists and has content (idempotency)
if [[ -f "${COLLAPSED}" ]] && [[ -s "${COLLAPSED}" ]]; then
    exit 0
fi

# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "${COLLAPSED}"

# Process each VCF
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    
    # Extract variants: sample, chrom, pos, ref, alt, af
    # AF is in INFO field, typically as AF=x.xxxx
    # We use bcftools query to extract fields
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%INFO/AF]\n' "${VCF_GZ}" | \
    while IFS=$'\t' read -r CHROM POS REF ALT AF; do
        # Handle multi-allelic sites if necessary, but lofreq usually splits them
        # AF might be missing or multiple values; take the first one
        # If AF is missing, we might need to calculate it or skip
        # For simplicity, we assume AF is present in INFO
        if [[ -z "${AF}" ]]; then
            AF="NA"
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE}" "${CHROM}" "${POS}" "${REF}" "${ALT}" "${AF}"
    done >> "${COLLAPSED}"
done

exit 0