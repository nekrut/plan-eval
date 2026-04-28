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
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" | \
            samtools view -bS -o "${BAM}"
    fi

    # Sort and index BAM if not already done
    if [[ ! -f "${BAI}" ]]; then
        samtools sort -o "${BAM}" "${BAM}"
        samtools index "${BAM}"
    fi

    # Call variants using lofreq
    if [[ ! -f "${VCF_GZ}" ]]; then
        lofreq call --bcf -f "${REF}" -o "${VCF_GZ}" "${BAM}"
    fi

    # Index VCF if not already done
    if [[ ! -f "${TBI}" ]]; then
        tabix -p vcf "${VCF_GZ}"
    fi
done

# Create collapsed TSV
COLLAPSED="${OUT_DIR}/collapsed.tsv"

# Check if collapsed file already exists and has content (idempotency)
if [[ -f "${COLLAPSED}" ]] && [[ -s "${COLLAPSED}" ]]; then
    exit 0
fi

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "${COLLAPSED}"

# Process each VCF and append to collapsed TSV
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    
    # Extract variants using bcftools query and format as TSV
    # Columns: CHROM, POS, REF, ALT, AF (from FORMAT/AF or INFO/AF)
    # lofreq typically outputs AF in INFO field or FORMAT field
    # We'll use bcftools query to extract the necessary fields
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF_GZ}" | \
    while IFS=$'\t' read -r CHROM POS REF ALT AF; do
        # Skip lines where AF is missing or not a number
        if [[ -z "${AF}" ]] || [[ "${AF}" == "." ]]; then
            continue
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE}" "${CHROM}" "${POS}" "${REF}" "${ALT}" "${AF}"
    done >> "${COLLAPSED}"
done