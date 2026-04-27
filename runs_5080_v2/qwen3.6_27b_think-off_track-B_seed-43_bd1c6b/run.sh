#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4

# Sample list
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory
mkdir -p "${OUT_DIR}"

# 1. Index Reference (idempotent)
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

# 2. Per-sample processing
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Check if final VCF is already present; if so, skip sample
    if [[ -f "${VCF_GZ}" && -f "${TBI}" ]]; then
        continue
    fi

    # Alignment
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" "${REF}" "${R1}" "${R2}" | \
            samtools view -b -@ "${THREADS}" -o "${BAM}" -
    fi

    # Sorting and Indexing
    if [[ ! -f "${BAI}" ]]; then
        samtools sort -@ "${THREADS}" -o "${BAM}" "${BAM}"
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # Variant Calling with Lofreq
    # Using --call-indels to be thorough, though mtDNA is mostly SNVs.
    # --min-af 0.01 to capture low frequency variants typical in mtDNA heteroplasmy.
    if [[ ! -f "${VCF_GZ}" ]]; then
        lofreq call --call-indels -f "${REF}" -o "${VCF_GZ}" --min-af 0.01 "${BAM}"
        # Lofreq outputs uncompressed VCF by default if not gzipped extension, 
        # but we specified .vcf.gz so it should be bgzipped. 
        # Ensure it is bgzipped and indexed.
        if [[ ! -f "${TBI}" ]]; then
            tabix -p vcf "${VCF_GZ}"
        fi
    fi
done

# 3. Collapse VCFs into TSV
COLLAPSED="${OUT_DIR}/collapsed.tsv"

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "${COLLAPSED}"

for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    
    # Extract INFO field for AF (Allele Frequency)
    # bcftools query format:
    # %CHROM %POS %REF %ALT %INFO/AF
    # Note: Lofreq AF is in INFO field. 
    # We filter for PASS or . and valid variants.
    
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${VCF_GZ}" | \
    while IFS=$'\t' read -r CHROM POS REF ALT AF; do
        # Handle multi-allelic sites if any (though rare in mtDNA amplicons)
        # bcftools query splits ALT and AF by comma if multi-allelic.
        # We need to expand them.
        
        IFS=',' read -ra ALT_ARR <<< "${ALT}"
        IFS=',' read -ra AF_ARR <<< "${AF}"
        
        for i in "${!ALT_ARR[@]}"; do
            ALT_VAL="${ALT_ARR[$i]}"
            AF_VAL="${AF_ARR[$i]:-0}"
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE}" "${CHROM}" "${POS}" "${REF}" "${ALT_VAL}" "${AF_VAL}"
        done
    done >> "${COLLAPSED}"
done