#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# Index reference
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS_DIR}/${SAMPLE}.bam"
    VCF_UNCOMPRESSED="${RESULTS_DIR}/${SAMPLE}.vcf"
    VCF="${RESULTS_DIR}/${SAMPLE}.vcf.gz"
    BAI="${RESULTS_DIR}/${SAMPLE}.bam.bai"
    
    # Alignment
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | \
        samtools view -bS - | \
        samtools sort -o "$BAM" -
    fi
    
    # Index BAM
    if [[ ! -f "$BAI" ]]; then
        samtools index "$BAM"
    fi
    
    # Variant calling
    if [[ ! -f "$VCF" ]]; then
        lofreq call -f "$REF" -o "$VCF_UNCOMPRESSED" "$BAM"
        bgzip -f "$VCF_UNCOMPRESSED"
        tabix -p vcf "$VCF"
    fi
done

# Create collapsed table
if [[ ! -f "${RESULTS_DIR}/collapsed.tsv" ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        
        for SAMPLE in "${SAMPLES[@]}"; do
            VCF="${RESULTS_DIR}/${SAMPLE}.vcf.gz"
            bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%AF\n" "$VCF" | \
            awk -v sample="$SAMPLE" '{print sample "\t" $0}'
        done
    } > "${RESULTS_DIR}/collapsed.tsv"
fi