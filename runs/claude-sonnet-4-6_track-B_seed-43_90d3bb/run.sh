#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"

mkdir -p "${RESULTS}"

# Index reference if needed
[[ -f "${REF}.bwt" ]] || bwa index "${REF}"
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW}/${SAMPLE}_1.fq.gz"
    R2="${RAW}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS}/${SAMPLE}.bam"
    IQ_BAM="${RESULTS}/${SAMPLE}.iq.bam"
    VCF="${RESULTS}/${SAMPLE}.vcf.gz"

    # Align and sort
    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "${REF}" "${R1}" "${R2}" \
          | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    # Index sorted BAM
    [[ -f "${BAM}.bai" ]] || samtools index "${BAM}"

    # Add indel quality scores (required for lofreq --call-indels)
    if [[ ! -f "${IQ_BAM}" ]]; then
        lofreq indelqual --dindel -f "${REF}" -o "${IQ_BAM}" "${BAM}"
        samtools index "${IQ_BAM}"
    fi

    # Call variants
    if [[ ! -f "${VCF}" ]]; then
        lofreq call-parallel --pp-threads "${THREADS}" \
            -f "${REF}" \
            --call-indels \
            -o "${RESULTS}/${SAMPLE}.vcf" \
            "${IQ_BAM}"
        bgzip "${RESULTS}/${SAMPLE}.vcf"
        tabix -p vcf "${VCF}"
    fi
done

# Build collapsed table
COLLAPSED="${RESULTS}/collapsed.tsv"
if [[ ! -f "${COLLAPSED}" ]]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${RESULTS}/${SAMPLE}.vcf.gz" >> "${COLLAPSED}"
    done
fi