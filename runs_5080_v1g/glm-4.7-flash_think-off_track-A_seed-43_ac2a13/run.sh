#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# 2. Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa
fi

# 3-5. Per-sample alignment, sorting, and indexing
for SAMPLE in "${SAMPLES[@]}"; do
    RAW_1="data/raw/${SAMPLE}_1.fq.gz"
    RAW_2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"

    if [[ ! -f "${BAI}" ]]; then
        bwa mem -t "${THREADS}" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" "${RAW_1}" "${RAW_2}" | samtools sort -@ "${THREADS}" -o "${BAM}"
        samtools index -@ "${THREADS}" "${BAM}"
    fi
done

# 6. Variant calling with lofreq
for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="results/${SAMPLE}.vcf.gz.tbi"

    if [[ ! -f "${VCF_TBI}" ]]; then
        lofreq call-parallel --pp-threads "${THREADS}" --verbose --ref data/ref/chrM.fa --out "${VCF}" --sig --bonf "results/${SAMPLE}.bam"
        bgzip -c "${VCF}" > "${VCF_GZ}"
        tabix -p vcf "${VCF_GZ}"
        rm "${VCF}"
    fi
done

# 8. Collapse step
COLLAPSED="results/collapsed.tsv"
HEADER="sample\tchrom\tpos\tref\talt\taf"

if [[ ! -f "${COLLAPSED}" ]]; then
    {
        echo -e "${HEADER}"
        for SAMPLE in "${SAMPLES[@]}"; do
            bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${SAMPLE}.vcf.gz"
        done
    } > "${COLLAPSED}"
fi