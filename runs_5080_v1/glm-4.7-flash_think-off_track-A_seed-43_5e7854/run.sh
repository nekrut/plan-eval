#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RAW="data/raw"

mkdir -p results

# 2. Reference indexing (once)
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
fi

# 3-5. Per-sample alignment, sorting, and indexing
for SAMPLE in "${SAMPLES[@]}"; do
    INPUT_1="${RAW}/${SAMPLE}_1.fq.gz"
    INPUT_2="${RAW}/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"
    VCF="results/${SAMPLE}.vcf"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="results/${SAMPLE}.vcf.gz.tbi"

    if [[ -f "${VCF_TBI}" ]]; then
        continue
    fi

    bwa mem -t "${THREADS}" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" "${REF}" "${INPUT_1}" "${INPUT_2}" | \
        samtools sort -@ "${THREADS}" -o "${BAM}"
    samtools index -@ "${THREADS}" "${BAM}"

    lofreq call-parallel --pp-threads "${THREADS}" -f "${REF}" "${BAM}" -o "${VCF}"
    bgzip -c "${VCF}" > "${VCF_GZ}"
    tabix -p vcf "${VCF_GZ}"
    rm -f "${VCF}"
done

# 8. Collapse step
COLLAPSED="results/collapsed.tsv"
HEADER="sample\tchrom\tpos\tref\talt\taf"

if [[ -f "${COLLAPSED}" ]]; then
    # Check if any input VCF is newer than the TSV
    NEED_REBUILD=0
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="results/${SAMPLE}.vcf.gz"
        if [[ "${VCF_GZ}" -nt "${COLLAPSED}" ]]; then
            NEED_REBUILD=1
            break
        fi
    done
    if [[ "${NEED_REBUILD}" -eq 0 ]]; then
        exit 0
    fi
fi

{
    echo -e "${HEADER}"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${SAMPLE}.vcf.gz"
    done
} > "${COLLAPSED}"