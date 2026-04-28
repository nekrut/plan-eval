#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF_SRC="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"
REF_DIR="${OUT_DIR}/ref"
REF="${REF_DIR}/chrM.fa"

mkdir -p "${OUT_DIR}" "${REF_DIR}"

if [ ! -s "${REF}" ]; then
    cp "${REF_SRC}" "${REF}"
fi

if [ ! -s "${REF}.bwt" ]; then
    bwa index "${REF}"
fi

if [ ! -s "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for sample in "${SAMPLES[@]}"; do
    bam="${OUT_DIR}/${sample}.bam"
    bai="${OUT_DIR}/${sample}.bam.bai"
    vcf="${OUT_DIR}/${sample}.vcf.gz"
    tbi="${OUT_DIR}/${sample}.vcf.gz.tbi"
    r1="${RAW_DIR}/${sample}_1.fq.gz"
    r2="${RAW_DIR}/${sample}_2.fq.gz"

    if [ ! -s "${bam}" ]; then
        tmp_bam="${OUT_DIR}/${sample}.unsorted.bam"
        bwa mem -t "${THREADS}" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "${REF}" "${r1}" "${r2}" \
            | samtools view -@ "${THREADS}" -b -F 4 -o "${tmp_bam}" -
        samtools sort -@ "${THREADS}" -o "${bam}" "${tmp_bam}"
        rm -f "${tmp_bam}"
    fi

    if [ ! -s "${bai}" ]; then
        samtools index -@ "${THREADS}" "${bam}"
    fi

    if [ ! -s "${vcf}" ]; then
        raw_vcf="${OUT_DIR}/${sample}.lofreq.vcf"
        lofreq call-parallel --pp-threads "${THREADS}" -f "${REF}" -o "${raw_vcf}" "${bam}"
        bcftools view -Oz -o "${vcf}" "${raw_vcf}"
        rm -f "${raw_vcf}"
    fi

    if [ ! -s "${tbi}" ]; then
        tabix -p vcf "${vcf}"
    fi
done

collapsed="${OUT_DIR}/collapsed.tsv"
if [ ! -s "${collapsed}" ]; then
    tmp_collapsed="${collapsed}.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${tmp_collapsed}"
    for sample in "${SAMPLES[@]}"; do
        vcf="${OUT_DIR}/${sample}.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${vcf}" \
            | awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$1,$2,$3,$4,$5}' \
            >> "${tmp_collapsed}"
    done
    mv "${tmp_collapsed}" "${collapsed}"
fi