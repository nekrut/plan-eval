#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"

mkdir -p results

if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi

for SAMPLE in ${SAMPLES}; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="results/${SAMPLE}.bam.bai"
    VCF="results/${SAMPLE}.vcf"
    VCFGZ="results/${SAMPLE}.vcf.gz"
    TBI="results/${SAMPLE}.vcf.gz.tbi"

    if [ ! -f "${BAM}" ] || [ "${R1}" -nt "${BAM}" ] || [ "${R2}" -nt "${BAM}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "${REF}" "${R1}" "${R2}" \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    if [ ! -f "${BAI}" ] || [ "${BAM}" -nt "${BAI}" ]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    if [ ! -f "${TBI}" ] || [ "${BAM}" -nt "${TBI}" ]; then
        lofreq call-parallel --pp-threads "${THREADS}" --verbose \
            --ref "${REF}" \
            --out "${VCF}" \
            "${BAM}"
        bgzip -f "${VCF}"
        tabix -p vcf "${VCFGZ}"
    fi
done

COLLAPSED="results/collapsed.tsv"
REBUILD=0
if [ ! -f "${COLLAPSED}" ]; then
    REBUILD=1
else
    for SAMPLE in ${SAMPLES}; do
        if [ "results/${SAMPLE}.vcf.gz" -nt "${COLLAPSED}" ]; then
            REBUILD=1
            break
        fi
    done
fi

if [ "${REBUILD}" -eq 1 ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in ${SAMPLES}; do
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "results/${SAMPLE}.vcf.gz" >> "${COLLAPSED}"
    done
fi