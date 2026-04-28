#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF=data/ref/chrM.fa

mkdir -p results

if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

if [ ! -f "${REF}.bwt" ] || [ ! -f "${REF}.amb" ] || [ ! -f "${REF}.ann" ] || [ ! -f "${REF}.pac" ] || [ ! -f "${REF}.sa" ]; then
    bwa index "${REF}"
fi

for sample in "${SAMPLES[@]}"; do
    fq1="data/raw/${sample}_1.fq.gz"
    fq2="data/raw/${sample}_2.fq.gz"
    bam="results/${sample}.bam"
    bai="results/${sample}.bam.bai"
    vcf="results/${sample}.vcf"
    vcfgz="results/${sample}.vcf.gz"
    tbi="results/${sample}.vcf.gz.tbi"

    if [ ! -s "${bam}" ] || [ "${fq1}" -nt "${bam}" ] || [ "${fq2}" -nt "${bam}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "${REF}" "${fq1}" "${fq2}" \
            | samtools sort -@ "${THREADS}" -o "${bam}" -
    fi

    if [ ! -s "${bai}" ] || [ "${bam}" -nt "${bai}" ]; then
        samtools index -@ "${THREADS}" "${bam}"
    fi

    if [ ! -s "${vcfgz}" ] || [ "${bam}" -nt "${vcfgz}" ]; then
        lofreq call-parallel --pp-threads "${THREADS}" -f "${REF}" -o "${vcf}" "${bam}"
        bgzip -f "${vcf}"
    fi

    if [ ! -s "${tbi}" ] || [ "${vcfgz}" -nt "${tbi}" ]; then
        tabix -p vcf "${vcfgz}"
    fi

    if [ -f "${vcf}" ]; then
        rm -f "${vcf}"
    fi
done

collapsed=results/collapsed.tsv
need_rebuild=0
if [ ! -s "${collapsed}" ]; then
    need_rebuild=1
else
    for sample in "${SAMPLES[@]}"; do
        if [ "results/${sample}.vcf.gz" -nt "${collapsed}" ]; then
            need_rebuild=1
            break
        fi
    done
fi

if [ "${need_rebuild}" -eq 1 ]; then
    tmp="${collapsed}.tmp"
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz"
        done
    } > "${tmp}"
    mv "${tmp}" "${collapsed}"
fi