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

    if [ ! -s "${tbi}" ] || [ "${bam}" -nt "${tbi}" ]; then
        rm -f "${vcf}" "${vcfgz}" "${tbi}"
        lofreq call-parallel --pp-threads "${THREADS}" -f "${REF}" -o "${vcf}" "${bam}"
        bgzip "${vcf}"
        tabix -p vcf "${vcfgz}"
    fi
done

tsv=results/collapsed.tsv
need_tsv=0
if [ ! -s "${tsv}" ]; then
    need_tsv=1
else
    for sample in "${SAMPLES[@]}"; do
        if [ "results/${sample}.vcf.gz" -nt "${tsv}" ]; then
            need_tsv=1
            break
        fi
    done
fi

if [ "${need_tsv}" -eq 1 ]; then
    tmp="${tsv}.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${tmp}"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz" >> "${tmp}"
    done
    mv "${tmp}" "${tsv}"
fi