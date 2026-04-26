#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"

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
    vcf_tmp="results/${sample}.vcf"
    vcf_gz="results/${sample}.vcf.gz"
    vcf_tbi="results/${sample}.vcf.gz.tbi"

    if [ ! -f "${bam}" ] || [ "${fq1}" -nt "${bam}" ] || [ "${fq2}" -nt "${bam}" ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            "${REF}" "${fq1}" "${fq2}" \
            | samtools sort -@ "${THREADS}" -o "${bam}"
    fi

    if [ ! -f "${bai}" ] || [ "${bam}" -nt "${bai}" ]; then
        samtools index -@ "${THREADS}" "${bam}"
    fi

    if [ ! -f "${vcf_gz}" ] || [ "${bam}" -nt "${vcf_gz}" ]; then
        rm -f "${vcf_tmp}" "${vcf_gz}" "${vcf_tbi}"
        lofreq call-parallel --pp-threads "${THREADS}" -f "${REF}" -o "${vcf_tmp}" "${bam}"
        bgzip -f "${vcf_tmp}"
    fi

    if [ ! -f "${vcf_tbi}" ] || [ "${vcf_gz}" -nt "${vcf_tbi}" ]; then
        tabix -f -p vcf "${vcf_gz}"
    fi
done

collapsed="results/collapsed.tsv"
need_rebuild=0
if [ ! -f "${collapsed}" ]; then
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
    tmp_collapsed="${collapsed}.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${tmp_collapsed}"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz" >> "${tmp_collapsed}"
    done
    mv "${tmp_collapsed}" "${collapsed}"
fi