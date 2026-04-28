#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW="data/raw"
RESULTS="results"

SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p "${RESULTS}"

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "${REF}"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW}/${SAMPLE}_1.fq.gz"
    R2="${RAW}/${SAMPLE}_2.fq.gz"
    BAM="${RESULTS}/${SAMPLE}.bam"
    VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"

    if [[ ! -f "${BAM}" ]] || [[ ! -f "${BAM}.bai" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "${REF}" "${R1}" "${R2}" \
            | samtools sort -@ "${THREADS}" -o "${RESULTS}/${SAMPLE}.sorted.bam"

        lofreq indelqual --dindel \
            -f "${REF}" \
            -o "${BAM}" \
            "${RESULTS}/${SAMPLE}.sorted.bam"

        rm -f "${RESULTS}/${SAMPLE}.sorted.bam"
        samtools index "${BAM}"
    fi

    if [[ ! -f "${VCF_GZ}" ]] || [[ ! -f "${VCF_GZ}.tbi" ]]; then
        lofreq call-parallel \
            --pp-threads "${THREADS}" \
            -f "${REF}" \
            --call-indels \
            -o "${RESULTS}/${SAMPLE}.vcf" \
            "${BAM}"

        bgzip -f "${RESULTS}/${SAMPLE}.vcf"
        tabix -p vcf "${VCF_GZ}"
    fi
done

COLLAPSED="${RESULTS}/collapsed.tsv"
if [[ ! -f "${COLLAPSED}" ]]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${COLLAPSED}"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${RESULTS}/${SAMPLE}.vcf.gz"
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${VCF_GZ}" >> "${COLLAPSED}"
    done
fi