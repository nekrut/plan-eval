#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAW="data/raw"
OUT="results"

mkdir -p "${OUT}"

[[ -f "${REF}.bwt" ]] || bwa index "${REF}"
[[ -f "${REF}.fai" ]] || samtools faidx "${REF}"

SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${OUT}/${SAMPLE}.bam"
    IQ_BAM="${OUT}/${SAMPLE}.iq.bam"
    VCF_GZ="${OUT}/${SAMPLE}.vcf.gz"

    if [[ ! -f "${BAM}" ]]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
            "${REF}" \
            "${RAW}/${SAMPLE}_1.fq.gz" \
            "${RAW}/${SAMPLE}_2.fq.gz" \
        | samtools sort -@ "${THREADS}" -o "${BAM}" -
    fi

    [[ -f "${BAM}.bai" ]] || samtools index -@ "${THREADS}" "${BAM}"

    if [[ ! -f "${IQ_BAM}" ]]; then
        lofreq indelqual --dindel -f "${REF}" -o "${IQ_BAM}" "${BAM}"
        samtools index "${IQ_BAM}"
    fi

    if [[ ! -f "${VCF_GZ}" ]]; then
        lofreq call \
            -f "${REF}" \
            --call-indels \
            -o "${OUT}/${SAMPLE}.vcf" \
            "${IQ_BAM}"
        bcftools view -Oz -o "${VCF_GZ}" "${OUT}/${SAMPLE}.vcf"
        tabix -p vcf "${VCF_GZ}"
        rm -f "${OUT}/${SAMPLE}.vcf"
    fi
done

if [[ ! -f "${OUT}/collapsed.tsv" ]]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query \
            -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
            "${OUT}/${SAMPLE}.vcf.gz" \
        >> "${OUT}/collapsed.tsv"
    done
fi