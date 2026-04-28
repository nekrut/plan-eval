#!/usr/bin/env bash
set -euo pipefail

THREADS=4
RAW=data/raw
REF_DIR=data/ref
REF=${REF_DIR}/chrM.fa
OUT=results
mkdir -p "${OUT}"

SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

# Index reference (idempotent)
if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

for s in "${SAMPLES[@]}"; do
    BAM="${OUT}/${s}.bam"
    BAI="${OUT}/${s}.bam.bai"
    VCF="${OUT}/${s}.vcf.gz"
    TBI="${OUT}/${s}.vcf.gz.tbi"

    if [ ! -s "${BAM}" ] || [ ! -s "${BAI}" ]; then
        TMP_SORT="${OUT}/${s}.sort.tmp.bam"
        TMP_VIT="${OUT}/${s}.viterbi.tmp.bam"
        TMP_INDEL="${OUT}/${s}.indel.tmp.bam"
        bwa mem -t "${THREADS}" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
            "${REF}" "${RAW}/${s}_1.fq.gz" "${RAW}/${s}_2.fq.gz" \
            | samtools sort -@ "${THREADS}" -o "${TMP_SORT}" -
        samtools index -@ "${THREADS}" "${TMP_SORT}"
        lofreq viterbi -f "${REF}" "${TMP_SORT}" \
            | samtools sort -@ "${THREADS}" -o "${TMP_VIT}" -
        samtools index -@ "${THREADS}" "${TMP_VIT}"
        lofreq indelqual --dindel -f "${REF}" -o "${TMP_INDEL}" "${TMP_VIT}"
        samtools sort -@ "${THREADS}" -o "${BAM}" "${TMP_INDEL}"
        samtools index -@ "${THREADS}" "${BAM}"
        rm -f "${TMP_SORT}" "${TMP_SORT}.bai" "${TMP_VIT}" "${TMP_VIT}.bai" "${TMP_INDEL}"
    fi

    if [ ! -s "${VCF}" ] || [ ! -s "${TBI}" ]; then
        TMP_VCF="${OUT}/${s}.lofreq.tmp.vcf"
        lofreq call --call-indels -f "${REF}" -o "${TMP_VCF}" "${BAM}"
        bcftools sort -O z -o "${VCF}" "${TMP_VCF}"
        tabix -f -p vcf "${VCF}"
        rm -f "${TMP_VCF}"
    fi
done

COLLAPSED="${OUT}/collapsed.tsv"
if [ ! -s "${COLLAPSED}" ]; then
    TMP_COL="${OUT}/collapsed.tmp.tsv"
    {
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for s in "${SAMPLES[@]}"; do
            bcftools norm -m- -f "${REF}" "${OUT}/${s}.vcf.gz" 2>/dev/null \
                | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
                | awk -v sm="${s}" 'BEGIN{OFS="\t"}{print sm,$1,$2,$3,$4,$5}'
        done
    } > "${TMP_COL}"
    mv "${TMP_COL}" "${COLLAPSED}"
fi