#!/usr/bin/env bash
set -euo pipefail

THREADS=4
RAW=data/raw
REF_DIR=data/ref
REF=${REF_DIR}/chrM.fa
OUT=results
mkdir -p "${OUT}"

SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

if [ ! -f "${REF}.bwt" ]; then
    bwa index "${REF}"
fi
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

for sample in "${SAMPLES[@]}"; do
    bam="${OUT}/${sample}.bam"
    bai="${OUT}/${sample}.bam.bai"
    vcf="${OUT}/${sample}.vcf.gz"
    tbi="${OUT}/${sample}.vcf.gz.tbi"
    r1="${RAW}/${sample}_1.fq.gz"
    r2="${RAW}/${sample}_2.fq.gz"

    if [ ! -f "${bam}" ] || [ ! -f "${bai}" ]; then
        tmp_bam="${OUT}/${sample}.unindel.bam"
        bwa mem -t "${THREADS}" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "${REF}" "${r1}" "${r2}" \
            | samtools sort -@ "${THREADS}" -o "${tmp_bam}" -
        samtools index -@ "${THREADS}" "${tmp_bam}"
        lofreq indelqual --dindel -f "${REF}" -o "${bam}" "${tmp_bam}"
        samtools index -@ "${THREADS}" "${bam}"
        rm -f "${tmp_bam}" "${tmp_bam}.bai"
    fi

    if [ ! -f "${vcf}" ] || [ ! -f "${tbi}" ]; then
        raw_vcf="${OUT}/${sample}.lofreq.vcf"
        lofreq call-parallel --pp-threads "${THREADS}" --call-indels -f "${REF}" -o "${raw_vcf}" "${bam}"
        bcftools view "${raw_vcf}" -Oz -o "${vcf}"
        tabix -p vcf "${vcf}"
        rm -f "${raw_vcf}"
    fi
done

collapsed="${OUT}/collapsed.tsv"
if [ ! -f "${collapsed}" ]; then
    tmp_collapsed="${OUT}/collapsed.tsv.tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${tmp_collapsed}"
    for sample in "${SAMPLES[@]}"; do
        vcf="${OUT}/${sample}.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${vcf}" \
            | awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$1,$2,$3,$4,$5}' >> "${tmp_collapsed}"
    done
    mv "${tmp_collapsed}" "${collapsed}"
fi