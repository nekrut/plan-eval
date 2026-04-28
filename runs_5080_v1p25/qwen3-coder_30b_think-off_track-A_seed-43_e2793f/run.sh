#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"

mkdir -p results
cd data/ref
REF=chrM.fa
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "${REF}"
    bwa index "${REF}"
fi
cd ../..

for sample in $SAMPLES; do
    BAM="results/${sample}.bam"
    BAI="results/${sample}.bam.bai"
    VCF="results/${sample}.vcf"
    VCF_GZ="results/${sample}.vcf.gz"
    VCF_GZ_TBI="results/${sample}.vcf.gz.tbi"

    if [[ -f "${VCF_GZ_TBI}" ]] && [[ "${VCF_GZ_TBI}" -nt "results/${sample}.bam" ]]; then
        continue
    fi

    if [[ ! -f "${BAM}" ]] || [[ "${BAM}" -nt "data/raw/${sample}_1.fq.gz" ]] || [[ "${BAM}" -nt "data/raw/${sample}_2.fq.gz" ]]; then
        bwa mem -t "${THREADS}" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa \
            data/raw/"${sample}"_1.fq.gz \
            data/raw/"${sample}"_2.fq.gz \
            | samtools sort -@ "${THREADS}" -o "${BAM}"
    fi

    if [[ ! -f "${BAI}" ]] || [[ "${BAI}" -nt "${BAM}" ]]; then
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    if [[ ! -f "${VCF}" ]] || [[ "${VCF}" -nt "${BAM}" ]]; then
        lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o "${VCF}" "${BAM}"
    fi

    if [[ ! -f "${VCF_GZ}" ]] || [[ "${VCF_GZ}" -nt "${VCF}" ]]; then
        bgzip -c "${VCF}" > "${VCF_GZ}"
        tabix -p vcf "${VCF_GZ}"
        rm -f "${VCF}"
    fi
done

if [[ ! -f results/collapsed.tsv ]] || [[ results/collapsed.tsv -nt results/M117-bl.vcf.gz.tbi ]] || [[ results/collapsed.tsv -nt results/M117-ch.vcf.gz.tbi ]] || [[ results/collapsed.tsv -nt results/M117C1-bl.vcf.gz.tbi ]] || [[ results/collapsed.tsv -nt results/M117C1-ch.vcf.gz.tbi ]]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
    for sample in $SAMPLES; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz" >> results/collapsed.tsv
    done
fi