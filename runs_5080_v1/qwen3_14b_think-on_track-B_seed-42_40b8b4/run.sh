#!/usr/bin/env bash
set -euo pipefail

if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for sample in "${samples[@]}"; do
    bam="results/${sample}.bam"
    bai="results/${sample}.bam.bai"
    vcf="results/${sample}.vcf.gz"
    tbi="results/${sample}.vcf.gz.tbi"

    if [ ! -f "$bam" ]; then
        bwa mem -t 4 data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools view -Sb - | \
        samtools sort -o "$bam" - && \
        samtools index "$bam" -@ 4
    fi

    if [ ! -f "$vcf" ]; then
        bcftools mpileup -f data/ref/chrM.fa "$bam" | \
        bcftools call -mv -O z -o "$vcf" --ploidy 1 -@ 4
        tabix -p vcf "$vcf"
    fi
done

if [ ! -f results/collapsed.tsv ]; then
    > results/collapsed.tsv
    for sample in "${samples[@]}"; do
        vcf="results/${sample}.vcf.gz"
        bcftools query -f '%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$vcf" >> results/collapsed.tsv
    done
fi