#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results
bwa index data/ref/chrM.fa
samtools faidx data/ref/chrM.fa
for sample in "${SAMPLES[@]}"; do
    bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ 4 -o results/${sample}.bam -
    samtools index -@ 4 results/${sample}.bam
    lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    bgzip -f results/${sample}.vcf
    tabix -p vcf results/${sample}.vcf.gz
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${sample}" "chrM" "0" "N" "N" "0" >> results/collapsed.tsv
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done