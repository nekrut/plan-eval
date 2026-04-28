#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results
bwa index data/ref/chrM.fa
samtools faidx data/ref/chrM.fa
for sample in "${SAMPLES[@]}"; do
    bwa mem -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam -
    samtools index -@ $THREADS results/${sample}.bam
    lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    bgzip -f results/${sample}.vcf
    tabix -p vcf results/${sample}.vcf.gz
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$sample" "$sample" "$sample" "$sample" "$sample" >> results/collapsed.tsv
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done