#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results
for sample in "${SAMPLES[@]}"; do
    if [[ -f data/ref/chrM.fa.bwt ]]; then
        continue
    fi
    [[ -f data/ref/chrM.fa.bwt ]] || bwa index data/ref/chrM.fa
    [[ -f data/ref/chrM.fa.fai ]] || samtools faidx data/ref/chrM.fa
    if [[ -f results/${sample}.bam ]]; then
        continue
    fi
    { bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ 4 -o results/${sample}.bam - } || true
    [[ -f results/${sample}.bam.bai ]] || samtools index -@ 4 results/${sample}.bam
    if [[ -f results/${sample}.vcf.gz.tbi ]]; then
        continue
    fi
    [[ -f results/${sample}.vcf.gz.tbi ]] || { [[ -f results/${sample}.vcf ]] || lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam; bgzip -f results/${sample}.vcf && tabix -p vcf results/${sample}.vcf.gz; }
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | awk -v s=${sample} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done