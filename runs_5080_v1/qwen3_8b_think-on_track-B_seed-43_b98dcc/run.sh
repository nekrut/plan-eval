#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
samtools faidx data/ref/chrM.fa

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    if [ -f "results/${sample}.bam" ]; then
        continue
    fi
    bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools view -S -b | \
        samtools sort -o results/${sample}.bam
    samtools index results/${sample}.bam
    lofreq call -f data/ref/chrM.fa -o results/${sample}.vcf.gz results/${sample}.bam
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' results/${sample}.vcf.gz | \
        awk '{print "'${sample}'\t"$0}' > results/${sample}_collapsed.tsv
done

cat results/*.collapsed.tsv | sort -k1,1 -k2,2n | uniq > results/collapsed.tsv