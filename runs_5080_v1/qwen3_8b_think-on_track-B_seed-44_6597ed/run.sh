#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
samtools faidx data/ref/chrM.fa

if [ ! -f results/collapsed.tsv ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    if [ -f results/${sample}.bam ]; then
        echo "Skipping ${sample} as output already exists."
        continue
    fi
    bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -o results/${sample}.bam -T tmp
    samtools index results/${sample}.bam
    lofreq call -f data/ref/chrM.fa -o results/${sample}.vcf.gz results/${sample}.bam
    bgzip results/${sample}.vcf.gz
    tabix -p vcf results/${sample}.vcf.gz
done

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    if [ ! -f results/${sample}.vcf.gz ]; then
        continue
    fi
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | awk -F'\t' '{print "'${sample}'"\t$1\t$2\t$3\t$4\t$5}' >> results/collapsed.tsv
done