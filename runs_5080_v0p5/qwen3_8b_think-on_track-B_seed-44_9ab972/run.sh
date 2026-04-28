#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

if [ ! -f data/ref/chrM.fa.bai ]; then
    samtools index data/ref/chrM.fa
fi

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    if [ -f results/${sample}.bam ]; then
        continue
    fi

    bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools sort -o results/${sample}.bam

    if [ ! -f results/${sample}.bam.bai ]; then
        samtools index results/${sample}.bam
    fi

    if [ ! -f results/${sample}.vcf.gz ]; then
        lofreq call -f data/ref/chrM.fa -o results/${sample}.vcf.gz results/${sample}.bam
    fi

    if [ ! -f results/${sample}.tsv ]; then
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
            awk -v sample="$sample" '{print sample, $1, $2, $3, $4, $5}' > results/${sample}.tsv
    fi
done

if [ ! -f results/collapsed.tsv ]; then
    cat results/*.tsv | awk 'NR==1 {print "sample\tchrom\tpos\tref\talt\taf"} {print}' > results/collapsed.tsv
fi