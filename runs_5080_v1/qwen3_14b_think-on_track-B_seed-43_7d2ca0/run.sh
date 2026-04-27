#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

if [ -f "results/collapsed.tsv" ]; then
    exit 0
fi

if [ ! -f "data/ref/chrM.fa.fai" ]; then
    samtools faidx data/ref/chrM.fa
fi
if [ ! -f "data/ref/chrM.fa.bwt" ]; then
    bwa index data/ref/chrM.fa
fi

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    if [ -f "results/${sample}.bam" ] && [ -f "results/${sample}.bam.bai" ] && [ -f "results/${sample}.vcf.gz" ] && [ -f "results/${sample}.vcf.gz.tbi" ]; then
        continue
    fi

    if [ ! -f "results/${sample}.bam" ] || [ ! -f "results/${sample}.bam.bai" ]; then
        bwa mem -R "@RG\tID:${sample}\tSM:${sample}" data/ref/chrM.fa \
            data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
            samtools view -bS - | \
            samtools sort -o results/${sample}.bam
        samtools index results/${sample}.bam
    fi

    if [ ! -f "results/${sample}.vcf.gz" ] || [ ! -f "results/${sample}.vcf.gz.tbi" ]; then
        lofreq call -f data/ref/chrM.fa results/${sample}.bam > results/${sample}.vcf
        bgzip results/${sample}.vcf
        tabix -p vcf results/${sample}.vcf.gz
    fi

    bcftools query -s $sample -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
        awk -v sample="$sample" '{print sample, $1, $2, $3, $4, $5}' >> results/collapsed.tmp
done

sort -k2,2n -k3,3n results/collapsed.tmp | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6}' > results/collapsed.tsv
rm results/collapsed.tmp