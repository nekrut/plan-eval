#!/usr/bin/env bash
set -euo pipefail

# Index reference
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

# Process each sample
for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    # Alignment and sorting
    if [ ! -f results/${sample}.bam ]; then
        mkdir -p results
        bwa mem data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
        samtools view -Sb - | \
        samtools sort -o results/${sample}.bam
        samtools index results/${sample}.bam
    fi

    # Variant calling
    if [ ! -f results/${sample}.vcf.gz ]; then
        lofreq call -f data/ref/chrM.fa -o results/${sample}.vcf.gz results/${sample}.bam
        tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Generate collapsed TSV
if [ ! -f results/collapsed.tsv ]; then
    > results/collapsed.tsv
    for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
        awk -v sample=${sample} '{print sample, $1, $2, $3, $4, $5}' >> results/collapsed.tsv
    done
fi