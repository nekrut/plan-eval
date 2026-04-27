#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

# Index reference
if [ ! -f results/chrM.fa.bai ]; then
    samtools faidx data/ref/chrM.fa -o results/chrM.fa.fai
    samtools index data/ref/chrM.fa -o results/chrM.fa.bai
fi

# Function to process a sample
process_sample() {
    local sample=$1
    local fq1="data/raw/${sample}_1.fq.gz"
    local fq2="data/raw/${sample}_2.fq.gz"
    local bam="results/${sample}.bam"
    local bai="results/${sample}.bam.bai"
    local vcf="results/${sample}.vcf.gz"
    local tbi="results/${sample}.vcf.gz.tbi"

    # Mapping
    if [ ! -f "$bai" ]; then
        bwa mem -t 4 data/ref/chrM.fa "$fq1" "$fq2" | \
            samtools view -b - | \
            samtools sort -@ 2 -o "$bam"
        samtools index "$bam" -o "$bai"
    fi

    # Variant Calling
    if [ ! -f "$tbi" ]; then
        lofreq call-parallel -f data/ref/chrM.fa -m lofreq -t 4 "$bam" -o "$vcf"
        tabix -p vcf "$vcf"
    fi
}

# Process all samples
process_sample "M117-bl"
process_sample "M117-ch"
process_sample "M117C1-bl"
process_sample "M117C1-ch"

# Collapsing VCFs
if [ ! -f results/collapsed.tsv ]; then
    bcftools concat -Oz -o results/merged.vcf.gz results/M117-bl.vcf.gz results/M117-ch.vcf.gz results/M117C1-bl.vcf.gz results/M117C1-ch.vcf.gz
    tabix -p vcf results/merged.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' results/merged.vcf.gz > results/collapsed.tsv
    rm results/merged.vcf.gz results/merged.vcf.gz.tbi
fi