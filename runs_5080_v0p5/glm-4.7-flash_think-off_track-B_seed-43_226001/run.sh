#!/usr/bin/env bash
set -euo pipefail

THREADS=4
mkdir -p results

# Index reference
if [ ! -f data/ref/chrM.fa.bwt ]; then
    bwa index data/ref/chrM.fa
fi

# Function to process a sample
process_sample() {
    local sample=$1
    local fq1="data/raw/${sample}_1.fq.gz"
    local fq2="data/raw/${sample}_2.fq.gz"
    local bam="results/${sample}.bam"
    local vcf="results/${sample}.vcf.gz"

    # Mapping
    if [ ! -f "$bam" ]; then
        bwa mem -t "$THREADS" data/ref/chrM.fa "$fq1" "$fq2" | \
        samtools view -b - > "$bam"
        samtools index "$bam"
    fi

    # Variant calling
    if [ ! -f "$vcf" ]; then
        lofreq call-parallel -f data/ref/chrM.fa -o "$vcf" -m lofreq -t "$THREADS" "$bam"
        tabix -p vcf "$vcf"
    fi
}

# Process all samples
process_sample "M117-bl"
process_sample "M117-ch"
process_sample "M117C1-bl"
process_sample "M117C1-ch"

# Collapse VCFs
if [ ! -f results/collapsed.tsv ]; then
    bcftools concat -Oz -o results/collapsed.vcf.gz \
        results/M117-bl.vcf.gz results/M117-ch.vcf.gz \
        results/M117C1-bl.vcf.gz results/M117C1-ch.vcf.gz
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' results/collapsed.vcf.gz > results/collapsed.tsv
    tabix -p vcf results/collapsed.vcf.gz
fi

exit 0