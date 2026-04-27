#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing (once)
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
    bwa index data/ref/chrM.fa
fi

# Per-sample alignment and processing
for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".bam ]; then
        continue
    fi

    # Alignment
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
    samtools sort -@ "$THREADS" -o results/"${sample}".bam

    # Indexing
    samtools index -@ "$THREADS" results/"${sample}".bam

    # Variant calling
    lofreq call-parallel --pp-threads "$THREADS" \
        -f data/ref/chrM.fa \
        -o results/"${sample}".vcf \
        results/"${sample}".bam

    # VCF compression and indexing
    bgzip results/"${sample}".vcf
    tabix -p vcf results/"${sample}".vcf.gz
    rm results/"${sample}".vcf
done

# Collapse variants
if [ ! -f results/collapsed.tsv ]; then
    for sample in "${samples[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
            results/"${sample}".vcf.gz
    done | \
    awk 'NR == 1 {print "sample\tchrom\tpos\tref\talt\taf"} NR > 1 {print}' > results/collapsed.tsv
fi

# Idempotency check for final outputs
for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".vcf.gz.tbi ]; then
        if [ -n "$(find results/"${sample}".vcf.gz -newer results/collapsed.tsv)" ]; then
            bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
                results/"${sample}".vcf.gz >> results/collapsed.tsv
        fi
    fi
done

# Final check for collapsed.tsv
if [ -f results/collapsed.tsv ]; then
    awk 'NR == 1 {print "sample\tchrom\tpos\tref\talt\taf"} NR > 1 {print}' > temp.tsv
    cat results/collapsed.tsv >> temp.tsv
    mv temp.tsv results/collapsed.tsv
fi