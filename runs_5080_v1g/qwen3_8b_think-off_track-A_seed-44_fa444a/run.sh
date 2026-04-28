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

    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
    samtools sort -@ "$THREADS" -o results/"${sample}".bam

    samtools index -@ "$THREADS" results/"${sample}".bam
done

# Variant calling with lofreq
for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".vcf.gz ]; then
        continue
    fi

    lofreq call-parallel --pp-threads "$THREADS" --verbose \
        --ref data/ref/chrM.fa --out results/"${sample}".vcf \
        --sig --bonf results/"${sample}".bam
done

# VCF compression and indexing
for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".vcf.gz ]; then
        continue
    fi

    bgzip results/"${sample}".vcf
    tabix -p vcf results/"${sample}".vcf.gz
    rm results/"${sample}".vcf
done

# Collapse variants into tsv
if [ ! -f results/collapsed.tsv ]; then
    cat <(echo "sample\tchrom\tpos\tref\talt\taf") \
        $(for sample in "${samples[@]}"; do
            bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
                results/"${sample}".vcf.gz
        done) > results/collapsed.tsv
fi

# Idempotency check
for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".bam ]; then
        if [ -n "$(find data/ref/ -newer results/"${sample}".bam)" ]; then
            samtools index -@ "$THREADS" results/"${sample}".bam
        fi
    fi
done

for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".vcf.gz ]; then
        if [ -n "$(find data/ref/ -newer results/"${sample}".vcf.gz)" ]; then
            lofreq call-parallel --pp-threads "$THREADS" --verbose \
                --ref data/ref/chrM.fa --out results/"${sample}".vcf \
                --sig --bonf results/"${sample}".bam
        fi
    fi
done

for sample in "${samples[@]}"; do
    if [ -f results/"${sample}".vcf.gz ]; then
        if [ -n "$(find data/ref/ -newer results/"${sample}".vcf.gz)" ]; then
            bgzip results/"${sample}".vcf
            tabix -p vcf results/"${sample}".vcf.gz
            rm results/"${sample}".vcf
        fi
    fi
done

if [ -f results/collapsed.tsv ]; then
    if [ -n "$(find data/ref/ -newer results/collapsed.tsv)" ]; then
        cat <(echo "sample\tchrom\tpos\tref\talt\taf") \
            $(for sample in "${samples[@]}"; do
                bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
                    results/"${sample}".vcf.gz
            done) > results/collapsed.tsv
    fi
fi