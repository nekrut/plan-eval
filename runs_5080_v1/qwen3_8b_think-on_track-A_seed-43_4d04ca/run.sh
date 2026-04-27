#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
    bwa index data/ref/chrM.fa
fi

# Per-sample alignment
for sample in "${samples[@]}"; do
    if [ ! -f results/"$sample".bam ]; then
        bwa mem -t "$THREADS" -R "@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA" data/raw/"$sample"_1.fq.gz data/raw/"$sample"_2.fq.gz | samtools sort -@ "$THREAD
    fi
done

# BAM indexing
for sample in "${samples[@]}"; do
    if [ ! -f results/"$sample".bam.bai ]; then
        samtools index -@ "$THREADS" results/"$sample".bam
    fi
done

# Variant calling
for sample in "${samples[@]}"; do
    if [ ! -f results/"$sample".vcf ]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/"$sample".vcf results/"$sample".bam
    fi
done

# VCF compression and indexing
for sample in "${samples[@]}"; do
    if [ ! -f results/"$sample".vcf.gz ]; then
        bgzip results/"$sample".vcf
        tabix -p vcf results/"$sample".vcf.gz
        rm results/"$sample".vcf
    fi
done

# Collapse step
if [ ! -f results/collapsed.tsv ]; then
    for sample in "${samples[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/"$sample".vcf.gz
    done | awk 'NR == 1 {print "sample\tchrom\tpos\tref\talt\taf"}1' > results/collapsed.tsv
else
    for sample in "${samples[@]}"; do
        if [ "$(stat -c %Y results/"$sample".vcf.gz)" -gt "$(stat -c %Y results/collapsed.tsv)" ]; then
            for sample in "${samples[@]}"; do
                bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/"$sample".vcf.gz
            done | awk 'NR == 1 {print "sample\tchrom\tpos\tref\talt\taf"}1' > results/collapsed.tsv
            break
        fi
    done
fi