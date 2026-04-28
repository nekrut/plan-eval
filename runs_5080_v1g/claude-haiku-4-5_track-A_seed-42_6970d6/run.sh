#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Prepare results directory
mkdir -p results

# Reference indexing (once, in data/ref/)
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

if [ ! -f data/ref/chrM.fa.amb ]; then
    bwa index data/ref/chrM.fa
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    # Alignment: bwa mem → samtools sort → BAM
    if [ ! -f results/"${sample}".bam ]; then
        bwa mem -t "${THREADS}" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa \
            data/raw/"${sample}"_1.fq.gz \
            data/raw/"${sample}"_2.fq.gz | \
            samtools sort -@ "${THREADS}" -o results/"${sample}".bam
    fi

    # BAM indexing
    if [ ! -f results/"${sample}".bam.bai ]; then
        samtools index -@ "${THREADS}" results/"${sample}".bam
    fi

    # Variant calling with lofreq
    if [ ! -f results/"${sample}".vcf.gz.tbi ]; then
        rm -f results/"${sample}".vcf results/"${sample}".vcf.gz results/"${sample}".vcf.gz.tbi

        lofreq call-parallel --pp-threads "${THREADS}" --verbose \
            --ref data/ref/chrM.fa \
            --out results/"${sample}".vcf \
            --sig \
            --bonf \
            results/"${sample}".bam

        # Compress VCF
        bgzip results/"${sample}".vcf

        # Index compressed VCF
        tabix -p vcf results/"${sample}".vcf.gz
    fi
done

# Collapse step
collapse_needed=false
if [ ! -f results/collapsed.tsv ]; then
    collapse_needed=true
else
    for sample in "${SAMPLES[@]}"; do
        if [ -f results/"${sample}".vcf.gz ] && [ results/"${sample}".vcf.gz -nt results/collapsed.tsv ]; then
            collapse_needed=true
            break
        fi
    done
fi

if [ "$collapse_needed" = true ]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/"${sample}".vcf.gz
        done
    } > results/collapsed.tsv
fi