#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
    bwa index data/ref/chrM.fa
fi

# Per-sample alignment and processing
for sample in "${SAMPLES[@]}"; do
    # Skip if outputs already exist
    if [[ -f results/"${sample}".vcf.gz.tbi ]]; then
        continue
    fi

    # Alignment
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
        data/ref/chrM.fa \
        data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
    samtools sort -@ "$THREADS" -o results/"${sample}".bam

    # Indexing
    samtools index -@ "$THREADS" results/"${sample}".bam

    # Variant calling
    lofreq call-parallel --pp-threads "$THREADS" --verbose \
        --ref data/ref/chrM.fa \
        --out results/"${sample}".vcf \
        --sig --bonf \
        results/"${sample}".bam

    # Compression and indexing
    bgzip results/"${sample}".vcf
    tabix -p vcf results/"${sample}".vcf.gz
    rm results/"${sample}".vcf
done

# Collapse variants
if [[ ! -f results/collapsed.tsv ]]; then
    cat /dev/null > results/collapsed.tsv
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
            results/"${sample}".vcf.gz >> results/collapsed.tsv
    done
    echo -e "sample\tchrom\tpos\tref\talt\taf" | tee results/collapsed.tsv
fi

# Final idempotency check
for sample in "${SAMPLES[@]}"; do
    if [[ -f results/"${sample}".vcf.gz.tbi ]]; then
        touch results/"${sample}".vcf.gz.tbi
    fi
done
touch results/collapsed.tsv