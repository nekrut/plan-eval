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

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    # Skip if outputs already exist and are up-to-date
    if [[ -f results/"${sample}".vcf.gz.tbi && -f results/"${sample}".bam && -f results/"${sample}".bam.bai ]]; then
        continue
    fi

    # Alignment
    if [[ ! -f results/"${sample}".bam ]]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o results/"${sample}".bam
    fi

    # Indexing
    if [[ ! -f results/"${sample}".bam.bai ]]; then
        samtools index -@ "$THREADS" results/"${sample}".bam
    fi

    # Variant calling
    if [[ ! -f results/"${sample}".vcf.gz ]]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref data/ref/chrM.fa \
            --out results/"${sample}".vcf \
            --sig --bonf \
            results/"${sample}".bam
        bgzip results/"${sample}".vcf
        tabix -p vcf results/"${sample}".vcf.gz
        rm results/"${sample}".vcf
    fi
done

# Collapse variants
if [[ ! -f results/collapsed.tsv ]]; then
    cat /dev/null > results/collapsed.tsv
    echo -e "sample\tchrom\tpos\tref\talt\taf" >> results/collapsed.tsv
    for sample in "${SAMPLES[@]}"; do
        if [[ -f results/"${sample}".vcf.gz ]]; then
            bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
                --sample "${sample}" results/"${sample}".vcf.gz >> results/collapsed.tsv
        fi
    done
fi

# Final idempotency check
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/"${sample}".vcf.gz.tbi ]]; then
        tabix -p vcf results/"${sample}".vcf.gz
    fi
done

exit 0