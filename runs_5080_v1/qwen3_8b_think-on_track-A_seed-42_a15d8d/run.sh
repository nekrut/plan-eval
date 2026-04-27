#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi
if [[ ! -f data/ref/chrM.fa.amb ]]; then
    bwa index data/ref/chrM.fa
fi

# Per-sample processing
for sample in "${samples[@]}"; do
    # Skip if outputs already exist and are up-to-date
    if [[ -f results/"${sample}".vcf.gz.tbi && -f results/"${sample}".bam && -f results/"${sample}".bam.bai ]]; then
        continue
    fi

    # Alignment
    if [[ ! -f results/"${sample}".bam ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o results/"${sample}".bam
    fi

    # BAM indexing
    if [[ ! -f results/"${sample}".bam.bai ]]; then
        samtools index -@ "$THREADS" results/"${sample}".bam
    fi

    # Variant calling
    if [[ ! -f results/"${sample}".vcf.gz ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa \
            results/"${sample}".bam > results/"${sample}".vcf
        bgzip results/"${sample}".vcf
        tabix -p vcf results/"${sample}".vcf.gz
        rm results/"${sample}".vcf
    fi
done

# Collapsed TSV
if [[ ! -f results/collapsed.tsv || -n $(find results/*.vcf.gz -newer results/collapsed.tsv) ]]; then
    cat /dev/null > results/collapsed.tsv
    for sample in "${samples[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
            results/"${sample}".vcf.gz >> results/collapsed.tsv
    done
    echo -e "sample\tchrom\tpos\tref\talt\taf" | tee results/collapsed.tsv > /dev/null
fi

# Final idempotency check
for sample in "${samples[@]}"; do
    if [[ ! -f results/"${sample}".vcf.gz.tbi ]]; then
        exit 1
    fi
done
if [[ ! -f results/collapsed.tsv ]]; then
    exit 1
fi

exit 0