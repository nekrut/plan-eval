#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.sai ]; then
    bwa index data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

# Process samples
for sample in "${SAMPLES[@]}"; do
    # Alignment and sorting
    if [ ! -f results/"${sample}".bam ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/"${sample}"_1.fq.gz data/raw/"${sample}"_2.fq.gz | samtools sort -@ "$THREADS" -o results/"${sample}".bam -
    fi

    # BAM index
    if [ ! -f results/"${sample}".bam.bai ]; then
        samtools index -@ "$THREADS" results/"${sample}".bam
    fi

    # Variant calling
    if [ ! -f results/"${sample}".vcf; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/"${sample}".vcf results/"${sample}".bam
    fi

    # VCF compression and tabix
    if [ ! -f results/"${sample}".vcf.gz.tbi ]; then
        # Use bgzip if available, otherwise use gzip (bgzip is standard for vcf)
        if command -v bgzip > /dev/null; then
            bgzip -f results/"${sample}".vcf
        else
            gzip -c results/"${sample}".vcf > results/"${sample}".vcf.gz
            mv results/"${sample}".vcf.gz.tmp results/"${sample}".vcf.gz 2>/dev/null || mv results/"${sample}".vcf.gz_original results/"${sample}".vcf.gz
        fi
        # Note: Since bgzip is standard in the provided toolset context for vcf, we assume it's available or handled by samtools logic
        # To ensure strict adherence to tools:
        if [ ! -f results/"${sample}".vcf.gz ]; then
            gzip -c results/"${sample}".vcf > results/"${sample}".vcf.gz
        fi
        tabix -p vcf results/"${sample}".vcf.gz
    fi
done

# Create collapsed table header if it doesn't exist or is empty
if [ ! -s results/collapsed.tsv ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
fi

# Append data to collapsed table (only for new entries)
# To ensure idempotency and avoid duplicates, we check if the specific sample entry exists or just rebuild the file
# Given the requirement "be idempotent", we will construct the final table from all vcf files.
> results/collapsed.tsv
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

for sample in "${SAMPLES[@]}"; do
    if [ -f results/"${sample}".vcf.gz ]; then
        bcftools query -l results/"${sample}".vcf.gz | xargs -I {} bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/"${sample}".vcf.gz | awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
    fi
done