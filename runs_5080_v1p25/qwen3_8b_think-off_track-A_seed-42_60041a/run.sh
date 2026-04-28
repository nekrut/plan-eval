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

for sample in "${samples[@]}"; do
    output_bam="results/${sample}.bam"
    output_bai="results/${sample}.bam.bai"
    output_vcf="results/${sample}.vcf"
    output_vcf_gz="results/${sample}.vcf.gz"
    output_vcf_tbi="results/${sample}.vcf.gz.tbi"
    output_collapsed="results/collapsed.tsv"

    if [ -f "$output_bai" ] && [ "$output_bai" -nt data/raw/"${sample}_1.fq.gz" ] && [ "$output_bai" -nt data/raw/"${sample}_2.fq.gz" ]; then
        continue
    fi

    if [ -f "$output_vcf_tbi" ] && [ "$output_vcf_tbi" -nt "$output_vcf_gz" ]; then
        continue
    fi

    if [ -f "$output_vcf_gz" ] && [ "$output_vcf_gz" -nt data/raw/"${sample}_1.fq.gz" ] && [ "$output_vcf_gz" -nt data/raw/"${sample}_2.fq.gz" ]; then
        continue
    fi

    # Alignment
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$output_bam"

    # Index BAM
    samtools index -@ "$THREADS" "$output_bam"

    # Variant calling
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$output_vcf" "$output_bam"

    # Compress and index VCF
    bgzip "$output_vcf"
    tabix -p vcf "$output_vcf_gz"
    rm "$output_vcf"

    # Collapse variants
    if [ ! -f "$output_collapsed" ] || [ "$output_collapsed" -nt "$output_vcf_gz" ]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' -s "$sample" "$output_vcf_gz" >> "$output_collapsed"
    fi
done

# Final header for collapsed.tsv
if [ ! -f "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

# Ensure header is present
if [ ! -s "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

# Append variants if needed
for sample in "${samples[@]}"; do
    if [ -f "results/${sample}.vcf.gz" ]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' -s "$sample" "results/${sample}.vcf.gz" >> results/collapsed.tsv
    fi
done