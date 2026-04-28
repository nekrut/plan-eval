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

for sample in "${samples[@]}"; do
    output_bam="results/${sample}.bam"
    output_bai="results/${sample}.bam.bai"
    output_vcf="results/${sample}.vcf"
    output_vcf_gz="results/${sample}.vcf.gz"
    output_vcf_tbi="results/${sample}.vcf.gz.tbi"
    output_collapsed="results/collapsed.tsv"

    if [[ -f "$output_bai" && -f "$output_vcf_gz" && -f "$output_vcf_tbi" && -f "$output_collapsed" ]]; then
        continue
    fi

    # Alignment
    if [ ! -f "$output_bam" ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$output_bam"
    fi

    # Index BAM
    if [ ! -f "$output_bai" ]; then
        samtools index -@ "$THREADS" "$output_bam"
    fi

    # Variant calling
    if [ ! -f "$output_vcf" ]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$output_vcf" "$output_bam"
    fi

    # Compress and index VCF
    if [ ! -f "$output_vcf_gz" ]; then
        bgzip "$output_vcf"
        tabix -p vcf "$output_vcf_gz"
        rm "$output_vcf"
    fi

    # Collapse variants
    if [ ! -f "$output_collapsed" ]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$output_vcf_gz" > "$output_collapsed"
    fi
done

# Final collapse step
if [ ! -f "results/collapsed.tsv" ]; then
    cat results/*.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > results/collapsed.tsv
fi

# Add header to collapsed.tsv
if [ ! -f "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

# Ensure header is present
if ! grep -q '^sample\tchrom\tpos\tref\talt\taf$' results/collapsed.tsv; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" >> results/collapsed.tsv
fi

# Final check for idempotency
if [[ -f "$output_bai" && -f "$output_vcf_gz" && -f "$output_vcf_tbi" && -f "$output_collapsed" ]]; then
    exit 0
fi