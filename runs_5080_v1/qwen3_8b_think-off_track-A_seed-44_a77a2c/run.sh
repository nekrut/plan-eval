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
    output_bam="results/${sample}.bam"
    output_bai="results/${sample}.bam.bai"
    output_vcf="results/${sample}.vcf"
    output_vcf_gz="results/${sample}.vcf.gz"
    output_vcf_tbi="results/${sample}.vcf.gz.tbi"
    output_collapsed="results/collapsed.tsv"

    if [ -f "$output_bai" ] && [ "$output_bai" -nt data/raw/"${sample}_1.fq.gz" ] && [ "$output_bai" -nt data/raw/"${sample}_2.fq.gz" ]; then
        continue
    fi

    if [ ! -f "$output_bam" ] || [ "$output_bam" -nt data/raw/"${sample}_1.fq.gz" ] || [ "$output_bam" -nt data/raw/"${sample}_2.fq.gz" ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$output_bam"
    fi

    if [ ! -f "$output_bai" ] || [ "$output_bai" -nt "$output_bam" ]; then
        samtools index -@ "$THREADS" "$output_bam"
    fi

    if [ ! -f "$output_vcf" ] || [ "$output_vcf" -nt "$output_bam" ]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$output_vcf" "$output_bam"
    fi

    if [ ! -f "$output_vcf_gz" ] || [ "$output_vcf_gz" -nt "$output_vcf" ]; then
        bgzip "$output_vcf"
    fi

    if [ ! -f "$output_vcf_tbi" ] || [ "$output_vcf_tbi" -nt "$output_vcf_gz" ]; then
        tabix -p vcf "$output_vcf_gz"
    fi

    if [ ! -f "$output_collapsed" ] || [ "$output_collapsed" -nt "$output_vcf_gz" ]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' -s "$output_vcf_gz" > "$output_collapsed"
    fi
done

# Final header for collapsed.tsv
if [ ! -f "results/collapsed.tsv" ] || [ "results/collapsed.tsv" -nt "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "results/collapsed.tsv"
fi

# Append all samples' variant data
for sample in "${samples[@]}"; do
    if [ -f "results/${sample}.vcf.gz" ]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' -s "results/${sample}.vcf.gz" >> "results/collapsed.tsv"
    fi
done

# Ensure final header is present
if [ ! -s "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "results/collapsed.tsv"
fi

# Final check for collapsed.tsv
if [ -f "results/collapsed.tsv" ]; then
    if [ "results/collapsed.tsv" -nt "results/collapsed.tsv" ]; then
        echo "Collapsed TSV is up to date."
    fi
fi