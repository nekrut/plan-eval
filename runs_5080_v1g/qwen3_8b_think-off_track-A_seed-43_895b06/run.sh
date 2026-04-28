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

    if [[ -f "$output_bai" && -f "$output_vcf_gz" && -f "$output_vcf_tbi" ]]; then
        continue
    fi

    # Alignment
    if [ ! -f "$output_bam" ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$output_bam"
    fi

    # Index BAM
    if [ ! -f "$output_bai" ]; then
        samtools index -@ "$THREADS" "$output_bam"
    fi

    # Variant calling
    if [ ! -f "$output_vcf" ]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref data/ref/chrM.fa --out "$output_vcf" \
            --sig --bonf "$output_bam"
    fi

    # Compress and index VCF
    if [ ! -f "$output_vcf_gz" ]; then
        bgzip "$output_vcf"
        tabix -p vcf "$output_vcf_gz"
    fi

    # Collapse variants
    if [ ! -f "$output_collapsed" ]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
            "$output_vcf_gz" > "$output_collapsed"
    fi
done

# Final collapse step
if [ ! -f "results/collapsed.tsv" ]; then
    cat results/*.vcf.gz | bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > results/collapsed.tsv
fi

# Ensure header is present
if [ ! -s "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi

# Final check for idempotency
for sample in "${samples[@]}"; do
    if [ ! -f "results/${sample}.bam.bai" ]; then
        samtools index -@ "$THREADS" "results/${sample}.bam"
    fi
    if [ ! -f "results/${sample}.vcf.gz.tbi" ]; then
        tabix -p vcf "results/${sample}.vcf.gz"
    fi
done

# Final check for collapsed.tsv
if [ ! -s "results/collapsed.tsv" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
fi