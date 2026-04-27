#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai || ! -f data/ref/chrM.fa.amb ]]; then
    samtools faidx data/ref/chrM.fa
    bwa index data/ref/chrM.fa
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam_file="results/${sample}.bam"
    vcf_file="results/${sample}.vcf.gz"
    tbi_file="results/${sample}.vcf.gz.tbi"

    # Skip if output exists and is up-to-date
    if [[ -f "$bam_file" && -f "$vcf_file" && -f "$tbi_file" ]]; then
        if [[ "$bam_file" -nt "$vcf_file" || "$tbi_file" -nt "$vcf_file" ]]; then
            continue
        fi
    fi

    # Alignment
    if [[ ! -f "$bam_file" ]]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$bam_file"
    fi

    # Indexing
    if [[ ! -f "results/${sample}.bam.bai" ]]; then
        samtools index -@ "$THREADS" "$bam_file"
    fi

    # Variant calling
    if [[ ! -f "$vcf_file" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" \
            -f data/ref/chrM.fa \
            -i "$bam_file" \
            -o "results/${sample}.vcf"
    fi

    # Compression and indexing
    if [[ ! -f "$vcf_file" ]]; then
        bgzip -f "results/${sample}.vcf"
        tabix -p vcf "results/${sample}.vcf.gz"
        rm -f "results/${sample}.vcf"
    fi
done

# Collapsed TSV
collapsed_tsv="results/collapsed.tsv"
if [[ ! -f "$collapsed_tsv" || "$(find results/*.vcf.gz -newer "$collapsed_tsv" 2>/dev/null)" ]]; then
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
            -s "results/${sample}.vcf.gz" \
            -o "results/${sample}.tmp" || true
        cat "results/${sample}.tmp" >> "$collapsed_tsv"
        rm -f "results/${sample}.tmp"
    done

    # Add header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed_tsv"
    tail -n +2 results/collapsed.tsv >> "$collapsed_tsv"
    rm -f results/collapsed.tsv
    mv results/collapsed.tsv.tmp "$collapsed_tsv"
fi