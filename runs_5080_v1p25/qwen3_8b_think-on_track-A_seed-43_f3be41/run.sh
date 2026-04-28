#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi
if [[ ! -f data/ref/chrM.fa.amb ]]; then
    bwa index data/ref/chrM.fa
fi

# Per-sample alignment and processing
for sample in "${SAMPLES[@]}"; do
    bam_file="results/${sample}.bam"
    vcf_file="results/${sample}.vcf"
    vcf_gz="results/${sample}.vcf.gz"
    vcf_tbi="results/${sample}.vcf.gz.tbi"

    # Skip if output exists and is up-to-date
    if [[ -f "$vcf_tbi" && -nt "$vcf_tbi" "$bam_file" ]]; then
        continue
    fi

    # Alignment
    if [[ ! -f "$bam_file" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$bam_file"
    fi

    # Indexing
    if [[ ! -f "results/${sample}.bam.bai" ]]; then
        samtools index -@ "$THREADS" "$bam_file"
    fi

    # Variant calling
    if [[ ! -f "$vcf_file" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$vcf_file" "$bam_file"
    fi

    # Compression and indexing
    if [[ ! -f "$vcf_gz" ]]; then
        bgzip "$vcf_file"
        tabix -p vcf "$vcf_gz"
        rm "$vcf_file"
    fi
done

# Collapsed TSV
collapsed_tsv="results/collapsed.tsv"
if [[ ! -f "$collapsed_tsv" || -nt "$collapsed_tsv" "${SAMPLES[@]/#results/%}" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed_tsv"
    for sample in "${SAMPLES[@]}"; do
        vcf_gz="results/${sample}.vcf.gz"
        if [[ -f "$vcf_gz" ]]; then
            bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$vcf_gz" >> "$collapsed_tsv"
        fi
    done
fi