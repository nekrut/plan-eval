#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="data/raw"
REF_DIR="data/ref"
RESULTS_DIR="results"

REF_FASTA="$REF_DIR/chrM.fa"
THREADS=4

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p "$RESULTS_DIR"

# Check if all outputs exist
ALL_EXIST=true
for sample in "${SAMPLES[@]}"; do
    bam="$RESULTS_DIR/${sample}.bam"
    bai="$RESULTS_DIR/${sample}.bam.bai"
    vcf="$RESULTS_DIR/${sample}.vcf.gz"
    tbi="$RESULTS_DIR/${sample}.vcf.gz.tbi"
    if [[ ! -f "$bam" || ! -f "$bai" || ! -f "$vcf" || ! -f "$tbi" ]]; then
        ALL_EXIST=false
        break
    fi
done
if [[ "$ALL_EXIST" = true && -f "$RESULTS_DIR/collapsed.tsv" ]]; then
    exit 0
fi

# Index reference if needed
if [[ ! -f "$REF_FASTA.bwt" ]]; then
    bwa index "$REF_FASTA"
fi
if [[ ! -f "$REF_FASTA.fai" ]]; then
    samtools faidx "$REF_FASTA"
fi

# Process each sample
for sample in "${SAMPLES[@]}"; do
    bam="$RESULTS_DIR/${sample}.bam"
    bai="$RESULTS_DIR/${sample}.bam.bai"
    vcf="$RESULTS_DIR/${sample}.vcf.gz"
    tbi="$RESULTS_DIR/${sample}.vcf.gz.tbi"

    if [[ ! -f "$bam" ]]; then
        fq1="$RAW_DIR/${sample}_1.fq.gz"
        fq2="$RAW_DIR/${sample}_2.fq.gz"
        bwa mem -t "$THREADS" -M -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "$REF_FASTA" "$fq1" "$fq2" |
            samtools view -Sb -@ "$THREADS" - |
            samtools sort -@ "$THREADS" -o "$bam" -
    fi

    if [[ ! -f "$bai" ]]; then
        samtools index "$bam"
    fi

    if [[ ! -f "$vcf" ]]; then
        lofreq call -f "$REF_FASTA" -b "$bam" -o "$vcf" -q 20 -Q 20 -p "$THREADS" -s "$sample"
    fi

    if [[ ! -f "$tbi" ]]; then
        tabix -p vcf "$vcf"
    fi
done

# Create collapsed.tsv
collapsed="$RESULTS_DIR/collapsed.tsv"
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed"
bcftools query -f '%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$RESULTS_DIR"/*.vcf.gz >> "$collapsed"

exit 0