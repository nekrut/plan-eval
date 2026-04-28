#!/usr/bin/env bash
set -euo pipefail

# Define paths relative to current working directory
RAW_DIR="data/raw"
REF_DIR="data/ref"
RES_DIR="results"

# Create results directory if it doesn't exist
mkdir -p "$RES_DIR"

# Define reference genome path
REF_FAI="${REF_DIR}/chrM.fa.fai"
REF_IDX="${REF_DIR}/chrM.fa.fai"
REF_FASTA="${REF_DIR}/chrM.fa"

# Create index for reference if not present
if [[ ! -f "$REF_FAI" ]]; then
    samtools faidx "$REF_FASTA"
fi

# Define samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Function to process a single sample
process_sample() {
    local sample="$1"
    local sample_dir="${RES_DIR}/${sample}"
    local sample_bam="${sample_dir}/${sample}.bam"
    local sample_bai="${sample_dir}/${sample}.bam.bai"
    local sample_vcf="${sample_dir}/${sample}.vcf.gz"
    local sample_vcf_tbi="${sample_dir}/${sample}.vcf.gz.tbi"
    
    # Create sample directory
    mkdir -p "$sample_dir"
    
    # Check if BAM already exists (idempotency)
    if [[ -f "$sample_bam" ]]; then
        return 0
    fi
    
    # Align reads
    bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}\tPL:MiSeq" \
        "$REF_FASTA" \
        "${RAW_DIR}/${sample}_1.fq.gz" \
        "${RAW_DIR}/${sample}_2.fq.gz" \
        | samtools view -b - > "$sample_bam"
    
    # Sort and index BAM
    samtools sort -o "${sample_bam}.tmp" "$sample_bam"
    mv "${sample_bam}.tmp" "$sample_bam"
    samtools index "$sample_bam"
    
    # Call variants with LoFreq
    lofreq call \
        -b "$sample_bam" \
        -o "${sample_dir}/${sample}.vcf" \
        --min-depth 1 \
        --min-fraction 0.01 \
        --min-qual 20 \
        --min-aae 0.5 \
        --min-aae-std 0.5 \
        --min-aae-std-frac 0.5 \
        --min-aae-frac 0.5 \
        --min-aae-frac-std 0.5 \
        --min-aae-frac-std-frac 0.5 \
        --min-aae-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac 0.5 \
        --min-aae-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std-frac-std 0.5 \
        --min-aae-frac-std-