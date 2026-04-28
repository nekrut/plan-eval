#!/usr/bin/env bash
set -euo pipefail

# Define paths relative to current working directory
RAW_DIR="data/raw"
REF_DIR="data/ref"
RES_DIR="results"

# Create results directory if it doesn't exist
mkdir -p "$RES_DIR"

# Define sample names and their paired-end files
declare -A SAMPLES
SAMPLES["M117-bl"]="M117-bl"
SAMPLES["M117-ch"]="M117-ch"
SAMPLES["M117C1-bl"]="M117C1-bl"
SAMPLES["M117C1-ch"]="M117C1-ch"

# Function to process a single sample
process_sample() {
    local sample_name="$1"
    local sample_dir="$RES_DIR/${sample_name}"
    
    # Create sample-specific directory
    mkdir -p "$sample_dir"
    
    # Define input files
    local fq1="${RAW_DIR}/${sample_name}_1.fq.gz"
    local fq2="${RAW_DIR}/${sample_name}_2.fq.gz"
    local ref="${REF_DIR}/chrM.fa"
    
    # Define output files
    local bam="${sample_dir}/${sample_name}.bam"
    local bai="${sample_dir}/${sample_name}.bam.bai"
    local vcf="${sample_dir}/${sample_name}.vcf.gz"
    local vcf_tbi="${sample_dir}/${sample_name}.vcf.gz.tbi"
    
    # Check if BAM already exists (idempotency)
    if [[ -f "$bam" ]]; then
        return 0
    fi
    
    # Align reads with BWA-MEM
    bwa mem -t 4 -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:MiSeq" \
        "$ref" "$fq1" "$fq2" | \
        samtools view -bS - > "$bam"
    
    # Sort and index BAM
    samtools sort -o "$bam" "$bam"
    samtools index "$bam"
    
    # Call variants with LoFreq
    lofreq call \
        -b "$bam" \
        -o "${sample_dir}/${sample_name}.vcf" \
        --min-depth 1 \
        --min-fraction 0.01 \
        --min-qual 20 \
        --min-var-freq 0.01 \
        --p-value 0.05 \
        --min-var-freq 0.01
    
    # Compress VCF and create index
    tabix -p vcf "${sample_dir}/${sample_name}.vcf"
    mv "${sample_dir}/${sample_name}.vcf" "$vcf"
    mv "${sample_dir}/${sample_name}.vcf.tbi" "$vcf_tbi"
}

# Function to collapse variants across all samples
collapse_variants() {
    local collapsed_file="${RES_DIR}/collapsed.tsv"
    
    # Check if any VCF exists (idempotency)
    if [[ ! -f "${RES_DIR}/M117-bl/M117-bl.vcf.gz" ]]; then
        echo "Error: No VCF files found. Run processing first." >&2
        exit 1
    fi
    
    # Create header
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$collapsed_file"
    
    # Process each sample's VCF
    for sample in "${!SAMPLES[@]}"; do
        local vcf="${RES_DIR}/${sample}/${sample}.vcf.gz"
        if [[ -f "$vcf" ]]; then
            # Extract variants and format for collapse
            bcftools query -f '%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$vcf" >> "$collapsed_file"
        fi
    done
}

# Main execution
# Process each sample
for sample in "${!SAMPLES[@]}"; do
    process_sample "$sample"
done

# Collapse variants
collapse_variants

exit 0