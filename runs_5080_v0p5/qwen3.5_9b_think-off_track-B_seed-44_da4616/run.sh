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
    local r1="${RAW_DIR}/${sample_name}_1.fq.gz"
    local r2="${RAW_DIR}/${sample_name}_2.fq.gz"
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
    bwa mem -t 4 -R "@RG\tID:${sample_name}\tSM:${sample_name}\tPL:MiSeq\tPN:MiSeq\tCL:MiSeq" \
        "$ref" "$r1" "$r2" | \
    samtools view -b - > "$bam"
    
    # Sort and index BAM
    samtools sort -o "$bam" "$bam"
    samtools index "$bam"
    
    # Call variants with LoFreq
    lofreq call \
        -o "${sample_dir}/${sample_name}.vcf" \
        -b "$bam" \
        -f "$ref" \
        --min-coverage 1 \
        --min-fraction 0.01 \
        --min-qual 20 \
        --min-depth 1
    
    # Convert to VCF and index
    bcftools view -Oz -o "$vcf" "${sample_dir}/${sample_name}.vcf"
    tabix "$vcf"
    
    # Cleanup intermediate VCF
    rm -f "${sample_dir}/${sample_name}.vcf"
}

# Function to generate collapsed table
generate_collapsed_table() {
    local output_file="${RES_DIR}/collapsed.tsv"
    
    # Check if output already exists (idempotency)
    if [[ -f "$output_file" ]]; then
        return 0
    fi
    
    # Write header
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$output_file"
    
    # Collect all VCF files
    local vcf_files=()
    for sample in "${!SAMPLES[@]}"; do
        local vcf="${RES_DIR}/${sample}/${sample}.vcf.gz"
        if [[ -f "$vcf" ]]; then
            vcf_files+=("$vcf")
        fi
    done
    
    # If no VCFs found, exit
    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        return 0
    fi
    
    # Concatenate all VCFs and filter for mitochondrial variants (chrM)
    # Use bcftools to merge and then awk to extract AF
    bcftools concat -Oz -o "${RES_DIR}/merged.vcf.gz" "${vcf_files[@]}"
    
    # Filter for chrM and extract required fields
    bcftools view -f PASS -i 'CHROM=="chrM"' "${RES_DIR}/merged.vcf.gz" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {
        sample = "'${SAMPLES[$1]}'"
        # AF is in INFO field, extract it
        match($0, /AF=([0-9.]+)/, arr)
        af = arr[1]
        if (af == "") af = "NA"
        print sample, $1, $2, $4, $5, af
    }' >> "$output_file"
    
    # Cleanup merged VCF
    rm -f "${RES_DIR}/merged.vcf.gz"
}

# Process each sample
for sample in "${!SAMPLES[@]}"; do
    process_sample "$sample"
done

# Generate collapsed table
generate_collapsed_table

exit 0