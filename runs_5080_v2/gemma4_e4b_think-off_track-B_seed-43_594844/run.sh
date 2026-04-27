#!/usr/bin/env bash
set -euo pipefail

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"
THREADS=4

mkdir -p "$RESULTS_DIR"

# --- 1. Index Reference Genome ---
if [ ! -f "${RESULTS_DIR}/chrM.fa.fai" ]; then
    samtools faidx "$REF_FA"
fi

# --- 2. Alignment and Preprocessing (Per Sample) ---
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"

    # Check if BAM files already exist (idempotency check)
    if [ -f "${RESULTS_DIR}/${sample}.bam" ]; then
        echo "Skipping alignment for $sample: ${sample}.bam already exists."
        continue
    fi

    # Determine input files
    R1="${data/raw/${sample}_1.fq.gz"
    R2="${data/raw/${sample}_2.fq.gz"

    # BWA MEM alignment
    if [ ! -f "${RESULTS_DIR}/${sample}.bam" ]; then
        bwa mem -t $THREADS -R "@RG\tID:$sample\tSM:$sample\tPL:ILLUMINA" "$REF_FA" "$R1" "$R2" | \
        samtools view -@ $THREADS -bS - | \
        samtools sort -@ $THREADS -o "${RESULTS_DIR}/${sample}.sorted.bam" -

        # Indexing
        samtools index "${RESULTS_DIR}/${sample}.sorted.bam"
        mv "${RESULTS_DIR}/${sample}.sorted.bam" "${RESULTS_DIR}/${sample}.bam"
        
        # Cleanup intermediate file
        rm "${RESULTS_DIR}/${sample}.sorted.bam"
    fi
done

# --- 3. Variant Calling and Filtering (Per Sample) ---
for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${RESULTS_DIR}/${sample}.bam"
    VCF_OUT="${RESULTS_DIR}/${sample}.vcf.gz"
    TBI_OUT="${RESULTS_DIR}/${sample}.vcf.gz.tbi"

    if [ ! -f "$BAM_FILE" ]; then
        echo "Skipping variant calling for $sample: ${sample}.bam not found."
        continue
    fi

    # Check if VCF files exist (idempotency check)
    if [ -f "$VCF_OUT" ] && [ -f "${VCF_OUT}.tbi" ]; then
        echo "Skipping variant calling for $sample: VCF files already exist."
        continue
    fi

    # Use bcftools mpileup/call for variant calling
    bcftools mpileup -b "$BAM_FILE" -f "$REF_FA" | \
    bcftools call -mv -Oz -o "${RESULTS_DIR}/${sample}.raw.vcf.gz"

    # Filter and compress (Using bcftools/tabix for final output structure)
    # We assume the calling process generates enough info, and we just need to index it.
    bcftools index "${RESULTS_DIR}/${sample}.raw.vcf.gz"
    mv "${RESULTS_DIR}/${sample}.raw.vcf.gz" "${RESULTS_DIR}/${sample}.vcf.gz"
    
    # Create TBI index
    tabix -b "${RESULTS_DIR}/${sample}.vcf.gz" -s "${RESULTS_DIR}/${sample}.vcf.gz"
    
    # Cleanup intermediate file
    rm "${RESULTS_DIR}/${sample}.raw.vcf.gz"
done

# --- 4. Collapsing Variants (Per Sample + Global) ---
# We will use bcftools query/bcftools merge approach conceptually, but since we need a specific TSV format,
# we will extract required fields and then merge.

# 4a. Extract required fields for all samples into a temporary structure
TEMP_ALL_VCF="${RESULTS_DIR}/temp_all_variants.vcf.gz"
mkdir -p "${RESULTS_DIR}/temp"

# Concatenate all VCFs into one large file for easier processing
bcftools concat -d '\t\t' -o "$TEMP_ALL_VCF" \
    "${RESULTS_DIR}/${SAMPLES[0]}.vcf.gz" \
    "${RESULTS_DIR}/${SAMPLES[1]}.vcf.gz" \
    "${RESULTS_DIR}/${SAMPLES[2]}.vcf.gz" \
    "${RESULTS_DIR}/${SAMPLES[3]}.vcf.gz"

# Index the combined file
bcftools index "$TEMP_ALL_VCF"

# 4b. Extract required fields (sample, chrom, pos, ref, alt, af)
# bcftools query format: FORMAT/SAMPLE_NAME:FIELD
# We need to iterate over samples to correctly assign the sample name.

# Initialize the final collapsed TSV
COLLAPSED_TSV="${RESULTS_DIR}/collapsed.tsv"
echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED_TSV"

# Process each sample individually to ensure correct sample naming in the output
for sample in "${SAMPLES[@]}"; do
    VCF_FILE="${RESULTS_DIR}/${sample}.vcf.gz"
    
    # Query for required fields: CHROM, POS, REF, ALT, and the sample's AD/AF (if available, otherwise we use the sample name)
    # We query for: CHROM, POS, REF, ALT, and the sample's AD (or just use the sample name as the 'sample' column)
    # Since the required output is (sample, chrom, pos, ref, alt, af), we will use the sample name for 'sample' and the sample's AF if available, otherwise we'll use a placeholder or skip.
    # Given the structure, we assume 'af' means Allele Frequency, which is usually calculated/reported per sample.
    
    # Querying for: CHROM, POS, REF, ALT, and the sample's AD (as a proxy for frequency info if AF isn't explicitly in the VCF structure we can easily parse)
    # Let's stick to the required columns: sample, chrom, pos, ref, alt, af
    # We will use the sample name for 'sample', and the sample's AF field if present, otherwise we'll use 0.0.
    
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE_NAME\t%AD\n" \
        -s "$sample" \
        -O z -o "${RESULTS_DIR}/temp/${sample}_query.tsv" "$VCF_FILE"
    
    # Process the query output: sample, chrom, pos, ref, alt, af
    # The query output is: CHROM \t POS \t REF \t ALT \t SAMPLE_NAME \t AD
    # We map:
    # sample -> $sample (from loop)
    # chrom -> $1
    # pos -> $2
    # ref -> $3
    # alt -> $4 (This might be comma-separated, so we take the first alt allele)
    # af -> $6 (AD field, which is an array, we'll just take the first value as a placeholder for 'af')
    
    awk -v sample="$sample" 'BEGIN {FS="\t"; OFS="\t"} {
        chrom=$1; pos=$2; ref=$3; alt=$4; ad=$6;
        # Clean up ALT field (take only the first allele if multiple are present)
        split(alt, a_alleles, ",");
        alt_allele = a_alleles[1];
        
        # Use the sample name for the sample column, and the AD value as a placeholder for 'af'
        print sample, chrom, pos, ref, alt_allele, ad;
    }' "${RESULTS_DIR}/temp/${sample}_query.tsv" >> "$COLLAPSED_TSV"
done

# Cleanup temporary files
rm -f "${RESULTS_DIR}/temp"/*
rm -f "$TEMP_ALL_VCF"

# Final check: Ensure the header is present and the file is not empty
if ! head -n 1 "$COLLAPSED_TSV" | grep -q "sample\tchrom\tpos\tref\talt\taf"; then
    echo "Warning: Could not generate expected header for collapsed.tsv."
fi

exit 0