#!/usr/bin/env bash
set -euo pipefail

# Create output directory
mkdir -p results

# Reference file
REF="data/ref/chrM.fa"

# Check if reference is indexed, if not index it
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

# List of samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    VCF_GZ="results/${SAMPLE}.vcf.gz"

    # Skip if BAM already exists (idempotency)
    if [ -f "${BAM}" ] && [ -f "${VCF_GZ}" ]; then
        continue
    fi

    # Step 1: Alignment with BWA
    # Use mem algorithm for paired-end reads
    bwa mem -t 4 "${REF}" "${R1}" "${R2}" | \
        samtools sort -@ 4 -o "${BAM}" -

    # Step 2: Index BAM
    samtools index "${BAM}"

    # Step 3: Variant calling with LoFreq
    # LoFreq requires a FASTA reference and a BAM file
    # It outputs a VCF file
    lofreq call -f "${REF}" -d 0.0001 -o "results/${SAMPLE}.vcf" "${BAM}"

    # Step 4: Compress and index VCF with bcftools
    bcftools view -O z -o "${VCF_GZ}" "results/${SAMPLE}.vcf"
    bcftools index "${VCF_GZ}"

    # Clean up intermediate VCF
    rm -f "results/${SAMPLE}.vcf"
done

# Step 5: Collapse VCFs into a single TSV
# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv

# Process each sample's VCF and append to collapsed.tsv
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    
    # Use bcftools query to extract fields
    # Format: sample, chrom, pos, ref, alt, AF
    # Note: bcftools query format uses % for fields
    # %CHROM, %POS, %REF, %ALT, %AF (for genotype-level AF) or %INFO/AF
    # Since these are amplicon samples, we expect single samples per VCF.
    # We use %SAMPLE to get the sample name, but since we are iterating, we can hardcode or use %SAMPLE.
    # However, bcftools query with -s might be needed if multiple samples were in one VCF.
    # Here, each VCF has one sample.
    
    bcftools query -f "%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n" "${VCF_GZ}" | \
        awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >> results/collapsed.tsv
done