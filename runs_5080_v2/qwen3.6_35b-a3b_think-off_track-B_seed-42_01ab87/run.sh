#!/usr/bin/env bash
set -euo pipefail

# Create output directory
mkdir -p results

# Reference file paths
REF=data/ref/chrM.fa
REF_BWT="${REF}.bwt"
REF_PAC="${REF}.pac"
REF_ANI="${REF}.ann"
REF_AMB="${REF}.amb"
REF_SA="${REF}.sa"

# Check if reference is indexed, if not index it
if [[ ! -f "$REF_BWT" ]]; then
    bwa index "$REF"
fi

# Define samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Function to process a single sample
process_sample() {
    local sample=$1
    local fq1="data/raw/${sample}_1.fq.gz"
    local fq2="data/raw/${sample}_2.fq.gz"
    local bam="results/${sample}.bam"
    local bai="${bam}.bai"
    local vcf_gz="results/${sample}.vcf.gz"
    local vcf_tbi="${vcf_gz}.tbi"

    # Skip if all outputs exist (idempotency)
    if [[ -f "$bam" && -f "$bai" && -f "$vcf_gz" && -f "$vcf_tbi" ]]; then
        return 0
    fi

    # 1. Mapping with BWA-MEM
    # Use -M for mark shorter split hits as secondary (compatible with GATK/samtools)
    # Use -R for read group info
    bwa mem -t 4 -M -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" "$REF" "$fq1" "$fq2" | \
        samtools sort -@ 4 -o "$bam" -

    # 2. Index BAM
    samtools index "$bam"

    # 3. Variant Calling with LoFreq
    # LoFreq requires a FASTA reference and a BAM file
    # --nproc 4 for parallel processing
    lofreq call -f "$REF" -r 1-16569 -i "$bam" -o "results/${sample}.vcf" --nproc 4

    # 4. Compress and Index VCF
    bgzip -c "results/${sample}.vcf" > "$vcf_gz"
    rm "results/${sample}.vcf" # Clean up uncompressed VCF
    tabix -p vcf "$vcf_gz"
}

# Process all samples in parallel
pids=()
for sample in "${SAMPLES[@]}"; do
    process_sample "$sample" &
    pids+=($!)
done

# Wait for all samples to finish
for pid in "${pids[@]}"; do
    wait "$pid"
done

# 5. Collapse VCFs into a single TSV
# We need to extract: sample, chrom, pos, ref, alt, af
# LoFreq VCF format typically has:
# CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1 ...
# AF is usually in the INFO field or derived from AD/DP in FORMAT.
# LoFreq outputs AF in the INFO field as "AF=..." or in the FORMAT as "AF".
# Let's use bcftools to normalize and extract fields.

# First, create a temporary directory for intermediate files
mkdir -p results/tmp

# Process each VCF to extract the required columns
# Using bcftools query to get specific fields
# LoFreq VCFs have the sample name in the last column header.
# We want: sample, chrom, pos, ref, alt, af
# Note: LoFreq might output multiple ALT alleles or complex variants. We'll handle simple SNPs/Indels.
# The AF field in LoFreq VCF is often in the INFO field.

for sample in "${SAMPLES[@]}"; do
    vcf="results/${sample}.vcf.gz"
    # Extract CHROM, POS, REF, ALT, and AF from INFO
    # LoFreq INFO field contains AF. Format: ...;AF=0.5;...
    # We use bcftools query to parse this.
    # If AF is not in INFO, we might need to calculate from AD/DP in FORMAT.
    # However, LoFreq explicitly outputs AF in INFO.
    
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf" | \
    while IFS=$'\t' read -r chrom pos ref alt af; do
        # Handle missing AF (shouldn't happen for called variants but good practice)
        if [[ -z "$af" ]]; then
            af="."
        fi
        # Handle multi-allelic sites: bcftools query might output comma-separated values for ALT and AF
        # We need to explode these into separate rows
        IFS=',' read -ra alts <<< "$alt"
        IFS=',' read -ra afs <<< "$af"
        
        for i in "${!alts[@]}"; do
            local_alt="${alts[$i]}"
            local_af="${afs[$i]}"
            if [[ -z "$local_af" ]]; then
                local_af="."
            fi
            # Skip if alt is empty or AF is empty after split
            if [[ -n "$local_alt" && -n "$local_af" ]]; then
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$chrom" "$pos" "$ref" "$local_alt" "$local_af"
            fi
        done
    done >> results/tmp/all_variants.tsv
done

# Create the final collapsed.tsv with header
printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv
cat results/tmp/all_variants.tsv >> results/collapsed.tsv

# Clean up temp directory
rm -rf results/tmp