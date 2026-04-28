#!/usr/bin/env bash
set -euo pipefail

# Create output directory
mkdir -p results

# Reference files
REF=data/ref/chrM.fa
REF_BWT="${REF}.bwt"
REF_PAC="${REF}.pac"
REF_ANI="${REF}.ann"
REF_AMB="${REF}.amb"
REF_SA="${REF}.sa"
REF_FAI="${REF}.fai"

# Index reference if not present
if [ ! -f "$REF_BWT" ]; then
    bwa index "$REF"
fi

# List of samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="${BAM}.bai"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="${VCF_GZ}.tbi"

    # Step 1: Alignment with BWA-MEM
    if [ ! -f "$BAM" ]; then
        bwa mem -t 4 "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "$BAM" -
    fi

    # Step 2: Index BAM
    if [ ! -f "$BAI" ]; then
        samtools index -@ 4 "$BAM"
    fi

    # Step 3: Variant calling with LoFreq
    # LoFreq requires a FASTA reference and a BAM file
    # It performs base quality recalibration and variant calling
    if [ ! -f "$VCF_GZ" ]; then
        lofreq call-parallel \
            -f "$REF" \
            -r "$BAM" \
            -o "$VCF_GZ" \
            --nproc 4 \
            --force
    fi

    # Step 4: Index VCF
    if [ ! -f "$VCF_TBI" ]; then
        tabix -p vcf "$VCF_GZ"
    fi
done

# Step 5: Collapse all VCFs into a single TSV
# Columns: sample chrom pos ref alt af
# We need to extract these fields from each VCF file

# Create header
printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv

# Process each VCF
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    
    # Use bcftools query to extract fields
    # FORMAT: %CHROM\t%POS\t%REF\t%ALT\t%AF
    # Note: AF is in the INFO field, usually as AF=... or we can use %AF if available in FORMAT
    # For LoFreq, the AF is typically in the INFO field. Let's use bcftools to extract it.
    # LoFreq VCF format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
    # INFO contains AF=...
    
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" | \
    while IFS=$'\t' read -r chrom pos ref alt af; do
        # Handle multi-allelic sites where AF might be a comma-separated list
        # We'll take the first AF value for simplicity, or split if needed
        # LoFreq typically reports single alleles per line or uses comma separation
        # If AF is empty or missing, skip or set to 0
        if [ -z "$af" ]; then
            af="0"
        fi
        # If AF is a list (e.g., "0.5,0.5"), take the first one
        af=$(echo "$af" | cut -d',' -f1)
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "$alt" "$af"
    done >> results/collapsed.tsv
done