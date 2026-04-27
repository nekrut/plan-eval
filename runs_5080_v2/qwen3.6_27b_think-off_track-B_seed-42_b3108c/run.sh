#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4

# Sample list
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory
mkdir -p "$OUT_DIR"

# 1. Index Reference (Idempotent: check if .fai exists)
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# 2. Per-sample processing
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Check if final VCF index exists to skip sample if already done
    if [[ -f "$TBI" ]]; then
        continue
    fi

    # Alignment with BWA
    # Use bwa mem, pipe to samtools view to convert to BAM, sort, and output
    # We use a temporary unsorted BAM in memory/pipe to avoid disk I/O for intermediate files
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | \
        samtools view -bS -@ "$THREADS" - | \
        samtools sort -@ "$THREADS" -o "$BAM" -

    # Index BAM
    samtools index "$BAM"

    # Variant Calling with Lofreq
    # lofreq call -> filter -> mpileup -> call -> vcf
    # lofreq call --call-indels is good for amplicons
    # We use --min-af 0.01 to catch low frequency variants typical in mtDNA
    # Output is VCF
    
    # lofreq call produces VCF. We pipe to bcftools to ensure proper formatting and bgzip.
    # Note: lofreq call output is VCF. We need to bgzip it.
    
    # Step 1: Call variants
    # lofreq call -f <ref> -o <vcf> <bam>
    # We use --min-af 0.001 to be sensitive, but filter later if needed. 
    # For mtDNA, we often want high sensitivity.
    
    # Create a temporary VCF file for lofreq output
    TMP_VCF="${OUT_DIR}/${SAMPLE}.vcf"
    
    lofreq call \
        --call-indels \
        --min-af 0.001 \
        --min-qual 20 \
        -f "$REF" \
        -o "$TMP_VCF" \
        "$BAM"

    # Step 2: Normalize and filter VCF
    # Remove duplicates, sort, and bgzip
    # bcftools norm -f <ref> is good practice
    # bcftools view -e 'FILTER="low"' to remove low quality calls if any
    
    # Sort and normalize
    bcftools sort -Oz -o "$VCF_GZ" "$TMP_VCF"
    
    # Index VCF
    tabix -p vcf "$VCF_GZ"
    
    # Clean up temporary VCF
    rm -f "$TMP_VCF"
done

# 3. Collapse VCFs into a single TSV
# Columns: sample  chrom  pos  ref  alt  af

COLLAPSED_TSV="${OUT_DIR}/collapsed.tsv"

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED_TSV"

for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    
    # Extract variants from VCF
    # bcftools query to get CHROM, POS, REF, ALT, AF
    # AF is in INFO field. For lofreq, it's usually AF or VAF.
    # lofreq uses AF in INFO.
    # Format: %CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
    
    # Handle case where AF might be missing or multi-allelic
    # We assume biallelic for simplicity in TSV, or split if multi-allelic
    
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" | \
    while IFS=$'\t' read -r CHROM POS REF ALT AF; do
        # If AF is missing, use .
        if [[ -z "$AF" || "$AF" == "." ]]; then
            AF="."
        fi
        
        # Handle multi-allelic sites by splitting ALT and AF
        # ALT is comma-separated, AF is comma-separated
        IFS=',' read -ra ALT_ARR <<< "$ALT"
        IFS=',' read -ra AF_ARR <<< "$AF"
        
        for i in "${!ALT_ARR[@]}"; do
            ALT_SINGLE="${ALT_ARR[$i]}"
            AF_SINGLE="${AF_ARR[$i]:-}"
            if [[ -z "$AF_SINGLE" ]]; then
                AF_SINGLE="."
            fi
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$CHROM" "$POS" "$REF" "$ALT_SINGLE" "$AF_SINGLE"
        done
    done >> "$COLLAPSED_TSV"
done

# Sort the collapsed TSV by sample, chrom, pos for consistency (optional but good practice)
# We keep the header at the top
head -n 1 "$COLLAPSED_TSV" > "${COLLAPSED_TSV}.tmp"
tail -n +2 "$COLLAPSED_TSV" | sort -t$'\t' -k1,1 -k2,2 -k3,3n >> "${COLLAPSED_TSV}.tmp"
mv "${COLLAPSED_TSV}.tmp" "$COLLAPSED_TSV"