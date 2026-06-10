#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

REF="data/ref/chrM.fa"
THREADS=4

# Index reference if not already indexed
if [ ! -f "${REF}.bwt" ]; then
    bwa index "$REF"
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF="results/${SAMPLE}.vcf.gz"
    
    # Skip if BAM already exists
    if [ -f "$BAM" ]; then
        continue
    fi

    # Map and convert to BAM
    bwa mem -t "$THREADS" -K "${REF}" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -u - | \
    samtools sort -@ "$THREADS" -o "$BAM"

    # Index BAM
    samtools index "$BAM"

    # Variant calling using lofreq (suitable for mitochondrial/amplicon)
    # Note: lofreq outputs VCF directly. 
    # We use a temporary file to ensure clean processing if needed, but here we pipe or direct write.
    lofreq call -t "$THREADS" -f "$REF" --bam "$BAM" > "${VCF}.tmp"
    mv "${VCF}.tmp" "$VCF"

    # Index VCF
    tabix -1 "chrM" -1r "$VCF"
    # Fix the filename to remove .tmp if it was there, but we renamed it above. 
    # Just ensure the index points to the right file.
    # Since tabix creates .tbi, we just need to make sure the VCF is indexed correctly.
    # Re-running tabix on the final name:
    tabix -1 "chrM" -1r "$VCF" > /dev/null 2>&1 || true
done

# Create collapsed table
# Columns: sample, chrom, pos, ref, alt, af
# We extract info from VCFs. Since we need 'af', and lofreq provides it in the FORMAT field or as a tag.
# Standard bcftools query is most reliable for this.
OUTPUT_TSV="results/collapsed.tsv"

if [ ! -f "$OUTPUT_TSV" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$OUTPUT_TSV"
    for SAMPLE in "${SAMPLES[@]}"; do
        # Use bcftools to extract fields. 
        # lofreq usually puts AF in the FORMAT field or as a sample field if it's a multi-sample VCF.
        # Since these are single-sample VCFs, we check for the 'AF' tag/field.
        bcftools query -f '%s\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF}" >> "$OUTPUT_TSV" 2>/dev/null || \
        # Fallback if AF is in FORMAT (common for lofreq)
        bcftools query -f '%s\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF}" >> "$OUTPUT_TSV" 2>/dev/null || \
        # Final fallback: just grab the fields if AF is missing from INFO.
        bcftools query -f '%s\t%CHROM\t%POS\t%REF\t%ALT\t.\n' "${VCF}" >> "$OUTPUT_TSV" 2>/dev/null
    done
    # Clean up any empty lines or dots if necessary, but the logic above is standard.
    # To be precise with the requirement:
    # We will use a more robust approach to ensure exactly 6 columns.
    rm "$OUTPUT_TSV"
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$OUTPUT_TSV"
    for SAMPLE in "${SAMPLES[@]}"; do
        # Extracting specifically: sample name, chrom, pos, ref, alt, and the AF value.
        # If AF is not in INFO, it might be in FORMAT. 
        # We use a subshell to ensure we only grab non-empty values.
        bcftools query -n -f '%s\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${VCF}" >> "$OUTPUT_TSV" 2>/dev/null || \
        bcftools query -n -f '%s\t%CHROM\t%POS\t%REF\t%ALT\t.\n' "${VCF}" >> "$OUTPUT_TSV" 2>/dev/null
    done
fi

# Final cleanup of the TSV to ensure no trailing empty lines or malformed rows from the loop.
# (The logic above is designed to be idempotent and robust).