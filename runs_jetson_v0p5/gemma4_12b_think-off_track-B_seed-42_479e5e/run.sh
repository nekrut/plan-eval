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

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    
    # Skip if already exists (Idempotency)
    if [[ -f "$BAM" && -f "$VCF_GZ" ]]; then
        continue
    fi

    # Map and convert to BAM
    bwa mem -t "$THREADS" -K "${REF}" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -Sb - > "$BAM"
    
    # Index BAM
    samtools index "$BAM"

    # Variant calling with lofreq (handles low coverage and mitochondrial noise well)
    # lofreq outputs VCF directly
    lofreq call "$BAM" > "$VCF_GZ"
    
    # Index VCF
    tabix -1 "chrM" -1r "$VCF_GZ" > "${VCF_GZ}.tbi"
done

# Create collapsed table
# Columns: sample, chrom, pos, ref, alt, af
# We use bcftools query to extract fields and awk to format. 
# Note: lofreq VCFs usually contain the allele frequency in the FORMAT field or as a separate info tag depending on version;
# however, standard VCF parsing for 'af' is safest via bcftools.
echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # Extract fields: CHROM, POS, REF, ALT and the AF info field if present
    # If AF is not a standard tag in the specific VCF version, we default to 1.0 for variants found.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF" | \
    awk -v s="$SAMPLE" 'BEGIN{OFS="\t"} {if($5==".") $5="1.0"; print s,$1,$2,$3,$4,$5}' >> results/collapsed.tsv
done

# Clean up the trailing empty lines or potential formatting issues from multiple appends if necessary, 
# but the logic above ensures one line per variant per sample.
# The header is already written once at the start of the loop's first iteration (or we can just write it once before the loop).
# To ensure exactly one header:
# (The logic above actually appends the data; let's refine to ensure only one header exists)

# Refined Collapsed Logic:
# 1. Create a temp file for all data
# 2. Prepend header
# 3. Move to final
TEMP_COLLAPSE=$(mktemp)
echo -e "sample\tchrom\tpos\tref\talt\taf" > "$TEMP_COLLAPSE"
for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF" | \
    awk -v s="$SAMPLE" 'BEGIN{OFS="\t"} {if($5==".") $5="1.0"; print s,$1,$2,$3,$4,$5}' >> "$TEMP_COLLAPSE"
done
mv "$TEMP_COLLAPSE" results/collapsed.tsv