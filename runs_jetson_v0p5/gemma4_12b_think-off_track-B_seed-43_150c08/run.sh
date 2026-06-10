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
    VCF="results/${SAMPLE}.vcf.gz"
    
    # Skip if already exists
    if [[ -f "$BAM" && -f "$VCF" ]]; then
        continue
    fi

    # Map and sort/index BAM
    # Using bwa mem for alignment, samtools to produce a sorted, indexed bam
    bwa mem -t "$THREADS" -K "${REF}" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -1S - > "${BAM}.tmp"
    
    samtools sort -@ "$THREADS" -o "$BAM" "${BAM}.tmp"
    rm "${BAM}.tmp"
    
    samtools index "$BAM"
    
    # Variant calling with lofreq (preferred for mitochondrial/short reads)
    # lofreq outputs a VCF directly. 
    # We pipe to bcftools to ensure standard format and then compress/index.
    lofreq call -t "$THREADS" --ref "$REF" "$BAM" > "${VCF}.tmp"
    
    # Ensure the VCF is valid for tabix (bcftools filter/sort if needed)
    # lofreq output usually needs a header check, but we can just pipe to bcftools.
    bcftools view "${VCF}.tmp" > "${VCF}.raw"
    bgzip -c "${VCF}.raw" > "$VCF"
    tabix -1 "$VCF"
    rm "${VCF}.raw" "${VCF}.tmp"
done

# Create collapsed table
# Columns: sample, chrom, pos, ref, alt, af
# We iterate through the VCFs to extract data. 
# Note: lofreq provides AF in the FORMAT field or as a separate field depending on version.
# Standard bcftools query can extract these.

echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # Extract fields: SampleName, Chrom, Pos, Ref, Alt, and Allele Frequency (AF)
    # bcftools query -1 is the sample name if it exists in the header, 
    # but since we have a loop, we can just prepend the variable.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF" | \
    awk -v s="$SAMPLE" 'BEGIN{FS=OFS="\t"} {if($6==".") $6="0"; print s, $1, $2, $3, $4, $5}' >> results/collapsed.tsv
done

# Clean up any trailing empty lines or fix formatting if necessary (optional)
# The awk logic above handles the sample name injection.
# Re-running the loop to ensure exact column count and format:
# We'll overwrite the file with a clean pass to ensure no extra newlines from previous steps.
> results/collapsed.tsv
echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    bcftools query -n -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF" | \
    awk -v s="$SAMPLE" 'BEGIN{FS=OFS="\t"} {print s, $1, $2, $3, $4, $5}' >> results/collapsed.tsv
done