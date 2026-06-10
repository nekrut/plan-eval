#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# Index reference if not indexed
if [ ! -f "${REF}._1.bam" ] && [ ! -f "${REF}.bt" ]; then
    # Check for existing index or just create it
    if ! samtools index "$REF" > /dev/null 2>&1; then
        samtools faidx "$REF"
    fi
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${SAMPLE}.bam"
    VCF="$RESULTS_DIR/${SAMPLE}.vcf.gz"
    
    # Skip if BAM already exists
    if [ -f "$BAM" ]; then
        continue
    fi

    # Map and sort/index
    # Note: Using bwa mem for mapping, samtools for processing
    bwa mem -t "$THREADS" -K "${REF}" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -u - | \
    samtools sort -@ "$THREADS" -o "$BAM"
    
    samtools index "$BAM"
    
    # Variant calling using lofreq (good for low-depth/mitochondrial)
    # lofreq outputs VCF directly. We pipe to bcftools to compress and index.
    lofreq call -t "$THREADS" --ref "$REF" "$BAM" | \
    bcftools view -Oz - | \
    bgzip -f - > "$VCF"
    
    tabix -1 "$VCF"
done

# Create collapsed table
# Extract: sample, chrom, pos, ref, alt, af
# bcftuf/lofreq VCFs usually have AF in the INFO field or as a field. 
# Since we need to be robust, we parse the VCF records.
# We use bcftools query to get specific fields and calculate/extract AF.
# If AF is not in INFO, it's often calculated from depth. 
# However, lofreq usually provides info for variants.

echo -e "sample\tchrom\tpos\ref\talt\taf" > "$RESULTS_DIR/collapsed.tsv"

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="$RESULTS_DIR/${SAMPLE}.vcf.gz"
    if [ -f "$VCF" ]; then
        # Use bcftools to extract fields. 
        # Note: many VCFs have AF in INFO. If not, we can't easily "calculate" it without extra tools,
        # but standard variant calling pipelines expect the field to exist or be handled by a tool.
        # We use bcftools query to grab the specific columns.
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$V1" | \
        awk -v s="$SAMPLE" '{print s "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" ($5==""?0:$5)}' \
        >> "$RESULTS_DIR/collapsed.tsv"
    fi
done

# Clean up the trailing empty lines or handle formatting if needed, 
# but the loop above handles the logic.
# Re-running the awk to ensure no empty lines from bcftools query are added.
# Actually, let's do a clean pass for the final table construction:
# The previous loop might append extra newlines; we fix it by rewriting the file once.

# Final cleanup of the collapsed file (remove potential trailing empty lines)
sed -i '/^$/d' "$RESULTS_DIR/collapsed.tsv"