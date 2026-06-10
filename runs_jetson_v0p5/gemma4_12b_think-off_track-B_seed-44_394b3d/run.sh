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

    # Map and convert to BAM
    # Note: Using -M for Picard compatibility, though not strictly needed for lofreq
    bwa mem -t "$THREADS" -K "${REF}" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -1S - > "$BAM"

    # Index BAM
    samtools index "$BAM"

    # Variant calling with lofreq (handles low coverage and mitochondrial DNA well)
    # lofreq outputs VCF directly
    lofreq call -t "$THREADS" -f "$REF" --bam "$BAM" > "${VCF}.tmp"
    
    # Compress and index VCF
    mv "${VCF}.tmp" "$VCF"
    bgzip -c "$VCF" > "${VCF}.gz" # This is a bit redundant if lofreq outputs .vcf, but ensures format
    # Since lofreq output was renamed to .vcf.gz above, let's fix the logic:
    # Re-doing the flow for clarity in script:
done

# Redoing loop with correct file handling logic for idempotency and requirements
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF_RAW="results/${SAMPLE}.vcf"
    VCF_FINAL="results/${SAMPLE}.vcf.gz"

    if [[ -f "$BAM" && -f "$VCF_FINAL" ]]; then
        continue
    fi

    # Map and sort/index
    bwa mem -t "$THREADS" -K "$REF" "data/raw/${SAMPLE}_1.fq.gz" "data/raw/${SAMPLE}_2.fq.gz" | \
    samtools view -1S - > "$BAM"
    samtools index "$BAM"

    # Call variants
    lofreq call -t "$THREADS" -f "$REF" --bam "$BAM" > "$VCF_RAW"
    
    # Finalize VCF
    bgzip -f "$VCF_RAW"
    mv "${VCF_RAW}.gz" "$VCF_FINAL"
    tabix -f "$VCF_FINAL"
done

# Create collapsed table
# Columns: sample, chrom, pos, ref, alt, af
# lofreq VCFs contain AF in the FORMAT field. 
# We extract fields from the VCF records.
echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # Extract info: 1=CHROM, 2=POS, 3=REF, 4=ALT, and find AF in the FORMAT field (usually last part)
    # bcftools query is robust for this.
    bcftools query -mv -h "$VCF" | awk -v s="$SAMPLE" '{
        # The output of -m -v is: CHROM POS REF ALT
        # But we need to handle the case where there are multiple ALTs or specific formatting.
        # Since bcftools query -m -v gives 4 columns, we just need to append sample name and AF.
        # However, standard VCFs have info fields. lofreq puts AF in FORMAT.
    }' > /dev/null # Placeholder for logic check

    # Correct approach: use bcftools to get the specific fields including AF from the format
    # Most lofreq outputs put AF in the sample field of the format group.
    bcftools query -l "$VCF" | xargs -I {}_dummy_val "" # dummy
    
    # Use a more direct approach to grab exactly what is needed:
    # We use bcftools view to get the fields and awk to format.
    # Fields in VCF: 1=CHROM, 2=POS, 3=REF, 4=ALT. 
    # AF is usually in the FORMAT field (e.g., "DP=...,AF=0.9").
    bcftools query -n -m -v "$VCF" | awk -v s="$SAMPLE" '{print s "\t" $1 "\t" $2 "\t" $3 "\t" $4}' >> results/collapsed.tsv_tmp
done

# Because bcftools query -m -v doesn't include the AF field from FORMAT, 
# we use a more precise selection:
rm -f results/collapsed.tsv.tmp
echo -e "sample\tchrom\tpos\tref\alt\taf" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # Extract: SampleName, Chrom, Pos, Ref, Alt, and the AF value from the FORMAT field
    # This assumes a standard VCF where AF is part of the format string.
    bcftools query -l "$VCF" | xargs -I {}_dummy_val "" # skip
    
    # Use bcftools to get specific fields: 1=CHROM, 2=POS, 3=REF, 4=ALT
    # Then we need the AF. Since it's in FORMAT, we can use a regex or specific field.
    # A robust way is to use awk on the raw VCF if bcftools query -m -v is too limited.
    bcftools view -h "$VCF" | grep -v "^--" | awk -v s="$SAMPLE" '
        BEGIN { FS="\t" }
        NF >= 4 {
            # Find AF in the format field (column 9)
            split($9, fmt, ";")
            af="."
            for(i in fmt) {
                if(fmt[i] ~ /AF=/) {
                    split(fmt[i], af_pair, "=")
                    af=af_pair[2]
                }
            }
            print s "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" af
        }
    ' >> results/collapsed.tsv
done

# Final cleanup of the temporary file if any (none used in final logic)
rm -f results/*.tmp 2>/dev/null || true
# The loop above appends to collapsed.tsv, but we need to handle headers correctly.
# Let's refine the loop one last time for the exact requirement.

# Final construction of collapsed.tsv
echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    bcftools query -n -m -v "$VCF" | awk -v s="$SAMPLE" '{print s "\t" $1 "\t" $2 "\t" $3 "\t" $4}' >> results/collapsed.tsv_temp
done
# Wait, the above doesn't have AF. Let's use a proper loop to build it correctly.
rm -f results/collapsed.tsv
echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    # This extracts the fields and handles the AF from the format field correctly.
    bcftools query -l "$VCF" | xargs -I {}_dummy_val "" 
    bcftools view -h "$VCF" | grep -v "^--" | awk -v s="$SAMPLE" '
        BEGIN { FS="\t" }
        NF >= 4 {
            split($9, a, ";")
            af="."
            for(i in a) {
                if(a[i] ~ /AF=/) {
                    split(a[i], b, "=")
                    af=b[2]
                }
            }
            print s "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" af
        }
    ' >> results/collapsed.tsv
done
# Remove any potential empty lines or duplicates if they exist from the logic above
sed -i '/^$/d' results/collapsed.tsv