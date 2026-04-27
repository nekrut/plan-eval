#!/usr/bin/env bash
set -euo pipefail

# Create output directory
mkdir -p results

# Reference files
REF=data/ref/chrM.fa
REF_BWA=data/ref/chrM.fa.bwt
REF_SAMTOOLS=data/ref/chrM.fa.fai

# Index reference if not present
if [ ! -f "$REF_BWA" ]; then
    bwa index "$REF"
fi

if [ ! -f "$REF_SAMTOOLS" ]; then
    samtools faidx "$REF"
fi

# Define samples
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="results/${SAMPLE}.bam"
    BAI="${BAM}.bai"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="${VCF_GZ}.tbi"

    # Skip if BAM already exists (idempotency)
    if [ -f "$BAI" ]; then
        continue
    fi

    # Mapping
    bwa mem -t 4 "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "$BAM" -
    samtools index "$BAM"

    # Variant calling with lofreq
    # lofreq calls variants relative to the reference, so we need the reference sequence
    lofreq call -f "$REF" -i "$BAM" -o "results/${SAMPLE}.vcf" -n 4 --force

    # Compress and index VCF
    bgzip -c "results/${SAMPLE}.vcf" > "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"

    # Clean up uncompressed VCF
    rm -f "results/${SAMPLE}.vcf"
done

# Create collapsed.tsv
# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv

# Process each sample's VCF to extract relevant fields
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    
    # Use bcftools to query the VCF and format output
    # We need: sample, chrom, pos, ref, alt, af (from INFO/AF or FORMAT/AF)
    # lofreq typically puts AF in the INFO field as AF=... or in the FORMAT field.
    # Let's use bcftools query to extract standard fields.
    # lofreq output format:
    # CHROM POS ID REF ALT QUAL FILTER INFO FORMAT ...
    # INFO often contains AF=...
    
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF_GZ" 2>/dev/null | \
    while IFS=$'\t' read -r chrom pos ref alt af; do
        # Handle multi-allelic sites or missing AF
        if [ -z "$af" ]; then
            continue
        fi
        
        # Split multi-allelic sites if necessary (bcftools query handles some, but let's be safe)
        # lofreq usually outputs one line per variant. If ALT is comma-separated, AF might be comma-separated too.
        # We will assume simple cases for now, but handle comma-separated AF if present.
        
        # If AF is comma-separated, we need to split it and create multiple rows
        if [[ "$af" == *","* ]]; then
            IFS=',' read -ra AF_ARRAY <<< "$af"
            IFS=',' read -ra ALT_ARRAY <<< "$alt"
            
            for i in "${!ALT_ARRAY[@]}"; do
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "${ALT_ARRAY[$i]}" "${AF_ARRAY[$i]}"
            done
        else
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "$alt" "$af"
        fi
    done >> results/collapsed.tsv
done