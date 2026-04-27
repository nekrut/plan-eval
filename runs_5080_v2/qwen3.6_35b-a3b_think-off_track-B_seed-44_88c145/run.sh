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
    BAI="results/${SAMPLE}.bam.bai"
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    VCF_TBI="results/${SAMPLE}.vcf.gz.tbi"

    # Skip if BAM already exists (idempotency)
    if [ -f "$BAI" ]; then
        continue
    fi

    # Mapping
    bwa mem -t 4 "$REF" "$R1" "$R2" | \
        samtools view -bS - | \
        samtools sort -o "$BAM" -

    # Index BAM
    samtools index "$BAM"

    # Variant calling with lofreq
    # lofreq call requires a reference and a BAM
    lofreq call -f "$REF" -r "$REF" -o "${SAMPLE}.vcf" -t 4 "$BAM" 2>/dev/null || \
    lofreq call -f "$REF" -o "${SAMPLE}.vcf" -t 4 "$BAM" 2>/dev/null || \
    lofreq call -f "$REF" -o "${SAMPLE}.vcf" "$BAM"

    # Compress and index VCF
    bgzip -c "${SAMPLE}.vcf" > "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"

    # Clean up uncompressed VCF
    rm -f "${SAMPLE}.vcf"
done

# Create collapsed.tsv
# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv

# Extract variants from each VCF and format
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="results/${SAMPLE}.vcf.gz"
    
    # Skip if VCF doesn't exist
    if [ ! -f "$VCF_GZ" ]; then
        continue
    fi

    # Use bcftools to extract fields: CHROM, POS, REF, ALT, INFO/AF
    # bcftools query format: %sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n
    # Note: AF might be missing for some variants, so we handle that
    bcftools query -f "%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" 2>/dev/null | \
    while IFS=$'\t' read -r sample chrom pos ref alt af; do
        # Handle missing AF
        if [ -z "$af" ]; then
            af="."
        fi
        # Handle multi-allelic sites: bcftools query might output comma-separated values
        # We need to split them if necessary, but for simple collapsed table, 
        # we often just take the first or format appropriately.
        # Standard practice for collapsed tables with multi-allelic sites is to split rows.
        # However, bcftools query with %ALT gives comma separated.
        # Let's split ALT and corresponding AF if AF is also comma separated.
        
        # Simple approach: just output as is if single allele, or split if multi
        # For robustness, let's assume single allele per row for simplicity unless specified otherwise.
        # But lofreq can produce multi-allelic.
        # Let's use a more robust parsing if needed, but for this task, simple output is likely expected.
        
        # If ALT has commas, we need to split.
        if [[ "$alt" == *","* ]]; then
            IFS=',' read -ra ALT_ARR <<< "$alt"
            if [[ "$af" == *","* ]]; then
                IFS=',' read -ra AF_ARR <<< "$af"
                for i in "${!ALT_ARR[@]}"; do
                    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "${ALT_ARR[$i]}" "${AF_ARR[$i]}"
                done
            else
                # If AF is not comma separated but ALT is, it's ambiguous. 
                # Usually AF is also comma separated for multi-allelic.
                # If not, we might just output the first ALT or the whole string.
                # Let's output the first ALT and the first AF (if available) or just the whole ALT.
                # For this script, we'll just output the first ALT and corresponding AF if available.
                printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "${ALT_ARR[0]}" "$af"
            fi
        else
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "$alt" "$af"
        fi
    done >> results/collapsed.tsv
done