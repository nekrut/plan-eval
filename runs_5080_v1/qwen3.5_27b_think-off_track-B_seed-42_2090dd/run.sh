#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
REF_IDX="${REF}.bwt"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Index reference if not already done
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]] || [[ ! -f "${REF}.pac" ]] || [[ ! -f "${REF}.sa" ]]; then
    bwa index "$REF"
    samtools faidx "$REF"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${BAM}.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${VCF_GZ}.tbi"

    # Skip if all outputs exist
    if [[ -f "$BAM" ]] && [[ -f "$BAI" ]] && [[ -f "$VCF_GZ" ]] && [[ -f "$TBI" ]]; then
        continue
    fi

    # Align with BWA-MEM
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools view -bS - > "$BAM"

    # Sort and index BAM
    samtools sort -o "$BAM" "$BAM"
    samtools index "$BAM"

    # Call variants with LoFreq
    lofreq call -f "$REF" -o "${OUT_DIR}/${SAMPLE}.vcf" -t "$THREADS" "$BAM"

    # Compress and index VCF
    bgzip -c "${OUT_DIR}/${SAMPLE}.vcf" > "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"

    # Cleanup temporary uncompressed VCF
    rm -f "${OUT_DIR}/${SAMPLE}.vcf"
done

# Generate collapsed TSV
COLLAPSED="${OUT_DIR}/collapsed.tsv"
if [[ -f "$COLLAPSED" ]]; then
    exit 0
fi

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"

# Extract variants from all VCFs and append to collapsed file
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    # Use bcftools to extract required fields, handling multi-allelic sites if present
    # Format: sample, chrom, pos, ref, alt, af (AF info field)
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF_GZ" | \
    awk -v sample="$SAMPLE" '{print sample "\t" $0}' >> "$COLLAPSED"
done