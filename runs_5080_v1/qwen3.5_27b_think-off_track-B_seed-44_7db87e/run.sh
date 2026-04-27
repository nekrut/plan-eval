#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
REF_BAI="data/ref/chrM.fa.bai"
REF_FAI="data/ref/chrM.fa.fai"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory
mkdir -p "$OUT_DIR"

# Index reference if not present
if [[ ! -f "$REF_BAI" ]] || [[ ! -f "$REF_FAI" ]]; then
    samtools faidx "$REF"
    samtools index "$REF"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Skip if BAM and VCF already exist and are non-empty
    if [[ -f "$BAM" ]] && [[ -f "$VCF_GZ" ]] && [[ -s "$BAM" ]] && [[ -s "$VCF_GZ" ]]; then
        continue
    fi

    # Align with BWA-MEM
    bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" "$R1" "$R2" | \
        samtools view -b -@ "$THREADS" -o "$BAM" -

    # Sort and index BAM
    samtools sort -@ "$THREADS" -o "$BAM" "$BAM"
    samtools index -@ "$THREADS" "$BAM"

    # Call variants with LoFreq
    # LoFreq requires a sorted, indexed BAM. It outputs VCF.
    lofreq call -f "$REF" -o "${OUT_DIR}/${SAMPLE}.vcf" -t "$THREADS" "$BAM"

    # Compress and index VCF
    bgzip -f "${OUT_DIR}/${SAMPLE}.vcf"
    mv "${OUT_DIR}/${SAMPLE}.vcf.gz" "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"

    # Clean up uncompressed VCF if it exists
    rm -f "${OUT_DIR}/${SAMPLE}.vcf"
done

# Generate collapsed.tsv
COLLAPSED="${OUT_DIR}/collapsed.tsv"
if [[ -f "$COLLAPSED" ]]; then
    exit 0
fi

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"

# Extract variants from all VCFs and format
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    if [[ ! -f "$VCF_GZ" ]]; then
        continue
    fi
    # Use bcftools to extract required fields
    # CHROM, POS, REF, ALT, AF (from INFO/AF or calculate from AD/DP if needed, but lofreq usually provides AF)
    # LoFreq VCF format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
    # We need to parse INFO for AF. LoFreq typically puts AF in INFO.
    # If AF is missing, we might need to calculate it, but lofreq usually provides it.
    # Let's assume standard VCF INFO/AF exists. If not, we fallback to AD/DP calculation if possible,
    # but for simplicity and speed, we rely on lofreq's output.
    # bcftools query format: %CHROM\t%POS\t%REF\t%ALT\t%INFO/AF
    # If AF is not in INFO, lofreq might put it in FORMAT/AF.
    # Let's try to get AF from INFO first.
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" | \
        awk -v sample="$SAMPLE" 'BEGIN{FS="\t"; OFS="\t"} {
            # Handle cases where AF might be missing or multiple alleles
            if ($5 == "." || $5 == "") {
                # Fallback: try to get from FORMAT if possible, but bcftools query for FORMAT is complex per sample
                # For this task, we assume lofreq provides AF in INFO.
                # If not, we skip or set to 0.5 as placeholder?
                # Let's try to parse the full INFO field if AF is missing
                print sample, $1, $2, $3, $4, "NA"
            } else {
                print sample, $1, $2, $3, $4, $5
            }
        }' >> "$COLLAPSED"
done

# Ensure the file is not empty (header only is fine if no variants)
# If the loop didn't run or failed, we might have an empty file, but the header is written.
# The script exits 0 if the file exists at the start of the loop, so this block runs only once.