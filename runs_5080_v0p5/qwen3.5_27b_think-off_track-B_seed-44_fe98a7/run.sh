#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
REF_BASENAME="chrM"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory
mkdir -p "$OUT_DIR"

# Index reference if not already done
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi
if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Skip if all outputs exist
    if [[ -f "$BAM" && -f "$BAI" && -f "$VCF_GZ" && -f "$TBI" ]]; then
        continue
    fi

    # Alignment with BWA-MEM
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" "$REF" "$R1" "$R2" | \
        samtools view -b -@ "$THREADS" - > "$BAM"
    fi

    # Sort and index BAM
    if [[ ! -f "$BAI" ]]; then
        samtools sort -@ "$THREADS" -o "${BAM}.sorted" "$BAM"
        mv "${BAM}.sorted" "$BAM"
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Variant calling with LoFreq
    if [[ ! -f "$VCF_GZ" ]]; then
        lofreq call -f "$REF" -o "${OUT_DIR}/${SAMPLE}.vcf" -t "$THREADS" "$BAM"
        bgzip -c "${OUT_DIR}/${SAMPLE}.vcf" > "$VCF_GZ"
        rm -f "${OUT_DIR}/${SAMPLE}.vcf"
    fi

    # Index VCF
    if [[ ! -f "$TBI" ]]; then
        tabix -p vcf "$VCF_GZ"
    fi
done

# Generate collapsed.tsv
COLLAPSED="${OUT_DIR}/collapsed.tsv"
if [[ -f "$COLLAPSED" ]]; then
    exit 0
fi

# Header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"

# Process each VCF
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    # Extract sample, chrom, pos, ref, alt, AF (INFO/AF or calculate from AD/DP)
    # LoFreq VCF format: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
    # We need to parse INFO field for AF or calculate from FORMAT
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" | \
    awk -v sample="$SAMPLE" 'BEGIN{OFS="\t"} {print sample, $1, $2, $3, $4, $5}' >> "$COLLAPSED"
done