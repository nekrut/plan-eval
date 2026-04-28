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

    # Alignment with BWA
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools view -b - > "$BAM"
    fi

    # Sort and index BAM
    if [[ ! -f "$BAI" ]]; then
        samtools sort -o "${BAM}.sorted" "$BAM"
        mv "${BAM}.sorted" "$BAM"
        samtools index "$BAM"
    fi

    # Variant calling with LoFreq
    if [[ ! -f "$VCF_GZ" ]]; then
        lofreq call -f "$REF" -o "${OUT_DIR}/${SAMPLE}.vcf" "$BAM"
        bgzip "${OUT_DIR}/${SAMPLE}.vcf"
        mv "${OUT_DIR}/${SAMPLE}.vcf.gz" "$VCF_GZ"
    fi

    # Index VCF
    if [[ ! -f "$TBI" ]]; then
        tabix -p vcf "$VCF_GZ"
    fi
done

# Generate collapsed table
COLLAPSED="${OUT_DIR}/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        # Extract sample, chrom, pos, ref, alt, AF (INFO/AF)
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF_GZ" | \
        awk -v sample="$SAMPLE" '{print sample "\t" $0}'
    done | sort >> "$COLLAPSED"
fi