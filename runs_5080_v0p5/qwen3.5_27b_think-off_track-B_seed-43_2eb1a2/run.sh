#!/usr/bin/env bash
set -euo pipefail

# Configuration
REF="data/ref/chrM.fa"
REF_BASENAME="chrM"
RAW_DIR="data/raw"
OUT_DIR="results"
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Index reference if not already indexed
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# Build BWA index if not already present
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
        bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools view -bS - > "$BAM"
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

# Generate collapsed.tsv
COLLAPSED_TSV="${OUT_DIR}/collapsed.tsv"
if [[ -f "$COLLAPSED_TSV" ]]; then
    exit 0
fi

# Write header
printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED_TSV"

# Process each VCF to extract variants and calculate AF
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    if [[ ! -f "$VCF_GZ" ]]; then
        continue
    fi

    # Extract INFO fields and calculate AF
    # LoFreq VCF format: INFO contains DP (total depth) and DP4 (ref_fwd, ref_rev, alt_fwd, alt_rev)
    # AF = (alt_fwd + alt_rev) / DP
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%INFO/DP\t%INFO/DP4]\n' "$VCF_GZ" | \
    awk -v sample="$SAMPLE" '
    BEGIN { OFS="\t" }
    {
        chrom = $1
        pos = $2
        ref = $3
        alt = $4
        dp = $5
        dp4 = $6
        
        # Split DP4: ref_fwd,ref_rev,alt_fwd,alt_rev
        split(dp4, dp4_arr, ",")
        alt_fwd = dp4_arr[3]
        alt_rev = dp4_arr[4]
        
        if (dp > 0) {
            af = (alt_fwd + alt_rev) / dp
        } else {
            af = 0
        }
        
        print sample, chrom, pos, ref, alt, af
    }
    ' >> "$COLLAPSED_TSV"
done