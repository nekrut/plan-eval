#!/usr/bin/env bash
set -euo pipefail

# Configuration
THREADS=4
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create results directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Step 2: Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Step 3-7: Per-sample processing
for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${OUT_DIR}/${sample}.bam"
    BAI_FILE="${BAM_FILE}.bai"
    VCF_GZ_FILE="${OUT_DIR}/${sample}.vcf.gz"
    TBI_FILE="${VCF_GZ_FILE}.tbi"
    
    # Idempotency check: if final index exists and is newer than BAM, skip
    if [[ -f "$TBI_FILE" ]] && [[ "$TBI_FILE" -nt "$BAM_FILE" ]] && [[ "$TBI_FILE" -nt "${RAW_DIR}/${sample}_1.fq.gz" ]]; then
        continue
    fi
    
    # Step 3 & 4: Alignment and sorting
    # Note: -R argument uses literal \t as required by bwa mem
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "${RAW_DIR}/${sample}_1.fq.gz" "${RAW_DIR}/${sample}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$BAM_FILE" -
    
    # Step 5: BAM indexing
    samtools index -@ "$THREADS" "$BAM_FILE"
    
    # Step 6: Variant calling with lofreq
    # Output to temporary uncompressed VCF
    VCF_TMP="${OUT_DIR}/${sample}.vcf"
    lofreq call-parallel --pp-threads "$THREADS" --ref "$REF" --out "$VCF_TMP" "$BAM_FILE"
    
    # Step 7: VCF compression and indexing
    bgzip -c "$VCF_TMP" > "$VCF_GZ_FILE"
    rm -f "$VCF_TMP"
    tabix -p vcf "$VCF_GZ_FILE"
done

# Step 8: Collapse step
COLLAPSED_FILE="${OUT_DIR}/collapsed.tsv"
# Check if any VCF is newer than the collapsed file to decide if rebuild is needed
NEED_REBUILD=false
if [[ ! -f "$COLLAPSED_FILE" ]]; then
    NEED_REBUILD=true
else
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ_FILE="${OUT_DIR}/${sample}.vcf.gz"
        if [[ "$VCF_GZ_FILE" -nt "$COLLAPSED_FILE" ]]; then
            NEED_REBUILD=true
            break
        fi
    done
fi

if [[ "$NEED_REBUILD" == true ]]; then
    # Create header
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED_FILE"
    
    # Append data from all samples
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ_FILE="${OUT_DIR}/${sample}.vcf.gz"
        # bcftools query format: sample name (literal) + tab + CHROM + tab + POS + tab + REF + tab + ALT + tab + INFO/AF
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ_FILE" >> "$COLLAPSED_FILE"
    done
fi