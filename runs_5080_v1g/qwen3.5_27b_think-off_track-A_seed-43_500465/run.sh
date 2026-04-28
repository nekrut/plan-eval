#!/usr/bin/env bash
set -euo pipefail

# Configuration
THREADS=4
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
OUT_DIR="results"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Create results directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Step 2: Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Step 3-7: Per-sample processing loop
for sample in "${SAMPLES[@]}"; do
    BAM_FILE="${OUT_DIR}/${sample}.bam"
    BAI_FILE="${BAM_FILE}.bai"
    VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
    TBI_FILE="${VCF_GZ}.tbi"
    
    # Idempotency check: if final index exists, skip sample
    if [[ -f "$TBI_FILE" ]]; then
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
    lofreq call-parallel --pp-threads "$THREADS" \
        --ref "$REF" \
        --out "$VCF_TMP" \
        "$BAM_FILE"

    # Step 7: VCF compression and indexing
    bgzip -c "$VCF_TMP" > "$VCF_GZ"
    rm -f "$VCF_TMP"
    tabix -p vcf "$VCF_GZ"
done

# Step 8: Collapse variants into TSV
COLLAPSED_FILE="${OUT_DIR}/collapsed.tsv"
# Check if all VCFs exist and are newer than the TSV (or TSV doesn't exist)
NEED_REBUILD=false
if [[ ! -f "$COLLAPSED_FILE" ]]; then
    NEED_REBUILD=true
else
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
        if [[ "$VCF_GZ" -nt "$COLLAPSED_FILE" ]]; then
            NEED_REBUILD=true
            break
        fi
    done
fi

if [[ "$NEED_REBUILD" == true ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
            # bcftools query format: sample name literal, then fields
            # %INFO/AF might be missing for some variants; bcftools handles this gracefully
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED_FILE"
fi