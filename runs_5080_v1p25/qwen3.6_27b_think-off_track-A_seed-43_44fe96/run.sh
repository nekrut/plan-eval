#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"

mkdir -p "$OUT_DIR"

# 2. Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# 3-7. Per-sample alignment, sorting, indexing, calling, compression
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Idempotency guard: if final VCF index exists, skip sample
    if [[ -f "$TBI" ]]; then
        continue
    fi

    # Alignment and sorting
    if [[ ! -f "$BAM" ]] || [[ ! -f "$BAI" ]]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "$REF" \
            "${RAW_DIR}/${SAMPLE}_1.fq.gz" \
            "${RAW_DIR}/${SAMPLE}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$BAM"

        samtools index -@ "$THREADS" "$BAM"
    fi

    # Variant calling
    VCF="${OUT_DIR}/${SAMPLE}.vcf"
    if [[ ! -f "$VCF_GZ" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF" "$BAM"

        # Compression and indexing
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"

        # Remove intermediate uncompressed VCF
        rm -f "$VCF"
    fi
done

# 8. Collapse step
COLLAPSED="${OUT_DIR}/collapsed.tsv"
NEED_COLLAPSE=false

# Check if collapsed.tsv exists and is newer than all VCFs
if [[ -f "$COLLAPSED" ]]; then
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        if [[ "$VCF_GZ" -nt "$COLLAPSED" ]]; then
            NEED_COLLAPSE=true
            break
        fi
    done
else
    NEED_COLLAPSE=true
fi

if $NEED_COLLAPSE; then
    # Create temporary file for body
    TMP_BODY=$(mktemp)
    trap 'rm -f "$TMP_BODY"' EXIT

    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" >> "$TMP_BODY"
    done

    # Write header and body
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        cat "$TMP_BODY"
    } > "$COLLAPSED"
fi