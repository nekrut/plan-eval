#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
OUT_DIR="results"

mkdir -p "$OUT_DIR"

# Step 2: Reference indexing (idempotent)
if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

# Step 3-7: Per-sample processing
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    BAI="${OUT_DIR}/${SAMPLE}.bam.bai"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    # Check if final VCF index exists; if so, skip entire sample pipeline
    if [[ -f "$TBI" ]]; then
        continue
    fi

    # Step 3 & 4: Alignment and sorting
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "$REF" \
            "${RAW_DIR}/${SAMPLE}_1.fq.gz" \
            "${RAW_DIR}/${SAMPLE}_2.fq.gz" | \
        samtools sort -@ "$THREADS" -o "$BAM" -
    fi

    # Step 5: BAM indexing
    if [[ ! -f "$BAI" ]]; then
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Step 6: Variant calling
    VCF="${OUT_DIR}/${SAMPLE}.vcf"
    if [[ ! -f "$VCF_GZ" ]]; then
        # If VCF.gz doesn't exist, we need to call variants.
        # We check if the uncompressed VCF exists to avoid re-calling if only compression failed previously.
        if [[ ! -f "$VCF" ]]; then
            lofreq call-parallel \
                --pp-threads "$THREADS" \
                --verbose \
                --ref "$REF" \
                --out "$VCF" \
                "$BAM"
        fi

        # Step 7: Compress and index VCF
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm -f "$VCF"
    fi
done

# Step 8: Collapse step
COLLAPSED="${OUT_DIR}/collapsed.tsv"
# Check if collapsed.tsv needs rebuilding
NEED_REBUILD=0
if [[ ! -f "$COLLAPSED" ]]; then
    NEED_REBUILD=1
else
    # Check if any VCF.gz is newer than collapsed.tsv
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        if [[ "$VCF_GZ" -nt "$COLLAPSED" ]]; then
            NEED_REBUILD=1
            break
        fi
    done
fi

if [[ "$NEED_REBUILD" -eq 1 ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for SAMPLE in "${SAMPLES[@]}"; do
            VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
            bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED"
fi