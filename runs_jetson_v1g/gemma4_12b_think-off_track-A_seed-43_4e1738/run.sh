#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF_FA="data/ref/chrM.fa"
RESULTS_DIR="results"

mkdir -p "$RESULTS_DIR"

# 2. Reference indexing
if [ ! -f "data/ref/chrM.fa.fai" ]; then
    samtools faidx "$REF_FA"
fi

if [ ! -f "data/ref/chrM.fa.amb" ]; then
    bwa index "$REF_FA"
fi

# Process samples
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${SAMPLE}.bam"
    VCF_RAW="$RESULTS_DIR/${SAMPLE}.vcf"
    VCF_GZ="$RESULTS_DIR/${SAMPLE}.vcf.gz"
    VCF_TBI="$RESULTS_DIR/${SAMPLE}.vcf.gz.tbi"

    # 3 & 4. Alignment and Sorting
    if [ ! -f "$BAM" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "$REF_FA" \
            "data/raw/${SAMPLE}_1.fq.gz" \
            "data/raw/${SAMPLE}_2.fq.gz" | \
            samtools sort -@ "$THREADS" -o "$BAM"
    fi

    # 5. BAM indexing
    if [ ! -f "${BAM}.bai" ]; then
        samtools index -@ "$THREADS" "$BAM"
    fi

    # 6. Variant calling
    if [ ! -f "$VCF_RAW" ]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref "$REF_FA" \
            --out "$VCF_RAW" \
            --sig \
            --bonf \
            "$BAM"
    fi

    # 7. VCF compression and indexing
    if [ ! -f "$VCF_TBI" ]; then
        bgzip -c "$VCF_RAW" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF_RAW"
    fi
done

# 8. Collapse step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
# Check if any input VCF is newer than the TSV or if TSV doesn't exist
NEEDS_REBUILD=false
if [ ! -f "$COLLAPSED" ]; then
    NEEDS_REBUILD=true
else
    for SAMPLE in "${SAMPLES[@]}"; do
        if [ "$(find "$RESULTS_DIR/${SAMPLE}.vcf.gz" -newer "$COLLAPSED" 2>/dev/null)" ]; then
            NEEDS_REBUILD=true
            break
        fi
    done
fi

if [ "$NEEDS_REBUILD" = true ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RESULTS_DIR/${SAMPLE}.vcf.gz" >> "$COLLAPSED"
    done
fi