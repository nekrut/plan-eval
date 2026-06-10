#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_FA="data/ref/chrM.fa"

mkdir -p "$RESULTS_DIR"

# Reference indexing
if [ ! -f "data/ref/chrM.fa.fai" ]; then
    samtools faidx "$REF_FA"
fi
if [ ! -f "data/ref/chrM.fa.amb" ]; then
    bwa index "$REF_FA"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${SAMPLE}.bam"
    BAI="$RESULTS_DIR/${SAMPLE}.bam.bai"
    VCF_RAW="$RESULTS_DIR/${SAMPLE}.vcf"
    VCF_GZ="$RESULTS_DIR/${SAMPLE}.vcf.gz"
    TBI="$RESULTS_DIR/${SAMPLE}.vcf.gz.tbi"

    # Alignment and Sorting
    if [ ! -f "$BAI" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
            "$REF_FA" \
            "data/raw/${SAMPLE}_1.fq.gz" \
            "data/raw/${SAMPLE}_2.fq.gz" | \
            samtools sort -@ "$THREADS" -o "$BAM"
        
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Variant Calling
    if [ ! -f "$TBI" ]; then
        lofreq call-parallel --pp-threads "$THREADS" \
            --ref "$REF_FA" \
            "$BAM" > "$VCF_RAW"
        
        bgzip -c "$VCF_RAW" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF_RAW"
    fi
done

# Collapsed Table
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
# Check if any input VCF is newer than the TSV or if TSV doesn't exist
NEEDS_COLLAPSE=false
if [ ! -f "$COLLAPSED" ]; then
    NEEDS_COLLAPSE=true
else
    for SAMPLE in "${SAMPLES[@]}"; do
        if [ "$(find -newer "$RESULTS_DIR/${SAMPLE}.vcf.gz" "$COLLAPSED" 2>/dev/null | grep -q "$COLLAPSED")" ]; then
            # This logic is slightly complex in bash, simpler to check if any VCF is newer than TSV
            # We'll use a simple timestamp comparison for each sample
            if [ "$(1000*$(date +%s))" -gt "$(1000*$(date -r "$RESULTS_DIR/${SAMPLE}.vcf.gz" +%s))" ]; then # This is not quite right, let's use a loop check
                : 
            fi
        fi
    done
    # Re-evaluating logic: if any VCF is newer than the TSV file, we must rebuild.
    for SAMPLE in "${SAMPLES[@]}"; do
        if [ "$(ls -tr "$RESULTS_DIR/${SAMPLE}.vcf.gz" "$COLLAPSED" 2>/dev/null | grep -q "$COLLAPSED")" ]; then
             # If the TSV is older than any of the VCFs, we need to rebuild.
             # Since we can't easily do multi-file comparison in one go without complex logic:
             if [ "$(date -r "$RESULTS_DIR/${SAMPLE}.vcf.gz" +%s)" -gt "$(date -r "$COLLAPSED" +%s)" ]; then
                NEEDS_COLLAPSE=true
                break
             fi
        fi
    done
fi

if [ "$NEEDS_COLLAPSE" = true ] || [ ! -f "$COLLAPSED" ]; then
    echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RESULTS_DIR/${SAMPLE}.vcf.gz" >> "$COLLAPSED"
    done
fi