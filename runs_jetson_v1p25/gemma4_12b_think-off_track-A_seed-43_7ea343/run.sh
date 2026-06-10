#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_FA="data/ref/chrM.fa"

mkdir -p "$RESULTS_DIR"

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx "$REF_FA"
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
    bwa index "$REF_FA"
fi

for sample in "${SAMPLES[@]}"; do
    BAM="$RESULTS_DIR/${sample}.bam"
    BAI="$RESULTS_DIR/${sample}.bam.bai"
    VCF_RAW="$RESULTS_DIR/${sample}.vcf"
    VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
    VCF_TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"

    # Alignment and Sorting
    if [ ! -f "$BAI" ]; then
        bwa mem -t "$THREADS" \
            -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
            samtools sort -@ "$THREADS" -o "$BAM"
        samtools index -@ "$THREADS" "$BAM"
    fi

    # Variant Calling
    if [ ! -f "$VCF_TBI" ]; then
        lofreq call-parallel --pp-threads "$THREADS" \
            -f "$REF_FA" \
            -o "$VCF_RAW" \
            "$BAM"
        
        # Compression and Indexing (using bgzip via samtools or directly if available)
        # Since bgzip is standard for VCF, we use it to ensure compatibility with tabix.
        # If bgzip isn't explicitly in inventory but tabix/samtools are, 
        # usually bcftools view -O z is the fallback, but instructions specify bgzip.
        # However, since only samtools/bcftools/tabix are listed, and standard 
        # bioinformatics envs include bgzip as part of htslib (which tabix uses).
        bgzip -c "$VCF_RAW" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF_RAW"
    fi
done

# Collapse Step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
# Check if any input VCF is newer than the TSV
NEEDS_COLLAPSE=false
for sample in "${SAMPLES[@]}"; do
    if [ "$(find "$RESULTS_DIR/${sample}.vcf.gz" -newer "$COLLAPSED" 2>/dev/null || true)" ]; then
        # This logic is slightly flawed for "any", but we can check if file exists or is older
        :
    fi
done

# Re-evaluating idempotency for collapsed.tsv: 
# If the file doesn't exist OR any of the .vcf.gz files are newer than it, rebuild.
if [ ! -f "$COLLAPSED" ] || [ "$(ls -1_ "$RESULTS_DIR/"*.vcf.gz | xargs -I{} stat -c "%Y" {} | sort -n | tail -1)" -gt "$(stat -c "%Y" "$COLLAPSED" 2>/dev/null || echo 0)" ]; then
    # Note: The logic above is a bit complex for bash, let's simplify to "if any vcf.gz > collapsed.tsv"
    REBUILD=true
    for sample in "${SAMPLES[@]}"; do
        if [ "$(find "$RESULTS_DIR/${sample}.vcf.gz" -newer "$COLLAPSED" 2>/dev/null || true)" ]; then
            # This is a bit tricky with find, let's just check if the file exists and if it's newer than any vcf
            :
        fi
    done
    
    # Simple robust check: If file doesn't exist or any source is newer, rebuild.
    if [ ! -f "$COLLAPSED" ]; then REBUILD=true; fi
    for sample in "${SAMPLES[@]}"; do
        if [ "$(find "$RESULTS_DIR/${sample}.vcf.gz" -newer "$COLLAPSED" 2>/dev/null | wc -l)" -gt 0 ]; then
            REBUILD=true
        fi
    done

    if [ "$REBUILD" = true ]; then
        echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RESULTS_DIR/${sample}.vcf.gz" >> "$COLLAPSED"
        done
    fi
fi