#!/usr/bin/env bash
set -euo pipefail

# Globals
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
RESULTS_DIR="results"

# Create results directory
mkdir -p "$RESULTS_DIR"

# Reference indexing (once)
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${sample}_1.fq.gz"
    R2="${RAW_DIR}/${sample}_2.fq.gz"
    
    BAM="${RESULTS_DIR}/${sample}.bam"
    BAI="${RESULTS_DIR}/${sample}.bam.bai"
    VCF="${RESULTS_DIR}/${sample}.vcf"
    VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
    VCF_GZ_TBI="${RESULTS_DIR}/${sample}.vcf.gz.tbi"
    
    # Alignment and BAM sorting
    if [[ ! -f "$BAM" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "$R1" "$R2" \
            | samtools sort -@ "$THREADS" -o "$BAM"
    fi
    
    # BAM indexing
    if [[ ! -f "$BAI" ]]; then
        samtools index -@ "$THREADS" "$BAM"
    fi
    
    # Variant calling
    if [[ ! -f "$VCF_GZ" ]]; then
        if [[ -f "$VCF" ]]; then
            rm -f "$VCF"
        fi
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF" "$BAM"
        
        # VCF compression and indexing
        bgzip -c "$VCF" > "$VCF_GZ"
        rm -f "$VCF"
        tabix -p vcf "$VCF_GZ"
    else
        # Ensure index exists
        if [[ ! -f "$VCF_GZ_TBI" ]]; then
            tabix -p vcf "$VCF_GZ"
        fi
    fi
done

# Collapse step
COLLAPSED="${RESULTS_DIR}/collapsed.tsv"
REBUILD_COLLAPSED=false

if [[ ! -f "$COLLAPSED" ]]; then
    REBUILD_COLLAPSED=true
else
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
        if [[ "$VCF_GZ" -nt "$COLLAPSED" ]]; then
            REBUILD_COLLAPSED=true
            break
        fi
    done
fi

if [[ "$REBUILD_COLLAPSED" == true ]]; then
    {
        printf $'sample\tchrom\tpos\tref\talt\taf\n'
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="${RESULTS_DIR}/${sample}.vcf.gz"
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED"
fi