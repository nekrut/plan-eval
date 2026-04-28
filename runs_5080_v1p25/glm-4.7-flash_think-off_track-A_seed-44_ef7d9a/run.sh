#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
RAW_DIR="data/raw"
REF_DIR="data/ref"
REF_FA="$REF_DIR/chrM.fa"

mkdir -p "$RESULTS_DIR"

# 2. Reference indexing (once)
if [[ ! -f "$REF_FA.fai" ]]; then
    samtools faidx "$REF_FA"
fi
if [[ ! -f "$REF_FA.bwt" ]]; then
    bwa index "$REF_FA"
fi

# 3. Per-sample alignment
for SAMPLE in "${SAMPLES[@]}"; do
    INPUT_1="$RAW_DIR/${SAMPLE}_1.fq.gz"
    INPUT_2="$RAW_DIR/${SAMPLE}_2.fq.gz"
    BAM="$RESULTS_DIR/${SAMPLE}.bam"
    BAM_BAI="$RESULTS_DIR/${SAMPLE}.bam.bai"

    if [[ ! -f "$BAM_BAI" ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" "$INPUT_1" "$INPUT_2" | samtools sort -@ "$THREADS" -o "$BAM"
        samtools index -@ "$THREADS" "$BAM"
    fi
done

# 6. Variant calling
for SAMPLE in "${SAMPLES[@]}"; do
    VCF="$RESULTS_DIR/${SAMPLE}.vcf"
    VCF_GZ="$RESULTS_DIR/${SAMPLE}.vcf.gz"
    VCF_TBI="$RESULTS_DIR/${SAMPLE}.vcf.gz.tbi"

    if [[ ! -f "$VCF_TBI" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF_FA" -o "$VCF" "$RESULTS_DIR/${SAMPLE}.bam"
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF"
    fi
done

# 8. Collapse step
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for SAMPLE in "${SAMPLES[@]}"; do
            bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$RESULTS_DIR/${SAMPLE}.vcf.gz"
        done
    } > "$COLLAPSED"
fi