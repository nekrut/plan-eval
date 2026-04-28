#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
OUT_DIR="results"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p "$OUT_DIR"

if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
    BAM="${OUT_DIR}/${sample}.bam"
    BAI="${OUT_DIR}/${sample}.bam.bai"
    VCF_GZ="${OUT_DIR}/${sample}.vcf.gz"
    TBI="${OUT_DIR}/${sample}.vcf.gz.tbi"
    FQ1="${RAW_DIR}/${sample}_1.fq.gz"
    FQ2="${RAW_DIR}/${sample}_2.fq.gz"

    if [[ -f "$TBI" ]] && [[ "$TBI" -nt "$FQ1" ]] && [[ "$TBI" -nt "$FQ2" ]] && [[ "$TBI" -nt "$REF" ]]; then
        continue
    fi

    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -R "$RG" "$REF" "$FQ1" "$FQ2" | samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index -@ "$THREADS" "$BAM"

    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "${OUT_DIR}/${sample}.vcf" "$BAM"

    bgzip -c "${OUT_DIR}/${sample}.vcf" > "$VCF_GZ"
    rm -f "${OUT_DIR}/${sample}.vcf"
    tabix -p vcf "$VCF_GZ"
done

COLLAPSED="${OUT_DIR}/collapsed.tsv"
VCF_LIST=()
for sample in "${SAMPLES[@]}"; do
    VCF_LIST+=("${OUT_DIR}/${sample}.vcf.gz")
done

NEWEST_VCF=0
for vcf in "${VCF_LIST[@]}"; do
    if [[ -f "$vcf" ]]; then
        ts=$(stat -c %Y "$vcf" 2>/dev/null || stat -f %m "$vcf" 2>/dev/null)
        if [[ $ts -gt $NEWEST_VCF ]]; then
            NEWEST_VCF=$ts
        fi
    fi
done

if [[ ! -f "$COLLAPSED" ]] || [[ $NEWEST_VCF -gt $(stat -c %Y "$COLLAPSED" 2>/dev/null || stat -f %m "$COLLAPSED" 2>/dev/null || echo 0) ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${OUT_DIR}/${sample}.vcf.gz"
        done
    } > "$COLLAPSED"
fi