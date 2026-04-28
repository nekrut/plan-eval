#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RAW_DIR="data/raw"
RES_DIR="results"

mkdir -p "$RES_DIR"

if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.bwt" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
    BAM="${RES_DIR}/${sample}.bam"
    BAI="${RES_DIR}/${sample}.bam.bai"
    VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
    TBI="${RES_DIR}/${sample}.vcf.gz.tbi"
    FQ1="${RAW_DIR}/${sample}_1.fq.gz"
    FQ2="${RAW_DIR}/${sample}_2.fq.gz"

    if [[ -f "$TBI" ]] && [[ "$TBI" -nt "$FQ1" ]] && [[ "$TBI" -nt "$FQ2" ]] && [[ "$TBI" -nt "$REF" ]]; then
        continue
    fi

    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -R "$RG" "$REF" "$FQ1" "$FQ2" | samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index -@ "$THREADS" "$BAM"

    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "${RES_DIR}/${sample}.vcf" "$BAM"

    bgzip -c "${RES_DIR}/${sample}.vcf" > "$VCF_GZ"
    rm -f "${RES_DIR}/${sample}.vcf"
    tabix -p vcf "$VCF_GZ"
done

COLLAPSED="${RES_DIR}/collapsed.tsv"
NEWEST_VCF=0
for sample in "${SAMPLES[@]}"; do
    VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
    if [[ -f "$VCF_GZ" ]]; then
        VCF_MTIME=$(stat -c %Y "$VCF_GZ")
        if [[ "$VCF_MTIME" -gt "$NEWEST_VCF" ]]; then
            NEWEST_VCF=$VCF_MTIME
        fi
    fi
done

if [[ -f "$COLLAPSED" ]] && [[ "$NEWEST_VCF" -le $(stat -c %Y "$COLLAPSED") ]]; then
    exit 0
fi

{
    printf "sample\tchrom\tpos\tref\talt\taf\n"
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
    done
} > "$COLLAPSED"