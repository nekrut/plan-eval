#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
REF_DIR="data/ref"
RAW_DIR="data/raw"
RES_DIR="results"

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p "$RES_DIR"

if [[ ! -f "${REF}.fai" ]] || [[ ! -f "${REF}.amb" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
    BAM="${RES_DIR}/${sample}.bam"
    BAI="${BAM}.bai"
    VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
    TBI="${VCF_GZ}.tbi"
    FQ1="${RAW_DIR}/${sample}_1.fq.gz"
    FQ2="${RAW_DIR}/${sample}_2.fq.gz"

    if [[ -f "$TBI" ]] && [[ "$TBI" -nt "$FQ1" ]] && [[ "$TBI" -nt "$FQ2" ]] && [[ "$TBI" -nt "$REF" ]]; then
        continue
    fi

    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"

    bwa mem -t "$THREADS" -R "$RG" "$REF" "$FQ1" "$FQ2" | samtools sort -@ "$THREADS" -o "$BAM" -

    samtools index -@ "$THREADS" "$BAM"

    VCF_TMP="${RES_DIR}/${sample}.vcf"
    lofreq call-parallel -f "$REF" -o "$VCF_TMP" -t "$THREADS" "$BAM"

    bgzip -c "$VCF_TMP" > "$VCF_GZ"
    rm -f "$VCF_TMP"

    tabix -p vcf "$VCF_GZ"
done

COLLAPSED="${RES_DIR}/collapsed.tsv"
NEED_REBUILD=false

if [[ ! -f "$COLLAPSED" ]]; then
    NEED_REBUILD=true
else
    for sample in "${SAMPLES[@]}"; do
        VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
        if [[ -f "$VCF_GZ" ]] && [[ "$VCF_GZ" -nt "$COLLAPSED" ]]; then
            NEED_REBUILD=true
            break
        fi
    done
fi

if [[ "$NEED_REBUILD" == true ]]; then
    {
        printf "sample\tchrom\tpos\tref\talt\taf\n"
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED"
fi