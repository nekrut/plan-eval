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
    BAI="${BAM}.bai"
    VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
    TBI="${VCF_GZ}.tbi"

    if [[ -f "$TBI" ]] && [[ "$TBI" -nt "${RAW_DIR}/${sample}_1.fq.gz" ]] && [[ "$TBI" -nt "${RAW_DIR}/${sample}_2.fq.gz" ]]; then
        continue
    fi

    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$REF" "${RAW_DIR}/${sample}_1.fq.gz" "${RAW_DIR}/${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$BAM" -

    samtools index -@ "$THREADS" "$BAM"

    lofreq call-parallel --pp-threads "$THREADS" --ref "$REF" --out "${RES_DIR}/${sample}.vcf" "$BAM"

    bgzip -c "${RES_DIR}/${sample}.vcf" > "$VCF_GZ"
    rm -f "${RES_DIR}/${sample}.vcf"

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
        printf 'sample\tchrom\tpos\tref\talt\taf\n'
        for sample in "${SAMPLES[@]}"; do
            VCF_GZ="${RES_DIR}/${sample}.vcf.gz"
            bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ"
        done
    } > "$COLLAPSED"
fi