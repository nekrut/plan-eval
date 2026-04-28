#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
OUT_DIR="results"

mkdir -p "$OUT_DIR"

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    FQ1="data/raw/${SAMPLE}_1.fq.gz"
    FQ2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="$OUT_DIR/${SAMPLE}.bam"
    BAM_BAI="$OUT_DIR/${SAMPLE}.bam.bai"
    VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
    VCF_TBI="$OUT_DIR/${SAMPLE}.vcf.gz.tbi"

    if [[ ! -f "$VCF_TBI" ]]; then
        if [[ ! -f "$BAM" || "$FQ1" -nt "$BAM" || "$FQ2" -nt "$BAM" ]]; then
            RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
            bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$FQ1" "$FQ2" | \
            samtools sort -@ "$THREADS" -o "$BAM"
        fi

        if [[ ! -f "$BAM_BAI" || "$BAM" -nt "$BAM_BAI" ]]; then
            samtools index -@ "$THREADS" "$BAM"
        fi

        VCF_UNCOMPRESSED="$OUT_DIR/${SAMPLE}.vcf"
        if [[ ! -f "$VCF_UNCOMPRESSED" || "$BAM" -nt "$VCF_UNCOMPRESSED" ]]; then
            lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF_UNCOMPRESSED" "$BAM"
        fi

        if [[ ! -f "$VCF_GZ" || "$VCF_UNCOMPRESSED" -nt "$VCF_GZ" ]]; then
            bgzip -c "$VCF_UNCOMPRESSED" > "$VCF_GZ"
            rm "$VCF_UNCOMPRESSED"
        fi

        if [[ ! -f "$VCF_TBI" || "$VCF_GZ" -nt "$VCF_TBI" ]]; then
            tabix -p vcf "$VCF_GZ"
        fi
    fi
done

TSV="$OUT_DIR/collapsed.tsv"
NEEDS_TSV=false
if [[ ! -f "$TSV" ]]; then
    NEEDS_TSV=true
else
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
        if [[ "$VCF_GZ" -nt "$TSV" ]]; then
            NEEDS_TSV=true
            break
        fi
    done
fi

if [ "$NEEDS_TSV" = true ]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$TSV"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" >> "$TSV"
    done
fi

exit 0