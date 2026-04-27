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
    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"
    BAM="${OUT_DIR}/${SAMPLE}.bam"
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    VCF_TBI="${OUT_DIR}/${SAMPLE}.vcf.gz.tbi"

    if [[ ! -f "$VCF_TBI" ]]; then
        RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
        
        if [[ ! -f "$BAM" ]]; then
            bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2" | \
            samtools sort -@ "$THREADS" -o "$BAM"
        fi

        if [[ ! -f "${BAM}.bai" ]]; then
            samtools index -@ "$THREADS" "$BAM"
        fi

        VCF_TMP="${OUT_DIR}/${SAMPLE}.vcf"
        if [[ ! -f "$VCF_GZ" ]]; then
            lofreq call-parallel --pp-threads "$THREADS" -f "$REF" "$BAM" > "$VCF_TMP"
            bgzip -f "$VCF_TMP"
        fi
    fi
done

TSV_OUT="${OUT_DIR}/collapsed.tsv"
NEEDS_TSV=true

if [[ -f "$TSV_OUT" ]]; then
    for SAMPLE in "${SAMPLES[@]}"; do
        if [[ "${SAMPLES[0]}.vcf.gz" -nt "$TSV_OUT" ]] || [[ "${SAMPLES[0]}.vcf.gz.tbi" -nt "$TSV_OUT" ]]; then
            # Check if any sample VCF is newer than the TSV
            if [[ "${OUT_DIR}/${SAMPLE}.vcf.gz" -nt "$TSV_OUT" ]]; then
                NEEDS_TSV=false
                break
            fi
        fi
    done
fi

if [[ "$NEEDS_TSV" = true ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$TSV_OUT"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "${OUT_DIR}/${SAMPLE}.vcf.gz" >> "$TSV_OUT"
    done
fi

exit 0