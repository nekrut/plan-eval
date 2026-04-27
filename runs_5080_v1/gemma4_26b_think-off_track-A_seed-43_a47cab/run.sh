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
        if [[ ! -f "$BAM" ]]; then
            RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
            bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2" | \
            samtools sort -@ "$THREADS" -o "$BAM"
        fi

        if [[ ! -f "${BAM}.bai" ]]; then
            samtools index -@ "$THREADS" "$BAM"
        fi

        VCF_TMP="${OUT_DIR}/${SAMPLE}.vcf"
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF_TMP" "$BAM"
        
        bgzip -c "$VCF_TMP" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
        rm "$VCF_TMP"
    fi
done

TSV_OUT="${OUT_DIR}/collapsed.tsv"
NEEDS_TSV=true

for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
    if [[ -f "$TSV_OUT" ]]; then
        if [[ "$VCF_GZ" -nt "$TSV_OUT" ]]; then
            NEEDS_TSV=true
            break
        fi
    fi
done

if [ "$NEEDS_TSV" = true ]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$TSV_OUT"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="${OUT_DIR}/${SAMPLE}.vcf.gz"
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" >> "$TSV_OUT"
    done
fi

exit 0