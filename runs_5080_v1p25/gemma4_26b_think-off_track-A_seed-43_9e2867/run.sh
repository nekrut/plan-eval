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
    BAI="$OUT_DIR/${SAMPLE}.bam.bai"
    VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
    VCF_TBI="$OUT_DIR/${SAMPLE}.vcf.gz.tbi"

    if [[ ! -f "$VCF_TBI" ]] || [[ "$FQ1" -nt "$VCF_TBI" ]] || [[ "$FQ2" -nt "$VCF_TBI" ]]; then
        RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
        
        bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$FQ1" "$FQ2" | \
        samtools sort -@ "$THREADS" -o "$BAM"

        samtools index -@ "$THREADS" "$BAM"

        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$OUT_DIR/${SAMPLE}.vcf" "$BAM"

        bgzip -f "$OUT_DIR/${SAMPLE}.vcf"
        tabix -p vcf "$VCF_GZ"
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