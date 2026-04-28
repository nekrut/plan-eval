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
    BAM="$OUT_DIR/${SAMPLE}.bam"
    VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
    VCF_TBI="$OUT_DIR/${SAMPLE}.vcf.gz.tbi"

    if [[ -f "$VCF_TBI" && "$VCF_TBI" -nt "$R1" && "$VCF_TBI" -nt "$R2" ]]; then
        continue
    fi

    if [[ ! -f "$BAM" || "$BAM" -nt "$R1" || "$BAM" -nt "$R2" ]]; then
        RG_STR="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
        bwa mem -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2" | \
        samtools sort -@ "$THREADS" -o "$BAM"
        samtools index -@ "$THREADS" "$BAM"
    fi

    if [[ ! -f "$VCF_GZ" || "$VCF_GZ" -nt "$BAM" ]]; then
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref "$REF" --out "$OUT_DIR/${SAMPLE}.vcf" \
            --sig --bonf "$BAM"
        
        bgzip -f "$OUT_DIR/${SAMPLE}.vcf"
        tabix -p vcf "$VCF_GZ"
    fi
done

TSV="$OUT_DIR/collapsed.tsv"
NEEDS_TSV=true
for SAMPLE in "${SAMPLES[@]}"; do
    VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
    if [[ ! -f "$VCF_GZ" || ( -f "$TSV" && "$VCF_GZ" -nt "$TSV" ) ]]; then
        NEEDS_TSV=true
    fi
done

if [ "$NEEDS_TSV" = true ]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$TSV"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF_GZ="$OUT_DIR/${SAMPLE}.vcf.gz"
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCF_GZ" >> "$TSV"
    done
fi

exit 0