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
        VCF_UNCOMPRESSED="$OUT_DIR/${SAMPLE}.vcf"
        lofreq call-parallel --pp-threads "$THREADS" --verbose \
            --ref "$REF" --out "$VCF_UNCOMPRESSED" \
            --sig --bonf \
            "$BAM"
        
        bgzip -f "$VCF_UNCOMPRESSED"
        tabix -p vcf "$VCF_GZ"
    fi
done

TSV_FINAL="$OUT_DIR/collapsed.tsv"
NEEDS_TSV=false
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ ! -f "$OUT_DIR/${SAMPLE}.vcf.gz.tbi" ]]; then
        NEEDS_TSV=true
        break
    fi
    if [[ ! -f "$TSV_FINAL" || "$TSV_FINAL" -nt "$OUT_DIR/${SAMPLE}.vcf.gz.tbi" ]]; then
        NEEDS_TSV=true
        break
    fi
done

if [ "$NEEDS_TSV" = true ]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$TSV_FINAL"
    for SAMPLE in "${SAMPLES[@]}"; do
        bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$OUT_DIR/${SAMPLE}.vcf.gz" >> "$TSV_FINAL"
    done
fi

exit 0