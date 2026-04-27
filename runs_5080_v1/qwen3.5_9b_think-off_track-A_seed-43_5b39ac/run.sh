#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"
RAW="data/raw"
RES="results"

mkdir -p "$RES"

REF_IDX="$REF.fai"
REF_BWA="chrM.amb"
REF_BWA_ANN="chrM.ann"
REF_BWA_BWT="chrM.bwt"
REF_BWA_PAC="chrM.pac"
REF_BWA_SA="chrM.sa"

REF_INDEXED=0
if [[ ! -f "$REF_IDX" ]] || [[ ! -f "$REF_BWA" ]]; then
    samtools faidx "$REF"
    bwa index "$REF"
    REF_INDEXED=1
fi

for SAMPLE in $SAMPLES; do
    SAMPLE1="$RAW/${SAMPLE}_1.fq.gz"
    SAMPLE2="$RAW/${SAMPLE}_2.fq.gz"
    OUT_BAM="$RES/${SAMPLE}.bam"
    OUT_BAI="$RES/${SAMPLE}.bam.bai"
    OUT_VCF="$RES/${SAMPLE}.vcf"
    OUT_VCF_GZ="$RES/${SAMPLE}.vcf.gz"
    OUT_VCF_GZ_TBI="$RES/${SAMPLE}.vcf.gz.tbi"

    if [[ -f "$OUT_VCF_GZ_TBI" ]]; then
        continue
    fi

    if [[ -f "$OUT_BAM" ]]; then
        rm -f "$OUT_BAM" "$OUT_BAI" "$OUT_VCF" "$OUT_VCF_GZ" "$OUT_VCF_GZ_TBI"
    fi

    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -R "$RG" "$REF" "$SAMPLE1" "$SAMPLE2" 2>/dev/null | samtools sort -@ "$THREADS" -o "$OUT_BAM"
    samtools index -@ "$THREADS" "$OUT_BAM"

    lofreq call-parallel -f "$REF" -o "$OUT_VCF" --pp-threads "$THREADS" "$OUT_BAM" 2>/dev/null

    bgzip -f "$OUT_VCF"
    tabix -p vcf "$OUT_VCF_GZ"
    rm -f "$OUT_VCF"
done

{
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for SAMPLE in $SAMPLES; do
        bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$RES/${SAMPLE}.vcf.gz" 2>/dev/null
    done
} > "$RES/collapsed.tsv"