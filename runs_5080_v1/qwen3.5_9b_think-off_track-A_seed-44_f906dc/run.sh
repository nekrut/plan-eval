#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
REF="data/ref/chrM.fa"
RAW="data/raw"
RES="results"

mkdir -p "$RES"

REF_IDX="$REF.fai"
REF_BWA="chrM"

if [[ ! -f "$REF_IDX" ]]; then
    samtools faidx "$REF"
fi

if [[ ! -f "$REF_BWA.amb" ]]; then
    bwa index "$REF"
fi

for SAMPLE in $SAMPLES; do
    F1="$RAW/${SAMPLE}_1.fq.gz"
    F2="$RAW/${SAMPLE}_2.fq.gz"
    BAM="$RES/${SAMPLE}.bam"
    BAI="$RES/${SAMPLE}.bam.bai"
    VCF="$RES/${SAMPLE}.vcf"
    VCF_GZ="$RES/${SAMPLE}.vcf.gz"
    VCF_TBI="$RES/${SAMPLE}.vcf.gz.tbi"

    if [[ -f "$VCF_TBI" ]]; then
        continue
    fi

    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
    bwa mem -t $THREADS "$REF" "$RG" "$F1" "$F2" | samtools sort -@ $THREADS -o "$BAM"
    samtools index -@ $THREADS "$BAM"
    lofreq call-parallel -f "$REF" -o "$VCF" --pp-threads $THREADS "$BAM"
    bgzip "$VCF"
    tabix -p vcf "$VCF_GZ"
    rm -f "$VCF"
done

TSV="$RES/collapsed.tsv"
if [[ -f "$TSV" ]]; then
    exit 0
fi

{
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for SAMPLE in $SAMPLES; do
        bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$RES/${SAMPLE}.vcf.gz"
    done
} > "$TSV"