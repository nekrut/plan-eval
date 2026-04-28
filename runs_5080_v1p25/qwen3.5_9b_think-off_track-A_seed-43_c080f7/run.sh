#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES="M117-bl M117-ch M117C1-bl M117C1-ch"
mkdir -p results
REF=data/ref/chrM.fa
REF_IDX=chrM.fa.fai
BWA_IDX=data/ref/chrM.fa.amb
REF_INDEXED=1
if [[ ! -f "$REF_IDX" ]]; then
    samtools faidx "$REF"
    REF_INDEXED=0
fi
if [[ ! -f "$BWA_IDX" ]]; then
    bwa index "$REF"
    REF_INDEXED=0
fi
for SAMPLE in $SAMPLES; do
    F1="data/raw/${SAMPLE}_1.fq.gz"
    F2="data/raw/${SAMPLE}_2.fq.gz"
    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA"
    if [[ ! -f "results/${SAMPLE}.vcf.gz.tbi" ]] || [[ $(stat -c %y "results/${SAMPLE}.vcf.gz.tbi" | cut -d. -f1) != $(stat -c %y "$F1" | cut -d. -f1) ]]; then
        bwa mem -t $THREADS -R "$RG" "$REF" "$F1" "$F2" | samtools sort -@ $THREADS -o "results/${SAMPLE}.bam"
        samtools index -@ $THREADS "results/${SAMPLE}.bam"
        lofreq call-parallel --pp-threads 4 -f "$REF" -o "results/${SAMPLE}.vcf" "results/${SAMPLE}.bam"
        bgzip "results/${SAMPLE}.vcf"
        tabix -p vcf "results/${SAMPLE}.vcf.gz"
        rm "results/${SAMPLE}.vcf"
    fi
done
HEADER="sample\tchrom\tpos\tref\talt\taf"
TMPFILE=$(mktemp)
for SAMPLE in $SAMPLES; do
    if [[ -f "results/${SAMPLE}.vcf.gz.tbi" ]]; then
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${SAMPLE}.vcf.gz" >> "$TMPFILE"
    fi
done
echo "$HEADER" > results/collapsed.tsv
cat "$TMPFILE" >> results/collapsed.tsv
rm "$TMPFILE"