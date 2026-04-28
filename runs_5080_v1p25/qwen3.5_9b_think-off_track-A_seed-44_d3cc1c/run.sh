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

for s in $SAMPLES; do
    [ -f "$RES/${s}.vcf.gz.tbi" ] && continue
    [ -f "$REF_IDX" ] || samtools faidx "$REF"
    [ -f "$REF_BWA" ] || bwa index "$REF"
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$RAW/${s}_1.fq.gz" "$RAW/${s}_2.fq.gz" | \
        samtools sort -@ $THREADS -o "$RES/${s}.bam"
    samtools index -@ $THREADS "$RES/${s}.bam"
    lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "$RES/${s}.vcf" "$RES/${s}.bam"
    bgzip -f "$RES/${s}.vcf"
    tabix -p vcf "$RES/${s}.vcf.gz"
    rm -f "$RES/${s}.vcf"
done

{
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for s in $SAMPLES; do
        bcftools query -f '%sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$RES/${s}.vcf.gz"
    done
} > "$RES/collapsed.tsv"