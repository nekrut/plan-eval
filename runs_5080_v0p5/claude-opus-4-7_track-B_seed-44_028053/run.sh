#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF=data/ref/chrM.fa
RAW=data/raw
OUT=results
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p "$OUT"

if [ ! -f "${REF}.bwt" ]; then
    bwa index "$REF"
fi
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "$REF"
fi

for s in "${SAMPLES[@]}"; do
    bam="$OUT/$s.bam"
    bai="$OUT/$s.bam.bai"
    vcf="$OUT/$s.vcf.gz"
    tbi="$OUT/$s.vcf.gz.tbi"

    if [ ! -s "$bam" ] || [ ! -s "$bai" ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:$s\tSM:$s\tLB:$s\tPL:ILLUMINA" \
            "$REF" "$RAW/${s}_1.fq.gz" "$RAW/${s}_2.fq.gz" \
            | samtools sort -@ "$THREADS" -o "$bam" -
        samtools index -@ "$THREADS" "$bam"
    fi

    if [ ! -s "$vcf" ] || [ ! -s "$tbi" ]; then
        tmpvcf="$OUT/$s.lofreq.vcf"
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$tmpvcf" "$bam"
        bcftools view -Oz -o "$vcf" "$tmpvcf"
        tabix -p vcf "$vcf"
        rm -f "$tmpvcf"
    fi
done

collapsed="$OUT/collapsed.tsv"
if [ ! -s "$collapsed" ]; then
    tmp="$OUT/.collapsed.tsv.tmp"
    : > "$tmp"
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$tmp"
    for s in "${SAMPLES[@]}"; do
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/$s.vcf.gz" \
            | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$1,$2,$3,$4,$5}' >> "$tmp"
    done
    mv "$tmp" "$collapsed"
fi