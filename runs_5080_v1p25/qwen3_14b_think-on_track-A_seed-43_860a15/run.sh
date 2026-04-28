#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

if [ ! -f data/ref/chrM.fa.amb ]; then
    bwa index data/ref/chrM.fa
fi

for sample in "${samples[@]}"; do
    bam="results/${sample}.bam"
    if [ ! -f "$bam" ]; then
        bwa mem -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o "$bam"
    fi

    bai="results/${sample}.bam.bai"
    if [ ! -f "$bai" ]; then
        samtools index -@ $THREADS "$bam"
    fi

    vcf="results/${sample}.vcf"
    if [ ! -f "$vcf" ]; then
        lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o "$vcf" "$bam"
    fi

    vcf_gz="results/${sample}.vcf.gz"
    if [ ! -f "$vcf_gz" ]; then
        bgzip "$vcf"
        rm "$vcf"
    fi

    tbi="results/${sample}.vcf.gz.tbi"
    if [ ! -f "$tbi" ]; then
        tabix -p vcf "$vcf_gz"
    fi
done

collapsed="results/collapsed.tsv"
if [ ! -f "$collapsed" ] || ( [ -n "$(find results -name "*.vcf.gz" -newer "$collapsed")" ] ); then
    > "$collapsed"
    for sample in "${samples[@]}"; do
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/${sample}.vcf.gz >> "$collapsed"
    done
    sed -i '1s/.*/sample\tchrom\tpos\tref\talt\taf/' "$collapsed"
fi