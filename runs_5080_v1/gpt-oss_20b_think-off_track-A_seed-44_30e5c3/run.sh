#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

REF="data/ref/chrM.fa"
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi
if [[ ! -f "${REF}.amb" ]]; then
    bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
    fq1="data/raw/${sample}_1.fq.gz"
    fq2="data/raw/${sample}_2.fq.gz"
    bam="results/${sample}.bam"
    bai="results/${sample}.bam.bai"
    vcf="results/${sample}.vcf"
    vcfz="results/${sample}.vcf.gz"
    tbz="results/${sample}.vcf.gz.tbi"

    if [[ ! -f "$bam" || "$bam" -ot "$fq1" || "$bam" -ot "$fq2" ]]; then
        bwa mem -t "$THREADS" -R "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA" "$REF" "$fq1" "$fq2" | \
        samtools sort -@ "$THREADS" -o "$bam"
    fi

    if [[ ! -f "$bai" || "$bai" -ot "$bam" ]]; then
        samtools index -@ "$THREADS" "$bam"
    fi

    if [[ ! -f "$vcfz" || "$vcfz" -ot "$bam" ]]; then
        lofreq call-parallel -f "$REF" -i "$bam" -o "$vcf" --pp-threads "$THREADS"
        bcftools view -O z -o "$vcfz" "$vcf"
        rm -f "$vcf"
    fi

    if [[ ! -f "$tbz" || "$tbz" -ot "$vcfz" ]]; then
        tabix -p vcf "$vcfz"
    fi
done

collapsed="results/collapsed.tsv"
need=false
for sample in "${SAMPLES[@]}"; do
    tbz="results/${sample}.vcf.gz.tbi"
    if [[ ! -f "$collapsed" || "$tbz" -nt "$collapsed" ]]; then
        need=true
        break
    fi
done

if $need; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed"
    for sample in "${SAMPLES[@]}"; do
        vcfz="results/${sample}.vcf.gz"
        bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcfz" >> "$collapsed"
    done
fi

exit 0