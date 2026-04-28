#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for sample in "${samples[@]}"; do
    bam_file="results/${sample}.bam"
    bai_file="results/${sample}.bam.bai"
    vcf_file="results/${sample}.vcf.gz"
    tbi_file="results/${sample}.vcf.gz.tbi"

    if [ ! -f "$bam_file" ]; then
        bwa mem -t 4 data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -o "$bam_file"
        samtools index "$bam_file"
    fi

    if [ ! -f "$vcf_file" ]; then
        lofreq call -f data/ref/chrM.fa "$bam_file" | bcftools sort -o "$vcf_file"
        tabix -p vcf "$vcf_file"
    fi
done

collapsed_tsv="results/collapsed.tsv"
if [ ! -f "$collapsed_tsv" ]; then
    > "$collapsed_tsv"
    for sample in "${samples[@]}"; do
        vcf_file="results/${sample}.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%GT\n' "$vcf_file" | awk -v sample="$sample" 'NR==1 {next} {split($5, gt, "/"); alt_count = (gt[2] == "." ? 0 : gt[2]); af = alt_count / 2; print sample, $1, $2, $3, $4, af}' >> "$collapsed_tsv"
    done
fi