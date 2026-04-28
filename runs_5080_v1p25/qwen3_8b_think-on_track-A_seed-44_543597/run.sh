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
    bam_file="results/${sample}.bam"
    if [ -f "$bam_file" ]; then
        if [ "$bam_file" -nt data/raw/"${sample}_1.fq.gz" ] && [ "$bam_file" -nt data/raw/"${sample}_2.fq.gz" ]; then
            continue
        fi
    fi

    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$bam_file"
    samtools index -@ "$THREADS" "$bam_file"
done

for sample in "${samples[@]}"; do
    vcf_file="results/${sample}.vcf"
    vcf_gz="results/${sample}.vcf.gz"
    vcf_tbi="results/${sample}.vcf.gz.tbi"
    if [ -f "$vcf_gz" ]; then
        if [ "$vcf_gz" -nt "$vcf_file" ] && [ "$vcf_gz" -nt "$vcf_tbi" ]; then
            continue
        fi
    fi

    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$vcf_file" "results/${sample}.bam"
    bgzip "$vcf_file"
    tabix -p vcf "$vcf_gz"
    rm "$vcf_file"
done

collapsed_tsv="results/collapsed.tsv"
if [ -f "$collapsed_tsv" ]; then
    if [ "$collapsed_tsv" -nt results/*.vcf.gz ]; then
        exit 0
    fi
fi

for sample in "${samples[@]}"; do
    bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/"${sample}.vcf.gz" >> "$collapsed_tsv"
done

echo -e "sample\tchrom\tpos\tref\talt\taf" > "$collapsed_tsv"
cat results/*.vcf.gz | bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' >> "$collapsed_tsv"