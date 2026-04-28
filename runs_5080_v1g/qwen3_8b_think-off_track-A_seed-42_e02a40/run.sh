#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${samples[@]}"; do
    if [[ -f results/"${sample}".vcf.gz.tbi && -f results/"${sample}".vcf.gz && -f results/"${sample}".bam.bai && -f results/"${sample}".bam ]]; then
        continue
    fi

    if [[ ! -f data/ref/chrM.fa.fai || ! -f data/ref/chrM.fa.amb || ! -f data/ref/chrM.fa.ann || ! -f data/ref/chrM.fa.bwt || ! -f data/ref/chrM.fa.pac || ! -f data/ref/chrM.fa.sa ]]; then
        samtools faidx data/ref/chrM.fa
        bwa index data/ref/chrM.fa
    fi

    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/raw/"${sample}_1.fq.gz" data/raw/"${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o results/"${sample}".bam
    samtools index -@ "$THREADS" results/"${sample}".bam

    lofreq call-parallel --pp-threads "$THREADS" --verbose --ref data/ref/chrM.fa --out results/"${sample}".vcf --sig --bonf results/"${sample}".bam

    bgzip results/"${sample}".vcf
    tabix -p vcf results/"${sample}".vcf.gz
    rm results/"${sample}".vcf
done

for sample in "${samples[@]}"; do
    if [[ -f results/collapsed.tsv ]]; then
        if [[ -f results/"${sample}".vcf.gz.tbi && -f results/"${sample}".vcf.gz ]]; then
            continue
        fi
    fi

    bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/"${sample}".vcf.gz
done | awk 'NR == 1 {print "sample\tchrom\tpos\tref\talt\taf"} NR > 1 {print}' > results/collapsed.tsv

if [[ -f results/collapsed.tsv ]]; then
    if [[ -f results/"${samples[0]}"vcf.gz.tbi && -f results/"${samples[0]}"vcf.gz && -f results/"${samples[1]}"vcf.gz.tbi && -f results/"${samples[1]}"vcf.gz && -f results/"${samples[2]}"vcf.gz.tbi && -f results/"${samples[2]}"vcf.gz && -f results/"${samples[3]}"vcf.gz.tbi && -f results/"${samples[3]}"vcf.gz ]]; then
        exit 0
    fi
fi

exit 0