#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# 2. Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa
fi
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa
fi

# 3-5. Per-sample alignment, sorting, and indexing
for sample in "${SAMPLES[@]}"; do
    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    if [[ ! -f results/${sample}.bam.bai ]]; then
        bwa mem -t "$THREADS" -R "$RG" data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz \
            | samtools sort -@ "$THREADS" -o results/${sample}.bam
        samtools index -@ "$THREADS" results/${sample}.bam
    fi
done

# 6-7. Variant calling, compression, and indexing
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
        bgzip -c results/${sample}.vcf > results/${sample}.vcf.gz
        tabix -p vcf results/${sample}.vcf.gz
        rm results/${sample}.vcf
    fi
done

# 8. Collapse step
if [[ ! -f results/collapsed.tsv ]]; then
    {
        echo -e "sample\tchrom\tpos\tref\talt\taf"
        for sample in "${SAMPLES[@]}"; do
            bcftools query -f '${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz
        done
    } > results/collapsed.tsv
fi