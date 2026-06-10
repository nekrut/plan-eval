#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.1.h1.1.1.1.1.1 ]; then
    bwa index data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

for sample in "${SAMPLES[@]}"; do
    # Alignment and sorting
    if [ ! -f results/${sample}.bam ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
            samtools sort -@ "$THREADS" -o results/${sample}.bam -
    fi

    # BAM indexing
    if [ ! -f results/${sample}.bam.bai ]; then
        samtools index -@ "$THREADS" results/${sample}.bam
    fi

    # Variant calling
    if [ ! -f results/${sample}.vcf ]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    fi

    # VCF compression and tabix indexing
    if [ ! -f results/${sample}.vcf.gz ]; then
        bgzip -f results/${sample}.vcf
        tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Collapsed TSV generation
if [ ! -f results/collapsed.tsv ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | \
        awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$1,$2,$3,$4,$5}' >> results/collapsed.tsv
    done
fi