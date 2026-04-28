#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.bwt ]; then
    bwa index data/ref/chrM.fa
fi

if [ ! -f data/ref/chrM.fa.fai ]; then
    samtools faidx data/ref/chrM.fa
fi

# Check if all outputs already exist
all_exist=true
for sample in "${SAMPLES[@]}"; do
    bam=results/${sample}.bam
    bai=${bam}.bai
    vcf=results/${sample}.vcf.gz
    tbi=${vcf}.tbi
    if ! [ -f "$bam" ] || ! [ -f "$bai" ] || ! [ -f "$vcf" ] || ! [ -f "$tbi" ]; then
        all_exist=false
        break
    fi
done
if $all_exist && [ -f results/collapsed.tsv ]; then
    exit 0
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam=results/${sample}.bam
    bai=${bam}.bai
    vcf=results/${sample}.vcf.gz
    tbi=${vcf}.tbi

    # Alignment and sorting
    if ! [ -f "$bam" ]; then
        bwa mem -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
            samtools sort -@ $THREADS -o "$bam" -
    fi

    # BAM index
    if ! [ -f "$bai" ]; then
        samtools index -@ $THREADS "$bam"
    fi

    # Variant calling
    if ! [ -f "$vcf" ]; then
        lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o - "$bam" | \
            bcftools view -O z -o "$vcf" -
    fi

    # VCF index
    if ! [ -f "$tbi" ]; then
        tabix -p vcf "$vcf"
    fi
done

# Collapsed table
rm -f results/collapsed.tsv
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
    vcf=results/${sample}.vcf.gz
    if [ -f "$vcf" ]; then
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$vcf" | \
            awk -v s=${sample} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
    fi
done
exit 0