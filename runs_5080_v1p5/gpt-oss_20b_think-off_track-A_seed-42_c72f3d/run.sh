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

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
    bam=results/${sample}.bam
    bai=results/${sample}.bam.bai
    vcf=results/${sample}.vcf
    vcfgz=results/${sample}.vcf.gz
    tbi=results/${sample}.vcf.gz.tbi

    if [ ! -f "$bam" ]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ "$THREADS" -o "$bam" -
    fi

    if [ ! -f "$bai" ]; then
        samtools index -@ "$THREADS" "$bam"
    fi

    if [ ! -f "$vcf" ]; then
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$vcf" "$bam"
    fi

    if [ ! -f "$vcfgz" ]; then
        bcftools view -Oz -o "$vcfgz" "$vcf"
    fi

    if [ ! -f "$tbi" ]; then
        tabix -p vcf "$vcfgz"
    fi
done

# Collapsed TSV
collapsed=results/collapsed.tsv
if [ ! -f "$collapsed" ]; then
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$collapsed"
    for sample in "${SAMPLES[@]}"; do
        vcfgz=results/${sample}.vcf.gz
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$vcfgz" | awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> "$collapsed"
    done
fi

exit 0