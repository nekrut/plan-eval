#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"

mkdir -p results

# 2. Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx "$REF"
fi
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index "$REF"
fi

# 3. Per-sample alignment with bwa mem
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/${sample}.bam ]]; then
        bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
            data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz \
            | samtools sort -@ "$THREADS" -o results/${sample}.bam
    fi
done

# 5. BAM indexing
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ "$THREADS" results/${sample}.bam
    fi
done

# 6. Variant calling with lofreq call-parallel
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f results/${sample}.vcf.gz ]]; then
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" results/${sample}.bam \
            -o results/${sample}.vcf
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
            bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz
        done
    } > results/collapsed.tsv
fi