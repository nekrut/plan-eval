#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
REF=data/ref/chrM.fa
if [ ! -f "${REF}.fai" ]; then
    samtools faidx "${REF}"
fi

SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for sample in "${SAMPLES[@]}"; do
    BAM="results/${sample}.bam"
    BAI="results/${sample}.bam.bai"
    VCF="results/${sample}.vcf"
    VCFGZ="results/${sample}.vcf.gz"
    VTBI="results/${sample}.vcf.gz.tbi"
    INFO="results/${sample}.info"

    if [ ! -f "$BAM" ]; then
        bwa mem -t 4 "$REF" "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | \
            samtools view -Sb - | \
            samtools sort -o "$BAM"
        samtools index "$BAM"
    fi

    if [ ! -f "$VTBI" ]; then
        if [ ! -f "$VCF" ]; then
            lofreq call -f "$REF" "$BAM" -o "$VCF"
        fi
        if [ ! -f "$VTBI" ]; then
            bgzip "$VCF"
            tabix -p vcf "$VCFGZ"
        fi
    fi

    if [ ! -f "$INFO" ]; then
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCFGZ" > "$INFO"
    fi
done

COLLAPSED="results/collapsed.tsv"
if [ ! -f "$COLLAPSED" ]; then
    > "$COLLAPSED"
    for sample in "${SAMPLES[@]}"; do
        INFO="results/${sample}.info"
        awk -v sample="$sample" '{print sample "\t" $0}' "$INFO" >> "$COLLAPSED"
    done
fi