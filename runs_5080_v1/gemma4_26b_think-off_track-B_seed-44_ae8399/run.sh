#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

REF="data/ref/chrM.fa"
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

if [[ ! -f "${REF}.bwt" ]]; then
    bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

for SAMPLE in "${SAMPLES[@]}"; do
    BAM="results/${SAMPLE}.bam"
    VCF="results/${SAMPLE}.vcf.gz"

    if [[ -f "$VCF" ]]; then
        continue
    fi

    R1="data/raw/${SAMPLE}_1.fq.gz"
    R2="data/raw/${SAMPLE}_2.fq.gz"

    bwa mem -t 4 "$REF" "$R1" "$R2" | \
    samtools view -u - | \
    samtools sort -@ 4 -o "$BAM"

    samtools index "$BAM"

    lofreq call -f "$REF" -o "$VCF" "$BAM"
    
    if [[ ! -f "${VCF}.tbi" ]]; then
        tabix -p vcf "$VCF"
    fi
    
    # Ensure vcf is compressed and indexed correctly for downstream
    if [[ ! -f "$VCF" ]]; then
        exit 1
    fi
done

COLLAPSED="results/collapsed.tsv"
if [[ ! -f "$COLLAPSED" ]]; then
    printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="results/${SAMPLE}.vcf.gz"
        if [[ -f "$VCF" ]]; then
            bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n" "$VCF" >> "$COLLAPSED"
        fi
    done
fi

exit 0