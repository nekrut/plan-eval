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
    samtools view -@ 4 -bS - | \
    samtools sort -@ 4 -o "$BAM"

    samtools index "$BAM"

    lofreq call --region chrM:1-16569 -f "$BAM" > "${VCF}.tmp"
    
    bgzip -c "${VCF}.tmp" > "$VCF"
    tabix -p vcf "$VCF"
    rm "${VCF}.tmp"
done

HEADER="sample	chrom	pos	ref	alt	af"
echo -e "$HEADER" > results/collapsed.tsv

for SAMPLE in "${SAMPLES[@]}"; do
    VCF="results/${SAMPLE}.vcf.gz"
    if [[ -f "$VCF" ]]; then
        bcftools query -f "$SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n" "$VCF" >> results/collapsed.tsv
    fi
done

# Remove duplicates if any (though bcftools query on unique VCFs shouldn't)
sort -u -o results/collapsed.tsv results/collapsed.tsv 2>/dev/null || true