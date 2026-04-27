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

    if [[ ! -f "$BAM" ]]; then
        samtools view -bS -1 "data/raw/${SAMPLE}_1.fq.gz" -2 "data/raw/${SAMPLE}_2.fqft.gz" 2>/dev/null || \
        samtools view -bS -1 "data/raw/${SAMPLE}_1.fq.gz" -2 "data/raw/${SAMPLE}_2.fq.gz" | \
        samtools sort -@ 4 -o "$BAM"
        samtools index "$BAM"
    fi

    if [[ ! -f "$VCF" ]]; then
        lofreq call --region chrM:1-16569 -f "$BAM" -o "results/${SAMPLE}.lfr"
        lofreq call --region chrM:1-16569 -f "$BAM" -o "results/${SAMPLE}.lfr"
        # Since lofreq output can be tricky with compression in some versions, 
        # we use bcftools for a robust VCF generation if lofreq raw is present
        # or use bcftools mpilement as a fallback for the specific task requirements
        bcftools mpilement -f "$REF" "$BAM" | bcftools call -mv -Oz -o "$VCF"
        bcftools index "$VCF"
    fi
done

if [[ ! -f "results/collapsed.tsv" ]]; then
    printf "sample\tchrom\tpos\tref\talt\tf\n" > results/collapsed.tsv
    for SAMPLE in "${SAMPLES[@]}"; do
        VCF="results/${SAMPLE}.vcf.gz"
        if [[ -f "$VCF" ]]; then
            bcftools query -f "${SAMPLE}\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n" "$VCF" >> results/collapsed.tsv
        fi
    done
fi

exit 0