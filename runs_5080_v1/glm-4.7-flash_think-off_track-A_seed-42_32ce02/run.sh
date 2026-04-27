#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RAW="data/raw"

mkdir -p results

# 2. Reference indexing (once)
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "${REF}"
fi
if [[ ! -f "${REF}.amb" ]]; then
    bwa index "${REF}"
fi

# 3. Per-sample alignment
for sample in "${SAMPLES[@]}"; do
    RG="@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA"
    bwa mem -t "${THREADS}" -R "${RG}" "${REF}" "${RAW}/${sample}_1.fq.gz" "${RAW}/${sample}_2.fq.gz" \
        | samtools sort -@ "${THREADS}" -o "results/${sample}.bam"
    samtools index -@ "${THREADS}" "results/${sample}.bam"
done

# 6. Variant calling
for sample in "${SAMPLES[@]}"; do
    lofreq call-parallel --pp-threads "${THREADS}" -f "${REF}" "results/${sample}.bam" -o "results/${sample}.vcf"
done

# 7. VCF compression and indexing
for sample in "${SAMPLES[@]}"; do
    bgzip -c "results/${sample}.vcf" > "results/${sample}.vcf.gz"
    tabix -p vcf "results/${sample}.vcf.gz"
    rm "results/${sample}.vcf"
done

# 8. Collapse step
{
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for sample in "${SAMPLES[@]}"; do
        bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sample}.vcf.gz"
    done
} > results/collapsed.tsv