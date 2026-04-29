#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Step 1: BWA index
[[ -f data/ref/chrM.fa.bwt ]] || bwa index data/ref/chrM.fa

# Step 2: samtools faidx
[[ -f data/ref/chrM.fa.fai ]] || samtools faidx data/ref/chrM.fa

# Steps 3-6: Per-sample alignment, indexing, calling, compression
for sample in "${SAMPLES[@]}"; do
  # Step 3: Alignment and sort
  [[ -f results/"${sample}".bam ]] || {
    bwa mem -t "${THREADS}" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/"${sample}"_1.fq.gz data/raw/"${sample}"_2.fq.gz | samtools sort -@ "${THREADS}" -o results/"${sample}".bam -
  }

  # Step 4: BAM index
  [[ -f results/"${sample}".bam.bai ]] || samtools index -@ "${THREADS}" results/"${sample}".bam

  # Step 5: Variant calling
  [[ -f results/"${sample}".vcf || -f results/"${sample}".vcf.gz ]] || lofreq call-parallel --pp-threads "${THREADS}" -f data/ref/chrM.fa -o results/"${sample}".vcf results/"${sample}".bam

  # Step 6: Compression and indexing
  [[ -f results/"${sample}".vcf.gz.tbi ]] || {
    bgzip -f results/"${sample}".vcf
    tabix -p vcf results/"${sample}".vcf.gz
  }
done

# Step 7: Collapsed TSV (rebuild every run)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/"${sample}".vcf.gz | awk -v s="${sample}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done