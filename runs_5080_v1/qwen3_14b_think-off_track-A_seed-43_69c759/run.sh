#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
  samtools faidx data/ref/chrM.fa
fi

if [ ! -f data/ref/chrM.fa.amb ]; then
  bwa index data/ref/chrM.fa
fi

for sample in "${SAMPLES[@]}"; do
  BAM="results/${sample}.bam"
  BAI="results/${sample}.bam.bai"
  VCF_GZ="results/${sample}.vcf.gz"
  VCF_TBI="results/${sample}.vcf.gz.tbi"

  # Skip if all outputs already exist and are up to date
  if [ -f "$BAM" ] && [ -f "$BAI" ] && [ -f "$VCF_GZ" ] && [ -f "$VCF_TBI" ]; then
    continue
  fi

  # Alignment
  if [ ! -f "$BAM" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # Index BAM
  if [ ! -f "$BAI" ]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  # Variant calling
  if [ ! -f "$VCF_GZ" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${sample}.vcf "$BAM"
    bgzip -f results/${sample}.vcf
    mv results/${sample}.vcf.gz "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"
    rm results/${sample}.vcf
  fi
done

# Collapse VCFs into TSV
if [ ! -f results/collapsed.tsv ]; then
  cat $(printf "results/%s.vcf.gz\n" "${SAMPLES[@]}") | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
    sed "s/^/${SAMPLES[0]}\t/" > results/collapsed.tsv
fi

# Re-check all samples to ensure all outputs are present
for sample in "${SAMPLES[@]}"; do
  BAM="results/${sample}.bam"
  BAI="results/${sample}.bam.bai"
  VCF_GZ="results/${sample}.vcf.gz"
  VCF_TBI="results/${sample}.vcf.gz.tbi"

  if [ ! -f "$BAM" ] || [ ! -f "$BAI" ] || [ ! -f "$VCF_GZ" ] || [ ! -f "$VCF_TBI" ]; then
    exit 1
  fi
done

exit 0