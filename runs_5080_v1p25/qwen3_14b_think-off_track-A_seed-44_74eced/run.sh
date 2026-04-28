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
  # Skip if all outputs already exist and are up to date
  if [ -f results/${sample}.bam.bai ] && [ -f results/${sample}.vcf.gz.tbi ]; then
    continue
  fi

  # Alignment
  if [ ! -f results/${sample}.bam ]; then
    bwa mem -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ $THREADS -o results/${sample}.bam
  fi

  # BAM indexing
  if [ ! -f results/${sample}.bam.bai ]; then
    samtools index -@ $THREADS results/${sample}.bam
  fi

  # Variant calling
  if [ ! -f results/${sample}.vcf.gz ]; then
    lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    bgzip -f results/${sample}.vcf
    tabix -p vcf results/${sample}.vcf.gz
    rm -f results/${sample}.vcf
  fi
done

# Collapse VCFs into TSV
if [ ! -f results/collapsed.tsv ] || \
   ([ -n "$(find results/ -name "*.vcf.gz.tbi" -newer results/collapsed.tsv)" ] && [ -f results/collapsed.tsv ]); then
  cat $(printf "results/%s.vcf.gz\n" "${SAMPLES[@]}") | \
    bcftools query -f '%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
    sed "s/SAMPLE/${sample}/" | \
    sed -e '1i sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
fi