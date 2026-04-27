#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p results
cd "$(dirname "$0")"

# Reference indexing
if [[ ! -f data/ref/chrM.fa.amb ]]; then
  samtools faidx data/ref/chrM.fa
  bwa index data/ref/chrM.fa
fi

# Per-sample workflow
for sample in "${SAMPLES[@]}"; do
  if [[ -f results/${sample}.vcf.gz.tbi ]] && [[ results/${sample}.vcf.gz.tbi -nt data/raw/${sample}_1.fq.gz ]]; then
    continue
  fi

  bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz \
    | samtools sort -@ "$THREADS" -o results/${sample}.bam

  samtools index -@ "$THREADS" results/${sample}.bam

  lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam

  bgzip results/${sample}.vcf
  tabix -p vcf results/${sample}.vcf.gz
  rm results/${sample}.vcf
done

# Collapse step
if [[ ! -f results/collapsed.tsv ]] || [[ results/collapsed.tsv -lt results/M117-bl.vcf.gz ]]; then
  > results/collapsed.tsv
  echo -e "sample\tchrom\tpos\tref\talt\taf" > results/collapsed.tsv
  for sample in "${SAMPLES[@]}"; do
    bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz >> results/collapsed.tsv
  done
fi