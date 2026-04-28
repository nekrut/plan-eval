#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
results="results"
ref="data/ref/chrM.fa"

mkdir -p "$results"

# Reference indexing
if [ ! -f "$ref.fai" ] || [ ! -f "${ref}.amb" ]; then
  samtools faidx "$ref"
  bwa index "$ref"
fi

for sample in "${samples[@]}"; do
  bam="$results/${sample}.bam"
  bai="$results/${sample}.bam.bai"
  vcf="$results/${sample}.vcf.gz"
  vcf_tbi="$results/${sample}.vcf.gz.tbi"
  fastq1="data/raw/${sample}_1.fq.gz"
  fastq2="data/raw/${sample}_2.fq.gz"

  # Alignment and sorting
  if [ ! -f "$bam" ] || [ "$bam" -nt "$fastq1" ] || [ "$bam" -nt "$fastq2" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" "$ref" "$fastq1" "$fastq2" | \
      samtools sort -@ "$THREADS" -o "$bam"
  fi

  # BAM indexing
  if [ ! -f "$bai" ] || [ "$bai" -nt "$bam" ]; then
    samtools index -@ "$THREADS" "$bam"
  fi

  # Variant calling
  if [ ! -f "$vcf" ] || [ "$vcf" -nt "$bam" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f "$ref" -o "$results/${sample}.vcf" "$bam"
    bgzip -f "$results/${sample}.vcf"
    tabix -p vcf "$vcf"
    rm -f "$results/${sample}.vcf"
  fi
done

# Collapsed TSV
collapsed="$results/collapsed.tsv"
if [ ! -f "$collapsed" ] || \
   ([ "$(find "$results" -name "*.vcf.gz" -newer "$collapsed" | wc -l)" -gt 0 ] && [ "$collapsed" -nt "$results"/*.vcf.gz 2>/dev/null | grep -q '0$' ]); then
  > "$collapsed"
  echo -e "sample\tchrom\tpos\tref\talt\taf" >> "$collapsed"
  for sample in "${samples[@]}"; do
    bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$results/${sample}.vcf.gz" >> "$collapsed"
  done
fi