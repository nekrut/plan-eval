#!/usr/bin/env bash
set -euo pipefail

THREADS=4
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
results_dir="results"
ref="data/ref/chrM.fa"

mkdir -p "$results_dir"

# Reference indexing
if [ ! -f "${ref}.fai" ] || [ ! -f "${ref}.amb" ]; then
  samtools faidx "$ref"
  bwa index "$ref"
fi

for sample in "${samples[@]}"; do
  bam="$results_dir/${sample}.bam"
  bai="$results_dir/${sample}.bam.bai"
  vcf_gz="$results_dir/${sample}.vcf.gz"
  vcf_tbi="$results_dir/${sample}.vcf.gz.tbi"
  vcf_tmp="$results_dir/${sample}.vcf"

  # Alignment and sorting
  if [ ! -f "$bam" ] || [ "$bam" -nt "$bam" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
      "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | \
      samtools sort -@ "$THREADS" -o "$bam"
  fi

  # BAM indexing
  if [ ! -f "$bai" ] || [ "$bai" -nt "$bam" ]; then
    samtools index -@ "$THREADS" "$bam"
  fi

  # Variant calling
  if [ ! -f "$vcf_gz" ] || [ "$vcf_gz" -nt "$bam" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f "$ref" -o "$vcf_tmp" "$bam"
    bgzip -f "$vcf_tmp"
    tabix -p vcf "$vcf_gz"
    rm "$vcf_tmp"
  fi
done

# Collapsed TSV
collapsed="$results_dir/collapsed.tsv"
if [ ! -f "$collapsed" ] || \
  ([ "$(find "$results_dir" -name "*.vcf.gz" -newer "$collapsed" | wc -l)" -gt 0 ]); then
  > "$collapsed"
  echo -e "sample\tchrom\tpos\tref\talt\taf" >> "$collapsed"
  for sample in "${samples[@]}"; do
    bcftools query -f "$sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$results_dir/${sample}.vcf.gz" >> "$collapsed"
  done
fi