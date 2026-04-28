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
  VCF="results/${sample}.vcf"
  VCF_GZ="results/${sample}.vcf.gz"
  VCF_TBI="results/${sample}.vcf.gz.tbi"

  # Skip if all outputs already exist and are up to date
  if [ -f "$BAI" ] && [ -f "$VCF_TBI" ]; then
    if [ "$BAI" -nt "$VCF_TBI" ]; then
      continue
    fi
  fi

  # Alignment
  if [ ! -f "$BAM" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # BAM indexing
  if [ ! -f "$BAI" ]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  # Variant calling
  if [ ! -f "$VCF_GZ" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "$VCF" "$BAM"
    bgzip -f "$VCF"
    tabix -p vcf "$VCF_GZ"
    rm -f "$VCF"
  fi
done

# Collapse VCFs into TSV
COLLAPSED="results/collapsed.tsv"
if [ ! -f "$COLLAPSED" ] || \
   ([ -f "$COLLAPSED" ] && ([ "$(find results -name "*.vcf.gz" -newer "$COLLAPSED")" != "" ] || [ "$COLLAPSED" -lt "$(find results -name "*.vcf.gz" -printf "%T+\n" | sort -r | head -n1)" ])); then
  > "$COLLAPSED"
  echo -e "sample\tchrom\tpos\tref\talt\taf" >> "$COLLAPSED"
  for sample in "${SAMPLES[@]}"; do
    bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/${sample}.vcf.gz >> "$COLLAPSED"
  done
fi