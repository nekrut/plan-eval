#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_DIR="data/ref"
RAW_DIR="data/raw"

mkdir -p "$RESULTS_DIR"

# Reference indexing
if [ ! -f "$REF_DIR/chrM.fa.fai" ] || [ ! -f "$REF_DIR/chrM.fa.amb" ]; then
  samtools faidx "$REF_DIR/chrM.fa"
  bwa index "$REF_DIR/chrM.fa"
fi

for sample in "${SAMPLES[@]}"; do
  BAM="$RESULTS_DIR/${sample}.bam"
  BAI="$RESULTS_DIR/${sample}.bam.bai"
  VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
  VCF_TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"
  VCF_TMP="$RESULTS_DIR/${sample}.vcf"

  # Skip if all outputs already exist and are up to date
  if [ -f "$BAI" ] && [ -f "$VCF_GZ" ] && [ -f "$VCF_TBI" ]; then
    # Check if any input is newer than outputs
    if [ "$RAW_DIR/${sample}_1.fq.gz" -nt "$BAI" ] || [ "$RAW_DIR/${sample}_2.fq.gz" -nt "$BAI" ] || [ "$REF_DIR/chrM.fa" -nt "$BAI" ]; then
      :
    else
      continue
    fi
  fi

  # Alignment
  if [ ! -f "$BAM" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" "$REF_DIR/chrM.fa" "$RAW_DIR/${sample}_1.fq.gz" "$RAW_DIR/${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # Index BAM
  if [ ! -f "$BAI" ]; then
    samtools index -@ "$THREADS" "$BAM" "$BAI"
  fi

  # Variant calling
  if [ ! -f "$VCF_GZ" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF_DIR/chrM.fa" -o "$VCF_TMP" "$BAM"
    bgzip -f "$VCF_TMP"
    mv "$VCF_TMP.gz" "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"
    rm "$VCF_TMP"
    tabix -p vcf "$VCF_GZ"
  fi
done

# Collapse VCFs into TSV
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
HEADER="sample\tchrom\tpos\tref\talt\taf"
TMP_COLLAPSED="$RESULTS_DIR/collapsed.tmp"

if [ ! -f "$COLLAPSED" ] || [ "$(find "$RESULTS_DIR" -name "*.vcf.gz" -newer "$COLLAPSED")" ]; then
  > "$TMP_COLLAPSED"
  echo "$HEADER" > "$TMP_COLLAPSED"
  for sample in "${SAMPLES[@]}"; do
    bcftools query -f "$sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$RESULTS_DIR/${sample}.vcf.gz" >> "$TMP_COLLAPSED"
  done
  mv "$TMP_COLLAPSED" "$COLLAPSED"
fi