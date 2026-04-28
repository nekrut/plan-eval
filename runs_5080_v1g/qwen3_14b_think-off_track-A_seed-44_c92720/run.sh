#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_DIR="data/ref"
RAW_DIR="data/raw"

mkdir -p "$RESULTS_DIR"

# Reference indexing
if [[ ! -f "$REF_DIR/chrM.fa.fai" || ! -f "$REF_DIR/chrM.fa.amb" ]]; then
  samtools faidx "$REF_DIR/chrM.fa"
  bwa index "$REF_DIR/chrM.fa"
fi

for sample in "${SAMPLES[@]}"; do
  BAM="$RESULTS_DIR/${sample}.bam"
  BAI="$RESULTS_DIR/${sample}.bam.bai"
  VCF="$RESULTS_DIR/${sample}.vcf"
  VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
  VCF_TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"

  # Skip if all outputs already exist and are up to date
  if [[ -f "$BAI" && -f "$VCF_TBI" && "$BAI" -nt "$RAW_DIR/${sample}_1.fq.gz" && "$VCF_TBI" -nt "$RAW_DIR/${sample}_2.fq.gz" ]]; then
    continue
  fi

  # Alignment
  if [[ ! -f "$BAM" || "$BAM" -nt "$RAW_DIR/${sample}_1.fq.gz" || "$BAM" -nt "$RAW_DIR/${sample}_2.fq.gz" ]]; then
    bwa mem -t "$THREADS" -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" "$REF_DIR/chrM.fa" "$RAW_DIR/${sample}_1.fq.gz" "$RAW_DIR/${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # BAM indexing
  if [[ ! -f "$BAI" || "$BAI" -nt "$BAM" ]]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  # Variant calling
  if [[ ! -f "$VCF_TBI" || "$VCF_TBI" -nt "$BAM" ]]; then
    lofreq call-parallel --pp-threads "$THREADS" --verbose --ref "$REF_DIR/chrM.fa" --out "$VCF" --sig --bonf "$BAM"
    bgzip -f "$VCF"
    tabix -p vcf "$VCF_GZ"
    rm -f "$VCF"
  fi
done

# Collapse VCFs into TSV
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
HEADER="sample\tchrom\tpos\tref\talt\taf"
VCF_GZ_FILES=("$RESULTS_DIR/"*.vcf.gz)

if [[ ! -f "$COLLAPSED" || "$COLLAPSED" -nt "${VCF_GZ_FILES[-1]}" ]]; then
  > "$COLLAPSED"
  echo "$HEADER" > "$COLLAPSED"
  for vcf_gz in "${VCF_GZ_FILES[@]}"; do
    sample=$(basename "$vcf_gz" .vcf.gz)
    bcftools query -f "$sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_gz" >> "$COLLAPSED"
  done
fi