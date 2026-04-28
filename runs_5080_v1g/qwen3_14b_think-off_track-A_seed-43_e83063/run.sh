#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS_DIR="results"
REF_DIR="data/ref"
RAW_DIR="data/raw"

mkdir -p "$RESULTS_DIR"

# Step 2: Reference indexing (once)
REF_FA="$REF_DIR/chrM.fa"
if [ ! -f "$REF_FA.fai" ] || [ ! -f "$REF_FA.amb" ]; then
  samtools faidx "$REF_FA"
  bwa index "$REF_FA"
fi

# Step 3-5: Per-sample alignment, sort, index
for sample in "${SAMPLES[@]}"; do
  BAM="$RESULTS_DIR/${sample}.bam"
  BAI="$RESULTS_DIR/${sample}.bam.bai"
  VCF="$RESULTS_DIR/${sample}.vcf"
  VCF_GZ="$RESULTS_DIR/${sample}.vcf.gz"
  VCF_TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"

  # Skip if all outputs already exist and are up to date
  if [ -f "$BAI" ] && [ -f "$VCF_TBI" ]; then
    if [ "$BAI" -nt "$RAW_DIR/${sample}_1.fq.gz" ] && [ "$BAI" -nt "$RAW_DIR/${sample}_2.fq.gz" ]; then
      continue
    fi
  fi

  # Alignment
  if [ ! -f "$BAM" ] || [ "$BAM" -nt "$RAW_DIR/${sample}_1.fq.gz" ] || [ "$BAM" -nt "$RAW_DIR/${sample}_2.fq.gz" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" "$REF_FA" "$RAW_DIR/${sample}_1.fq.gz" "$RAW_DIR/${sample}_2.fq.gz" | samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # Index BAM
  if [ ! -f "$BAI" ] || [ "$BAI" -nt "$BAM" ]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  # Variant calling
  if [ ! -f "$VCF_TBI" ] || [ "$VCF_TBI" -nt "$BAM" ]; then
    lofreq call-parallel --pp-threads "$THREADS" --verbose --ref "$REF_FA" --out "$VCF" --sig --bonf "$BAM"
    bgzip -f "$VCF"
    tabix -p vcf "$VCF_GZ"
    rm -f "$VCF"
  fi
done

# Step 8: Collapse variants into TSV
COLLAPSED_TSV="$RESULTS_DIR/collapsed.tsv"
VCF_GZ_FILES=("$RESULTS_DIR/"*.vcf.gz)

if [ ! -f "$COLLAPSED_TSV" ] || [ "$COLLAPSED_TSV" -nt "${VCF_GZ_FILES[-1]}" ]; then
  > "$COLLAPSED_TSV"
  echo -e "sample\tchrom\tpos\tref\talt\taf" >> "$COLLAPSED_TSV"
  for vcf_gz in "${VCF_GZ_FILES[@]}"; do
    sample=$(basename "$vcf_gz" .vcf.gz)
    bcftools query -f "$sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_gz" >> "$COLLAPSED_TSV"
  done
fi