#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
RESULTS="results"

mkdir -p "$RESULTS"

# Reference indexing
REF="data/ref/chrM.fa"
if [ ! -f "${REF}.fai" ]; then
  samtools faidx "$REF"
fi
if [ ! -f "${REF}.amb" ]; then
  bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
  BAM="$RESULTS/${sample}.bam"
  BAI="$RESULTS/${sample}.bam.bai"
  VCF="$RESULTS/${sample}.vcf"
  VCFGZ="$RESULTS/${sample}.vcf.gz"
  VTBI="$RESULTS/${sample}.vcf.gz.tbi"

  # Alignment
  if [ ! -f "$BAM" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
      "data/raw/${sample}_1.fq.gz" "data/raw/${sample}_2.fq.gz" | \
      samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # BAM indexing
  if [ ! -f "$BAI" ] || [ "$BAI" -nt "$BAM" ]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  # Variant calling
  if [ ! -f "$VCFGZ" ] || [ "$VCFGZ" -nt "$BAM" ]; then
    lofreq call-parallel --pp-threads "$THREADS" --verbose --ref "$REF" --out "$VCF" --sig --bonf "$BAM"
    bgzip -f "$VCF"
    tabix -p vcf "$VCFGZ"
    rm -f "$VCF"
  fi
done

# Collapsed TSV
COLLAPSED="$RESULTS/collapsed.tsv"
TMPFILE=$(mktemp)
if [ ! -f "$COLLAPSED" ] || \
  ([ "$(find "$RESULTS" -name "*.vcf.gz" -newer "$COLLAPSED" | wc -l)" -gt 0 ]); then
  > "$TMPFILE"
  echo -e "sample\tchrom\tpos\tref\talt\taf" > "$COLLAPSED"
  for sample in "${SAMPLES[@]}"; do
    VCFGZ="$RESULTS/${sample}.vcf.gz"
    bcftools query -f "$sample\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$VCFGZ" >> "$TMPFILE"
  done
  cat "$TMPFILE" >> "$COLLAPSED"
  rm -f "$TMPFILE"
fi