#!/usr/bin/env bash
set -euo pipefail

THREADS=4
RESULTS_DIR="results"
DATA_RAW="data/raw"
DATA_REF="data/ref"

mkdir -p "$RESULTS_DIR"

# Check if reference is indexed, index if not
if [ ! -f "$DATA_REF/chrM.fa.fai" ]; then
  samtools faidx "$DATA_REF/chrM.fa"
fi

# Process each sample
for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
  # Check if output files already exist
  BAM="$RESULTS_DIR/${sample}.bam"
  BAI="$RESULTS_DIR/${sample}.bam.bai"
  VCF="$RESULTS_DIR/${sample}.vcf.gz"
  TBI="$RESULTS_DIR/${sample}.vcf.gz.tbi"

  if [ -f "$BAI" ] && [ -f "$TBI" ]; then
    continue
  fi

  # Align with bwa
  if [ ! -f "$BAM" ]; then
    bwa mem -t "$THREADS" "$DATA_REF/chrM.fa" "$DATA_RAW/${sample}_1.fq.gz" "$DATA_RAW/${sample}_2.fq.gz" | \
      samtools view -S -b - > "$BAM"
  fi

  # Sort and index BAM
  if [ ! -f "$BAI" ]; then
    samtools sort -o "$BAM" "$BAM" && \
      samtools index "$BAM" && \
      mv "${BAM}.bai" "$BAI"
  fi

  # Variant calling with lofreq
  if [ ! -f "$VCF" ]; then
    lofreq call -t "$THREADS" -f "$DATA_REF/chrM.fa" "$BAM" | \
      bcftools view -O z -o "$VCF" -
  fi

  # Index VCF
  if [ ! -f "$TBI" ]; then
    tabix -p vcf "$VCF"
  fi
done

# Collapse VCFs into a single TSV
if [ ! -f "$RESULTS_DIR/collapsed.tsv" ]; then
  cat "$RESULTS_DIR"/*.vcf.gz | \
    bcftools concat -o - | \
    bcftools query -f '%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%AF\n' | \
    awk 'BEGIN { OFS="\t"; print "sample", "chrom", "pos", "ref", "alt", "af" } { $1 = substr($1, 1, length($1)-1); print }' > "$RESULTS_DIR/collapsed.tsv"
fi