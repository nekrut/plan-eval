#!/usr/bin/env bash
set -euo pipefail

THREADS=4
RESULTS_DIR="results"
DATA_RAW="data/raw"
DATA_REF="data/ref"

mkdir -p "$RESULTS_DIR"

# Check if reference is indexed, if not, index it
if [ ! -f "$DATA_REF/chrM.fa.fai" ]; then
  samtools faidx "$DATA_REF/chrM.fa"
fi

# Create a list of samples
samples=()
for file in "$DATA_RAW"/*.fq.gz; do
  sample="${file##*/}"
  sample="${sample%_1.fq.gz}"
  if [[ ! " ${samples[@]} " =~ " $sample " ]]; then
    samples+=("$sample")
  fi
done

# Process each sample
for sample in "${samples[@]}"; do
  # Check if output files already exist
  if [ -f "$RESULTS_DIR/${sample}.bam" ] && [ -f "$RESULTS_DIR/${sample}.vcf.gz" ]; then
    continue
  fi

  # Align reads with bwa
  if [ ! -f "$RESULTS_DIR/${sample}.bam" ]; then
    bwa mem -t "$THREADS" "$DATA_REF/chrM.fa" "$DATA_RAW/${sample}_1.fq.gz" "$DATA_RAW/${sample}_2.fq.gz" | \
      samtools view -Sb - | \
      samtools sort -@ "$THREADS" -o "$RESULTS_DIR/${sample}.bam"
    samtools index "$RESULTS_DIR/${sample}.bam"
  fi

  # Variant calling with lofreq
  if [ ! -f "$RESULTS_DIR/${sample}.vcf.gz" ]; then
    lofreq mpileup -f "$DATA_REF/chrM.fa" "$RESULTS_DIR/${sample}.bam" | \
      lofreq call -f "$DATA_REF/chrM.fa" -o "$RESULTS_DIR/${sample}.vcf" -
    bgzip -f "$RESULTS_DIR/${sample}.vcf"
    tabix -p vcf "$RESULTS_DIR/${sample}.vcf.gz"
  fi
done

# Collapse VCFs into a single TSV
if [ ! -f "$RESULTS_DIR/collapsed.tsv" ]; then
  # Use bcftools to merge VCFs
  bcftools merge -o "$RESULTS_DIR/merged.vcf" -O v "${RESULTS_DIR}/*.vcf.gz"
  bgzip -f "$RESULTS_DIR/merged.vcf"
  tabix -p vcf "$RESULTS_DIR/merged.vcf.gz"

  # Use SnpSift to extract AF and format as TSV
  java -jar /usr/local/bin/SnpSift.jar extractFields "$RESULTS_DIR/merged.vcf.gz" \
    "CHROM POS ID REF ALT QUAL FILTER INFO" \
    -o "$RESULTS_DIR/collapsed.tsv" \
    -s "sample" \
    -f "AF" \
    -t
fi