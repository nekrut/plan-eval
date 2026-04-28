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

# Define samples
samples=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Process each sample
for sample in "${samples[@]}"; do
  # Define output paths
  bam_file="$RESULTS_DIR/${sample}.bam"
  bai_file="$RESULTS_DIR/${sample}.bam.bai"
  vcf_file="$RESULTS_DIR/${sample}.vcf.gz"
  tbi_file="$RESULTS_DIR/${sample}.vcf.gz.tbi"

  # Skip if output files already exist
  if [ -f "$bam_file" ] && [ -f "$vcf_file" ]; then
    continue
  fi

  # Align with bwa
  if [ ! -f "$bam_file" ]; then
    bwa mem -t "$THREADS" "$DATA_REF/chrM.fa" "$DATA_RAW/${sample}_1.fq.gz" "$DATA_RAW/${sample}_2.fq.gz" | \
      samtools view -Sb - | \
      samtools sort -@ "$THREADS" -o "$bam_file" -T "$RESULTS_DIR/${sample}.tmp" -
  fi

  # Index BAM
  if [ ! -f "$bai_file" ]; then
    samtools index "$bam_file" "$bai_file"
  fi

  # Variant calling with lofreq
  if [ ! -f "$vcf_file" ]; then
    lofreq mpileup -f "$DATA_REF/chrM.fa" "$bam_file" | \
      lofreq call -f "$DATA_REF/chrM.fa" -o "$vcf_file" - -t "$THREADS"
  fi

  # Index VCF
  if [ ! -f "$tbi_file" ]; then
    tabix -p vcf "$vcf_file"
  fi
done

# Collapse VCFs into a single TSV
collapsed_tsv="$RESULTS_DIR/collapsed.tsv"

if [ ! -f "$collapsed_tsv" ]; then
  # Use bcftools to merge VCFs
  bcftools merge -o "$RESULTS_DIR/merged.vcf" -O v -m none "${samples[@]/%/.vcf.gz}" && \
    # Convert to TSV with required columns
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$RESULTS_DIR/merged.vcf" | \
    awk 'NR==1 {print "sample\tchrom\tpos\tref\talt\taf"} NR>1 {print $0}' > "$collapsed_tsv" && \
    # Clean up
    rm "$RESULTS_DIR/merged.vcf"
fi