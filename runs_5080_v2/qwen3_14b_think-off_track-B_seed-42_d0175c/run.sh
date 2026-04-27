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
  if [ -f "$RESULTS_DIR/${sample}.vcf.gz" ]; then
    continue
  fi

  # Align and sort
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
    bgzip -c "$RESULTS_DIR/${sample}.vcf" > "$RESULTS_DIR/${sample}.vcf.gz"
    tabix -p vcf "$RESULTS_DIR/${sample}.vcf.gz"
  fi
done

# Collapse variants across samples
if [ ! -f "$RESULTS_DIR/collapsed.tsv" ]; then
  # Use bcftools to merge VCFs
  bcftools merge -o "$RESULTS_DIR/merged.vcf" -O v "${RESULTS_DIR}/*.vcf.gz"
  bgzip -c "$RESULTS_DIR/merged.vcf" > "$RESULTS_DIR/merged.vcf.gz"
  tabix -p vcf "$RESULTS_DIR/merged.vcf.gz"

  # Use SnpSift to collapse variants
  java -jar /usr/local/bin/SnpSift.jar collapse -v "$RESULTS_DIR/merged.vcf.gz" > "$RESULTS_DIR/collapsed.vcf"
  # Extract AF and other info
  grep -v "^#" "$RESULTS_DIR/collapsed.vcf" | \
  awk -F'\t' '{
    chrom = $1; pos = $2; ref = $4; alt = $5;
    split($7, info, ";");
    for (i in info) {
      if (info[i] ~ /^AF=/) {
        split(info[i], af, "=");
        af_val = af[2];
        break;
      }
    }
    print $0, af_val
  }' | \
  cut -f1,2,4,5,10,11 | \
  awk -F'\t' '{
    if (NF == 6) {
      print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6
    }
  }' > "$RESULTS_DIR/collapsed.tsv"
fi