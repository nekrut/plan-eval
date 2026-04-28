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
      samtools view -Sb - | \
      samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index "$BAM" "$BAI"
  fi

  # Variant calling with lofreq
  if [ ! -f "$VCF" ]; then
    lofreq call -t "$THREADS" -f "$DATA_REF/chrM.fa" "$BAM" | \
      bcftools view -O z -o "$VCF" -
    tabix -p vcf "$VCF"
  fi
done

# Collapse variants into a single TSV
COLLAPSED="$RESULTS_DIR/collapsed.tsv"
if [ ! -f "$COLLAPSED" ]; then
  # Extract variants from all VCFs
  for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' "$RESULTS_DIR/${sample}.vcf.gz" | \
      awk -v sample="$sample" '{print sample, $0}' >> "$COLLAPSED"
  done

  # Collapse by chrom, pos, ref, alt
  sort -k2,2 -k3,3 -k4,4 -k5,5 "$COLLAPSED" | \
    awk '{
      key = $2 FS $3 FS $4 FS $5
      if (!(key in seen)) {
        seen[key] = $0
      }
    }
    END {
      for (key in seen) {
        print seen[key]
      }
    }' | \
    sort -k2,2n -k3,3n > "$COLLAPSED.tmp" && mv "$COLLAPSED.tmp" "$COLLAPSED"
fi