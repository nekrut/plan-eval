#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF=data/ref/chrM.fa
RAW=data/raw/
RESULTS=results

mkdir -p "$RESULTS"

# Index reference if not already indexed
if [[ ! -f "$REF.bwt" ]]; then
  bwa index "$REF"
fi

# Process each sample
for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
  bam="$RESULTS/${sample}.bam"
  vcf="$RESULTS/${sample}.vcf.gz"
  
  # Skip if already processed
  if [[ -f "$bam" && -f "$vcf" ]]; then
    continue
  fi
  
  # Align reads
  bwa mem -t "$THREADS" "$REF" \
    "$RAW/${sample}_1.fq.gz" \
    "$RAW/${sample}_2.fq.gz" | \
  samtools view -b -@ "$THREADS" - | \
  samtools sort -@ "$THREADS" -o "$bam" -
  
  # Index BAM
  samtools index "$bam"
  
  # Call variants
  bcftools mpileup -Ou -f "$REF" "$bam" | \
  bcftools call -c -v -Oz -o "$vcf"
  
  # Index VCF
  tabix "$vcf"
done

# Collate variants into collapsed table
echo -e "sample\tchrom\tpos\tref\talt\taf" > "$RESULTS/collapsed.tsv"

for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
  vcf="$RESULTS/${sample}.vcf.gz"
  bcftools view "$vcf" | \
  awk -v sample="$sample" '
    BEGIN { OFS = "\t" }
    /^#/ { next }
    {
      chrom = $1
      pos = $2
      ref = $4
      alt = $5
      af = 0
      if ($8 ~ /AF=([^;]+)/) {
        af = $8
        gsub(/.*AF=/, "", af)
        gsub(/;.*/, "", af)
      }
      print sample, chrom, pos, ref, alt, af
    }
  ' >> "$RESULTS/collapsed.tsv"
done