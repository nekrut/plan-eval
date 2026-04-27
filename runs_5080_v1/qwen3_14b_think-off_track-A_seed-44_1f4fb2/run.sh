#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
  samtools faidx data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
  bwa index data/ref/chrM.fa
fi

for sample in "${SAMPLES[@]}"; do
  BAM="results/${sample}.bam"
  BAI="results/${sample}.bam.bai"
  VCF_GZ="results/${sample}.vcf.gz"
  VCF_TBI="results/${sample}.vcf.gz.tbi"

  # Skip if all outputs already exist and are up to date
  if [ -f "$BAM" ] && [ -f "$BAI" ] && [ -f "$VCF_GZ" ] && [ -f "$VCF_TBI" ]; then
    continue
  fi

  # Alignment
  if [ ! -f "$BAM" ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "$THREADS" -o "$BAM"
  fi

  # Index BAM
  if [ ! -f "$BAI" ]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  # Variant calling
  if [ ! -f "$VCF_GZ" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${sample}.vcf "$BAM"
    bgzip -c results/${sample}.vcf > "$VCF_GZ"
    tabix -p vcf "$VCF_GZ"
    rm -f results/${sample}.vcf
  fi
done

# Collapse VCFs into TSV
if [ ! -f results/collapsed.tsv ]; then
  cat $(printf "results/%s.vcf.gz\n" "${SAMPLES[@]}") | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
    sed "s/\t/\t${sample}\t/" | \
    sed '1i sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
else
  # Check if any VCF is newer than collapsed.tsv
  updated=false
  for sample in "${SAMPLES[@]}"; do
    VCF_GZ="results/${sample}.vcf.gz"
    if [ "$(date -r "$VCF_GZ" +%s)" -gt "$(date -r results/collapsed.tsv +%s)" ]; then
      updated=true
      break
    fi
  done

  if $updated; then
    cat $(printf "results/%s.vcf.gz\n" "${SAMPLES[@]}") | \
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
      sed "s/\t/\t${sample}\t/" | \
      sed '1i sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
  fi
fi