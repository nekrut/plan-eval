#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"

# Create results directory
mkdir -p results/

# === Reference indexing ===
if [[ ! -f "$REF.fai" ]]; then
  samtools faidx "$REF"
fi

if [[ ! -f "$REF.bwt" ]]; then
  bwa index "$REF"
fi

# === Per-sample processing ===
for sample in "${SAMPLES[@]}"; do
  FASTQ1="data/raw/${sample}_1.fq.gz"
  FASTQ2="data/raw/${sample}_2.fq.gz"
  BAM="results/${sample}.bam"
  BAI="results/${sample}.bam.bai"
  VCF_GZ="results/${sample}.vcf.gz"
  VCF_GZ_TBI="results/${sample}.vcf.gz.tbi"
  
  # Alignment and BAM creation
  if [[ ! -f "$BAM" ]]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      "$REF" "$FASTQ1" "$FASTQ2" | samtools sort -@ "$THREADS" -o "$BAM"
  fi
  
  # BAM indexing
  if [[ ! -f "$BAI" ]]; then
    samtools index -@ "$THREADS" "$BAM"
  fi
  
  # Variant calling
  if [[ ! -f "$VCF_GZ" ]]; then
    VCF_TEMP="results/${sample}.vcf"
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF_TEMP" "$BAM"
    bgzip "$VCF_TEMP"
  fi
  
  # VCF indexing
  if [[ ! -f "$VCF_GZ_TBI" ]]; then
    tabix -p vcf "$VCF_GZ"
  fi
done

# === Collapse to TSV ===
COLLAPSED="results/collapsed.tsv"

# Check if we need to rebuild collapsed.tsv
should_rebuild=0
if [[ ! -f "$COLLAPSED" ]]; then
  should_rebuild=1
else
  # Check if any VCF is newer than collapsed.tsv
  for sample in "${SAMPLES[@]}"; do
    if [[ "results/${sample}.vcf.gz" -nt "$COLLAPSED" ]]; then
      should_rebuild=1
      break
    fi
  done
fi

if [[ $should_rebuild -eq 1 ]]; then
  {
    printf "sample\tchrom\tpos\tref\talt\taf\n"
    for sample in "${SAMPLES[@]}"; do
      bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "results/${sample}.vcf.gz"
    done
  } > "$COLLAPSED"
fi