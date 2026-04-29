#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# Helper: ensure reference files exist (index & FASTA index)
ensure_ref_files() {
  if [[ ! -f data/ref/chrM.fa.bwt ]]; then
    bwa index data/ref/chrM.fa || exit 1
  fi
  if [[ ! -f data/ref/chrM.fa.fai ]]; then
    samtools faidx data/ref/chrM.fa || exit 1
  fi
}

# Main per-sample processing loop
for sample in "${SAMPLES[@]}"; do
  ensure_ref_files

  # Alignment & sorting
  bwa mem -t "$THREADS" \
      -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
      data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz |
    samtools sort -@ "$THREADS" -o "${OUTPUT_DIR}/${sample}.bam"

  # Index the BAM
  samtools index -@ "$THREADS" "${OUTPUT_DIR}/${sample}.bam" || exit 1

  # Variant calling with LoFreq
  lofreq call-parallel --pp-threads "$THREADS" \
      -f data/ref/chrM.fa -o "${OUTPUT_DIR}/${sample}.vcf" \
      "${OUTPUT_DIR}/${sample}.bam" || exit 1

  # Compress & index VCF
  bgzip -f "${OUTPUT_DIR}/${sample}.vcf"
  tabix -p vcf "${OUTPUT_DIR}/${sample}.vcf.gz" || exit 1
done

# Collapsed TSV (rebuild every run)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
      results/${sample}.vcf.gz |
    awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done