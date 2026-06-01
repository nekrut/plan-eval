#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# Helper: ensure reference files exist (index & FASTA index)
ref_dir="data/ref"
for ref_file in chrM.fa; do
  if [[ ! -f "${ref_dir}/${ref_file}" ]]; then
    echo "Error: ${ref_file} not found at ${ref_dir}/" >&2
    exit 1
  fi

  # BWA index (amb, ann, bwt, pac, sa)
  if [[ ! -f "${ref_dir}/${ref_file}.amb" ]] ||
      [[ ! -f "${ref_dir}/${ref_file}.ann" ]] ||
      [[ ! -f "${ref_dir}/${ref_file}.bwt" ]] ||
      [[ ! -f "${ref_dir}/${ref_file}.pac" ]] ||
      [[ ! -f "${ref_dir}/${ref_file}.sa" ]]; then
    bwa index "${ref_dir}/${ref_file}"
  fi

  # FASTA index (FAI)
  if [[ ! -f "${ref_dir}/${ref_file}.fai" ]]; then
    samtools faidx "${ref_dir}/${ref_file}"
  fi
done

# Function to run per-sample alignment, calling & processing
run_sample() {
  local sample=$1
  BAM="${OUTPUT_DIR}/${sample}.bam"
  VCF="${OUTPUT_DIR}/${sample}.vcf.gz"

  # Alignment + sorting (bwa mem)
  bwa mem -t "$THREADS" \
    -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
    "${ref_dir}/chrM.fa" "${data}/raw/${sample}_1.fq.gz" "${data}/raw/${sample}_2.fq.gz" |
  samtools sort -@ "$THREADS" -o "$BAM" -

  # Index the BAM
  samtools index -@ "$THREADS" "$BAM"

  # Variant calling with LoFreq (call-parallel)
  lofreq call-parallel --pp-threads "$THREADS" \
    -f "${ref_dir}/chrM.fa" -o "$VCF" "$BAM"

  # Compress & index VCF
  bgzip -f "$VCF"
  tabix -p vcf "$VCF.gz"
}

# Execute per-sample pipeline
for sample in "${SAMPLES[@]}"; do
  run_sample "$sample"
done

# Collapsed TSV (always rebuilt)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
    "$OUTPUT_DIR/${sample}.vcf.gz" |
  awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done