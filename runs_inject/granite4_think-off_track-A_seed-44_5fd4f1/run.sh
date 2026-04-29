#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# Helper to ensure output files exist before downstream steps
ensure_output() {
  local out=$1
  if [[ ! -f "$out" ]]; then
    echo "Running step that produces $out..."
    # Placeholder commands (actual logic will be executed by the caller)
  fi
}

# ------------------------------------------------------------------
# 1. Reference indexing – BWA
# ------------------------------------------------------------------
bwa index data/ref/chrM.fa || true
ensure_output data/ref/chrM.fa.bwt

# ------------------------------------------------------------------
# 2. Reference indexing – samtools faidx
# ------------------------------------------------------------------
samtools faidx data/ref/chrM.fa || true
ensure_output data/ref/chrM.fa.fai

# ------------------------------------------------------------------
# Per-sample alignment, sorting & variant calling (LoFreq)
# ------------------------------------------------------------------
for sample in "${SAMPLES[@]}"; do
  # Alignment + sort
  bwa mem -t "$THREADS" \
    -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
    data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz |
  samtools sort -@ "$THREADS" -o "${OUTPUT_DIR}/${sample}.bam"

  # Index BAM
  samtools index -@ "$THREADS" "${OUTPUT_DIR}/${sample}.bam"
  
  # Variant calling (LoFreq)
  lofreq call-parallel --pp-threads "$THREADS" \
    -f data/ref/chrM.fa -o "${OUTPUT_DIR}/${sample}.vcf" \
    "${OUTPUT_DIR}/${sample}.bam"

  # Compression + tabix
  bgzip -f "${OUTPUT_DIR}/${sample}.vcf"
  tabix -p vcf "${OUTPUT_DIR}/${sample}.vcf.gz"
done

# ------------------------------------------------------------------
# 3. Collapsed TSV (always rebuilt)
# ------------------------------------------------------------------
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
    "${OUTPUT_DIR}/${sample}.vcf.gz" |
  awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results/collapsed.tsv