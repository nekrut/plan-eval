#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# Helper to ensure output files exist
ensure_output() {
  local out=$1
  if [[ ! -f "$out" ]]; then
    echo "Running step that produces $out..."
    # ... (commands for the step) ...
  fi
}

# 1. Index reference with BWA
bwa index data/ref/chrM.fa || true
ensure_output data/ref/chrM.fa.bwt

# 2. Fasta index with samtools
samtools faidx data/ref/chrM.fa || true
ensure_output data/ref/chrM.fa.fai

# 3. Alignment & sorting per sample
for sample in "${SAMPLES[@]}"; do
  bwa mem -t $THREADS \
    -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
    data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz |
  samtools sort -@ $THREADS -o "$OUTPUT_DIR/${sample}.bam"
done

# 4. BAM indexing
for sample in "${SAMPLES[@]}"; do
  ensure_output "$OUTPUT_DIR/${sample}.bam.bai"
done

# 5. Variant calling with LoFreq (parallel)
for sample in "${SAMPLES[@]}"; do
  lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o "$OUTPUT_DIR/${sample}.vcf" "$OUTPUT_DIR/${sample}.bam"
done

# 6. Compress & index VCFs (remove original .vcf, create .vcf.gz + tbi)
for sample in "${SAMPLES[@]}"; do
  bgzip -f "$OUTPUT_DIR/${sample}.vcf" && \
  tabix -p vcf "$OUTPUT_DIR/${sample}.vcf.gz"
done

# 7. Collapsed TSV (always rebuild)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUTPUT_DIR/${sample}.vcf.gz" |
  awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done