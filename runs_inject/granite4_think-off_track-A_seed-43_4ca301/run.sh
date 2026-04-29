#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# Helper: ensure reference files exist before proceeding
bwa index data/ref/chrM.fa && samtools faidx data/ref/chrM.fa

for sample in "${SAMPLES[@]}"; do
  # Alignment & sorting (per-sample)
  bwa mem -t $THREADS \
    -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
    data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz |
    samtools sort -@ $THREADS -o "$OUTPUT_DIR/$sample.bam"

  # Index the BAM
  samtools index -@ $THREADS "$OUTPUT_DIR/$sample.bam"
  
  # Variant calling with LoFreq (parallel)
  lofreq call-parallel --pp-threads $THREADS \
    -f data/ref/chrM.fa -o "$OUTPUT_DIR/${sample}.vcf" "$OUTPUT_DIR/$sample.bam"

  # Compress & index VCF
  bgzip -f "$OUTPUT_DIR/${sample}.vcf"
  tabix -p vcf "$OUTPUT_DIR/${sample}.vcf.gz"
done

# Collapsed TSV (always rebuilt)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUTPUT_DIR/collapsed.tsv"

for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUTPUT_DIR/${sample}.vcf.gz" |
    awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUTPUT_DIR/collapsed.tsv"
done