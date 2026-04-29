#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
OUTPUT_DIR="results"
mkdir -p "$OUTPUT_DIR"

# Helper: ensure reference files exist
ref_path="data/ref/chrM.fa"
[[ -f $ref_path.fai ]] || samtools faidx "$ref_path"
[[ -f $ref_path.amb ]] && : || bwa index "$ref_path"

for sample in "${SAMPLES[@]}"; do
  BAM="${OUTPUT_DIR}/${sample}.bam"
  VCF_GZ="${OUTPUT_DIR}/${sample}.vcf.gz"
  
  # Alignment & sorting (bwa mem + samtools sort)
  bwa mem -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
    "$ref_path" "${data}/raw/${sample}_1.fq.gz" "${data}/raw/${sample}_2.fq.gz" |
  samtools sort -@ $THREADS -o "$BAM"

  # Index BAM
  samtools index -@ $THREADS "$BAM"
  
  # Variant calling with LoFreq
  LO_FREQ_CALL="${OUTPUT_DIR}/${sample}.vcf"
  lofreq call-parallel --pp-threads $THREADS -f "$ref_path" -o "$LO_FREQ_CALL" "$BAM"

  # Compress & index VCF
  bgzip -f "$LO_FREQ_CALL"
  tabix -p vcf "${OUTPUT_DIR}/${sample}.vcf.gz"
  
  # Clean up intermediate unzipped VCF
  rm "$LO_FREQ_CALL"
done

# Collapsed TSV (always rebuilt)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for sample in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUTPUT_DIR/${sample}.vcf.gz" |
  awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done