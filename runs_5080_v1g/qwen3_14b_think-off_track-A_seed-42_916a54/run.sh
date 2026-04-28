#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing
if [[ ! -f data/ref/chrM.fa.fai ]]; then
  samtools faidx data/ref/chrM.fa
fi

if [[ ! -f data/ref/chrM.fa.amb ]]; then
  bwa index data/ref/chrM.fa
fi

for sample in "${SAMPLES[@]}"; do
  # Skip if all outputs already exist and are up to date
  if [[ -f results/${sample}.bam.bai && -f results/${sample}.vcf.gz.tbi ]]; then
    continue
  fi

  # Alignment
  if [[ ! -f results/${sample}.bam ]]; then
    bwa mem -t $THREADS -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ $THREADS -o results/${sample}.bam
  fi

  # Index BAM
  if [[ ! -f results/${sample}.bam.bai ]]; then
    samtools index -@ $THREADS results/${sample}.bam
  fi

  # Variant calling
  if [[ ! -f results/${sample}.vcf.gz ]]; then
    lofreq call-parallel --pp-threads $THREADS --verbose --ref data/ref/chrM.fa \
      --out results/${sample}.vcf --sig --bonf results/${sample}.bam
    bgzip -f results/${sample}.vcf
    tabix -p vcf results/${sample}.vcf.gz
    rm -f results/${sample}.vcf
  fi
done

# Collapse VCFs into TSV
if [[ ! -f results/collapsed.tsv || $(find results/ -name "*.vcf.gz.tbi" -newer results/collapsed.tsv | wc -l) -gt 0 ]]; then
  cat $(printf "results/%s.vcf.gz\n" "${SAMPLES[@]}" | xargs -I{} bcftools query -f '{sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' {}) | \
    sed '1i sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
fi