#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

if [ ! -f "data/ref/chrM.fa.fai" ]; then
  samtools faidx data/ref/chrM.fa
fi

if [ ! -f "data/ref/chrM.fa.bwt" ]; then
  bwa index data/ref/chrM.fa
fi

for sample in "${SAMPLES[@]}"; do
  if [ ! -f "results/${sample}.bam" ]; then
    bwa mem -t "${THREADS}" \
      -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/ref/chrM.fa \
      data/raw/${sample}_1.fq.gz \
      data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "${THREADS}" -o "results/${sample}.bam"
  fi
  
  if [ ! -f "results/${sample}.bam.bai" ]; then
    samtools index -@ "${THREADS}" "results/${sample}.bam"
  fi
  
  if [ ! -f "results/${sample}.vcf.gz" ]; then
    if [ ! -f "results/${sample}.vcf" ]; then
      lofreq call-parallel --pp-threads "${THREADS}" --verbose \
        --ref data/ref/chrM.fa --out "results/${sample}.vcf" \
        --sig --bonf \
        "results/${sample}.bam"
    fi
    bgzip -f "results/${sample}.vcf"
  fi
  
  if [ ! -f "results/${sample}.vcf.gz.tbi" ]; then
    tabix -p vcf "results/${sample}.vcf.gz"
  fi
done

if [ ! -f "results/collapsed.tsv" ] || \
   [ "results/collapsed.tsv" -ot "results/M117-bl.vcf.gz" ] || \
   [ "results/collapsed.tsv" -ot "results/M117-ch.vcf.gz" ] || \
   [ "results/collapsed.tsv" -ot "results/M117C1-bl.vcf.gz" ] || \
   [ "results/collapsed.tsv" -ot "results/M117C1-ch.vcf.gz" ]; then
  {
    printf "sample\tchrom\tpos\tref\talt\taf\n"
    for sample in "${SAMPLES[@]}"; do
      bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
        "results/${sample}.vcf.gz"
    done
  } > "results/collapsed.tsv"
fi