#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

# Reference indexing (once, in data/ref/)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
  samtools faidx data/ref/chrM.fa
fi

if [[ ! -f data/ref/chrM.fa.amb ]]; then
  bwa index data/ref/chrM.fa
fi

# Per-sample alignment, sorting, indexing, and variant calling
for sample in "${SAMPLES[@]}"; do
  # bwa mem piped to samtools sort
  if [[ ! -f "results/${sample}.bam" || \
        "results/${sample}.bam" -ot "data/raw/${sample}_1.fq.gz" || \
        "results/${sample}.bam" -ot "data/raw/${sample}_2.fq.gz" ]]; then
    bwa mem -t ${THREADS} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/ref/chrM.fa \
      data/raw/${sample}_1.fq.gz \
      data/raw/${sample}_2.fq.gz | \
      samtools sort -@ ${THREADS} -o results/${sample}.bam
  fi

  # BAM indexing
  if [[ ! -f "results/${sample}.bam.bai" || \
        "results/${sample}.bam.bai" -ot "results/${sample}.bam" ]]; then
    samtools index -@ ${THREADS} results/${sample}.bam
  fi

  # Variant calling with lofreq and VCF compression
  if [[ ! -f "results/${sample}.vcf.gz" || \
        "results/${sample}.vcf.gz" -ot "results/${sample}.bam" ]]; then
    lofreq call-parallel --pp-threads ${THREADS} -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    bgzip -f results/${sample}.vcf
  fi

  # VCF indexing
  if [[ ! -f "results/${sample}.vcf.gz.tbi" || \
        "results/${sample}.vcf.gz.tbi" -ot "results/${sample}.vcf.gz" ]]; then
    tabix -p vcf results/${sample}.vcf.gz
  fi
done

# Collapse step → results/collapsed.tsv
REBUILD_COLLAPSED=0
if [[ ! -f results/collapsed.tsv ]]; then
  REBUILD_COLLAPSED=1
else
  for sample in "${SAMPLES[@]}"; do
    if [[ "results/${sample}.vcf.gz" -nt "results/collapsed.tsv" ]]; then
      REBUILD_COLLAPSED=1
      break
    fi
  done
fi

if [[ ${REBUILD_COLLAPSED} -eq 1 ]]; then
  {
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for sample in "${SAMPLES[@]}"; do
      bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/${sample}.vcf.gz
    done
  } > results/collapsed.tsv
fi

exit 0