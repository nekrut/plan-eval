#!/usr/bin/env bash
set -euo pipefail

# Globals
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Prepare results directory
mkdir -p results

# Step 2: Reference indexing (once)
if [[ ! -f data/ref/chrM.fa.fai ]]; then
  samtools faidx data/ref/chrM.fa
fi

if [[ ! -f data/ref/chrM.fa.amb ]]; then
  bwa index data/ref/chrM.fa
fi

# Steps 3-7: Per-sample processing
for sample in "${SAMPLES[@]}"; do
  # Step 3-4: Alignment and sort to BAM
  if [[ ! -f results/"${sample}".bam ]]; then
    bwa mem -t ${THREADS} -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/ref/chrM.fa data/raw/"${sample}"_1.fq.gz data/raw/"${sample}"_2.fq.gz | \
      samtools sort -@ ${THREADS} -o results/"${sample}".bam
  fi

  # Step 5: BAM indexing
  if [[ ! -f results/"${sample}".bam.bai ]]; then
    samtools index -@ ${THREADS} results/"${sample}".bam
  fi

  # Step 6-7: Variant calling, compression, and indexing
  if [[ ! -f results/"${sample}".vcf.gz ]]; then
    rm -f results/"${sample}".vcf
    lofreq call-parallel --pp-threads ${THREADS} -f data/ref/chrM.fa -o results/"${sample}".vcf results/"${sample}".bam
    bgzip -c results/"${sample}".vcf > results/"${sample}".vcf.gz
    rm results/"${sample}".vcf
  fi

  # Index VCF
  if [[ ! -f results/"${sample}".vcf.gz.tbi ]]; then
    tabix -p vcf results/"${sample}".vcf.gz
  fi
done

# Step 8: Collapse to TSV
should_rebuild=false
if [[ ! -f results/collapsed.tsv ]]; then
  should_rebuild=true
else
  for sample in "${SAMPLES[@]}"; do
    if [[ results/"${sample}".vcf.gz -nt results/collapsed.tsv ]]; then
      should_rebuild=true
      break
    fi
  done
fi

if [[ "$should_rebuild" == "true" ]]; then
  {
    printf "sample\tchrom\tpos\tref\talt\taf\n"
    for sample in "${SAMPLES[@]}"; do
      bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/"${sample}".vcf.gz
    done
  } > results/collapsed.tsv
fi

exit 0