#!/usr/bin/env bash
set -euo pipefail
THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
  samtools faidx -n data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
  bwa index data/ref/chrM.fa
fi

for sample in "${SAMPLES[@]}"; do
  # Alignment
  if [ ! -f results/${sample}.bam ]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "$THREADS" -o results/${sample}.bam
  fi

  # BAM indexing
  if [ ! -f results/${sample}.bam.bai ]; then
    samtools index -@ "$THREADS" results/${sample}.bam
  fi

  # Variant calling
  if [ ! -f results/${sample}.vcf.gz ]; then
    if [ ! -f results/${sample}.vcf ]; then
      lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    fi
    bgzip -f results/${sample}.vcf
    tabix -p vcf results/${sample}.vcf.gz
    rm -f results/${sample}.vcf
  fi
done

# Collapsed TSV
if [ ! -f results/collapsed.tsv ]; then
  printf "sample\tchrom\tpos\tref\talt\taf\n" > results/collapsed.tsv
fi

for sample in "${SAMPLES[@]}"; do
  vcf_file=results/${sample}.vcf.gz
  tsv_file=results/${sample}.tsv
  if [ ! -f "$tsv_file" ] || [ "$vcf_file" -nt "$tsv_file" ]; then
    bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcf_file" > "$tsv_file"
  fi
  if [ -f results/collapsed.tsv ] && [ "$tsv_file" -nt results/collapsed.tsv ]; then
    tail -n +2 "$tsv_file" >> results/collapsed.tsv
  fi
done