#!/usr/bin/env bash
set -euo pipefail

mkdir -p results
chmod -R a+rwx results

# Index reference
if [ ! -f data/ref/chrM.fa.fai ]; then
  samtools faidx data/ref/chrM.fa
fi

# Process each sample
for sample in M117-bl M117-ch M117C1-bl M117C1-ch; do
  if [ -f results/${sample}.bam ]; then
    continue
  fi

  # Align with bwa
  bwa mem -t 4 data/ref/chrM.fa \
    data/raw/${sample}_1.fq.gz \
    data/raw/${sample}_2.fq.gz | \
  samtools view -Sb - > results/${sample}.bam

  # Index BAM
  samtools index results/${sample}.bam

  # Call variants with lofreq
  lofreq call -t 4 -f data/ref/chrM.fa results/${sample}.bam | \
  bcftools view -b - > results/${sample}.vcf.gz

  # Index VCF
  bcftools index results/${sample}.vcf.gz

  # Collapse variants and extract AF
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' results/${sample}.vcf.gz | \
  awk -v sample="$sample" '{print sample, $1, $2, $3, $4}' >> results/collapsed.tsv
done

# Finalize collapsed.tsv
sort -k1,1 -k2,2n results/collapsed.tsv | \
  uniq -f 0 | \
  cut -f1-5 | \
  sort -k1,1 -k2,2n > results/collapsed.tsv

# Clean up intermediate files
find results -type f -name "*.bam" -exec rm -f {} \;
find results -type f -name "*.bam.bai" -exec rm -f {} \;
find results -type f -name "*.vcf.gz" -exec rm -f {} \;
find results -type f -name "*.vcf.gz.tbi" -exec rm -f {} \;