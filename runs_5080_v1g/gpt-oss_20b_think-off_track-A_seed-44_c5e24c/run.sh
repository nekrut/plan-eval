#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

# Create results directory if needed
mkdir -p results

# Reference indexing
if [ ! -f data/ref/chrM.fa.fai ]; then
  samtools faidx data/ref/chrM.fa
fi
if [ ! -f data/ref/chrM.fa.amb ]; then
  bwa index data/ref/chrM.fa
fi

# Per-sample processing
for sample in "${SAMPLES[@]}"; do
  bam=results/${sample}.bam
  bai=results/${sample}.bam.bai
  vcf=results/${sample}.vcf
  vcfz=results/${sample}.vcf.gz
  vcfztbi=results/${sample}.vcf.gz.tbi

  # Alignment and sorting
  if [ ! -f "$bam" ] || \
     [ data/raw/${sample}_1.fq.gz -nt "$bam" ] || \
     [ data/raw/${sample}_2.fq.gz -nt "$bam" ]; then
    bwa mem -t "$THREADS" -R "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA" \
      data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "$THREADS" -o "$bam"
  fi

  # BAM indexing
  if [ ! -f "$bai" ] || [ "$bam" -nt "$bai" ]; then
    samtools index -@ "$THREADS" "$bam"
  fi

  # Variant calling
  if [ ! -f "$vcf" ] || [ "$bam" -nt "$vcf" ]; then
    lofreq call-parallel --pp-threads "$THREADS" --verbose \
      --ref data/ref/chrM.fa --out "$vcf" \
      --sig --bonf "$bam"
  fi

  # VCF compression
  if [ ! -f "$vcfz" ] || [ "$vcf" -nt "$vcfz" ]; then
    samtools view -O z -o "$vcfz" "$vcf"
    rm -f "$vcf"
  fi

  # VCF indexing
  if [ ! -f "$vcfztbi" ] || [ "$vcfz" -nt "$vcfztbi" ]; then
    tabix -p vcf "$vcfz"
  fi
done

# Collapsed table
collapsed=results/collapsed.tsv
needs_rebuild=false
for sample in "${SAMPLES[@]}"; do
  vcfz=results/${sample}.vcf.gz
  if [ ! -f "$collapsed" ] || [ "$vcfz" -nt "$collapsed" ]; then
    needs_rebuild=true
    break
  fi
done

if $needs_rebuild; then
  {
    echo -e "sample\tchrom\tpos\tref\talt\taf"
    for sample in "${SAMPLES[@]}"; do
      bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/${sample}.vcf.gz
    done
  } > "$collapsed"
fi

exit 0