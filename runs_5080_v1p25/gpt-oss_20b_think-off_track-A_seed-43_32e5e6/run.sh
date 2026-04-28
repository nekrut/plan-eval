#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p results

REF=data/ref/chrM.fa

# Reference indexing
if [ ! -f "${REF}.fai" ]; then
  samtools faidx "$REF"
fi
if [ ! -f "${REF}.amb" ]; then
  bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
  fq1=data/raw/${sample}_1.fq.gz
  fq2=data/raw/${sample}_2.fq.gz
  bam=results/${sample}.bam
  bai=results/${sample}.bam.bai
  vcf=results/${sample}.vcf
  vcfgz=results/${sample}.vcf.gz
  vcfgz_tbi=results/${sample}.vcf.gz.tbi

  # Alignment
  if [ ! -f "$bai" ] || [ ! -f "$bam" ] || [ "$fq1" -nt "$bam" ] || [ "$fq2" -nt "$bam" ]; then
    bwa mem -t "$THREADS" -R "@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA" "$REF" "$fq1" "$fq2" | \
      samtools sort -@ "$THREADS" -o "$bam"
    samtools index -@ "$THREADS" "$bam"
  fi

  # Variant calling
  if [ ! -f "$vcfgz_tbi" ] || [ ! -f "$vcfgz" ] || [ "$bam" -nt "$vcfgz" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$vcf" "$bam"
    bcftools view -O z -o "$vcfgz" "$vcf"
    tabix -p vcf "$vcfgz"
    rm -f "$vcf"
  fi
done

# Collapse VCFs
collapsed=results/collapsed.tsv
rebuild=false
for sample in "${SAMPLES[@]}"; do
  vcfgz=results/${sample}.vcf.gz
  if [ ! -f "$collapsed" ] || [ "$vcfgz" -nt "$collapsed" ]; then
    rebuild=true
    break
  fi
done

if $rebuild; then
  rm -f "$collapsed"
  printf "sample\tchrom\tpos\tref\talt\taf\n" > "$collapsed"
  for sample in "${SAMPLES[@]}"; do
    vcfgz=results/${sample}.vcf.gz
    bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" "$vcfgz" >> "$collapsed"
  done
fi

exit 0