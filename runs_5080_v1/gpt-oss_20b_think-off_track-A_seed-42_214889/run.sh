#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

mkdir -p results

REF=data/ref/chrM.fa
if [ ! -f "${REF}.fai" ]; then
  samtools faidx "$REF"
fi
if [ ! -f "${REF}.amb" ]; then
  bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
  BAM=results/${sample}.bam
  BAI=results/${sample}.bam.bai
  VCF=results/${sample}.vcf
  VCFZ=results/${sample}.vcf.gz
  VCFZTBI=results/${sample}.vcf.gz.tbi

  if [ ! -f "$BAI" ]; then
    RG="@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA"
    bwa mem -t "$THREADS" -R "$RG" "$REF" data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | \
      samtools sort -@ "$THREADS" -o "$BAM"
    samtools index -@ "$THREADS" "$BAM"
  fi

  if [ ! -f "$VCFZTBI" ]; then
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -b "$BAM" -o "$VCF"
    bcftools view -O z -o "$VCFZ" "$VCF"
    tabix -p vcf "$VCFZ"
    rm -f "$VCF"
  fi
done

COLLAPSED=results/collapsed.tsv
REBUILD=false
if [ ! -f "$COLLAPSED" ]; then
  REBUILD=true
else
  for sample in "${SAMPLES[@]}"; do
    if [ "results/${sample}.vcf.gz" -nt "$COLLAPSED" ]; then
      REBUILD=true
      break
    fi
  done
fi

if $REBUILD; then
  printf "sample\tchrom\tpos\tref\talt\taf\n" > "$COLLAPSED"
  for sample in "${SAMPLES[@]}"; do
    bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" results/${sample}.vcf.gz >> "$COLLAPSED"
  done
fi

exit 0