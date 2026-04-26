#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
RAWDIR="data/raw"
OUTDIR="results"

mkdir -p "$OUTDIR"

if [[ ! -f "$REF.fai" ]]; then
  samtools faidx "$REF"
fi

if [[ ! -f "$REF.bwt" ]]; then
  bwa index "$REF"
fi

for sample in "${SAMPLES[@]}"; do
  R1="$RAWDIR/${sample}_1.fq.gz"
  R2="$RAWDIR/${sample}_2.fq.gz"
  BAM="$OUTDIR/${sample}.bam"
  BAI="${BAM}.bai"
  
  if [[ ! -f "$BAM" ]]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" \
      "$REF" "$R1" "$R2" | \
      samtools sort -@ "$THREADS" -o "$BAM"
  fi
  
  if [[ ! -f "$BAI" ]]; then
    samtools index -@ "$THREADS" "$BAM"
  fi
done

for sample in "${SAMPLES[@]}"; do
  BAM="$OUTDIR/${sample}.bam"
  VCF="$OUTDIR/${sample}.vcf"
  VCF_GZ="${VCF}.gz"
  VCF_GZ_TBI="${VCF_GZ}.tbi"
  
  if [[ ! -f "$VCF_GZ_TBI" ]] || [[ "$BAM" -nt "$VCF_GZ_TBI" ]]; then
    lofreq call-parallel --pp-threads "$THREADS" \
      -f "$REF" \
      -o "$VCF" \
      "$BAM"
    
    bgzip -f "$VCF"
    
    tabix -p vcf "$VCF_GZ"
  fi
done

COLLAPSED="$OUTDIR/collapsed.tsv"

needs_rebuild=false
if [[ ! -f "$COLLAPSED" ]]; then
  needs_rebuild=true
else
  for sample in "${SAMPLES[@]}"; do
    if [[ "$OUTDIR/${sample}.vcf.gz" -nt "$COLLAPSED" ]]; then
      needs_rebuild=true
      break
    fi
  done
fi

if [[ "$needs_rebuild" == true ]]; then
  {
    printf "sample\tchrom\tpos\tref\talt\taf\n"
    for sample in "${SAMPLES[@]}"; do
      bcftools query -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n" \
        "$OUTDIR/${sample}.vcf.gz"
    done
  } > "$COLLAPSED"
fi

exit 0