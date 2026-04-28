#!/usr/bin/env bash
set -euo pipefail

THREADS=4
REF="data/ref/chrM.fa"
RAWDIR="data/raw"
OUTDIR="results"

# Create output directory
mkdir -p "$OUTDIR"

# Index reference if needed
if [[ ! -f "${REF}.bwt" ]]; then
  bwa index "$REF"
fi

if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

# Sample list
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

for SAMPLE in "${SAMPLES[@]}"; do
  R1="${RAWDIR}/${SAMPLE}_1.fq.gz"
  R2="${RAWDIR}/${SAMPLE}_2.fq.gz"
  
  BAM="${OUTDIR}/${SAMPLE}.bam"
  BAI="${BAM}.bai"
  VCF="${OUTDIR}/${SAMPLE}.vcf.gz"
  TBI="${VCF}.tbi"
  
  # Check if this sample is done
  if [[ -f "$BAM" && -f "$BAI" && -f "$VCF" && -f "$TBI" ]]; then
    continue
  fi
  
  # BWA mapping and convert to BAM
  if [[ ! -f "$BAM" ]]; then
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | \
    samtools view -b | \
    samtools sort -@ "$THREADS" -o "$BAM"
    samtools index -@ "$THREADS" "$BAM"
  elif [[ ! -f "$BAI" ]]; then
    samtools index -@ "$THREADS" "$BAM"
  fi
  
  # Variant calling with lofreq
  if [[ ! -f "$VCF" ]]; then
    TEMP_VCF="${OUTDIR}/${SAMPLE}.vcf"
    lofreq call --call-indels \
      -f "$REF" \
      -o "$TEMP_VCF" \
      "$BAM"
    
    # Compress and index VCF
    bcftools view -O z -o "$VCF" "$TEMP_VCF"
    rm "$TEMP_VCF"
    tabix -p vcf "$VCF"
  fi
done

# Create collapsed.tsv
if [[ ! -f "${OUTDIR}/collapsed.tsv" ]]; then
  TEMP_COLLAPSED="${OUTDIR}/collapsed.tsv.tmp"
  {
    printf "sample\tchrom\tpos\tref\talt\taf\n"
    
    for SAMPLE in "${SAMPLES[@]}"; do
      VCF="${OUTDIR}/${SAMPLE}.vcf.gz"
      
      # Extract variants and format as TSV
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF" | \
      while IFS=$'\t' read chrom pos ref alt af; do
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$SAMPLE" "$chrom" "$pos" "$ref" "$alt" "$af"
      done
    done
  } > "$TEMP_COLLAPSED"
  mv "$TEMP_COLLAPSED" "${OUTDIR}/collapsed.tsv"
fi