#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results

mkdir -p "$OUT"
: > "$OUT/failures.log"
SURVIVORS=()
OK=0

try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# Step 1: Reference preparation
if [[ -f data/ref/chrM.fa.bwt ]]; then
  : # skip
else
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

if [[ -f data/ref/chrM.fa.fai ]]; then
  : # skip
else
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# Step 2: Per-sample processing
for s in "${SAMPLES[@]}"; do
  # Step 2a: Align + sort
  align_sample() {
    bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
      | samtools sort -@ 4 -o results/${s}.bam -
  }
  
  if [[ -f results/${s}.bam ]] && samtools quickcheck results/${s}.bam; then
    : # skip
  else
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_sample; then
      continue
    fi
  fi
  
  # Step 2b: Index BAM
  if [[ -s results/${s}.bam.bai ]]; then
    : # skip
  else
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ 4 results/${s}.bam; then
      continue
    fi
  fi
  
  # Step 2c: Variant calling
  lofreq_skip=false
  if [[ -s results/${s}.vcf.gz ]] && [[ -s results/${s}.vcf.gz.tbi ]] && bcftools view -h results/${s}.vcf.gz > /dev/null; then
    lofreq_skip=true
  elif [[ -s results/${s}.vcf ]] && bcftools view -h results/${s}.vcf > /dev/null; then
    lofreq_skip=true
  fi
  
  if ! $lofreq_skip; then
    if ! try "$s" vcf '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${s}.vcf results/${s}.bam; then
      continue
    fi
  fi
  
  # Step 2d: Compress and index VCF
  if [[ -s results/${s}.vcf.gz ]] && [[ -s results/${s}.vcf.gz.tbi ]] && bcftools view -h results/${s}.vcf.gz > /dev/null; then
    : # skip both
  else
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip -f results/${s}.vcf; then
      continue
    fi
    
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf results/${s}.vcf.gz; then
      continue
    fi
  fi
  
  # Step 2e: Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Step 3: Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

for s in "${SURVIVORS[@]}"; do
  collapse_sample() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
      | awk -v s=${s} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
  }
  
  if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_sample; then
    : # continue even if one sample fails
  fi
done

# Step 4: Summary and exit
TOTAL=${#SAMPLES[@]}

if (( OK >= TOTAL )); then
  echo "[run.sh] $OK/$TOTAL samples completed; no failures" >&2
else
  if [[ -s results/failures.log ]]; then
    FIRST_FAIL=$(head -1 results/failures.log | awk -F'\t' '{print $1" failed at "$2}')
    echo "[run.sh] $OK/$TOTAL samples completed; $FIRST_FAIL — see results/failures.log" >&2
  else
    echo "[run.sh] $OK/$TOTAL samples completed; unknown failure — see results/failures.log" >&2
  fi
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi