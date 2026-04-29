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

# Helper function
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# 1a. bwa index
[[ -f data/ref/chrM.fa.bwt ]] || try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
  echo "[run.sh] reference index failed" >&2
  exit 1
fi

# 1b. samtools faidx
[[ -f data/ref/chrM.fa.fai ]] || try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa
if [[ ! -f data/ref/chrM.fa.fai ]]; then
  echo "[run.sh] reference faidx failed" >&2
  exit 1
fi

# Per-sample loop
for s in "${SAMPLES[@]}"; do
  
  # Step 2a: Align + sort
  align_one() {
    bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
      | samtools sort -@ 4 -o results/${s}.bam -
  }
  
  if [[ -f results/${s}.bam ]] && samtools quickcheck results/${s}.bam > /dev/null 2>&1; then
    : # skip
  else
    if ! try "$s" align "[[ -s results/${s}.bam ]] && samtools quickcheck results/${s}.bam" -- align_one; then
      continue
    fi
  fi
  
  # Step 2b: BAM index
  if [[ -s results/${s}.bam.bai ]]; then
    : # skip
  else
    if ! try "$s" index "[[ -s results/${s}.bam.bai ]]" -- samtools index -@ 4 results/${s}.bam; then
      continue
    fi
  fi
  
  # Step 2c: Variant calling
  if [[ -s results/${s}.vcf.gz && -s results/${s}.vcf.gz.tbi ]] && bcftools view -h results/${s}.vcf.gz > /dev/null 2>&1; then
    : # skip if already compressed and indexed
  elif [[ -s results/${s}.vcf ]] && bcftools view -h results/${s}.vcf > /dev/null 2>&1; then
    : # skip if VCF exists and is valid
  else
    if ! try "$s" lofreq "[[ -s results/${s}.vcf ]] && bcftools view -h results/${s}.vcf > /dev/null" -- lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${s}.vcf results/${s}.bam; then
      continue
    fi
  fi
  
  # Step 2d: Compress + tabix
  if [[ -s results/${s}.vcf.gz && -s results/${s}.vcf.gz.tbi ]] && bcftools view -h results/${s}.vcf.gz > /dev/null 2>&1; then
    : # skip both
  else
    if ! try "$s" bgzip "[[ -s results/${s}.vcf.gz ]] && bcftools view -h results/${s}.vcf.gz > /dev/null" -- bgzip -f results/${s}.vcf; then
      continue
    fi
    if ! try "$s" tabix "[[ -s results/${s}.vcf.gz.tbi ]]" -- tabix -p vcf results/${s}.vcf.gz; then
      continue
    fi
  fi
  
  # Step 2e: Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
      | awk -v s="${s}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
  }
  try "$s" collapse "[[ -s results/collapsed.tsv ]]" -- collapse_one
done

# Final summary
TOTAL=${#SAMPLES[@]}
if [[ -s "$OUT/failures.log" ]]; then
  first_fail_line=$(awk -F'\t' '!seen[$1]++' "$OUT/failures.log" | head -1)
  first_fail_sample=$(echo "$first_fail_line" | cut -f1)
  first_fail_step=$(echo "$first_fail_line" | cut -f2)
  echo "[run.sh] ${OK}/${TOTAL} samples completed; ${first_fail_sample} failed at step ${first_fail_step} — see results/failures.log" >&2
else
  echo "[run.sh] ${OK}/${TOTAL} samples completed; no failures" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi