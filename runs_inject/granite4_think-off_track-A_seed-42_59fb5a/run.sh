#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
FAIL_LOG="${OUT}/failures.log"
SUCCESS_LOG="${OUT}/success.log"

mkdir -p "$OUT"

# Idempotency tracking
SURVIVORS=()
declare -i OK=0

# Defensive try wrapper
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" || printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$FAIL_LOG"
  return 1
}

# Helper to check reference index existence and samtools faidx existence
check_ref() {
  local sample="$1"
  [[ -f data/ref/chrM.fa.bwt ]] && validate='[[ -s data/ref/chrM.fa.bwt ]]' || validate='[[ -n \$(samtools bwa) ]]'
  try "$sample" reference_index "$validate" -- bwa index "$REF"
}

# Step 1: Reference preparation
check_ref __ref__

# Step 2a: Align & sort per sample
for s in "${SAMPLES[@]}"; do
  # BAM generation and sorting (bwa mem + samtools sort)
  try "$s" align '[[ -f results/${s}.bam ]] && samtools quickcheck results/${s}.bam' \
    -- bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" | samtools sort -@ $THREADS -o results/${s}.bam -
  
  # Step 2b: BAM indexing
  try "$s" bam_index '[[ -f results/${s}.bam.bai ]]' -- samtools index results/${s}.bam
  
  # Step 2c: Variant calling (lofreq)
  try "$s" vcf_calling '[[ -f results/${s}.vcf.gz ]] && bcftools view -h results/${s}.vcf.gz > /dev/null' \
    -- lofreq call-parallel --pp-threads $THREADS -f "$REF" -o results/${s}.vcf results/${s}.bam
  
  # Step 2d: Compression + tabix
  bgzip -f results/${s}.vcf
  tabix -p vcf results/${s}.vcf.gz
  
  # Mark as survivor if all steps succeeded
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Step 3: Collapsed TSV (sample, chrom, pos, ref, alt, af)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done

# Step 4: Write success summary to stderr (only once script finishes)
if (( OK == ${#SAMPLES[@]} )); then
  echo "[run.sh] 4/4 samples completed; no failures" >&2
else
  # Extract failed sample names and first failing step per sample
  FAILED=$(cut -f1 "$FAIL_LOG" | sort -u | grep -v '^__ref__$' | paste -sd, -)
  FIRST_STEP=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$FAIL_LOG")
  echo "[run.sh] ${OK}/${#SAMPLES[@]} samples completed; ${FAILED} failed at step ${FIRST_STEP}" >&2
fi

exit 0