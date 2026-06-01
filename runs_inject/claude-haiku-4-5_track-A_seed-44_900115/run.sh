#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Paths
REF="data/ref/chrM.fa"
OUT="results"

# Initialize output directory
mkdir -p "$OUT"

# Initialize failure log (truncate, no header)
: > "$OUT/failures.log"

# Tracking
SURVIVORS=()
OK=0

# Helper function: try <sample> <step_label> <validation_cmd_string> -- <cmd...>
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# ============================================================================
# Step 1: Reference preparation
# ============================================================================

# 1a. bwa index
if [[ -f "$REF.bwt" ]]; then
  : # Skip, already indexed
else
  if ! try __ref__ bwa_index "[[ -s $REF.bwt ]]" -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ -f "$REF.fai" ]]; then
  : # Skip, already indexed
else
  if ! try __ref__ faidx "[[ -s $REF.fai ]]" -- samtools faidx "$REF"; then
    echo "[run.sh] reference faidx failed" >&2
    exit 1
  fi
fi

# ============================================================================
# Step 2: Per-sample processing
# ============================================================================

for s in "${SAMPLES[@]}"; do
  
  # 2a. Align + sort
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" > /dev/null 2>&1; then
    : # Skip, BAM already valid
  else
    align_one() {
      bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "$OUT/${s}.bam" -
    }
    if ! try "$s" align "[[ -s $OUT/${s}.bam ]] && samtools quickcheck $OUT/${s}.bam > /dev/null 2>&1" -- align_one; then
      continue
    fi
  fi
  
  # 2b. BAM index
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    : # Skip, index already exists
  else
    if ! try "$s" index "[[ -s $OUT/${s}.bam.bai ]]" -- samtools index -@ "$THREADS" "$OUT/${s}.bam"; then
      continue
    fi
  fi
  
  # 2c. Variant calling
  if [[ -s "$OUT/${s}.vcf.gz" ]] && bcftools view -h "$OUT/${s}.vcf.gz" > /dev/null 2>&1; then
    : # Skip, lofreq already done
  elif [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" > /dev/null 2>&1; then
    : # Skip, lofreq already done
  else
    if ! try "$s" lofreq "[[ -s $OUT/${s}.vcf ]] && bcftools view -h $OUT/${s}.vcf > /dev/null 2>&1" -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"; then
      continue
    fi
  fi
  
  # 2d. Compress + tabix
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" > /dev/null 2>&1; then
    : # Skip, both already exist and valid
  else
    # bgzip
    if ! try "$s" bgzip "[[ -s $OUT/${s}.vcf.gz ]] && bcftools view -h $OUT/${s}.vcf.gz > /dev/null 2>&1" -- bgzip -f "$OUT/${s}.vcf"; then
      continue
    fi
    # tabix
    if ! try "$s" tabix "[[ -s $OUT/${s}.vcf.gz.tbi ]]" -- tabix -p vcf "$OUT/${s}.vcf.gz"; then
      continue
    fi
  fi
  
  # Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
  
done

# ============================================================================
# Step 3: Collapsed TSV
# ============================================================================

# Write header (always)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

# Append rows from survivors
for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
  }
  if ! try "$s" collapse "[[ -s $OUT/collapsed.tsv ]]" -- collapse_one; then
    : # Failure logged
  fi
done

# ============================================================================
# Step 4: Final summary
# ============================================================================

TOTAL=${#SAMPLES[@]}

# Build summary message
if (( OK == TOTAL )); then
  SUMMARY="[run.sh] $OK/$TOTAL samples completed; no failures"
  EXIT_CODE=0
else
  FIRST_FAILURE=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" 2>/dev/null | head -1)
  FIRST_FAILURE="${FIRST_FAILURE:-unknown failure}"
  SUMMARY="[run.sh] $OK/$TOTAL samples completed; $FIRST_FAILURE — see results/failures.log"
  if (( OK >= 1 )); then
    EXIT_CODE=0
  else
    EXIT_CODE=1
  fi
fi

echo "$SUMMARY" >&2
exit "$EXIT_CODE"