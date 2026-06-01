#!/usr/bin/env bash
set -euo pipefail

# Constants and paths
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
OUT="results"

# Setup
mkdir -p "$OUT"
: > "$OUT/failures.log"
SURVIVORS=()
OK=0

# Defensive helper function
try() {  # try <sample> <step_label> <validation_cmd_string> -- <cmd...>
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

# 1a: bwa index
if [[ -f "$REF.bwt" ]]; then
  : # already indexed
else
  if ! try __ref__ bwa_index "[[ -s $REF.bwt ]]" -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# 1b: samtools faidx
if [[ -f "$REF.fai" ]]; then
  : # already indexed
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
  
  # Step 2a: Align + sort → BAM
  if [[ -f "$OUT/$s.bam" ]] && samtools quickcheck "$OUT/$s.bam" >/dev/null 2>&1; then
    : # BAM already exists and is valid
  else
    align_one() {
      bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "$OUT/${s}.bam" -
    }
    if ! try "$s" align "[[ -s $OUT/$s.bam ]] && samtools quickcheck $OUT/$s.bam" -- align_one; then
      continue
    fi
  fi

  # Step 2b: Index BAM
  if [[ -s "$OUT/$s.bam.bai" ]]; then
    : # Index already exists
  else
    if ! try "$s" index "[[ -s $OUT/$s.bam.bai ]]" -- samtools index -@ "$THREADS" "$OUT/$s.bam"; then
      continue
    fi
  fi

  # Step 2c: Call variants
  if [[ -s "$OUT/$s.vcf.gz" && -s "$OUT/$s.vcf.gz.tbi" ]] && bcftools view -h "$OUT/$s.vcf.gz" >/dev/null 2>&1; then
    : # VCF already compressed and indexed
  elif [[ -s "$OUT/$s.vcf" ]] && bcftools view -h "$OUT/$s.vcf" >/dev/null 2>&1; then
    : # VCF already called
  else
    if ! try "$s" call "[[ -s $OUT/$s.vcf ]] && bcftools view -h $OUT/$s.vcf > /dev/null" -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$OUT/$s.vcf" "$OUT/$s.bam"; then
      continue
    fi
  fi

  # Step 2d: Compress and index VCF
  if [[ -s "$OUT/$s.vcf.gz" && -s "$OUT/$s.vcf.gz.tbi" ]] && bcftools view -h "$OUT/$s.vcf.gz" >/dev/null 2>&1; then
    : # Already compressed and indexed
  else
    # bgzip (only if .vcf still exists)
    if [[ -f "$OUT/$s.vcf" ]]; then
      if ! try "$s" bgzip "[[ -s $OUT/$s.vcf.gz ]] && bcftools view -h $OUT/$s.vcf.gz > /dev/null" -- bgzip -f "$OUT/$s.vcf"; then
        continue
      fi
    fi
    
    # tabix
    if ! try "$s" tabix "[[ -s $OUT/$s.vcf.gz.tbi ]]" -- tabix -p vcf "$OUT/$s.vcf.gz"; then
      continue
    fi
  fi

  # Mark as survivor
  SURVIVORS+=("$s")
  OK=$((OK + 1))
done

# ============================================================================
# Step 3: Collapsed TSV
# ============================================================================

printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" \
      | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
  }
  if ! try "$s" collapse "[[ -s $OUT/collapsed.tsv ]]" -- collapse_one; then
    : # log the failure but don't stop
  fi
done

# ============================================================================
# Step 4: Final summary
# ============================================================================

TOTAL=${#SAMPLES[@]}
if [[ $OK -eq $TOTAL ]]; then
  echo "[run.sh] $OK/$TOTAL samples completed; no failures" >&2
else
  FIRST_FAIL=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" 2>/dev/null | head -1 || echo "unknown")
  echo "[run.sh] $OK/$TOTAL samples completed; $FIRST_FAIL — see results/failures.log" >&2
fi

# Exit policy
if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi