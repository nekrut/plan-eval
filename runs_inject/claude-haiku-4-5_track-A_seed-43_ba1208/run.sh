#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
OUT="results"

# Setup
mkdir -p "$OUT"
: > "$OUT/failures.log"
SURVIVORS=()
OK=0

# Helper function: try
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# Step 1: Reference preparation

# 1a. bwa index
if [[ -f "$REF.bwt" ]]; then
  : # already indexed
else
  if ! try __ref__ bwa_index "[[ -s $REF.bwt ]]" -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ -f "$REF.fai" ]]; then
  : # already indexed
else
  if ! try __ref__ faidx "[[ -s $REF.fai ]]" -- samtools faidx "$REF"; then
    echo "[run.sh] reference faidx failed" >&2
    exit 1
  fi
fi

# Step 2: Per-sample loop
for s in "${SAMPLES[@]}"; do
  
  # Step 2a: Align + sort → BAM
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" 2>/dev/null; then
    : # already aligned
  else
    align_one() {
      bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ 4 -o "$OUT/${s}.bam" -
    }
    if ! try "$s" align "[[ -s $OUT/${s}.bam ]] && samtools quickcheck $OUT/${s}.bam 2>/dev/null" -- align_one; then
      continue
    fi
  fi
  
  # Step 2b: Index BAM
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    : # already indexed
  else
    if ! try "$s" index "[[ -s $OUT/${s}.bam.bai ]]" -- samtools index -@ 4 "$OUT/${s}.bam"; then
      continue
    fi
  fi
  
  # Step 2c: Variant calling
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" &>/dev/null; then
    : # already complete
  elif [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" &>/dev/null; then
    : # already called, will compress below
  else
    call_vcf() {
      lofreq call-parallel --pp-threads 4 -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"
    }
    if ! try "$s" lofreq "[[ -s $OUT/${s}.vcf ]] && bcftools view -h $OUT/${s}.vcf > /dev/null 2>&1" -- call_vcf; then
      continue
    fi
  fi
  
  # Step 2d: Compress + tabix
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" &>/dev/null; then
    : # already compressed and indexed
  else
    compress_vcf() {
      bgzip -f "$OUT/${s}.vcf"
    }
    if ! try "$s" bgzip "[[ -s $OUT/${s}.vcf.gz ]] && bcftools view -h $OUT/${s}.vcf.gz > /dev/null 2>&1" -- compress_vcf; then
      continue
    fi
    
    index_vcf() {
      tabix -p vcf "$OUT/${s}.vcf.gz"
    }
    if ! try "$s" tabix "[[ -s $OUT/${s}.vcf.gz.tbi ]]" -- index_vcf; then
      continue
    fi
  fi
  
  # Step 2e: Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Step 3: Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" \
      | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
  }
  if ! try "$s" collapse "[[ -s $OUT/collapsed.tsv ]]" -- collapse_one; then
    : # continue
  fi
done

# Step 4: Final summary
FAILED_SAMPLES=$(cut -f1 "$OUT/failures.log" 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, || true)
if [[ -z "$FAILED_SAMPLES" ]]; then
  FAILED_INFO="no failures"
else
  FIRST_FAIL=$(awk -F'\t' '!seen[$1]++{print $1" failed at step "$2}' "$OUT/failures.log" 2>/dev/null | head -1)
  FAILED_INFO="$FIRST_FAIL — see results/failures.log"
fi

printf '[run.sh] %d/%d samples completed; %s\n' "$OK" "${#SAMPLES[@]}" "$FAILED_INFO" >&2

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi