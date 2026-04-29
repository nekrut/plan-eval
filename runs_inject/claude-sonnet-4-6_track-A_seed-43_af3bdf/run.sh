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

# --- 1a. bwa index ---
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# --- 1b. samtools faidx ---
if [[ ! -f data/ref/chrM.fa.fai ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    exit 1
  fi
fi

# --- 2. Per-sample loop ---
for s in "${SAMPLES[@]}"; do

  # 2a. Align + sort
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam" 2>/dev/null; then
    :
  else
    align_one() {
      bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        data/ref/chrM.fa "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
    }
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # 2b. BAM index
  if [[ ! -s "results/${s}.bam.bai" ]]; then
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- \
        samtools index -@ "$THREADS" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2c. Variant calling
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && \
      bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  elif [[ -f "results/${s}.vcf" ]] && \
      bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; then
    :
  else
    if ! try "$s" lofreq \
        '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- \
        lofreq call-parallel --pp-threads "$THREADS" \
        -f data/ref/chrM.fa -o "results/${s}.vcf" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2d. Compress + tabix
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && \
      bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  else
    if ! try "$s" bgzip \
        '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- \
        bgzip -f "results/${s}.vcf"; then
      continue
    fi
    if ! try "$s" tabix \
        '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- \
        tabix -p vcf "results/${s}.vcf.gz"; then
      continue
    fi
  fi

  # 2e. Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))

done

# --- 3. Collapsed TSV ---
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

if [[ ${#SURVIVORS[@]} -gt 0 ]]; then
  for s in "${SURVIVORS[@]}"; do
    collapse_one() {
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" \
        | awk -v s="${s}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
    }
    try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one || true
  done
fi

# --- 4. Final summary ---
TOTAL=${#SAMPLES[@]}

if (( OK == TOTAL )); then
  echo "[run.sh] ${OK}/${TOTAL} samples completed; no failures" >&2
else
  FAIL_SUMMARY=$(awk -F'\t' '$1 != "__ref__" && !seen[$1]++ {print $1" failed at step "$2}' \
    "$OUT/failures.log" | paste -sd'; ')
  echo "[run.sh] ${OK}/${TOTAL} samples completed; ${FAIL_SUMMARY} — see results/failures.log" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi