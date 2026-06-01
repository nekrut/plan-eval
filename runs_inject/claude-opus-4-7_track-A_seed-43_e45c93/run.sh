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

align_one() {
  local sm="$1"
  bwa mem -t "$THREADS" -R "@RG\tID:${sm}\tSM:${sm}\tLB:${sm}\tPL:ILLUMINA" \
    "$REF" "data/raw/${sm}_1.fq.gz" "data/raw/${sm}_2.fq.gz" \
    | samtools sort -@ "$THREADS" -o "$OUT/${sm}.bam" -
}

collapse_one() {
  local sm="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${sm}.vcf.gz" \
    | awk -v s="$sm" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
}

# 1. Reference preparation
if [[ ! -f "${REF}.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s '"$REF"'.bwt ]]' -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

if [[ ! -f "${REF}.fai" ]]; then
  if ! try __ref__ faidx '[[ -s '"$REF"'.fai ]]' -- samtools faidx "$REF"; then
    echo "[run.sh] reference faidx failed" >&2
    exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do

  # 2a. Align + sort
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck '"$OUT"'/'"$s"'.bam' -- align_one "$s"; then
      continue
    fi
  fi

  # 2b. BAM index
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" bam_index '[[ -s '"$OUT"'/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2c. Variant calling (lofreq)
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" >/dev/null 2>&1; then
    :
  elif [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" >/dev/null 2>&1; then
    :
  else
    if ! try "$s" lofreq '[[ -s '"$OUT"'/'"$s"'.vcf ]] && bcftools view -h '"$OUT"'/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2d. bgzip + tabix
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" >/dev/null 2>&1; then
    :
  else
    if [[ ! -s "$OUT/${s}.vcf.gz" ]]; then
      if ! try "$s" bgzip '[[ -s '"$OUT"'/'"$s"'.vcf.gz ]] && bcftools view -h '"$OUT"'/'"$s"'.vcf.gz > /dev/null' -- bgzip -f "$OUT/${s}.vcf"; then
        continue
      fi
    fi
    if [[ ! -s "$OUT/${s}.vcf.gz.tbi" ]]; then
      if ! try "$s" tabix '[[ -s '"$OUT"'/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "$OUT/${s}.vcf.gz"; then
        continue
      fi
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

if [[ ${#SURVIVORS[@]} -gt 0 ]]; then
  for s in "${SURVIVORS[@]}"; do
    if ! try "$s" collapse '[[ -s '"$OUT"'/collapsed.tsv ]]' -- collapse_one "$s"; then
      :
    fi
  done
fi

# 4. Final summary + exit code
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  echo "[run.sh] ${OK}/${TOTAL} samples completed; no failures" >&2
else
  FIRST_FAILS=$(awk -F'\t' '$1 != "__ref__" && !seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" | tr '\n' ';' | sed 's/;$//')
  echo "[run.sh] ${OK}/${TOTAL} samples completed; ${FIRST_FAILS} — see results/failures.log" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi