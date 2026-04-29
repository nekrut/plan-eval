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
  [[ "${1:-}" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

align_one() {
  local s="$1"
  bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
    "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
    | samtools sort -@ "$THREADS" -o "$OUT/${s}.bam" -
}

collapse_one() {
  local s="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" \
    | awk -v sm="$s" 'BEGIN{OFS="\t"} {print sm,$0}' >> "$OUT/collapsed.tsv"
}

if [[ ! -s "${REF}.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s "'"$REF"'.bwt" ]]' -- bwa index "$REF"; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

if [[ ! -s "${REF}.fai" ]]; then
  if ! try __ref__ faidx '[[ -s "'"$REF"'.fai" ]]' -- samtools faidx "$REF"; then
    printf '[run.sh] reference faidx failed\n' >&2
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do
  # 2a align + sort
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck "'"$OUT/${s}.bam"'"' -- align_one "$s"; then
      continue
    fi
  fi

  # 2b BAM index
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" index '[[ -s "'"$OUT/${s}.bam.bai"'" ]]' -- samtools index -@ "$THREADS" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2c+2d combined: produce vcf.gz + tbi
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" >/dev/null 2>&1; then
    SURVIVORS+=("$s")
    OK=$((OK+1))
    continue
  fi

  if [[ ! -s "$OUT/${s}.vcf.gz" ]]; then
    if ! { [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" >/dev/null 2>&1; }; then
      if ! try "$s" lofreq '[[ -s "'"$OUT/${s}.vcf"'" ]] && bcftools view -h "'"$OUT/${s}.vcf"'" > /dev/null' -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"; then
        continue
      fi
    fi
    if ! try "$s" bgzip '[[ -s "'"$OUT/${s}.vcf.gz"'" ]] && bcftools view -h "'"$OUT/${s}.vcf.gz"'" > /dev/null' -- bgzip -f "$OUT/${s}.vcf"; then
      continue
    fi
  fi

  if [[ ! -s "$OUT/${s}.vcf.gz.tbi" ]]; then
    if ! try "$s" tabix '[[ -s "'"$OUT/${s}.vcf.gz.tbi"'" ]]' -- tabix -p vcf "$OUT/${s}.vcf.gz"; then
      continue
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in ${SURVIVORS[@]+"${SURVIVORS[@]}"}; do
  try "$s" collapse '[[ -s "'"$OUT/collapsed.tsv"'" ]]' -- collapse_one "$s" || true
done

TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  failed=$(awk -F'\t' '$1!="__ref__" && !seen[$1]++ {printf "%s failed at step %s; ", $1, $2}' "$OUT/failures.log")
  printf '[run.sh] %d/%d samples completed; %s— see %s/failures.log\n' "$OK" "$TOTAL" "$failed" "$OUT" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi