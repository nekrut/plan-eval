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

# 1a. bwa index
if [[ ! -f "${REF}.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index "$REF"; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ ! -f "${REF}.fai" ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx "$REF"; then
    printf '[run.sh] reference faidx failed\n' >&2
    exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do

  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
      | samtools sort -@ "$THREADS" -o "$OUT/${s}.bam" -
  }

  # 2a. align + sort
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck '"$OUT"'/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # 2b. bam index
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" bam_index '[[ -s '"$OUT"'/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2c. variant calling (skip whole section if .vcf.gz + .tbi already valid)
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" >/dev/null 2>&1; then
    SURVIVORS+=("$s")
    OK=$((OK + 1))
    continue
  fi

  if [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" >/dev/null 2>&1; then
    :
  elif [[ -s "$OUT/${s}.vcf.gz" ]] && bcftools view -h "$OUT/${s}.vcf.gz" >/dev/null 2>&1; then
    :
  else
    if ! try "$s" lofreq '[[ -s '"$OUT"'/'"$s"'.vcf ]] && bcftools view -h '"$OUT"'/'"$s"'.vcf > /dev/null' -- \
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2d. bgzip + tabix
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

  SURVIVORS+=("$s")
  OK=$((OK + 1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

collapse_one() {
  local cs="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${cs}.vcf.gz" \
    | awk -v s="$cs" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
}

if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    try "$s" collapse '[[ -s '"$OUT"'/collapsed.tsv ]]' -- collapse_one "$s" || true
  done
fi

# 4. Summary + exit
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )) && [[ ! -s "$OUT/failures.log" ]]; then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  failed_list=$(cut -f1 "$OUT/failures.log" 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, - || true)
  first_steps=$(awk -F'\t' '!seen[$1]++{printf "%s%s failed at %s", (NR>1?"; ":""), $1, $2} END{print ""}' "$OUT/failures.log" 2>/dev/null || true)
  if [[ -z "${failed_list:-}" ]]; then
    printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
  else
    printf '[run.sh] %d/%d samples completed; %s — see results/failures.log\n' "$OK" "$TOTAL" "$first_steps" >&2
  fi
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi