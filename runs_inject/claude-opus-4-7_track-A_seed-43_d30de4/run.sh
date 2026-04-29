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
  local s="$1"
  bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
    "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
    | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
}

collapse_one() {
  local sm="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sm}.vcf.gz" \
    | awk -v s="$sm" 'BEGIN{OFS="\t"}{print s,$0}' >> "results/collapsed.tsv"
}

# --- Reference prep ---
if [[ ! -s data/ref/chrM.fa.bwt ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

if [[ ! -s data/ref/chrM.fa.fai ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    printf '[run.sh] reference faidx failed\n' >&2
    exit 1
  fi
fi

# --- Per-sample loop ---
for s in "${SAMPLES[@]}"; do
  sample_ok=1

  # 2a align
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one "$s"; then
      sample_ok=0
    fi
  fi
  (( sample_ok )) || continue

  # 2b index
  if [[ -s "results/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"; then
      sample_ok=0
    fi
  fi
  (( sample_ok )) || continue

  # 2c lofreq + 2d bgzip/tabix with fine-grained guards
  need_2c=1
  need_bgzip=1
  need_tabix=1

  if [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    need_2c=0
    need_bgzip=0
    if [[ -s "results/${s}.vcf.gz.tbi" ]]; then
      need_tabix=0
    fi
  fi
  if (( need_2c )) && [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; then
    need_2c=0
  fi

  if (( need_2c )); then
    if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- \
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"; then
      sample_ok=0
    fi
  fi
  (( sample_ok )) || continue

  if (( need_bgzip )); then
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- \
        bgzip -f "results/${s}.vcf"; then
      sample_ok=0
    fi
  fi
  (( sample_ok )) || continue

  if (( need_tabix )); then
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "results/${s}.vcf.gz"; then
      sample_ok=0
    fi
  fi
  (( sample_ok )) || continue

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# --- Collapsed TSV ---
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one "$s"; then
      :
    fi
  done
fi

# --- Summary ---
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  failed_list=$(cut -f1 "$OUT/failures.log" 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, - || true)
  first_steps=$(awk -F'\t' '$1!="__ref__" && !seen[$1]++ {printf "%s at %s ", $1, $2}' "$OUT/failures.log" 2>/dev/null || true)
  printf '[run.sh] %d/%d samples completed; %s failed (%s) — see results/failures.log\n' \
    "$OK" "$TOTAL" "${failed_list:-none}" "${first_steps:-none}" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi