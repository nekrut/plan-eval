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

if [[ ! -s "$REF.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s '"$REF"'.bwt ]]' -- bwa index "$REF"; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

if [[ ! -s "$REF.fai" ]]; then
  if ! try __ref__ faidx '[[ -s '"$REF"'.fai ]]' -- samtools faidx "$REF"; then
    printf '[run.sh] reference faidx failed\n' >&2
    exit 1
  fi
fi

align_one() {
  local s="$1"
  bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
    "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
    | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
}

collapse_one() {
  local s="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" \
    | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
}

for s in "${SAMPLES[@]}"; do
  if ! { [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam"; }; then
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one "$s"; then
      continue
    fi
  fi

  if [[ ! -s "results/${s}.bam.bai" ]]; then
    if ! try "$s" bam_index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"; then
      continue
    fi
  fi

  vcf_ok=0
  if [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; then
    vcf_ok=1
  fi
  if [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    vcf_ok=1
  fi
  if (( vcf_ok == 0 )); then
    if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"; then
      continue
    fi
  fi

  if ! { [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; }; then
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip -f "results/${s}.vcf"; then
      continue
    fi
  fi

  if [[ ! -s "results/${s}.vcf.gz.tbi" ]]; then
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "results/${s}.vcf.gz"; then
      continue
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one "$s" || true
  done
fi

TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  set +o pipefail
  first_fail=$(awk -F'\t' '!seen[$1]++ && $1 != "__ref__" {print $1" failed at "$2}' "$OUT/failures.log" | paste -sd';' -)
  set -o pipefail
  printf '[run.sh] %d/%d samples completed; %s — see results/failures.log\n' "$OK" "$TOTAL" "$first_fail" >&2
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi