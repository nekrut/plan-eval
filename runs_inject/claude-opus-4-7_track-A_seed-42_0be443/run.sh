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
CURRENT_SAMPLE=""

try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

align_one() {
  bwa mem -t "$THREADS" -R "@RG\tID:${CURRENT_SAMPLE}\tSM:${CURRENT_SAMPLE}\tLB:${CURRENT_SAMPLE}\tPL:ILLUMINA" \
    data/ref/chrM.fa "data/raw/${CURRENT_SAMPLE}_1.fq.gz" "data/raw/${CURRENT_SAMPLE}_2.fq.gz" \
    | samtools sort -@ "$THREADS" -o "results/${CURRENT_SAMPLE}.bam" -
}

collapse_one() {
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${CURRENT_SAMPLE}.vcf.gz" \
    | awk -v s="${CURRENT_SAMPLE}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
}

# 1a. bwa index
if [[ ! -s data/ref/chrM.fa.bwt ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ ! -s data/ref/chrM.fa.fai ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    printf '[run.sh] reference faidx failed\n' >&2
    exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do
  CURRENT_SAMPLE="$s"

  # 2a — align + sort
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck "results/'"$s"'.bam"' -- align_one; then
      continue
    fi
  fi

  # 2b — BAM index
  if [[ -s "results/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" index '[[ -s "results/'"$s"'.bam.bai" ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2c+2d — fully done short-circuit
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && bcftools view -h "results/${s}.vcf.gz" >/dev/null 2>&1; then
    SURVIVORS+=("$s")
    OK=$((OK+1))
    continue
  fi

  # 2c — lofreq (only if neither vcf nor vcf.gz is present/valid)
  if [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" >/dev/null 2>&1; then
    :
  elif [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" >/dev/null 2>&1; then
    :
  else
    if ! try "$s" lofreq '[[ -s "results/'"$s"'.vcf" ]] && bcftools view -h "results/'"$s"'.vcf" >/dev/null' -- \
        lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "results/${s}.vcf" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2d — bgzip (skip if vcf.gz already exists)
  if [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" >/dev/null 2>&1; then
    :
  else
    if ! try "$s" bgzip '[[ -s "results/'"$s"'.vcf.gz" ]] && bcftools view -h "results/'"$s"'.vcf.gz" >/dev/null' -- \
        bgzip -f "results/${s}.vcf"; then
      continue
    fi
  fi

  # 2d — tabix
  if [[ -s "results/${s}.vcf.gz.tbi" ]]; then
    :
  else
    if ! try "$s" tabix '[[ -s "results/'"$s"'.vcf.gz.tbi" ]]' -- tabix -p vcf "results/${s}.vcf.gz"; then
      continue
    fi
  fi

  # 2e — survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    CURRENT_SAMPLE="$s"
    if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one; then
      continue
    fi
  done
fi

# 4. Summary + exit
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )) && [[ ! -s "$OUT/failures.log" ]]; then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  FIRST_FAIL=$(awk -F'\t' '$1!="__ref__" && !seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" | paste -sd';' -)
  if [[ -z "$FIRST_FAIL" ]]; then
    printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
  else
    printf '[run.sh] %d/%d samples completed; %s — see results/failures.log\n' "$OK" "$TOTAL" "$FIRST_FAIL" >&2
  fi
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi