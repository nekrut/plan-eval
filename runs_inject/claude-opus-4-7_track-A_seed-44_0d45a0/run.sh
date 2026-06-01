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

# 1a. bwa index
if [[ ! -f data/ref/chrM.fa.bwt ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ ! -f data/ref/chrM.fa.fai ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do

  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      data/ref/chrM.fa "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
      | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
  }

  # 2a. align + sort
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam"; then
    :
  else
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # 2b. BAM index
  if [[ -s "results/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2c. variant calling
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  elif [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; then
    :
  else
    if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "results/${s}.vcf" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2d. compress + tabix
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  else
    if [[ ! -s "results/${s}.vcf.gz" ]]; then
      if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip -f "results/${s}.vcf"; then
        continue
      fi
    fi
    if [[ ! -s "results/${s}.vcf.gz.tbi" ]]; then
      if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "results/${s}.vcf.gz"; then
        continue
      fi
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

collapse_one() {
  local sample="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sample}.vcf.gz" \
    | awk -v s="$sample" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
}

if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one "$s"; then
      continue
    fi
  done
fi

# 4. summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  FAILED_LIST=$(cut -f1 "$OUT/failures.log" 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, -)
  FIRST_FAIL=$(awk -F'\t' '!seen[$1]++ && $1!="__ref__"{printf "%s failed at step %s",$1,$2; exit}' "$OUT/failures.log" 2>/dev/null)
  if [[ -z "${FIRST_FAIL:-}" ]]; then
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