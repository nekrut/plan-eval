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

do_align() {
  local s="$1"
  bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
    "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
    | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
}

do_collapse() {
  local s="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" \
    | awk -v s="${s}" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
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
  # 2a. Align + sort
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- do_align "$s"; then continue; fi
  fi

  # 2b. BAM index
  if [[ -s "results/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" bam_index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"; then continue; fi
  fi

  # 2c/2d. Variant calling + bgzip + tabix
  need_lofreq=1
  need_bgzip=1
  need_tabix=1

  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] \
     && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    need_lofreq=0; need_bgzip=0; need_tabix=0
  fi

  if (( need_lofreq )); then
    if [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; then
      need_lofreq=0
    fi
  fi

  if (( need_lofreq )); then
    if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null 2>&1' \
        -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"; then
      continue
    fi
  fi

  if (( need_bgzip )); then
    if [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
      need_bgzip=0
    fi
  fi

  if (( need_bgzip )); then
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null 2>&1' \
        -- bgzip -f "results/${s}.vcf"; then
      continue
    fi
  fi

  if (( need_tabix )); then
    if [[ -s "results/${s}.vcf.gz.tbi" ]]; then
      need_tabix=0
    fi
  fi

  if (( need_tabix )); then
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "results/${s}.vcf.gz"; then
      continue
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- do_collapse "$s" || true
  done
fi

# 4. Summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  msg=$(awk -F'\t' '$1 != "__ref__" && !seen[$1]++ {printf "%s failed at step %s; ", $1, $2}' "$OUT/failures.log")
  printf '[run.sh] %d/%d samples completed; %s— see results/failures.log\n' "$OK" "$TOTAL" "$msg" >&2
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi