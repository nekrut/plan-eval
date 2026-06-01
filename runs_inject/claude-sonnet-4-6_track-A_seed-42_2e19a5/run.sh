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
if [[ ! -f "${REF}.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ ! -f "${REF}.fai" ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx "$REF"; then
    exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do

  # Step 2a — Align + sort
  if ! { [[ -f "${OUT}/${s}.bam" ]] && samtools quickcheck "${OUT}/${s}.bam" 2>/dev/null; }; then
    align_one() {
      bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "${OUT}/${s}.bam" -
    }
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # Step 2b — BAM index
  if ! [[ -s "${OUT}/${s}.bam.bai" ]]; then
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- \
        samtools index -@ "$THREADS" "${OUT}/${s}.bam"; then
      continue
    fi
  fi

  # Step 2c — Variant calling
  if ! { { [[ -s "${OUT}/${s}.vcf.gz" ]] && bcftools view -h "${OUT}/${s}.vcf.gz" > /dev/null 2>&1; } || \
         { [[ -s "${OUT}/${s}.vcf" ]]    && bcftools view -h "${OUT}/${s}.vcf"    > /dev/null 2>&1; }; }; then
    if ! try "$s" lofreq \
        '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- \
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" \
        -o "${OUT}/${s}.vcf" "${OUT}/${s}.bam"; then
      continue
    fi
  fi

  # Step 2d — Compress + tabix
  if ! { [[ -s "${OUT}/${s}.vcf.gz" ]] && [[ -s "${OUT}/${s}.vcf.gz.tbi" ]] && \
         bcftools view -h "${OUT}/${s}.vcf.gz" > /dev/null 2>&1; }; then
    if ! [[ -s "${OUT}/${s}.vcf.gz" ]]; then
      if ! try "$s" bgzip \
          '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- \
          bgzip -f "${OUT}/${s}.vcf"; then
        continue
      fi
    fi
    if ! [[ -s "${OUT}/${s}.vcf.gz.tbi" ]]; then
      if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- \
          tabix -p vcf "${OUT}/${s}.vcf.gz"; then
        continue
      fi
    fi
  fi

  # Step 2e — Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK + 1))

done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"

if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    collapse_one() {
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${OUT}/${s}.vcf.gz" \
        | awk -v samp="$s" 'BEGIN{OFS="\t"}{print samp,$0}' >> "${OUT}/collapsed.tsv"
    }
    try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one || true
  done
fi

# 4. Final summary + exit
TOTAL=${#SAMPLES[@]}

if (( OK == TOTAL )); then
  echo "[run.sh] ${OK}/${TOTAL} samples completed; no failures" >&2
else
  failed_msg=$(awk -F'\t' '!seen[$1]++ && $1 != "__ref__" {print $1" failed at "$2}' \
    "${OUT}/failures.log" | paste -sd', ')
  echo "[run.sh] ${OK}/${TOTAL} samples completed; ${failed_msg} — see results/failures.log" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi