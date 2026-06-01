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

# 1. Reference preparation
if [[ -f "$REF.bwt" ]]; then
  :
else
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

if [[ -f "$REF.fai" ]]; then
  :
else
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    echo "[run.sh] reference faidx failed" >&2
    exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do
  # 2a. Align + sort
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" > /dev/null 2>&1; then
    :
  else
    align_one() {
      bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
        | samtools sort -@ 4 -o "$OUT/${s}.bam" -
    }
    if ! try "$s" align '[[ -s '"$OUT"'/'"$s"'.bam ]] && samtools quickcheck '"$OUT"'/'"$s"'.bam > /dev/null 2>&1' -- align_one; then
      continue
    fi
  fi

  # 2b. BAM index
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" index '[[ -s '"$OUT"'/'"$s"'.bam.bai ]]' -- samtools index -@ 4 "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2c. Variant calling
  if [[ -s "$OUT/${s}.vcf.gz" ]] && [[ -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  elif [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" > /dev/null 2>&1; then
    :
  else
    if ! try "$s" lofreq '[[ -s '"$OUT"'/'"$s"'.vcf ]] && bcftools view -h '"$OUT"'/'"$s"'.vcf > /dev/null 2>&1' -- lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o "$OUT/${s}.vcf" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # 2d. Compress + tabix
  if [[ -s "$OUT/${s}.vcf.gz" ]] && [[ -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  else
    if [[ -s "$OUT/${s}.vcf" ]]; then
      if ! try "$s" bgzip '[[ -s '"$OUT"'/'"$s"'.vcf.gz ]] && bcftools view -h '"$OUT"'/'"$s"'.vcf.gz > /dev/null 2>&1' -- bgzip -f "$OUT/${s}.vcf"; then
        continue
      fi
      if ! try "$s" tabix '[[ -s '"$OUT"'/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "$OUT/${s}.vcf.gz"; then
        continue
      fi
    fi
  fi

  # 2e. Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" \
      | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
  }
  if ! try "$s" collapse '[[ -s '"$OUT"'/collapsed.tsv ]]' -- collapse_one; then
    :
  fi
done

# 4. Final summary
TOTAL=${#SAMPLES[@]}
FAILED_LIST=$(cut -f1 "$OUT/failures.log" 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, - 2>/dev/null || echo "")

if [[ $OK -eq $TOTAL ]]; then
  echo "[run.sh] $OK/$TOTAL samples completed; no failures" >&2
  exit 0
else
  if [[ -n "$FAILED_LIST" ]]; then
    echo "[run.sh] $OK/$TOTAL samples completed; $FAILED_LIST failed — see $OUT/failures.log" >&2
  else
    echo "[run.sh] $OK/$TOTAL samples completed; no failures" >&2
  fi
  if (( OK >= 1 )); then
    exit 0
  else
    exit 1
  fi
fi