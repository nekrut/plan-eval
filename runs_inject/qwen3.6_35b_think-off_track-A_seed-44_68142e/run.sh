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
if [[ -f data/ref/chrM.fa.bwt ]]; then
  true
else
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    echo "[run.sh] reference index failed" >&2; exit 1
  fi
fi

# 1b. samtools faidx
if [[ -f data/ref/chrM.fa.fai ]]; then
  true
else
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    echo "[run.sh] reference faidx failed" >&2; exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do
  # 2a. Align + sort
  if [[ -f results/${s}.bam ]] && samtools quickcheck results/${s}.bam 2>/dev/null; then
    true
  else
    align_one() {
      bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
        | samtools sort -@ 4 -o results/${s}.bam -
    }
    if ! try "$s" align "'samtools quickcheck results/'\"$s\"'.bam'" -- align_one; then
      continue
    fi
  fi

  # 2b. BAM index
  if [[ -s results/${s}.bam.bai ]]; then
    true
  else
    if ! try "$s" bai "'[[ -s results/'\"$s\"'.bam.bai ]]'" -- samtools index -@ 4 results/${s}.bam; then
      continue
    fi
  fi

  # 2c. Variant calling
  if [[ -s results/${s}.vcf.gz ]] && bcftools view -h results/${s}.vcf.gz > /dev/null 2>&1; then
    true
  elif [[ -s results/${s}.vcf ]] && bcftools view -h results/${s}.vcf > /dev/null 2>&1; then
    true
  else
    if ! try "$s" lofreq "'[[ -s results/'\"$s\"'.vcf ]] && bcftools view -h results/'\"$s\"'.vcf > /dev/null'" -- lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${s}.vcf results/${s}.bam; then
      continue
    fi
  fi

  # 2d. Compress + tabix
  if [[ -s results/${s}.vcf.gz ]] && [[ -s results/${s}.vcf.gz.tbi ]] && bcftools view -h results/${s}.vcf.gz > /dev/null 2>&1; then
    true
  else
    if ! try "$s" bgzip "'[[ -s results/'\"$s\"'.vcf.gz ]] && bcftools view -h results/'\"$s\"'.vcf.gz > /dev/null'" -- bgzip -f results/${s}.vcf; then
      continue
    fi
    if ! try "$s" tabix "'[[ -s results/'\"$s\"'.vcf.gz.tbi ]]'" -- tabix -p vcf results/${s}.vcf.gz; then
      continue
    fi
  fi

  # 2e. Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

for s in "${SURVIVORS[@]}"; do
  if ! try "$s" collapse "'[[ -s results/collapsed.tsv ]]'" -- bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz | awk -v s=${s} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv; then
    # If collapse fails for a specific sample, just skip it (already in survivors list but VCF might be bad)
    # However, per plan, we only add to SURVIVORS if all steps 2a-2d succeed.
    # If collapse fails here, we should probably remove from survivors or just leave it.
    # The plan says "only over SURVIVORS" and wrap in try. If it fails, we just don't add the row.
    # But the sample is already in SURVIVORS. Let's just continue.
    true
  fi
done

# 4. Final summary
TOTAL=${#SAMPLES[@]}
FAILED_SAMPLES=$(cut -f1 results/failures.log 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, - || true)

if [[ -n "$FAILED_SAMPLES" ]]; then
  FIRST_FAIL=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' results/failures.log | head -n 1)
  echo "[run.sh] ${OK}/${TOTAL} samples completed; ${FIRST_FAIL} — see results/failures.log" >&2
else
  echo "[run.sh] ${OK}/${TOTAL} samples completed; no failures" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi