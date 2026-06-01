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
  : # exists
else
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa; then
    echo "[run.sh] reference index failed" >&2; exit 1
  fi
fi

# 1b. samtools faidx
if [[ -f data/ref/chrM.fa.fai ]]; then
  : # exists
else
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa; then
    echo "[run.sh] reference faidx failed" >&2; exit 1
  fi
fi

# 2. Per-sample loop
for s in "${SAMPLES[@]}"; do
  # 2a. Align + sort
  bam_file="$OUT/$s.bam"
  if [[ -f "$bam_file" ]] && samtools quickcheck "$bam_file" 2>/dev/null; then
    : # skip
  else
    align_one() {
      bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ 4 -o "$bam_file" -
    }
    if ! try "$s" align "'samtools quickcheck $bam_file'" -- align_one; then
      continue
    fi
  fi

  # 2b. BAM index
  bai_file="$bam_file.bai"
  if [[ -s "$bai_file" ]]; then
    : # skip
  else
    if ! try "$s" bai "'[[ -s $bai_file ]]' " -- samtools index -@ 4 "$bam_file"; then
      continue
    fi
  fi

  # 2c. Variant calling
  vcf_file="$OUT/$s.vcf"
  vcf_gz_file="$OUT/$s.vcf.gz"
  if [[ -s "$vcf_gz_file" ]] && bcftools view -h "$vcf_gz_file" > /dev/null 2>&1; then
    : # skip 2c and 2d
  elif [[ -s "$vcf_file" ]] && bcftools view -h "$vcf_file" > /dev/null 2>&1; then
    : # skip 2c, proceed to 2d check
  else
    if ! try "$s" lofreq "'[[ -s $vcf_file ]] && bcftools view -h $vcf_file > /dev/null'" -- lofreq call-parallel --pp-threads 4 -f "$REF" -o "$vcf_file" "$bam_file"; then
      continue
    fi
  fi

  # 2d. Compress + tabix
  if [[ -s "$vcf_gz_file" ]] && [[ -s "$vcf_gz_file.tbi" ]] && bcftools view -h "$vcf_gz_file" > /dev/null 2>&1; then
    : # skip
  else
    # bgzip
    if ! try "$s" bgzip "'[[ -s $vcf_gz_file ]] && bcftools view -h $vcf_gz_file > /dev/null'" -- bgzip -f "$vcf_file"; then
      continue
    fi
    # tabix
    if ! try "$s" tabix "'[[ -s $vcf_gz_file.tbi ]]' " -- tabix -p vcf "$vcf_gz_file"; then
      continue
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  vcf_gz_file="$OUT/$s.vcf.gz"
  if ! try "$s" collapse "'[[ -s $OUT/collapsed.tsv ]]' " -- bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$vcf_gz_file" | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"; then
    # If collapse fails for a specific sample, it's already in survivors list but row missing.
    # The plan says "Wrap in try... so a single broken VCF only loses that row."
    # We don't remove from SURVIVORS here, just note failure if needed, but plan doesn't specify removing.
    # However, if the row wasn't added, the sample is effectively missing from collapsed.tsv.
    # The plan implies we just continue.
    :
  fi
done

# 4. Final summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  echo "[run.sh] ${OK}/${TOTAL} samples completed; no failures" >&2
else
  FAILED_SAMPLES=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd, -)
  FIRST_FAIL_STEP=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" | head -1)
  # If there are failures, construct the message. If FAILED_SAMPLES is empty but OK < TOTAL, it might be ref failure.
  if [[ -z "$FAILED_SAMPLES" ]]; then
    echo "[run.sh] ${OK}/${TOTAL} samples completed; reference failure — see results/failures.log" >&2
  else
    echo "[run.sh] ${OK}/${TOTAL} samples completed; ${FAILED_SAMPLES} failed at step ${FIRST_FAIL_STEP} — see results/failures.log" >&2
  fi
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi