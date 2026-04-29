#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF="data/ref/chrM.fa"
OUT="results"

# Setup
mkdir -p "$OUT"
: > "$OUT/failures.log"

SURVIVORS=()
OK=0

# Defensive helper: try <sample> <step_label> <validation_cmd_string> -- <cmd...>
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# Step 1a: Reference index (bwa)
if [[ -f "$REF.bwt" ]]; then
  : # already indexed
else
  if ! try __ref__ bwa_index "[[ -s $REF.bwt ]]" -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

# Step 1b: Reference faidx (samtools)
if [[ -f "$REF.fai" ]]; then
  : # already indexed
else
  if ! try __ref__ faidx "[[ -s $REF.fai ]]" -- samtools faidx "$REF"; then
    echo "[run.sh] reference faidx failed" >&2
    exit 1
  fi
fi

# Step 2: Per-sample loop
for s in "${SAMPLES[@]}"; do

  # Step 2a: Align + sort → results/{s}.bam
  if [[ -f "$OUT/${s}.bam" ]] && samtools quickcheck "$OUT/${s}.bam" &>/dev/null; then
    : # already aligned
  else
    align_one() {
      bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ 4 -o "$OUT/${s}.bam" -
    }
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # Step 2b: BAM index → results/{s}.bam.bai
  if [[ -s "$OUT/${s}.bam.bai" ]]; then
    : # already indexed
  else
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ 4 "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # Step 2c: Variant calling → results/{s}.vcf
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" &>/dev/null; then
    : # already have compressed and indexed VCF
  elif [[ -s "$OUT/${s}.vcf" ]] && bcftools view -h "$OUT/${s}.vcf" &>/dev/null; then
    : # VCF already exists and is valid
  else
    if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null 2>&1' -- lofreq call-parallel --pp-threads 4 -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"; then
      continue
    fi
  fi

  # Step 2d: Compress + tabix → results/{s}.vcf.gz and results/{s}.vcf.gz.tbi
  if [[ -s "$OUT/${s}.vcf.gz" && -s "$OUT/${s}.vcf.gz.tbi" ]] && bcftools view -h "$OUT/${s}.vcf.gz" &>/dev/null; then
    : # already compressed and indexed
  elif [[ -s "$OUT/${s}.vcf" ]]; then
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null 2>&1' -- bgzip -f "$OUT/${s}.vcf"; then
      continue
    fi
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "$OUT/${s}.vcf.gz"; then
      continue
    fi
  else
    continue
  fi

  # Step 2e: Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))

done

# Step 3: Collapsed TSV (over survivors)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" \
      | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
  }
  try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one || true
done

# Step 4: Final summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  echo "[run.sh] $OK/$TOTAL samples completed; no failures" >&2
else
  FIRST_FAILURE=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" | head -1)
  echo "[run.sh] $OK/$TOTAL samples completed; $FIRST_FAILURE — see $OUT/failures.log" >&2
fi

# Exit policy
if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi