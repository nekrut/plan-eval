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

try() {  # try <sample> <step_label> <validation_cmd_string> -- <cmd...>
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# Reference preparation
bwa index "$REF"
[[ -f "$REF.bwt" ]] || try __ref__ bwa_index '[[ -s "$REF.bwt" ]]' -- bwa index "$REF"
if [ $? -ne 0 ]; then echo "[run.sh] reference index failed" >&2; exit 1; fi

samtools faidx "$REF"
[[ -f "$REF.fai" ]] || try __ref__ faidx '[[ -s "$REF.fai" ]]' -- samtools faidx "$REF"
if [ $? -ne 0 ]; then exit 1; fi

# Per-sample loop
for s in "${SAMPLES[@]}"; do
  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
  }
  if ! try "$s" align '[[ -s "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam"' -- align_one; then continue; fi

  samtools index -@ "$THREADS" "results/${s}.bam"
  [[ -s "results/${s}.bam.bai" ]] || try "$s" index '[[ -s "results/${s}.bam.bai" ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"

  lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"
  if [[ ! -f "results/${s}.vcf.gz" ]] && [[ -f "results/${s}.vcf" ]]; then
    if ! try "$s" vcf_check '[[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null'; then continue; fi
  else
    if ! try "$s" vcf_check '[[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null'; then continue; fi
  fi

  bgzip -f "results/${s}.vcf"
  try "$s" bgzip '[[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null' -- bgzip -f "results/${s}.vcf"

  tabix -p vcf "results/${s}.vcf.gz"
  try "$s" tabix '[[ -s "results/${s}.vcf.gz.tbi" ]]' -- tabix -p vcf "results/${s}.vcf.gz"

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "results/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "results/collapsed.tsv"
done

# Final summary
failed_samples=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
first_failures=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log")

if [[ "$OK" == "4" ]]; then
  echo "[run.sh] 4/4 samples completed; no failures"
else
  echo "[run.sh] $OK/$(${#SAMPLES[@]}) samples completed; $failed_samples failed — see results/failures.log"
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi