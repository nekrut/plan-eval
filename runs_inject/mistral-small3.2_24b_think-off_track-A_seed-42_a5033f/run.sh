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

if [[ ! -f "$REF.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

if [[ ! -f "$REF.fai" ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx "$REF"; then
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam"; then
    continue
  fi

  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
      | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
  }

  if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
    continue
  fi

  if [[ -s "results/${s}.bam.bai" ]]; then
    continue
  fi

  if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ "$THREADS" "results/${s}.bam"; then
    continue
  fi

  if [[ -f "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null; then
    continue
  fi

  if [[ -f "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null; then
    continue
  fi

  if ! try "$s" call '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"; then
    continue
  fi

  if [[ -s "results/${s}.vcf.gz" ]] && [[ -s "results/${s}.vcf.gz.tbi" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null; then
    continue
  fi

  if ! try "$s" compress '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip -f "results/${s}.vcf"; then
    continue
  fi

  if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "results/${s}.vcf.gz"; then
    continue
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"; then
    continue
  fi
done

failed_samples=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
failed_steps=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log")

if [[ -n "$failed_samples" ]]; then
  echo "[run.sh] $OK/${#SAMPLES[@]} samples completed; $failed_steps — see $OUT/failures.log" >&2
else
  echo "[run.sh] ${#SAMPLES[@]}/${#SAMPLES[@]} samples completed; no failures" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi