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

# 1a. bwa index
if [[ ! -s "${REF}.bwt" ]]; then
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index "$REF"; then
    printf '[run.sh] reference index failed\n' >&2
    exit 1
  fi
fi

# 1b. samtools faidx
if [[ ! -s "${REF}.fai" ]]; then
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx "$REF"; then
    printf '[run.sh] reference faidx failed\n' >&2
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do

  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
      | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
  }

  lofreq_one() {
    lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"
  }

  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" \
      | awk -v sn="$s" 'BEGIN{OFS="\t"}{print sn,$0}' >> "results/collapsed.tsv"
  }

  # 2a. align + sort
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam" 2>/dev/null; then
    :
  else
    if ! try "$s" align "samtools quickcheck results/${s}.bam" -- align_one; then continue; fi
  fi

  # 2b. bam index
  if [[ ! -s "results/${s}.bam.bai" ]]; then
    if ! try "$s" bam_index "[[ -s results/${s}.bam.bai ]]" -- samtools index -@ "$THREADS" "results/${s}.bam"; then continue; fi
  fi

  # 2c/2d. lofreq + bgzip + tabix
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && bcftools view -h "results/${s}.vcf.gz" >/dev/null 2>&1; then
    :
  else
    if [[ ! -s "results/${s}.vcf.gz" ]]; then
      if [[ ! -s "results/${s}.vcf" ]] || ! bcftools view -h "results/${s}.vcf" >/dev/null 2>&1; then
        if ! try "$s" lofreq "[[ -s results/${s}.vcf ]] && bcftools view -h results/${s}.vcf > /dev/null" -- lofreq_one; then continue; fi
      fi
      if ! try "$s" bgzip "[[ -s results/${s}.vcf.gz ]] && bcftools view -h results/${s}.vcf.gz > /dev/null" -- bgzip -f "results/${s}.vcf"; then continue; fi
    fi
    if [[ ! -s "results/${s}.vcf.gz.tbi" ]]; then
      if ! try "$s" tabix "[[ -s results/${s}.vcf.gz.tbi ]]" -- tabix -p vcf "results/${s}.vcf.gz"; then continue; fi
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. collapsed.tsv
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "results/collapsed.tsv"
if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    collapse_one() {
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" \
        | awk -v sn="$s" 'BEGIN{OFS="\t"}{print sn,$0}' >> "results/collapsed.tsv"
    }
    try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one || true
  done
fi

# 4. summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )) && [[ ! -s "$OUT/failures.log" ]]; then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  failed_list="$(cut -f1 "$OUT/failures.log" 2>/dev/null | sort -u | grep -v '^__ref__$' | paste -sd, - || true)"
  first_step="$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" | paste -sd';' - || true)"
  printf '[run.sh] %d/%d samples completed; %s failed at step %s — see results/failures.log\n' "$OK" "$TOTAL" "${failed_list:-none}" "${first_step:-none}" >&2
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi