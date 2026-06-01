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

  index_one() {
    samtools index -@ "$THREADS" "results/${s}.bam"
  }

  call_one() {
    lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o "results/${s}.vcf" "results/${s}.bam"
  }

  bgzip_one() {
    bgzip -f "results/${s}.vcf"
  }

  tabix_one() {
    tabix -p vcf "results/${s}.vcf.gz"
  }

  # 2a. align
  if [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam"; then
    :
  else
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # 2b. index BAM
  if [[ -s "results/${s}.bam.bai" ]]; then
    :
  else
    if ! try "$s" index_bam '[[ -s results/'"$s"'.bam.bai ]]' -- index_one; then
      continue
    fi
  fi

  # 2c + 2d. variant calling + compress + tabix
  if [[ -s "results/${s}.vcf.gz" && -s "results/${s}.vcf.gz.tbi" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
    :
  else
    # Step 2c — call
    if [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; then
      :
    elif [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
      :
    else
      if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- call_one; then
        continue
      fi
    fi

    # Step 2d — bgzip
    if [[ -s "results/${s}.vcf.gz" ]] && bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; then
      :
    else
      if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip_one; then
        continue
      fi
    fi

    # Step 2d — tabix
    if [[ -s "results/${s}.vcf.gz.tbi" ]]; then
      :
    else
      if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix_one; then
        continue
      fi
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv

collapse_one() {
  local sm="$1"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${sm}.vcf.gz" \
    | awk -v s="${sm}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
}

if (( ${#SURVIVORS[@]} > 0 )); then
  for s in "${SURVIVORS[@]}"; do
    if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one "$s"; then
      continue
    fi
  done
fi

# 4. Final summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  FAILED_LIST=$(cut -f1 results/failures.log | sort -u | grep -v '^__ref__$' | paste -sd, -)
  FIRST_FAIL=$(awk -F'\t' '!seen[$1]++ && $1!="__ref__"{print $1" failed at "$2}' results/failures.log | paste -sd'; ' -)
  printf '[run.sh] %d/%d samples completed; %s — see results/failures.log\n' "$OK" "$TOTAL" "${FIRST_FAIL:-$FAILED_LIST}" >&2
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi