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
  BAM="$OUT/${s}.bam"
  BAI="$OUT/${s}.bam.bai"
  VCF="$OUT/${s}.vcf"
  VCFGZ="$OUT/${s}.vcf.gz"
  TBI="$OUT/${s}.vcf.gz.tbi"

  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
      | samtools sort -@ "$THREADS" -o "$BAM" -
  }

  # 2a. align + sort
  if [[ -f "$BAM" ]] && samtools quickcheck "$BAM" 2>/dev/null; then
    :
  else
    if ! try "$s" align 'samtools quickcheck "'"$BAM"'"' -- align_one; then
      continue
    fi
  fi

  # 2b. BAM index
  if [[ -s "$BAI" ]]; then
    :
  else
    if ! try "$s" bam_index '[[ -s "'"$BAI"'" ]]' -- samtools index -@ "$THREADS" "$BAM"; then
      continue
    fi
  fi

  # 2c + 2d. variant calling, compress, tabix
  if [[ -s "$VCFGZ" && -s "$TBI" ]] && bcftools view -h "$VCFGZ" >/dev/null 2>&1; then
    :
  else
    if [[ ! ( -s "$VCF" && $(bcftools view -h "$VCF" 2>/dev/null | wc -l) -gt 0 ) ]]; then
      rm -f "$VCF"
      if ! try "$s" lofreq '[[ -s "'"$VCF"'" ]] && bcftools view -h "'"$VCF"'" > /dev/null' -- \
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "$VCF" "$BAM"; then
        continue
      fi
    fi
    if ! try "$s" bgzip '[[ -s "'"$VCFGZ"'" ]] && bcftools view -h "'"$VCFGZ"'" > /dev/null' -- \
      bgzip -f "$VCF"; then
      continue
    fi
    if ! try "$s" tabix '[[ -s "'"$TBI"'" ]]' -- tabix -p vcf "$VCFGZ"; then
      continue
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]+"${SURVIVORS[@]}"}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" \
      | awk -v sm="$s" 'BEGIN{OFS="\t"}{print sm,$0}' >> "$OUT/collapsed.tsv"
  }
  try "$s" collapse '[[ -s "'"$OUT"'/collapsed.tsv" ]]' -- collapse_one || true
done

# 4. Final summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  FAILED_LIST=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd, -)
  FIRST_FAIL=$(awk -F'\t' '!seen[$1]++ && $1!="__ref__"{printf "%s failed at step %s; ", $1, $2}' "$OUT/failures.log")
  printf '[run.sh] %d/%d samples completed; %s failed (%s) — see results/failures.log\n' \
    "$OK" "$TOTAL" "${FAILED_LIST:-none}" "${FIRST_FAIL%; }" >&2
fi

if (( OK >= 1 )); then exit 0; else exit 1; fi