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

# Per-sample loop
for s in "${SAMPLES[@]}"; do

  # 2a. Align + sort
  if ! { [[ -f "results/${s}.bam" ]] && samtools quickcheck "results/${s}.bam" 2>/dev/null; }; then
    align_one() {
      bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        data/ref/chrM.fa "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ "$THREADS" -o "results/${s}.bam" -
    }
    val_align='samtools quickcheck results/'"$s"'.bam'
    if ! try "$s" align "$val_align" -- align_one; then continue; fi
  fi

  # 2b. BAM index
  if [[ ! -s "results/${s}.bam.bai" ]]; then
    val_bai='[[ -s results/'"$s"'.bam.bai ]]'
    if ! try "$s" index "$val_bai" -- samtools index -@ "$THREADS" "results/${s}.bam"; then continue; fi
  fi

  # 2c. Variant calling
  vcf_skip=0
  if { [[ -s "results/${s}.vcf.gz" ]] && [[ -s "results/${s}.vcf.gz.tbi" ]] && \
       bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; }; then
    vcf_skip=1
  elif { [[ -s "results/${s}.vcf" ]] && bcftools view -h "results/${s}.vcf" > /dev/null 2>&1; }; then
    vcf_skip=1
  fi
  if (( vcf_skip == 0 )); then
    val_vcf='[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null'
    if ! try "$s" lofreq "$val_vcf" -- \
        lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "results/${s}.vcf" "results/${s}.bam"; then
      continue
    fi
  fi

  # 2d. Compress + tabix
  if ! { [[ -s "results/${s}.vcf.gz" ]] && [[ -s "results/${s}.vcf.gz.tbi" ]] && \
         bcftools view -h "results/${s}.vcf.gz" > /dev/null 2>&1; }; then
    val_bgzip='[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null'
    if ! try "$s" bgzip "$val_bgzip" -- bgzip -f "results/${s}.vcf"; then continue; fi
    val_tabix='[[ -s results/'"$s"'.vcf.gz.tbi ]]'
    if ! try "$s" tabix "$val_tabix" -- tabix -p vcf "results/${s}.vcf.gz"; then continue; fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))

done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "results/${s}.vcf.gz" \
      | awk -v s="${s}" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
  }
  val_collapse='[[ -s results/collapsed.tsv ]]'
  if ! try "$s" collapse "$val_collapse" -- collapse_one; then
    :
  fi
done

# 4. Final summary
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  fail_detail=$(awk -F'\t' '$1 != "__ref__" && !seen[$1]++ {print $1" failed at "$2; exit}' \
    "$OUT/failures.log" || true)
  printf '[run.sh] %d/%d samples completed; %s — see results/failures.log\n' \
    "$OK" "$TOTAL" "$fail_detail" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi