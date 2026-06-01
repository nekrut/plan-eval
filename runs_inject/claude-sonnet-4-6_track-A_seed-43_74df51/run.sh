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
  [[ "${1:-}" == "--" ]] && shift
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
    echo "[run.sh] reference faidx failed" >&2
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do

  # 2a. Align + sort
  align_one() {
    bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
      | samtools sort -@ 4 -o "${OUT}/${s}.bam" -
  }

  if ! { [[ -f "${OUT}/${s}.bam" ]] && samtools quickcheck "${OUT}/${s}.bam"; }; then
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # 2b. Index BAM
  if [[ ! -s "${OUT}/${s}.bam.bai" ]]; then
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- \
         samtools index -@ 4 "${OUT}/${s}.bam"; then
      continue
    fi
  fi

  # 2c. Variant calling
  if { [[ -s "${OUT}/${s}.vcf.gz" ]] && [[ -s "${OUT}/${s}.vcf.gz.tbi" ]] && \
       bcftools view -h "${OUT}/${s}.vcf.gz" > /dev/null; } || \
     { [[ -s "${OUT}/${s}.vcf" ]] && bcftools view -h "${OUT}/${s}.vcf" > /dev/null; }; then
    : # skip lofreq
  else
    if ! try "$s" lofreq \
         '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- \
         lofreq call-parallel --pp-threads 4 -f "$REF" -o "${OUT}/${s}.vcf" "${OUT}/${s}.bam"; then
      continue
    fi
  fi

  # 2d. Compress + tabix
  if ! { [[ -s "${OUT}/${s}.vcf.gz" ]] && [[ -s "${OUT}/${s}.vcf.gz.tbi" ]] && \
         bcftools view -h "${OUT}/${s}.vcf.gz" > /dev/null; }; then
    if ! try "$s" bgzip \
         '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- \
         bgzip -f "${OUT}/${s}.vcf"; then
      continue
    fi
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- \
         tabix -p vcf "${OUT}/${s}.vcf.gz"; then
      continue
    fi
  fi

  # 2e. Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK + 1))

done

# 3. Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"

for s in "${SURVIVORS[@]+"${SURVIVORS[@]}"}"; do
  collapse_one() {
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${OUT}/${s}.vcf.gz" \
      | awk -v samp="$s" 'BEGIN{OFS="\t"}{print samp,$0}' >> "${OUT}/collapsed.tsv"
  }
  if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- collapse_one; then
    : # row lost; continue to next survivor
  fi
done

# 4. Final summary + exit
TOTAL=${#SAMPLES[@]}

if (( OK == TOTAL )); then
  printf '[run.sh] %d/%d samples completed; no failures\n' "$OK" "$TOTAL" >&2
else
  FAIL_MSG=$(awk -F'\t' \
    '$1 != "__ref__" && !seen[$1]++ {msgs = msgs (msgs == "" ? "" : "; ") $1 " failed at " $2} END {print msgs}' \
    "${OUT}/failures.log")
  printf '[run.sh] %d/%d samples completed; %s — see results/failures.log\n' \
    "$OK" "$TOTAL" "${FAIL_MSG:-unknown failure}" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi