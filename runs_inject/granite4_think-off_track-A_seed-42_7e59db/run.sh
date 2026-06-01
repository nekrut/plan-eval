#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
FAIL_LOG="${OUT}/failures.log"
SUCCESS_LOG="${OUT}/success.log"

mkdir -p "$OUT"

# Idempotency tracking
SURVIVORS=()
declare -i OK=0

# Defensive try wrapper
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  if [[ "$1" == "--" ]]; then shift; fi
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" || printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$FAIL_LOG"
}

# Helper to write success log entry
success() {
  local sample="$1"
  echo "${sample}" >> "$SUCCESS_LOG"
}

### Reference preparation (once)
bwa index "$REF"
if ! try __ref__ bwa_index '[[ -s '"$REF".bwt ]]' -- bwa index "$REF"; then
  printf '%s\t%s\t%s\n' "__ref__" "bwa_index" "command_or_validation_failed" >> "$FAIL_LOG"
  exit 1
fi

samtools faidx "$REF"
if ! try __ref__ faidx '[[ -s '"$REF".fai ]]' -- samtools faidx "$REF"; then
  exit 1
fi

### Per-sample processing loop
for s in "${SAMPLES[@]}"; do
  # Alignment & sorting
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" | \
      samtools sort -@ $THREADS -o "${OUT}/${s}.bam"
  }
  if ! try "$s" align '[[ -f '"${OUT}/${s}".bam ]] && samtools quickcheck '"${OUT}/${s}".bam' -- align_one; then
    continue
  fi

  # BAM indexing
  samtools index -@ $THREADS "${OUT}/${s}.bam"
  if ! try "$s" bam_index '[[ -s '"${OUT}/${s}".bai ]]' -- samtools index -@ $THREADS "${OUT}/${s}.bam"; then
    continue
  fi

  # Variant calling (lofreq)
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "${OUT}/${s}.vcf" "${OUT}/${s}.bam"
  if ! try "$s" vcf_call '[[ -s '"${OUT}/${s}".vcf ]] && bcftools view -h '"${OUT}/${s}".vcf > /dev/null' -- lofreq call-parallel ...; then
    continue
  fi

  # Compression & tabix
  bgzip -f "${OUT}/${s}.vcf"
  if ! try "$s" vcf_bgzip '[[ -s '"${OUT}/${s}".vcf.gz ]]' -- bgzip -f "${OUT}/${s}.vcf"; then
    continue
  fi

  tabix -p vcf "${OUT}/${s}.vcf.gz"
  if ! try "$s" vcf_tabix '[[ -s '"${OUT}/${s}".vcf.gz.tbi ]]' -- tabix -p vcf "${OUT}/${s}.vcf.gz"; then
    continue
  fi

  # Mark success and increment counter
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (tab-delimited with header)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${OUT}/${s}.vcf.gz" | \
    awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "${OUT}/collapsed.tsv"
done

# Summary
if (( OK == ${#SAMPLES[@]} )); then
  echo "4/4 samples completed; no failures" >&2
else
  failed=$(cut -f1 "$FAIL_LOG" | sort -u | grep -v '^__ref__$' | paste -sd,)
  if [[ $failed ]]; then
    echo "${OK}/${#SAMPLES[@]} samples completed; ${failed} failed at step(s) — see $FAIL_LOG" >&2
  else
    echo "4/4 samples completed; no failures" >&2
  fi
fi

exit 0