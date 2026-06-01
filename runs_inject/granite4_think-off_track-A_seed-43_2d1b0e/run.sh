#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
REF=data/ref/chrM.fa
OUT=results
mkdir -p "$OUT"
: > "$OUT/failures.log"

# Samples list
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

# Helper function to retry commands with validation
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# Idempotent validation functions for each step
bwa_index_validate() { [[ -f "$REF".bwt ]]; }
faidx_validate()     { [[ -f "$REF".fai ]]; }

# Reference preparation (run once)
if ! try __ref__ bwa_index '[[ -s '$REF'.bwt ]]' -- bwa index "$REF"; then
  echo "[run.sh] reference index failed" >&2; exit 1
fi

if ! try __ref__ faidx '[[ -s '$REF'.fai ]]' -- samtools faidx "$REF"; then
  exit 1
fi


# Per-sample processing loop
for s in "${SAMPLES[@]}"; do
  # Align & sort BAM
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz |
      samtools sort -@ $THREADS -o "$OUT/${s}.bam" -
  }
  if ! try "$s" align '[[ -f "$OUT/${s}.bam" && samtools quickcheck "$OUT/${s}.bam" ]]'
     -- align_one; then
    continue
  fi

  # BAM index
  if ! try "$s" bam_index '[[ -s "$OUT/${s}.bam.bai" ]]' -- samtools index "$OUT/${s}.bam"; then
    continue
  fi

  # Variant calling (lofreq)
  vcf_call() {
    lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"
  }
  if ! try "$s" vcf_call '[[ -s "$OUT/${s}.vcf ]] && bcftools view -h "$OUT/${s}.vcf" > /dev/null' -- vcf_call; then
    continue
  fi

  # Compress & index VCF
  bgzip() { bgzip -f "$1"; }
  tabix() { tabix -p vcf "$1.gz"; }

  if ! try "$s" vcf_compress '[[ -s "$OUT/${s}.vcf.gz && -s "$OUT/${s}.vcf.gz.tbi" ]]'
     -- bgzip "$OUT/${s}.vcf" && tabix "$OUT/${s}.vcf.gz"; then
    continue
  fi

  # Mark as survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (only over surviving samples)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" |
    awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
done

# Summary and exit
FAIL_SAMPLES=$(cut -f1 -d, "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
if [ -s "$OUT/failures.log" ]; then
  echo "[run.sh] $(( ${#SAMPLES[@]} - FAIL_SAMPLES ))/4 samples completed; $FAIL_SAMPLES failed(s) at step(s) — see $OUT/failures.log"
else
  echo "[run.sh] 4/4 samples completed; no failures"
fi
exit 0