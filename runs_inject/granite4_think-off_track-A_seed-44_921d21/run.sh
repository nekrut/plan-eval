#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
FAILURES="${OUT}/failures.log"
SURVIVORS=()
OK=0

# Defensive try wrapper
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" || printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$FAILURES"
}

# Idempotency guard: ensure failures log is empty before each run
> "${FAILURES}"

# Reference preparation (once)
bwa index "$REF"
if ! try __ref__ bwa_index '[[ -s "$REF".bwt ]]' -- bwa index "$REF"; then exit 1; fi

samtools faidx "$REF"
if ! try __ref__ faidx '[[ -s "$REF".fai ]]' -- samtools faidx "$REF"; then exit 1; fi

# Per-sample processing
for s in "${SAMPLES[@]}"; do
  # Align & sort BAM
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
      | samtools sort -@ $THREADS -o "${OUT}/${s}.bam"
  }
  try "$s" align '[[ -f "${OUT}/${s}.bam" && samtools quickcheck ${OUT}/${s}.bam ]]' -- align_one

  # BAM index
  samtools index -@ $THREADS "${OUT}/${s}.bam"
  if ! try "$s" bam_index '[[ -s "${OUT}/${s}.bai ]]'; then continue; fi

  # Variant calling (lofreq)
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "${OUT}/${s}.vcf" "${OUT}/${s}.bam"
  if ! try "$s" vcf_call '[[ -s "${OUT}/${s}.vcf ]] && bcftools view -h ${OUT}/${s}.vcf > /dev/null' -- lofreq call-parallel ...

  # Compress & index VCF
  bgzip -f "${OUT}/${s}.vcf"
  tabix -p vcf "${OUT}/${s}.vcf.gz"
  if ! try "$s" vcf_compress '[[ -s "${OUT}/${s}.vcf.gz" && -s "${OUT}/${s}.vcf.gz.tbi ]]'; then continue; fi

  # Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (only over surviving samples)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "${OUT}/${s}.vcf.gz" \
    | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "${OUT}/collapsed.tsv"
done

# Summary & exit
FAILED=$(cut -f1 "$FAILURES" | sort -u | grep -v '^__ref__$' | paste -sd,)
if (( OK == ${#SAMPLES[@]} )); then
  echo "[run.sh] 4/4 samples completed; no failures"
else
  echo "[run.sh] $OK/$[${#SAMPLES[@]}] samples completed; $FAILED failed at step(s) — see results/failures.log" >&2
fi
exit 1