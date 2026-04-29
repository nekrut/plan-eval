#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
mkdir -p "$OUT"

# Failure tracking
: > "$OUT/failures.log"
SURVIVORS=()
OK=0

# Defensive try wrapper
try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  if [[ "$1" == "--" ]]; then shift; fi
  if eval "$validate"; then return 0; fi
  "$@" && eval "$validate" || printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
}

# Reference preparation (idempotent)
bwa index "$REF" && try __ref__ '[[ -s '$REF'.bwt ]]' -- bwa index "$REF"
if [[ ! -f $REF.bwt ]]; then echo "[run.sh] reference index failed"; exit 1; fi
samtools faidx "$REF" && try __ref__ '[[ -s '$REF'.fai ]]' -- samtools faidx "$REF"
if [[ ! -f $REF.fai ]]; then echo "[run.sh] reference indexing failed"; exit 1; fi

# Per-sample processing
for s in "${SAMPLES[@]}"; do
  # Align & sort BAM
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz | \
      samtools sort -@ $THREADS -o "$OUT/${s}.bam" -
  }
  try "$s" align '[[ -f "$OUT/${s}.bam" && samtools quickcheck "$OUT/${s}.bam" ]]' -- align_one

  # BAM indexing
  samtools index "$OUT/${s}.bam" && try "$s" index '[[ -s "$OUT/${s}.bam.bai" ]]' -- samtools index "$OUT/${s}.bam"
  if [[ ! -f $OUT/${s}.bam.bai ]]; then continue; fi

  # Variant calling
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "$OUT/${s}.vcf" "$OUT/${s}.bam"

  # Compression & tabix
  bgzip -f "$OUT/${s}.vcf"
  tabix -p vcf "$OUT/${s}.vcf.gz" && try "$s" vcf '[[ -s "$OUT/${s}.vcf.gz.tbi" ]]' -- tabix -p vcf "$OUT/${s}.vcf.gz"

  # Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (only over surviving samples)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/${s}.vcf.gz" | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
done

# Summary & exit
FAILED=$(cut -f1 -d, "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
if [[ $FAILED ]]; then
  echo "[run.sh] $OK/${#SAMPLES[@]} samples completed; $FAILED failed at step(s) — see results/failures.log"
else
  echo "[run.sh] 4/4 samples completed; no failures"
fi
exit 0