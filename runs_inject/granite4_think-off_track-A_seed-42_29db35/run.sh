#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
mkdir -p "$OUT"

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
if [[ ! -f $REF.bwt ]]; then printf '[run.sh] reference index failed\n' >&2; exit 1; fi
samtools faidx "$REF" && try __ref__ '[[ -s '$REF'.fai ]]' -- samtools faidx "$REF"
if [[ ! -f $REF.fai ]]; then exit 1; fi

# Per-sample processing loop
for s in "${SAMPLES[@]}"; do
  # Alignment & sorting
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
      | samtools sort -@ $THREADS -o results/${s}.bam -
  }
  try "$s" align '[[ -f results/'"$s".bam && samtools quickcheck results/'"$s".bam ]]' -- align_one

  # BAM indexing
  samtools index -@ $THREADS results/${s}.bam && try "$s" idx '[[ -s results/'"$s".bam.bai ]]' -- samtools index results/${s}.bam

  # Variant calling (lofreq)
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o results/${s}.vcf results/${s}.bam && try "$s" vcf '[[ -s results/'"$s".vcf ]]' -- \
    lofreq call-parallel --pp-threads $THREADS -f "$REF" -o results/${s}.vcf results/${s}.bam

  # Compression & tabix
  bgzip -f results/${s}.vcf && try "$s" bgzip '[[ -s results/'"$s".vcf.gz ]]' -- bgzip results/${s}.vcf
  tabix -p vcf results/${s}.vcf.gz && try "$s" tabix '[[ -s results/'"$s".vcf.gz.tbi ]]' -- tabix results/${s}.vcf.gz

  # Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (only over surviving samples)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
    | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
done

# Summary & exit
FAILED=$(cut -f1 -d, "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
if [[ $FAILED ]]; then
  echo "[run.sh] $OK/${#SAMPLES[@]} samples completed; $FAILED failed at step(s) — see results/failures.log" >&2
else
  echo "[run.sh] 4/4 samples completed; no failures" >&2
fi
exit 0