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
  if eval "$validate" && eval "$sample $@"; then return 0; else
    printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
    return 1
  fi
}

# Reference preparation (idempotent)
bwa index "$REF" || try __ref__ bwa_index '[[ -s '"$REF".bwt ]]' -- bwa index "$REF"
samtools faidx "$REF" || try __ref__ faidx '[[ -s '"$REF".fai ]]' -- samtools faidx "$REF"

# Per-sample processing
for s in "${SAMPLES[@]}"; do
  # Alignment & sorting
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" \
      data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz | samtools sort -@ $THREADS -o results/${s}.bam -
  }
  try "$s" align '[[ -f results/'"$s".bam && samtools quickcheck results/'"$s".bam ]]' -- align_one

  # BAM indexing
  samtools index -@ $THREADS results/${s}.bam || continue

  # Variant calling (lofreq)
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o results/${s}.vcf results/${s}.bam || continue

  # Compression & tabix
  bgzip -f results/${s}.vcf && tabix -p vcf results/${s}.vcf.gz || continue

  # Mark as survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (sample, chrom, pos, ref, alt, AF)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
    | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done

# Summary & exit
FAILED=$(cut -f1 -q "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
if [[ $FAILED ]]; then
  echo "[run.sh] $OK/$[${#SAMPLES[@]}] samples completed; $FAILED failed at step $(awk -F'\t' '!seen[$1]++{print $2}' "$OUT/failures.log") — see results/failures.log" >&2
else
  echo "[run.sh] 4/4 samples completed; no failures" >&2
fi
exit 0