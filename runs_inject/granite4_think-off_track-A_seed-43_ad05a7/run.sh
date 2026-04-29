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
  if eval "$validate" && "$@"; then return 0; fi
  if eval "$validate" && "$@"; then return 0; else printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"; fi
}

# Reference preparation (idempotent)
bwa index "$REF"
[[ -f "$REF".bwt ]] || try __ref__ bwa_index '[[ -s "$REF".bwt ]]' -- bwa index "$REF"

samtools faidx "$REF"
[[ -f "$REF".fai ]] || try __ref__ faidx '[[ -s "$REF".fai ]]' -- samtools faidx "$REF"

# Per-sample processing
for s in "${SAMPLES[@]}"; do
  # Align & sort BAM
  align_one() {
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" "$REF" data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz | \
      samtools sort -@ $THREADS -o results/${s}.bam -
  }
  try "$s" align '[[ -f results/'"$s".bam && samtools quickcheck results/'"$s".bam ]]' -- align_one

  # BAM index
  samtools index -@ $THREADS results/${s}.bam
  [[ -s results/${s}.bai ]] || continue

  # Variant calling (lofreq)
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o results/${s}.vcf results/${s}.bam
  [[ -s results/${s}.vcf.gz ]] || continue

  # Compression & tabix
  bgzip -f results/${s}.vcf
  tabix -p vcf results/${s}.vcf.gz
  [[ -s results/${s}.vcf.gz && -s results/${s}.vcf.gz.tbi ]] || continue

  SURVIVORS+=("$s")
  ((OK++))
done

# Collapsed TSV (only over surviving samples)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz | \
    awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"
done

# Summary & exit
failed_samples=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
if [[ $failed_samples ]]; then
  echo "[run.sh] $OK/$[${#SAMPLES[@]}] samples completed; $failed_samples failed at step(s) — see results/failures.log"
else
  echo "[run.sh] 4/4 samples completed; no failures"
fi
exit 0