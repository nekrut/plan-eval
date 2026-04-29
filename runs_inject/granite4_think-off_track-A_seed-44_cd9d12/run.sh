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
bwa index data/ref/chrM.fa
try __ref__ '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa

samtools faidx data/ref/chrM.fa
try __ref__ '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa

# Per-sample processing loop
for s in "${SAMPLES[@]}"; do
  # Alignment & sorting
  align_one() {
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
      | samtools sort -@ "$THREADS" -o results/${s}.bam -
  }
  try "$s" align '[[ -f results/'"$s".bam && samtools quickcheck results/'"$s".bam ]]' -- align_one

  # BAM indexing
  samtools index -@ "$THREADS" results/${s}.bam > /dev/null
  try "$s" idx '[[ -s results/'"$s".bam.bai ]]' -- samtools index -@ "$THREADS" results/${s}.bam

  # Variant calling (lofreq)
  lofreq call-parallel --pp-threads "$THREADS" -f data/ref/chrM.fa -o results/${s}.vcf results/${s}.bam
  try "$s" vcf '[[ -s results/'"$s".vcf ]] && bcftools view -h results/'"$s".vcf > /dev/null' -- lofreq call-parallel ...

  # Compression & tabix
  bgzip -f results/${s}.vcf
  tabix -p vcf results/${s}.vcf.gz
  try "$s" comp '[[ -s results/'"$s".vcf.gz && -s results/'"$s".vcf.gz.tbi ]] && bcftools view -h results/'"$s".vcf.gz > /dev/null' -- bgzip ...

  # Mark success
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Collapsed TSV (only over surviving samples)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
    | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done

# Summary & exit
failed_samples=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
if (( OK == ${#SAMPLES[@]} )); then
  echo "[run.sh] 4/4 samples completed; no failures"
else
  echo "[run.sh] $OK/$[${#SAMPLES[@]}] samples completed; $failed_samples failed" >&2
fi
exit 0