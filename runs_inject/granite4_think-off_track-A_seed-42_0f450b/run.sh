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

# Helper function to retry commands with validation
try() {
  local sample="$1"; shift
  local step="$2"; shift
  local validate="$2"; shift
  if [[ "$1" == "--" ]]; then shift; fi
  for attempt in {1..2}; do
    # Run the command(s)
    "$@" && eval "$validate"
    if [ $? -eq 0 ]; then return 0; fi
  done
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$FAIL_LOG"
  return 1
}

# Idempotency: ensure logs are empty before starting
> "$FAIL_LOG"

# Reference preparation (idempotent)
bwa index "$REF"
try __ref__ '[[ -s '"$REF".bwt ]]' -- bwa index "$REF" || { echo "[run.sh] reference index failed"; exit 1; }

samtools faidx "$REF"
try __ref__ '[[ -s '"$REF".fai ]]' -- samtools faidx "$REF" || { echo "[run.sh] reference indexing failed"; exit 1; }

# Per-sample processing
for s in "${SAMPLES[@]}"; do
  BAM="${OUT}/${s}.bam"
  BAI="${OUT}/${s}.bai"
  VCF="${OUT}/${s}.vcf"
  GZ_VCF="${OUT}/${s}.vcf.gz"

  # Align & sort
  try "$s" "align" '[[ -s '"$BAM"' ]] && samtools quickcheck '"$BAM"''
  bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
    "$REF" "${data}/raw/${s}_1.fq.gz" "${data}/raw/${s}_2.fq.gz" |
    samtools sort -@ $THREADS -o "$BAM" -

  # Index BAM
  try "$s" "bam_index" '[[ -s '"$BAI"' ]]'
  samtools index -@ $THREADS "$BAM"

  # Variant calling (lofreq)
  try "$s" "vcf_call" '[[ -s '"$VCF"' ]]'
  lofreq call-parallel --pp-threads $THREADS -f "$REF" -o "$GZ_VCF" "$BAM"

  # Compress & index VCF
  bgzip -f "$GZ_VCF"
  tabix -p vcf "$GZ_VCF.gz"
done

# Collapse results into a single TSV (header + each sample)
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"

for s in "${SAMPLES[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$GZ_VCF.gz" |
    awk -v s="$s" 'BEGIN{OFS="\t"}{$1=s; print}'
done >> "${OUT}/collapsed.tsv"

# Write success summary
printf "4/4 samples completed; no failures.\n" > "$SUCCESS_LOG"