#!/usr/bin/env bash
set -euo pipefail

# Constants
THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
FAILURES="${OUT}/failures.log"
SURVIVORS=()
mkdir -p "$OUT"

# Defensive try function (retries once on failure)
try() {
  local sample="$1"; shift
  local step="$1"; shift
  local validate="$1"; shift # validation command string
  if [[ "$1" == "--" ]]; then shift; fi
  if eval "$validate"; then return 0; fi
  if eval "$validate" && eval "$sample $step"; then return 0; fi
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$FAILURES"
  return 1
}

# Idempotency guard for reference index
if [[ ! -f "${REF}.bwt" ]]; then
  try __ref__ bwa_index '[[ -s '${REF}'.bwt ]]' -- bwa index ${REF}
else
  # Index already exists, skip validation
  : > "${FAILURES}"
fi

# Idempotency guard for reference FASTA indexing
if [[ ! -f "${REF}.fai" ]]; then
  try __ref__ faidx '[[ -s '${REF}'.fai ]]' -- samtools faidx ${REF}
else
  : > "${FAILURES}"
fi

# Process each sample
for s in "${SAMPLES[@]}"; do
  BAM="${OUT}/${s}.bam"
  BAI="${BAM}.bai"

  # Align & sort (bwa mem + samtools sort)
  try "$s" align '[[ -f '${BAM} ]] && samtools quickcheck '${BAM}'' -- \
    bwa mem -t $THREADS -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" ${REF} data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz | samtools sort -@ $THREADS -o ${BAM} -

  # BAM index
  try "$s" bam_index '[[ -f '${BAI}' ]]' -- samtools index -@ $THREADS ${BAM}

  # Variant calling (lofreq)
  VCF="${OUT}/${s}.vcf"
  try "$s" vcf_calling '[[ -f '${VCF} ]]' -- lofreq call-parallel --pp-threads $THREADS -f ${REF} -o ${VCF} ${BAM}

  # Compress & index VCF
  BGZIP="${OUT}/${s}.vcf.gz"
  TABIX="${OUT}/${s}.vcf.gz.tbi"
  try "$s" vcf_compress '[[ -f '${BGZIP}' && -f '${TABIX}' ]]' -- bgzip -f ${VCF}
  try "$s" vcf_tabix '[[ -f '${TABIX}']]' -- tabix -p vcf ${BGZIP}

  # Mark survivor
  SURVIVORS+=("$s")
done

# Collapsed TSV (sample, chrom, pos, ref, alt, af)
: > "${OUT}/collapsed.tsv"
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "${OUT}/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' ${OUT}/${s}.vcf.gz | \
    awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "${OUT}/collapsed.tsv"
done

# Summary
OK=$(grep -c '^' <(printf '%s\n' "${SURVIVORS[@]}"))
TOTAL=${#SAMPLES[@]}
if (( OK == TOTAL )); then
  echo "[run.sh] $TOTAL/Total samples completed; no failures" >&2
else
  # Extract failed sample(s)
  FAILED=$(cut -f1 -d, "$FAILURES" | sort -u | grep -v '^__ref__$' | paste -sd, -)
  if (( OK > 0 )); then echo "[run.sh] $OK/Total samples completed; $TOTAL-OK failures: $FAILED" >&2 else echo "[run.sh] $OK/Total samples completed; ALL failures: $FAILED" >&2 fi
  exit 1
fi