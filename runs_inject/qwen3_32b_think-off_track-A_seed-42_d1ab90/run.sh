#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")
REF=data/ref/chrM.fa
OUT=results
mkdir -p "$OUT"
: > "$OUT/failures.log"
SURVIVORS=()
OK=0

try() {
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

# Reference preparation
if [[ -f "$REF.bwt" && -s "$REF.bwt" ]]; then
  samtools quickcheck "$REF.bwt"
else
  if ! try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

if [[ -f "$REF.fai" && -s "$REF.fai" ]]; then
  samtools quickcheck "$REF.fai"
else
  if ! try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx "$REF"; then
    echo "[run.sh] reference index failed" >&2
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do
  # Step 2a: Align + sort
  if [[ -f "$OUT/$s.bam" && samtools quickcheck "$OUT/$s.bam" ]]; then
    :
  else
    align_one() {
      bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
        "$REF" "data/raw/${s}_1.fq.gz" "data/raw/${s}_2.fq.gz" \
        | samtools sort -@ 4 -o "$OUT/$s.bam" -
    }
    if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
      continue
    fi
  fi

  # Step 2b: BAM index
  if [[ -s "$OUT/$s.bam.bai" ]]; then
    :
  else
    if ! try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ 4 "$OUT/$s.bam"; then
      continue
    fi
  fi

  # Step 2c: Variant calling
  if [[ -s "$OUT/$s.vcf.gz" && -s "$OUT/$s.vcf.gz.tbi" && bcftools view -h "$OUT/$s.vcf.gz" > /dev/null ]]; then
    :
  elif [[ -s "$OUT/$s.vcf" && bcftools view -h "$OUT/$s.vcf" > /dev/null ]]; then
    :
  else
    if ! try "$s" vcf_call '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads 4 -f "$REF" -o "$OUT/$s.vcf" "$OUT/$s.bam"; then
      continue
    fi
  fi

  # Step 2d: Compress + tabix
  if [[ -s "$OUT/$s.vcf.gz" && -s "$OUT/$s.vcf.gz.tbi" && bcftools view -h "$OUT/$s.vcf.gz" > /dev/null ]]; then
    :
  else
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip -f "$OUT/$s.vcf"; then
      continue
    fi
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf "$OUT/$s.vcf.gz"; then
      continue
    fi
  fi

  # Mark survivor
  SURVIVORS+=("$s")
  OK=$((OK+1))
done

# Step 3: Collapsed TSV
printf 'sample\tchrom\tpos\tref\talt\taf\n' > "$OUT/collapsed.tsv"

for s in "${SURVIVORS[@]}"; do
  if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$OUT/$s.vcf.gz" \
    | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}' >> "$OUT/collapsed.tsv"; then
    continue
  fi
done

# Step 4: Final summary
if [[ -s "$OUT/failures.log" ]]; then
  failed_samples=$(cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,)
  first_failures=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' "$OUT/failures.log" | paste -sd'; ')
  echo "[run.sh] $OK/${#SAMPLES[@]} samples completed; $failed_samples failed — see $OUT/failures.log" >&2
else
  echo "[run.sh] ${#SAMPLES[@]}/${#SAMPLES[@]} samples completed; no failures" >&2
fi

if (( OK >= 1 )); then
  exit 0
else
  exit 1
fi