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

try() {  # try <sample> <step_label> <validation_cmd_string> -- <cmd...>
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

if [[ -f $REF.bwt ]]; then
  if ! try __ref__ bwa_index '[[ -s '$REF'.bwt ]]' -- bwa index "$REF"; then
    echo "[run.sh] reference index failed" >&2; exit 1
  fi
fi

if [[ -f $REF.fai ]]; then
  if ! try __ref__ faidx '[[ -s '$REF'.fai ]]' -- samtools faidx "$REF"; then
    exit 1
  fi
fi

for s in "${SAMPLES[@]}"; do
  align_one() { bwa mem -t $THREADS -R "@RG\tID:$s\tSM:$s\tLB:$s\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz | samtools sort -@ $THREADS -o results/$s.bam -; }
  if ! try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one; then
    continue
  fi

  if [[ -s results/$s.bam.bai ]]; then
    :
  else
    if ! try "$s" index_bam '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ $THREADS results/$s.bam; then
      continue
    fi
  fi

  if [[ -s results/$s.vcf.gz && -s results/$s.vcf.gz.tbi ]]; then
    :
  else
    lofreq_call() { lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o results/$s.vcf results/$s.bam; }
    if ! try "$s" lofreq '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq_call; then
      continue
    fi

    bgzip_one() { bgzip -f results/$s.vcf; }
    if ! try "$s" bgzip '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip_one; then
      continue
    fi

    tabix_one() { tabix -p vcf results/$s.vcf.gz; }
    if ! try "$s" tabix '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix_one; then
      continue
    fi
  fi

  SURVIVORS+=("$s")
  OK=$((OK+1))
done

printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for s in "${SURVIVORS[@]}"; do
  bcftools_query() { bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/$s.vcf.gz; }
  awk_one() { awk -v s=$s 'BEGIN{OFS="\t"}{print s,$0}'; }
  if ! try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- bcftools_query | awk_one >> results/collapsed.tsv; then
    :
  fi
done

failed_samples=$(cut -f1 results/failures.log | sort -u | grep -v '^__ref__$' | paste -sd,)
first_failures=$(awk -F'\t' '!seen[$1]++{print $1" failed at "$2}' results/failures.log)
echo "[run.sh] $OK/${#SAMPLES[@]} samples completed; $failed_samples failed at step $first_failures — see results/failures.log" >&2
if (( OK >= 1 )); then exit 0; else exit 1; fi