#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)
REF=data/ref/chrM.fa
OUT=results
mkdir -p "$OUT"
: > "$OUT/failures.log"
SURVIVORS=() OK=0

try() {  # try <sample> <step_label> <validation_cmd_string> -- <cmd...>
  local sample="$1" step="$2" validate="$3"; shift 3
  [[ "$1" == "--" ]] && shift
  if "$@" && eval "$validate"; then return 0; fi
  "$@" && eval "$validate" && return 0
  printf '%s\t%s\t%s\n' "$sample" "$step" "command_or_validation_failed" >> "$OUT/failures.log"
  return 1
}

bwa index data/ref/chrM.fa
try __ref__ bwa_index '[[ -s data/ref/chrM.fa.bwt ]]' -- bwa index data/ref/chrM.fa
samtools faidx data/ref/chrM.fa
try __ref__ faidx '[[ -s data/ref/chrM.fa.fai ]]' -- samtools faidx data/ref/chrM.fa

for s in "${SAMPLES[@]}"; do
  align_one() {
    bwa mem -t 4 -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      data/ref/chrM.fa data/raw/${s}_1.fq.gz data/raw/${s}_2.fq.gz \
      | samtools sort -@ 4 -o results/${s}.bam -
  }
  try "$s" align 'samtools quickcheck results/'"$s"'.bam' -- align_one

  samtools index -@ 3 results/${s}.bam
  try "$s" index '[[ -s results/'"$s"'.bam.bai ]]' -- samtools index -@ 3 results/${s}.bam

  lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${s}.vcf results/${s}.bam
  try "$s" call '[[ -s results/'"$s"'.vcf ]] && bcftools view -h results/'"$s"'.vcf > /dev/null' -- lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${s}.vcf results/${s}.bam

  bgzip -f results/${s}.vcf
  try "$s" bcf_gz '[[ -s results/'"$s"'.vcf.gz ]] && bcftools view -h results/'"$s"'.vcf.gz > /dev/null' -- bgzip -f results/${s}.vcf

  tabix -p vcf results/${s}.vcf.gz
  try "$s" tbi '[[ -s results/'"$s"'.vcf.gz.tbi ]]' -- tabix -p vcf results/${s}.vcf.gz

  SURVIVORS+=("$s") OK=$((OK+1))
done

printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
for s in "${SURVIVORS[@]}"; do
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
    | awk -v s=${s} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
  try "$s" collapse '[[ -s results/collapsed.tsv ]]' -- bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${s}.vcf.gz \
    | awk -v s=${s} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done

cut -f1 "$OUT/failures.log" | sort -u | grep -v '^__ref__$' | paste -sd,
if (( OK >= 3 )); then exit 2; else exit 4; fi