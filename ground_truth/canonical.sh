#!/usr/bin/env bash
set -euo pipefail

# Ground-truth chrM variant-calling workflow for Zenodo 5119008.
# This is the answer key against which all model-generated scripts are scored.
# Deterministic at -t 4 in the locked conda env.
#
# Inputs  (relative to cwd):
#   data/raw/{sample}_{1,2}.fq.gz   for sample in {M117-bl, M117-ch, M117C1-bl, M117C1-ch}
#   data/ref/chrM.fa                (decompressed)
# Outputs (relative to cwd, written under results/):
#   results/{sample}.bam, .bam.bai
#   results/{sample}.vcf.gz, .vcf.gz.tbi
#   results/collapsed.tsv           (sample\tchrom\tpos\tref\talt\taf)

THREADS=4
SAMPLES=(M117-bl M117-ch M117C1-bl M117C1-ch)

# Anchor on bench root regardless of caller's cwd
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$SCRIPT_DIR"   # ground_truth/ is our working dir

# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate bench

mkdir -p results

# Symlink data/ from bench root so canonical and model runs use the same cwd shape
[[ -e data ]] || ln -s "$BENCH_ROOT/data" data
REF=data/ref/chrM.fa

if [[ ! -f "${REF}.bwt" ]]; then
  bwa index "$REF"
fi
if [[ ! -f "${REF}.fai" ]]; then
  samtools faidx "$REF"
fi

for s in "${SAMPLES[@]}"; do
  R1="data/raw/${s}_1.fq.gz"
  R2="data/raw/${s}_2.fq.gz"
  BAM="results/${s}.bam"
  VCF="results/${s}.vcf.gz"

  if [[ ! -s "$BAM" ]]; then
    bwa mem -t "$THREADS" -R "@RG\tID:${s}\tSM:${s}\tLB:${s}\tPL:ILLUMINA" \
      "$REF" "$R1" "$R2" 2> "results/${s}.bwa.log" \
      | samtools sort -@ "$THREADS" -o "$BAM" -
    samtools index "$BAM"
  fi

  if [[ ! -s "$VCF" ]]; then
    lofreq call-parallel --pp-threads "$THREADS" \
      -f "$REF" -o "results/${s}.vcf" "$BAM" 2> "results/${s}.lofreq.log"
    bgzip -f "results/${s}.vcf"
    tabix -p vcf "$VCF"
  fi
done

# collapsed TSV across samples
{
  printf "sample\tchrom\tpos\tref\talt\taf\n"
  for s in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
      "results/${s}.vcf.gz" \
      | awk -v s="$s" 'BEGIN{OFS="\t"}{print s,$0}'
  done
} > results/collapsed.tsv

echo "canonical run complete: $(wc -l < results/collapsed.tsv) variant rows (incl. header)"
