You are designing an EXTREMELY DETAILED, ZERO-AMBIGUITY implementation recipe for a junior bash programmer who knows shell scripting fluently but has NEVER used bwa, samtools, lofreq, bcftools, or tabix before. The implementer must be able to translate your plan directly into a `run.sh` without making any independent decisions about flags, filenames, or argument ordering.

You will not write the bash code yourself — that's the implementer's job. But your plan must specify, for every command:
- exact tool name and subcommand
- every flag, with its exact value
- every input file path, with the exact form expected
- every output file path
- the exact relative ordering and what stdin/stdout pipes to where

Audience constraint: assume the implementer will copy-paste your invocations into a shell script. If you say "use threads", they will not pick a number. If you say "use the read group string", they will not know the format. Be concrete down to the literal string.

# Workflow goal
Per-sample variant calling on 4 paired-end MiSeq amplicon samples mapped to the human mitochondrial reference (chrM, 16,569 bp). Final outputs:

  results/{sample}.bam
  results/{sample}.bam.bai
  results/{sample}.vcf.gz
  results/{sample}.vcf.gz.tbi
  results/collapsed.tsv          # columns: sample, chrom, pos, ref, alt, af

# Inputs already on disk
- data/raw/{sample}_1.fq.gz, {sample}_2.fq.gz for sample in {M117-bl, M117-ch, M117C1-bl, M117C1-ch}
- data/ref/chrM.fa (decompressed, but NOT yet indexed)

# Tools available (only these; pinned in a conda env)
{TOOL_INVENTORY}

# Required content of your plan

For each numbered step, give:

1. The exact command line (or pipeline) the implementer must run. Use `{sample}` as a placeholder and tell them to expand it. Show the full invocation as a fenced one-liner of code, e.g. `bwa mem -t 4 -R "..." data/ref/chrM.fa data/raw/{sample}_1.fq.gz data/raw/{sample}_2.fq.gz`.
2. What that step's output file is (exact path).
3. An idempotency guard: the exact `[[ -f ... ]] || ...` test that would skip this step.
4. Any tool-specific gotchas (e.g. "bwa rejects real tab characters in -R; the literal text must be backslash-t"; "lofreq's `call-parallel` requires `--pp-threads`"; "bgzip operates in place, removing the source file").

Cover every step end-to-end:

a. Reference preparation: `bwa index data/ref/chrM.fa` and `samtools faidx data/ref/chrM.fa`. Specify what files each one produces (the BWT/SA/PAC/AMB/ANN set vs the .fai).
b. Per-sample alignment: `bwa mem -t 4 -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/{sample}_1.fq.gz data/raw/{sample}_2.fq.gz | samtools sort -@ 4 -o results/{sample}.bam -`. Note the exact RG string format (colons not equals; literal backslash-t separators that bwa expands internally; do NOT use printf, echo -e, or $'\t').
c. BAM indexing: `samtools index -@ 4 results/{sample}.bam`. Output: `results/{sample}.bam.bai`. Do NOT mark duplicates — this is amplicon data.
d. Variant calling: the EXACT command is `lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/{sample}.vcf results/{sample}.bam`. The BAM is positional (not behind a flag). The output is uncompressed `.vcf` at this stage.
e. VCF compression and indexing: `bgzip -f results/{sample}.vcf` (this removes `results/{sample}.vcf` and produces `results/{sample}.vcf.gz`), then `tabix -p vcf results/{sample}.vcf.gz` (produces `.vcf.gz.tbi`).
f. Collapse step: write a single header line `printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv`, then for each sample append `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/{sample}.vcf.gz | awk -v s={sample} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv`. Note: `%INFO/AF` (not `%AF`); awk's `OFS` must be `\t`; the bcftools format string uses real `\t` and `\n` (these are bcftools format codes, not shell escapes).

Required at the top of the plan:
- A `set -euo pipefail` directive expectation.
- A `THREADS=4` constant expectation.
- A `SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")` array expectation.
- A `mkdir -p results` step.

Required at the bottom:
- An idempotency note: every step must guard against re-doing work, and a second invocation on a populated `results/` must exit 0 with no work performed.
- A re-build note: rebuild `results/collapsed.tsv` from scratch on each run (cheap; avoids stale rows).

Style: numbered Markdown list with code-fenced invocations. ≤700 words. Do not include the bash boilerplate yourself — describe what the implementer must write.
