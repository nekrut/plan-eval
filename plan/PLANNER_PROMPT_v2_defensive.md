You are designing an EXTREMELY DETAILED, ZERO-AMBIGUITY implementation recipe for a junior bash programmer who knows shell scripting fluently but has NEVER used bwa, samtools, lofreq, bcftools, or tabix before. The implementer must be able to translate your plan directly into a `run.sh` without making any independent decisions about flags, filenames, or argument ordering.

You will not write the bash code yourself — that's the implementer's job. But your plan must specify, for every command:
- exact tool name and subcommand
- every flag, with its exact value
- every input file path, with the exact form expected
- every output file path
- the exact relative ordering and what stdin/stdout pipes to where

**Critical defensive-execution requirement.** Some tool invocations will fail at runtime — transient I/O hiccups, malformed input on a single sample, environment issues, occasional silent truncation. The implementer's `run.sh` MUST be defensive against these failures. The plan must specify, in addition to the happy-path commands above:

- **Output validation after every per-sample step.** After every step that writes a file (alignment, indexing, calling, compression, tabix), validate the output before treating that step as done. Specify the exact validation command (e.g. `samtools quickcheck results/{sample}.bam`, `bcftools view -h results/{sample}.vcf.gz > /dev/null`) and the exact failure response.
- **Retry-once on transient failure.** Specify that any tool invocation that exits non-zero MUST be retried exactly once with identical arguments before being treated as a permanent failure. (Many real failures are transient: dropped network, slow disk, momentary CPU starvation.) Give the exact bash construct.
- **Per-sample isolation.** A failure of one sample MUST NOT abort the whole workflow. Other samples must still be processed end-to-end. Specify how the script should keep going for the surviving samples.
- **Failure logging.** The script must record, in a structured way, which sample × step combinations failed and why. Specify the exact format and target file (e.g. `results/failures.log` with `sample\tstep\treason\n`).
- **Final summary line on stderr.** The very last line on stderr must be a single human-readable summary like `[run.sh] 3/4 samples completed; M117C1-ch failed at step 5 (lofreq) — see results/failures.log`. The implementer must compute and emit this.
- **Exit code policy.** The script exits 0 if AT LEAST ONE sample completes the full pipeline through `results/{sample}.vcf.gz.tbi`. Otherwise exit 1. (We deliberately do NOT bail on the first error — partial success is the correct outcome under tool flakiness.)
- **`set -euo pipefail` is still required**, but every command that may legitimately fail must be wrapped so the script does not exit on its first non-zero return. Specify the exact bash idiom (`if ! cmd; then ... fi` or a helper function).

Audience constraint: assume the implementer will copy-paste your invocations into a shell script. If you say "use threads", they will not pick a number. If you say "use the read group string", they will not know the format. If you say "validate the BAM", they will not know which command. Be concrete down to the literal string.

# Workflow goal
Per-sample variant calling on 4 paired-end MiSeq amplicon samples mapped to the human mitochondrial reference (chrM, 16,569 bp). Final outputs:

  results/{sample}.bam
  results/{sample}.bam.bai
  results/{sample}.vcf.gz
  results/{sample}.vcf.gz.tbi
  results/collapsed.tsv          # columns: sample, chrom, pos, ref, alt, af
  results/failures.log           # NEW: tab-separated sample\tstep\treason rows for any per-sample failures (zero-row file if no failures)

# Inputs already on disk
- data/raw/{sample}_1.fq.gz, {sample}_2.fq.gz for sample in {M117-bl, M117-ch, M117C1-bl, M117C1-ch}
- data/ref/chrM.fa (decompressed, but NOT yet indexed)

# Tools available (only these; pinned in a conda env)
{TOOL_INVENTORY}

# Required content of your plan

For each numbered step, give:

1. The exact command line (or pipeline) the implementer must run.
2. What that step's output file is (exact path).
3. An idempotency guard: the exact `[[ -f ... ]] || ...` test that would skip this step.
4. **An output-validation command** (after the step runs): the exact one-liner that returns 0 iff the output is structurally valid and non-empty. For BAM: `samtools quickcheck FILE`. For BAM index: `[[ -s FILE.bai ]]`. For VCF.GZ: `bcftools view -h FILE > /dev/null && [[ $(bcftools view -H FILE | wc -l) -ge 0 ]]`. (Zero variants is acceptable for our amplicon data; what we are checking is structural integrity, not biology.) For tabix index: `[[ -s FILE.tbi ]]`.
5. **A retry-once-then-skip pattern**: tell the implementer to retry the step exactly once on failure (re-running both the command and the validation) before logging the sample as failed at this step and continuing to the next sample.
6. Any tool-specific gotchas (e.g. "bwa rejects real tab characters in -R; the literal text must be backslash-t"; "lofreq's `call-parallel` requires `--pp-threads`"; "bgzip operates in place, removing the source file").

Cover every step end-to-end:

a. Reference preparation: `bwa index data/ref/chrM.fa` and `samtools faidx data/ref/chrM.fa`. These are once-only steps; if either fails after retry, the script must exit 1 (no samples can proceed without the reference index).
b. Per-sample alignment: `bwa mem -t 4 -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/{sample}_1.fq.gz data/raw/{sample}_2.fq.gz | samtools sort -@ 4 -o results/{sample}.bam -`. Note the exact RG string format (colons not equals; literal backslash-t separators that bwa expands internally; do NOT use printf, echo -e, or $'\t'). Validation: `samtools quickcheck results/{sample}.bam`.
c. BAM indexing: `samtools index -@ 4 results/{sample}.bam`. Output: `results/{sample}.bam.bai`. Validation: `[[ -s results/{sample}.bam.bai ]]`. Do NOT mark duplicates — this is amplicon data.
d. Variant calling: the EXACT command is `lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/{sample}.vcf results/{sample}.bam`. Validation: `[[ -s results/{sample}.vcf ]] && bcftools view -h results/{sample}.vcf > /dev/null`.
e. VCF compression and indexing: `bgzip -f results/{sample}.vcf` then `tabix -p vcf results/{sample}.vcf.gz`. Validation: `bcftools view -h results/{sample}.vcf.gz > /dev/null && [[ -s results/{sample}.vcf.gz.tbi ]]`.
f. Collapse step: only over the **surviving** samples (those that produced a valid `.vcf.gz.tbi`). Header line `printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv`, then for each surviving sample append `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/{sample}.vcf.gz | awk -v s={sample} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv`.

Required at the top of the plan:
- A `set -euo pipefail` directive expectation.
- A `THREADS=4` constant expectation.
- A `SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")` array expectation.
- A `mkdir -p results` step.
- A defensive helper or pattern (e.g. `try() { ... }`) the implementer can reuse. Specify its signature and behavior precisely.
- An empty-on-creation `results/failures.log` (header optional but consistent).

Required at the bottom:
- The exact final-summary stderr line format (see above).
- The exact exit-code policy (0 if any sample fully succeeded, 1 otherwise).
- An idempotency note: every step must guard against re-doing work; a second invocation on a populated `results/` (with no new failures) must exit 0 with no work performed.

Style: numbered Markdown list with code-fenced invocations. ≤900 words. Do not include the bash boilerplate yourself — describe what the implementer must write.
