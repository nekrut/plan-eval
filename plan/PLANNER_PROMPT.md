You are designing a per-step implementation plan for a junior bioinformatician. The implementer will read your plan and produce a single bash script. Your plan must be unambiguous about ordering, tool choices, and flag values, but it should not include the bash code itself — that is the implementer's job.

Audience: a competent shell programmer who knows generic bioinformatics tools but not the specific best-practice flags for amplicon mitochondrial variant calling.

# Workflow goal
Per-sample variant calling on 4 paired-end MiSeq amplicon samples mapped to the human mitochondrial reference (chrM, 16,569 bp). The final outputs are:

  results/{sample}.bam
  results/{sample}.bam.bai
  results/{sample}.vcf.gz
  results/{sample}.vcf.gz.tbi
  results/collapsed.tsv          # columns: sample, chrom, pos, ref, alt, af

# Inputs already on disk
- data/raw/{sample}_1.fq.gz, {sample}_2.fq.gz for sample in {M117-bl, M117-ch, M117C1-bl, M117C1-ch}
- data/ref/chrM.fa (decompressed but not indexed)

# Tools available (only these; pinned in a conda env)
{TOOL_INVENTORY}

# Your task
Write a numbered, ≤10-step plan covering:
1. Reference preparation (which indexes are needed, which tools, and in what order).
2. Per-sample alignment (specify the aligner and threading). The read-group string passed to `bwa mem -R` must follow this EXACT format, character-for-character: `@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA`. Tell the implementer:
   - The separators between key and value are colons (`:`), not equals signs (`=`).
   - The separators between fields are the two literal characters `\` and `t` (a backslash followed by t). Do NOT use `printf`, `echo -e`, `$'\t'`, or any other mechanism that converts those escapes into real tab characters; `bwa` rejects real tabs with "the read group line contained literal <tab> characters". The simplest correct form is to put the whole string inside a double-quoted shell argument with literal `\t` text inside, e.g. `-R "@RG\tID:..."`.
3. BAM post-processing (sort, index; whether duplicate marking is needed for amplicon data — note that amplicon data should NOT be duplicate-marked because PCR duplicates are expected).
4. Variant calling (specify caller, threading, parallelization choice, and the output format).
5. VCF post-processing (compression and indexing; what extension and indexer).
6. Collapse step (specify exactly how to derive the collapsed.tsv: which fields, which separator, header on/off, how to attach the sample name).
7. Idempotency expectations (re-running on a populated results/ dir must exit 0 without redoing finished work).

Constraints to mention to the implementer:
- THREADS=4 throughout.
- Use the "lofreq call-parallel" subcommand for variant calling.
- Use bgzip + tabix for VCF compression and indexing (not "bcftools view -O z" — be explicit).
- Do not perform duplicate marking; this is amplicon data.
- Use bcftools query for the collapse step; the AF tag is in INFO.

Output format: a single Markdown document with a numbered list and short descriptive bullets. No code blocks. No file paths outside data/, results/, and the cwd. ≤400 words.
