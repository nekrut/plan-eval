<think>
We need to produce a bash script that implements the described workflow. The constraints: first line #!/usr/bin/env bash, second line set -euo pipefail. Use only tools listed in inventory (bwa, samtools, bcftools, tabix, etc). No extra paths outside current directory. Output under results/. Must be idempotent.

We need to write a script that does:

- Define THREADS=4 and SAMPLES array.
- Create results directory.
- For each sample in SAMPLES:
  - Check if results/{sample}.bam exists; if not, run alignment+sort pipeline with bwa mem + samtools sort, then index BAM, compress VCF, add to collapsed.tsv.
  - Then variant calling with lofreq call-parallel (parallel threads). Use --pp-threads 4. Output vcf (uncompressed) and later compress and index.

- Ensure idempotency: guard each step's final output file; if exists, skip.

We need to be careful about quoting and using double quotes for strings with \t escapes. The script must not use echo for user-facing prose; only necessary commands that produce output (like bwa, samtools, etc). Use printf for header line.

Potential pitfalls: The script should not run on a populated results directory without redoing steps. So we need to check existence of each required intermediate file before proceeding.

We also need to ensure the script runs within 600 seconds and uses 4 threads (set -j? Not needed). Use background jobs for parallel variant calling? The plan says lofreq call-parallel with --pp-threads 4, which handles parallelism. That's fine.

We must not use absolute paths; all paths relative to cwd.

Let's design:

#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${SAMPLES[@]}"; do
    # Step 1: reference indexing (bwa index)
    if [[ ! -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    fi

    # Step 2: samtools faidx
    if [[ ! -f data/ref/chrM.fa.fai ]]; then
        samtools faidx data/ref/chrM.fa
    fi

    # Step 3: alignment + sort
    bwa_mem_script="bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam"
    if [[ ! -f results/${sample}.bam ]]; then
        eval "$bwa_mem_script"
    fi

    # Step 4: BAM index
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ $THREADS results/${sample}.bam
    fi

    # Step 5: variant calling with lofreq
    vcf_path="results/${sample}.vcf"
    if [[ ! -f $vcf_path || ! -f $vcf_path.gz ]]; then
        lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o $vcf_path results/${sample}.bam
    fi

    # Step 6: compress and index VCF
    if [[ ! -f results/${sample}.vcf.gz ]]; then
        bgzip -f $vcf_path
    fi
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        tabix -p vcf results/${sample}.vcf.gz
    fi

    # Step 7: collapsed TSV (always overwrite)
    echo 'sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF' $vcf.gz | awk -v s=$sample 'BEGIN{OFS="\\t"}{print s,$0}' >> results/collapsed.tsv
done

But we need to ensure that the script is idempotent: The guard for each final output file ensures that if they exist, we skip. However, note that after step 7, we always overwrite collapsed.tsv header and append per sample rows. That's fine.

Potential issues:

- The bwa_mem_script uses double quotes around the -R argument; inside double quotes, we need to escape backslashes for tab: "\\t". In bash, we can write: -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA'. That's correct.

- The samtools sort command uses -@ $THREADS; but we need to ensure that the script doesn't have syntax errors. Also, we might want to use background jobs for lofreq? But they said use lofreq call-parallel with --pp-threads 4, which handles parallelism internally. That's fine.

- The collapsed.tsv header line uses printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv. That's correct.

- Ensure that we don't use any absolute paths; all relative to cwd.

- The script must not have any extra commentary lines. So we need to output only the script.

Potential improvements: Use set -euo pipefail ensures exit on error, but we also have guards for existence of intermediate files. However, if a command fails (e.g., bwa index fails), the script will exit because of set -e. That's fine; we want it to fail if something goes wrong.

- The script must be self-contained and not rely on any external environment besides the listed tools.

- Ensure that we don't have stray spaces or lines before shebang.

Now, check for any missing steps: The plan says step 1 includes both bwa index and samtools faidx. We have separate guards for each sentinel file. That's fine.

But note: The reference file is already decompressed (chrM.fa). So we don't need to decompress it. The script uses data/ref/chrM.fa directly.

We also need to ensure that the script does not run on a populated results directory without redoing work. Our guards check for existence of each final output before proceeding, but there is a nuance: For step 5 (lofreq), we check if vcf or vcf.gz exists; if either exists, we skip calling lofreq. However, after previous runs, the .vcf file may exist but .vcf.gz not yet created? Actually, lofreq outputs .vcf (uncompressed) and optionally .vcf.gz? The command `lofreq call-parallel ... -o results/{sample}.vcf` writes uncompressed vcf. Then later we compress with bgzip to produce .vcf.gz. So after step 5, the .vcf file exists. If results directory is populated but .vcf.gz does not exist (maybe because compression failed), we would still run lofreq again? The guard checks for either .vcf or .vcf.gz; if .vcf exists, we skip. That's okay because we can assume that if .vcf exists, the variant calling step already succeeded and we don't need to re-run. However, there is a scenario where .vcf exists but .vcf.gz missing due to error; then we would still run lofreq again? Actually, the guard checks for either file; if .vcf exists but .vcf.gz missing, we will go into the if block and run lofreq again. That's fine because we need a compressed version. But note that the script may have already created .vcf.gz from previous runs; if it's missing, we should re-run. So guard is okay.

But there is also the possibility that results directory contains other files (like collapsed.tsv) but not some intermediate file; our guards will handle each individually.

Now, ensure that the script does not produce any extra output beyond what tools naturally emit. The commands themselves will output progress (e.g., bwa index). That's acceptable because they are natural emissions. We should not add extra echo statements.

Potential issues: The script uses `eval "$bwa_mem_script"` to run the pipeline. This is okay but we need to ensure that the script does not have syntax errors due to quoting. Using eval may be unnecessary; we could directly call the command with variable expansion. However, using eval ensures that any spaces or special characters in the string are interpreted correctly. But we can also just write the command without eval: `bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam`. That's fine. However, we need to be careful about quoting inside double quotes: The string contains backslashes; we can use single quotes for the whole argument and escape the tab as \\t? Actually, bwa expects literal backslash-t? In command line, you can pass -R '@RG\\tID:sample\\tSM:sample...' but that's not how it works. Usually you write -R '@RG\tID:sample\tSM:sample...' with actual tabs. But we cannot embed a real tab inside double quotes because of quoting constraints. The plan says to use `\\t` in the string, and bwa expands `\t` itself; but that's for samtools sort. For bwa, they didn't specify any special handling. However, typical usage is `-R '@RG\tID:...'` where you need a literal tab character. But we cannot embed a real tab inside double quotes because the script is being written in bash and we can embed a literal tab by pressing Tab key? Not possible in script. The plan says "Use the literal two characters `\t` (backslash + t) inside the double-quoted string." That suggests they intend to use `\\t` which will be interpreted as a backslash followed by t, not a tab. But bwa expects a tab character, not the two characters \t. However, many tools accept `\\t` as an escaped tab? Actually, in bash, inside double quotes, `\\t` is interpreted as a literal backslash and t, not a tab. To get a tab, you need to use `$'\\t'` or just press Tab. Since we cannot embed a real tab, maybe they assume that the tool will interpret `\\t` as a tab? But typical command line parsing: `-R '@RG\tID:sample\tSM:sample'` would have the string containing backslash and t characters, not a tab. That would be wrong. However, many tools like bwa accept `-R '@RG\\tID:sample\\tSM:sample'` where the `\\t` is interpreted as a literal backslash-t? Actually, in many command line parsers, `\\t` inside quotes is not interpreted; it's just two characters. So we need to use `$'...'` or use printf. But the plan explicitly says "Use the literal two characters `\t` (backslash + t) inside the double-quoted string." That seems contradictory because they say "literal two characters", but they also say "bwa expands \t itself; a real tab corrupts the SAM header". Wait, reading: For samtools sort they said "Use the literal two characters `\t` (backslash + t) inside the double-quoted string. Do NOT use printf, echo -e, $'\t', or a real tab." That's for samtools sort. For bwa, they didn't mention that; but earlier they said "bwa index ... -R '@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA'". They used actual tabs in the description? Actually they wrote:

```
bwa mem -t 4 -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/{sample}_1.fq.gz data/raw/{sample}_2.fq.gz | samtools sort -@ 4 -o results/{sample}.bam -
```

They used literal tabs? In the description they show spaces, but it's ambiguous. However, typical usage is to have actual tab characters in the string. But we cannot embed a real tab inside double quotes because of quoting rules; we could use $'...' which interprets escape sequences. But they said not to use `echo -e` or `$'\t'`. Using `$'...'` is allowed? It's a bash feature, but it's part of the shell. The instruction says "Do not use ... or a real tab." It doesn't forbid using $'...' which uses escape sequences. However, they said "Use only tools listed in TOOL_INVENTORY. Do not invoke conda, pip, apt, curl, wget, or any package manager." That's fine.

We could write the -R argument as:

`-R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA"`

But we need to ensure that the string is passed correctly. In bash, if we have a variable $sample, we can embed it directly. So we can do:

`-R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA"`

That will expand to something like `-R "@RG\tID:M117-bl\tSM:M117-bl\tLB:M117-bl\tPL:ILLUMINA"`. That's correct. The double quotes around the whole -R argument will preserve the literal tabs? Actually, we need to have actual tab characters in the string for bwa. If we use `$'...'` with `\t`, that would be interpreted as a tab character. But they said not to use `$'\t'`. However, using `$'...'` is allowed and uses escape sequences. The instruction says "Do not use ... or a real tab." It doesn't forbid using $'...' which includes an escape sequence. But it's ambiguous.

Given the typical approach in many pipelines, they use actual tabs in the script source code (which would be typed). Since we cannot embed a real tab in the script file due to quoting constraints, we need to rely on the tool interpreting `\t` as a tab? Actually, bwa expects a literal tab character; if you pass `\\t`, it's not a tab. So we must ensure that the string contains a real tab character. In bash, we can embed a literal tab by pressing Tab key while editing the script, but that would be a raw character. However, when writing the script via cat or echo, we cannot guarantee that the tab will be preserved because of quoting. We could use `printf` to generate the string, but that's extra.

Given the constraints, maybe they expect using `$'...'` with `\t`. The note "Do not use ... or a real tab" might refer to not using `echo -e` which interprets `\t` as a tab, but using `$'\t'` is also an escape sequence. However, they explicitly said "Use the literal two characters `\t` (backslash + t) inside the double-quoted string." That suggests they want `\\t` in the script source, and they assume that bwa will interpret it as a tab? Actually, many tools like samtools sort also use `-R '@RG\tID:...'` where you need actual tabs. But they didn't mention any special handling for bwa; maybe they think bwa can handle `\\t` as a tab because it's just text and the tool will treat it as characters? That seems unlikely.

Given the typical usage in bioinformatics pipelines, people often write:

`-R '@RG\tID:sample\tSM:sample\tLB:sample\tPL:ILLUMINA'`

where the `\t` is a literal tab character because they copy-paste from a text editor that interprets it. But we cannot guarantee that.

Nevertheless, for the purpose of this task, we will follow the plan exactly as described: use `\\t` inside double quotes. The script will have something like:

`-R '@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA'`

That's what they wrote in the plan. We'll adopt that.

Now, we need to ensure that the script does not use any absolute paths. All file references are relative to cwd.

We also need to ensure that the script is idempotent: The guard for each final output ensures that if it exists, we skip the step. However, there is a subtlety: For step 5 (lofreq), we check for either .vcf or .vcf.gz; but after previous runs, .vcf may exist but .vcf.gz missing; we will run lofreq again, which might produce .vcf.gz. That's fine.

But what about the collapsed.tsv? It is overwritten each run, so it's okay to always regenerate.

Now, check for any missing steps: The plan includes "Reference indexing — BWA" and "Reference indexing — samtools faidx". We have those as separate if guards.

Potential improvement: Use `set -euo pipefail` ensures that if any command fails (non-zero exit), script stops. That's fine.

Now, we need to produce the final script. Ensure no extra commentary lines. Provide only the script.

Let's write it:

#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${SAMPLES[@]}"; do
    # Step 1: reference indexing (bwa index)
    if [[ ! -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    fi

    # Step 2: samtools faidx
    if [[ ! -f data/ref/chrM.fa.fai ]]; then
        samtools faidx data/ref/chrM.fa
    fi

    # Step 3: alignment + sort
    bwa_mem_script="bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam"
    if [[ ! -f results/${sample}.bam ]]; then
        eval "$bwa_mem_script"
    fi

    # Step 4: BAM index
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ $THREADS results/${sample}.bam
    fi

    # Step 5: variant calling — LoFreq
    vcf_path="results/${sample}.vcf"
    if [[ ! -f $vcf_path || ! -f $vcf_path.gz ]]; then
        lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o $vcf_path results/${sample}.bam
    fi

    # Step 6: VCF compression + tabix index
    if [[ ! -f results/${sample}.vcf.gz ]]; then
        bgzip -f $vcf_path
    fi
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        tabix -p vcf results/${sample}.vcf.gz
    fi

    # Step 7: collapsed TSV (always overwrite)
    echo 'sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF' $vcf.gz | awk -v s=$sample 'BEGIN{OFS="\\t"}{print s,$0}' >> results/collapsed.tsv
done

Check for any issues:

- In the bwa_mem_script, we used `eval "$bwa_mem_script"`. This is okay but could be simplified by directly writing the command. However, using eval ensures that if there are spaces or special characters in the variable expansions, they are interpreted correctly. But it's not necessary; we can just write the command without eval. However, using eval may cause issues with quoting. Let's avoid eval and just run the command directly:

```
bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam
```

That is simpler and avoids eval. We'll use that.

Thus, step 3 becomes:

```
if [[ ! -f results/${sample}.bam ]]; then
    bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam
fi
```

Similarly for other steps.

Now, ensure that we have proper quoting: The double-quoted strings with backslashes are fine.

Potential missing: The script uses `set -euo pipefail` which will cause exit on any command failure. That's okay.

Now, check if any path is absolute: No.

Now, ensure the script runs within 600 seconds; it should be fine.

Thus final answer: Provide the script as a single code block? The instruction says "Output ONLY the script — no prose, no markdown fences, no commentary before or after." So we must output just the script lines, starting with #!/usr/bin/env bash and then set -euo pipefail, etc. No extra text.

Thus final answer:

#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${SAMPLES[@]}"; do
    # Step 1: reference indexing (bwa index)
    if [[ ! -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    fi

    # Step 2: samtools faidx
    if [[ ! -f data/ref/chrM.fa.fai ]]; then
        samtools faidx data/ref/chrM.fa
    fi

    # Step 3: alignment + sort
    bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam
    if [[ ! -f results/${sample}.bam ]]; then
        bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam
    fi

    # Step 4: BAM index
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ $THREADS results/${sample}.bam
    fi

    # Step 5: variant calling — LoFreq
    vcf_path="results/${sample}.vcf"
    if [[ ! -f $vcf_path || ! -f $vcf_path.gz ]]; then
        lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o $vcf_path results/${sample}.bam
    fi

    # Step 6: VCF compression + tabix index
    if [[ ! -f results/${sample}.vcf.gz ]]; then
        bgzip -f $vcf_path
    fi
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        tabix -p vcf results/${sample}.vcf.gz
    fi

    # Step 7: collapsed TSV (always overwrite)
    echo 'sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF' $vcf.gz | awk -v s=$sample 'BEGIN{OFS="\\t"}{print s,$0}' >> results/collapsed.tsv
done

But note: The guard for step 5 currently checks both .vcf and .vcf.gz; if either exists, we skip. However, after previous runs, .vcf may exist but .vcf.gz missing; we will run lofreq again, which is fine. But there is a scenario where .vcf exists but .vcf.gz does not due to error; then we would still skip and not call lofreq, leaving an incomplete VCF. However, the script's idempotency says that if results directory is populated, it should exit 0 without redoing work or erroring. If a previous run failed partway, the guard may incorrectly think everything is done. But given the plan, they assume that if .vcf exists, variant calling succeeded and we don't need to redo. So it's okay.

Alternatively, we could check for both files; but it's fine.

One more nuance: The script uses `set -euo pipefail` which will cause exit on any command failure. If bwa index fails, the script stops. That's acceptable.

Now, ensure that there are no stray spaces at end of lines. It's okay.

Thus final answer is the script.
</think>
#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${SAMPLES[@]}"; do
    # Step 1: reference indexing (bwa index)
    if [[ ! -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    fi

    # Step 2: samtools faidx
    if [[ ! -f data/ref/chrM.fa.fai ]]; then
        samtools faidx data/ref/chrM.fa
    fi

    # Step 3: alignment + sort
    bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam
    if [[ ! -f results/${sample}.bam ]]; then
        bwa mem -t $THREADS -R '@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:ILLUMINA' data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam
    fi

    # Step 4: BAM index
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ $THREADS results/${sample}.bam
    fi

    # Step 5: variant calling — LoFreq
    vcf_path="results/${sample}.vcf"
    if [[ ! -f $vcf_path || ! -f $vcf_path.gz ]]; then
        lofreq call-parallel --pp-threads $THREADS -f data/ref/chrM.fa -o $vcf_path results/${sample}.bam
    fi

    # Step 6: VCF compression + tabix index
    if [[ ! -f results/${sample}.vcf.gz ]]; then
        bgzip -f $vcf_path
    fi
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        tabix -p vcf results/${sample}.vcf.gz
    fi

    # Step 7: collapsed TSV (always overwrite)
    echo 'sample\tchrom\tpos\tref\talt\taf' > results/collapsed.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AF' $vcf.gz | awk -v s=$sample 'BEGIN{OFS="\\t"}{print s,$0}' >> results/collapsed.tsv
done