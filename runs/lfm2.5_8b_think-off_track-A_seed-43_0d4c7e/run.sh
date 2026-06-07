<think>
We need to produce a bash script that implements the described workflow. Must be idempotent, use only tools listed (bwa, samtools, bcftools, tabix, etc). The script must have shebang and set -euo pipefail as first two lines.

We need to ensure no absolute paths outside current working directory. All paths relative to cwd. The data/raw and data/ref directories are under the current directory (the script's location). So we can refer to them directly.

We need to define constants: THREADS=4, SAMPLES array.

Create output dir results if not exists.

Then for each sample, do steps 1-6 with guards. Then step 7 collapsed TSV.

We must be careful about quoting and escaping. Use double quotes for strings that contain variables but need to preserve literal backslashes? The script will be executed in bash; we can use single quotes for the -R string to avoid expansion, but inside double quotes we need to embed a literal backslash-t for tab. However the plan says: "Use the literal two characters \t (backslash + t) inside the double-quoted argument." That means we should write something like `bwa mem -t 4 -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA" ...`. In bash, if we use double quotes for the whole command line, the \t will be interpreted as a tab character unless we escape it. Actually inside double quotes, backslash escapes are not interpreted except for \", \\, and maybe \b, \f, etc. So to have literal \t we need to either use printf or use single quotes. But the plan says "Use the literal two characters \t (backslash + t) inside the double-quoted argument." That suggests they want us to write `\t` inside double quotes and rely on bwa's handling? Wait, bwa expects a tab character, not the two characters \t. The instruction: "Use the literal two characters \t (backslash + t) inside the double-quoted argument. Do NOT use printf, echo -e, $'\t', or a real tab." That is contradictory because if we write `\t` inside double quotes, it's just backslash-t, not a tab. However maybe they mean that when you pass to bwa, the string will be interpreted as a literal tab by the shell? Actually no, the shell does not interpret \t as tab unless you use $'...' or printf. But they explicitly say "Do NOT use ... $'\t', or a real tab." So they want us to embed the two characters backslash and t, which will be passed as such to bwa. However bwa will treat them as literal characters, not a tab. That would break alignment. But maybe they assume that bwa's -R argument can accept \t as a placeholder for a tab? I'm not sure. Let's read the plan again: "Use the literal two characters `\t` (backslash + t) inside the double-quoted argument." It might be that they want us to write `\t` in the script, and bwa will interpret it as a tab because of some internal handling? Actually bwa's -R expects a string; if you give it `\t`, it will treat it as two characters: backslash and t. That would not produce proper alignment. However maybe they intend that we use `printf` to generate the string, but they forbid using printf. Hmm.

Given typical usage, one would do `-R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"` with a literal tab character. But the instruction says not to use real tabs. So we need to embed `\t` as two characters. In many tools like bcftools, you can use `-F` etc. But for alignment, maybe they expect that the shell will interpret `\t` as a tab because of some environment? Actually in bash, inside double quotes, backslash escapes are not interpreted except for `\n`, `\t`, `\r`, `\b`, `\$`. Wait, I recall that in double-quoted strings, backslash escapes are interpreted: `\t` becomes a tab. Yes! In Bash, within double quotes, the escape sequences like `\t` are interpreted as a tab character. So if we write `-R "@RG\tID:sample"` inside double quotes, the shell will interpret `\t` as a tab. That's correct. So using double quotes is fine. The instruction says "Do NOT use printf, echo -e, $'\t', or a real tab." But they explicitly say to use literal two characters \t inside double-quoted argument. That seems contradictory but maybe they mean not to use `printf` to generate the string, but using double quotes with `\t` is okay because it's interpreted by the shell. So we can do that.

Thus we will write commands like:

bwa mem -t 4 -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA" ...

Now we need to ensure idempotency: each step has a guard checking for output file and if not exists, run the pipeline. For steps that produce multiple files (like bwa index), we need to check both outputs? The plan says: "Idempotency guard: `[[ -f data/ref/chrM.fa.bwt ]] || bwa index ...`". That's for step 1. For other steps, they have similar guards.

We'll implement each sample loop:

for sample in "${SAMPLES[@]}"; do
    # Step 1: index reference
    if [[ -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    else
        bwa index data/ref/chrM.fa
    fi

But we need to wrap the whole pipeline in braces? The plan says: "Guard: `[[ -f results/{sample}.bam ]] || { bwa mem ... | samtools sort ... ; }` — wrap the whole pipeline in braces so the guard covers both stages." So we can do:

if [[ ! -f results/{sample}.bam ]]; then
    bwa mem -t 4 -R "@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA" data/ref/chrM.fa data/raw/{sample}_1.fq.gz data/raw/{sample}_2.fq.gz | samtools sort -@ 4 -o results/{sample}.bam -
fi

But note that the guard uses `results/{sample}.bam` which may not exist if previous run failed partially. That's okay.

Similarly for step 2: check for .bwt file before indexing.

We need to ensure we use the correct paths relative to cwd. The data/raw and data/ref directories are under cwd.

Now step 5: variant calling. It says guard: `[[ -f results/{sample}.vcf || -f results/{sample}.vcf.gz ]] || lofreq call-parallel ...`. So we check if either .vcf or .vcf.gz exists, then run lofreq; else run lofreq (which will produce .vcf). But note that after step 6 they compress and index. However for variant calling we need the BAM file (already indexed) and reference. The guard ensures we don't re-run if vcf already present.

We also need to ensure that the output files are placed in results/. The script should not produce any extra output beyond what tools naturally emit. So we can omit echo statements.

Now step 6: compress and index. Guard: `[[ -f results/{sample}.vcf.gz.tbi ]] || { bgzip -f results/{sample}.vcf && tabix -p vcf results/{sample}.vcf.gz ; }`. This will compress the .vcf to .vcf.gz, then run tabix on that file. The guard ensures we only compress if not already present.

Now step 7: collapsed TSV. Must always overwrite header and append per sample rows.

We need to ensure that the script is idempotent for all steps. Also note that after step 6, .vcf.gz exists, but .vcf may have been deleted? Actually they say "the guard: ... || { bgzip -f results/{sample}.vcf && tabix -p vcf results/{sample}.vcf.gz ; }". This runs both commands sequentially; if bgzip fails, the whole thing fails. But we assume it works.

Now we need to produce a single bash script with all these steps.

We must not include any extra commentary. Provide only the script.

Let's construct:

#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${SAMPLES[@]}"; do
    # Step 1: reference indexing
    if [[ ! -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    fi

    # Step 2: samtools faidx
    if [[ ! -f data/ref/chrM.fa.fai ]]; then
        samtools faidx data/ref/chrM.fa
    fi

    # Step 3: alignment + sort
    bwa mem -t $THREADS -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam -
fi

    # Step 4: BAM index
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ $THREADS results/${sample}.bam
    fi

    # Step 5: variant calling
    if [[ ! -f results/${sample}.vcf || -f results/${sample}.vcf.gz ]]; then
        lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    fi

    # Step 6: compress and index
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        bgzip -f results/${sample}.vcf && tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Collapsed TSV
{
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
} > /dev/null 2>&1
for sample in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | awk -v s=${sample} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done

But we need to ensure that the script does not produce any extra output. The above includes a line `> results/collapsed.tsv > /dev/null 2>&1` which suppresses the header print. That's okay because it's internal; but they said "Do not write to stderr/stdout beyond what tools naturally emit". Printing the header is not from tools, but we can suppress it. However maybe better to just let it go to stdout (which would be captured). But we can keep it suppressed.

Alternatively we could do:

{
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
} && true

But that still prints header. The requirement: "Write all outputs under results/ (relative to cwd)." It doesn't forbid printing header. However they said "Do not write to stderr/stdout beyond what tools naturally emit; do not use echo for user-facing prose." Printing the header is not a tool output, but it's okay? Might be considered extra output. Safer to suppress.

We also need to ensure that the script is idempotent: after step 6, .vcf.gz exists, and we don't want to re-run steps 5-6 if they already exist. The guard for step 5 checks for .vcf or .vcf.gz; if .vcf exists, skip variant calling. However note that after step 6, .vcf may have been overwritten? Actually step 6 compresses the .vcf to .vcf.gz and leaves .vcf (the uncompressed) still there? The command `bgzip -f results/${sample}.vcf` deletes the original .vcf file. So after step 6, only .vcf.gz remains. Therefore step 5 guard will see .vcf.gz but not .vcf, so it will run lofreq again, which would reprocess and overwrite .vcf (uncompressed). That's fine because we want to regenerate the uncompressed VCF each time? The plan says "check both because step 6 will delete the .vcf and leave .vcf.gz." So they expect that after step 6, .vcf is gone, so variant calling must be re-run. So our guard for step 5 should check for .vcf or .vcf.gz; if either exists, we need to run lofreq. That's what we have.

Now ensure that the script uses only tools listed: bwa, samtools, bcftools, tabix, etc. All are present.

Check for any missing dependencies: The script uses `bgzip` and `tabix`. Both are in tool inventory.

Now check if any variable expansions need quoting: In the -R string we used double quotes around the whole argument; inside we have `{sample}` which is a variable, but it's inside double quotes, so it will be expanded. That's fine.

Potential issue: The script uses `printf` for header line. That's okay because it's just writing to stdout; but they said not to write beyond what tools emit. However the header is part of output format specification, and we need to produce it. It's acceptable.

Now ensure that the script does not have any syntax errors. Let's review:

- Shebang and set -euo pipefail: okay.
- THREADS=4; SAMPLES array defined.
- mkdir -p results
- for loop over samples.

Inside loop, each step has guard with `if [[ ! -f ... ]]; then command; fi`. That ensures if output exists, skip. However note that after step 6, .vcf.gz exists but .vcf does not; we want to run variant calling and compression again. The guard for step 5 checks for .vcf or .vcf.gz; it will be false because .vcf missing, so we go into lofreq call. That's correct.

But after step 5, we may have created .vcf (uncompressed). Then step 6 compresses it to .vcf.gz and deletes .vcf. So after step 6, .vcf is gone, .vcf.gz exists. The guard for step 7 will see .vcf.gz but not .vcf, so we run bcftools query again, which will overwrite .vcf (uncompressed). That's fine.

Now check that the script does not have any stray spaces or missing quotes. For example:

`bwa mem -t $THREADS -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam -`

The `-` after the pipe indicates end of input; that's fine.

Now check for any missing newline at end of script. Not required but okay.

Potential improvements: Use `set -euo pipefail` ensures exit on error. However if a command fails (like bwa index fails), the script will exit because of `set -e`. That's okay; we want to stop on errors.

Now ensure that the script does not use any absolute paths. All paths are relative to cwd, which is the directory containing data/raw and data/ref. That's fine.

Now check if any step uses a variable that might be undefined: `$THREADS` defined earlier; we used it in loops.

Now final script:

We'll produce exactly that.
</think>
#!/usr/bin/env bash
set -euo pipefail

THREADS=4
SAMPLES=("M117-bl" "M117-ch" "M117C1-bl" "M117C1-ch")

mkdir -p results

for sample in "${SAMPLES[@]}"; do
    # Step 1: reference indexing
    if [[ ! -f data/ref/chrM.fa.bwt ]]; then
        bwa index data/ref/chrM.fa
    fi

    # Step 2: samtools faidx
    if [[ ! -f data/ref/chrM.fa.fai ]]; then
        samtools faidx data/ref/chrM.fa
    fi

    # Step 3: alignment + sort
    bwa mem -t $THREADS -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" data/ref/chrM.fa data/raw/${sample}_1.fq.gz data/raw/${sample}_2.fq.gz | samtools sort -@ $THREADS -o results/${sample}.bam -
fi

    # Step 4: BAM index
    if [[ ! -f results/${sample}.bam.bai ]]; then
        samtools index -@ $THREADS results/${sample}.bam
    fi

    # Step 5: variant calling
    if [[ ! -f results/${sample}.vcf || -f results/${sample}.vcf.gz ]]; then
        lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/${sample}.vcf results/${sample}.bam
    fi

    # Step 6: compress and index
    if [[ ! -f results/${sample}.vcf.gz.tbi ]]; then
        bgzip -f results/${sample}.vcf && tabix -p vcf results/${sample}.vcf.gz
    fi
done

# Collapsed TSV
{
    printf 'sample\tchrom\tpos\tref\talt\taf\n' > results/collapsed.tsv
} && true
for sample in "${SAMPLES[@]}"; do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' results/${sample}.vcf.gz | awk -v s=${sample} 'BEGIN{OFS="\t"}{print s,$0}' >> results/collapsed.tsv
done