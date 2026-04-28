# plan-eval: a benchmark for recipe-driven bioinformatics workflow execution

A controlled comparison of frontier (Anthropic Claude 4.x) and open-weight (Ollama-hosted, 4–70 B) language models at executing an Opus-authored mtDNA variant-calling recipe, scored against a published canonical answer key. Two consumer hardware platforms (Jetson AGX Orin, RTX 5080), six plan-detail levels, two tracks (with-plan, no-plan), three seeds per cell, **n=403 scored runs**, **25 distinct models**.

## Abstract

We measure how faithfully a small or local language model converts a hyper-detailed Opus-authored recipe into an executable bioinformatics workflow, holding the substrate constant: per-sample variant calling on four paired-end mtDNA Illumina samples (Zenodo 5119008) against a published Galaxy Training Network reference workflow. Primary metric is M3, the macro-mean Jaccard index of `(chrom, pos, ref, alt)` PASS calls vs the canonical VCF, with an AF tolerance window of ±0.02. Three quantified findings: (i) plan-detail dominates model capability — going from a lean ~1200-token plan (v1) to a hyper-detailed ~1150-byte/command plan (v2) flips most local models from M3 ≈ 0.0 to 1.000 and leaves Anthropic models unchanged at 1.000; (ii) on Jetson MAXN, **13 of 14 free open-weight models score M3 = 1.000 ± 0.000 on v2 Track A**, with the smallest perfect model (`granite4`, 2.1 GB on disk) at 15 s/seed; (iii) the v1→v2 cliff for ≥27 B dense local models reduces to a single command line — adding the literal `lofreq call-parallel` invocation to v1 (variant v1.25) brings them all to 1.000, falsifying a "model size" hypothesis in favour of a "tool-CLI specificity" one. Tracks with no plan collapse all open-weight models to M3 = 0.000 ± 0.000 while leaving frontier Anthropic at M3 ≥ 0.667. Mechanically extracting the lofreq snippet from Galaxy's IUC tool registry (variant v1g) breaks Claude Haiku and small/MoE locals because the registry's Cheetah templates produce syntactically incomplete CLIs once their runtime parameter bindings are stripped. Implication: hand-authored, command-literal plans are the bottleneck for free-local recipe execution; the model class and the hardware are not.

## 1. Objective

### 1.1 Question

Can a small or local implementer model faithfully convert an Opus-authored recipe into an executable workflow, and how does this depend on recipe specificity?

### 1.2 Why mtDNA variant calling as substrate

The human mitochondrial chromosome (`chrM`, 16,569 bp) is a uniquely well-defined target for an LLM-execution benchmark:

- **Deterministic ground truth.** The dataset was published with a complete canonical workflow in the [Galaxy Training Network mitochondrial-short-variants tutorial](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/mitochondrial-short-variants/tutorial.html). The workflow author and the dataset author overlap with this paper's author, so the answer key is not contested.
- **Small, fast.** Four paired-end Illumina samples, ~838 KB compressed total. Each model's `run.sh` finishes in 8–17 s on either platform once the model emits a valid script.
- **Spans homoplasmic and heteroplasmic variants.** The four samples include both clean (AF ≈ 1.0) and low-frequency (AF down to ~0.04) calls, so an AF-tolerant scoring function is non-trivially exercised.
- **Toolchain is canonical and locked.** BWA-MEM → LoFreq → bcftools/SnpSift, all from a pinned bioconda environment, so model "creativity" in tool selection is constrained.

### 1.3 Why Jaccard on PASS 4-tuples as the primary metric

For each sample, the canonical and model-produced VCFs are compared as sets of 4-tuples `(chrom, pos, ref, alt)` over PASS-filtered records, with a per-call AF tolerance of ±0.02:

> M3 = mean over 4 samples of |M ∩ G| / |M ∪ G|

where M and G are the model's and the ground-truth tuple sets respectively. The metric is set-based, bounded `[0, 1]`, agnostic to caller-specific INFO/FORMAT fields, and averages cleanly across samples without weighting. It does not penalize a model for using a different but valid tool than the recipe specifies (e.g. `bcftools mpileup` instead of `lofreq`); only the *calls* are compared.

### 1.4 Recipe-then-implement framing

Workflow generation is split into two stages: a **plan** (Markdown recipe) authored once by Claude Opus 4.7, and an **implementation** (single self-contained `run.sh`) emitted by the model under test. This decomposition (a) separates planning capability from execution capability — a question that "ask Opus to do everything end-to-end" cannot answer; (b) mirrors how practitioners actually delegate (a senior writes the recipe, a junior writes the script); and (c) admits a clean no-plan control (Track B) that quantifies the marginal value of the plan.

## 2. Results

### 2.1 Plan-detail dominates model capability

![Figure 1. Mean M3 by model × plan variant.](figures/fig1_headline_heatmap.png)

The headline matrix (Figure 1) shows mean M3 across all (model × plan variant) cells, with both hardware platforms pooled. Plan variants are ordered from no-plan (Track B) on the left to most detailed (v2) on the right; intermediate columns v0.5 / v1 / v1g / v1.25 / v1.5 are defined in §4 and Table 2. The v2 column is uniformly green (M3 ≈ 1.0); the no-plan column is uniformly red except for Anthropic frontier models. Most non-frontier rows transition sharply between v1 and v2.

**Anthropic models are insensitive to plan detail on Track A.** Opus 4.7, Sonnet 4.6 and Haiku 4.5 all score M3 = 1.000 ± 0.000 (n = 3) on both v1 and v2 Track A. They have enough internal toolchain knowledge to fill in v1's ambiguities themselves.

**Local open-weight models are dominated by plan detail.** On the RTX 5080, only `qwen3.6:27b` (dense) reaches v2 levels on the lean v1 plan. All other tested open-weight models score M3 ∈ {0.000, 0.333} on v1 and 1.000 ± 0.000 on v2.

**On Jetson AGX Orin in MAXN power mode, 13 of 14 free open-weight models score M3 = 1.000 ± 0.000 on v2 Track A** (Table 1). The single failure (`nemotron-3-nano`, 24 B) is not a plan-detail or budget issue but a code-correctness bug: the model emits a script in 22 s that calls `bgzip results/${sample}.vcf` before `lofreq` has produced the file, exits 255 in 8 s. A second model (`olmo-3.1:32b`) is excluded from the n=14 count above because its three seeds time out at the 900 s wall budget — a real reasoning-mode failure that does not respect Ollama's `/no_think`.

### 2.2 The v1 cliff reduces to one missing command line

![Figure 2. The v1 cliff and its single-line repair.](figures/fig2_v1_cliff_repair.png)

The v1→v2 transition adds ≈1500 bytes of prose and code fences. Two intermediate plan variants isolate which part of that delta carries the load (Figure 2, RTX 5080, Track A, n = 3 per cell):

- **v1.25** = v1 + one extra code-fenced line, the literal `lofreq call-parallel --pp-threads $T -f $REF -o $OUT $BAM` invocation.
- **v1.5** = v2 with every prose paragraph and "Gotchas" subsection stripped — only numbered headings and code fences remain.

For ≥27 B dense models (qwen3.5:27b, qwen3.6:35b-a3b, qwen3-coder:30b, gemma4:26b) **v1.25 alone is sufficient** to bring M3 to 1.000 ± 0.000. The only tool whose CLI these models could not reconstruct from prose is `lofreq`, whose BAM is positional and not flagged. Smaller models (≤14 B, including qwen3:14b, qwen3.5:9b, qwen3:8b, gemma4:e4b) require **every** command line to be literalized; v1.5 brings them to 1.000 while v1.25 leaves them at 0.000–0.333.

The prose paragraphs in v2 — read-group escape-character warnings, in-place bgzip notes, format-string conventions — are not load-bearing for the implementer model. v1.5 omits all of them and produces equivalent M3 to v2 on the models that pass v2 at all. The variants are summarized in Table 2.

### 2.3 No-plan tracks isolate plan value from model capability

![Figure 3. Plan value: Track B vs Track A v1 vs Track A v2.](figures/fig3_plan_value.png)

Track B presents only the problem statement and a tool inventory — no recipe — to the model. Across all tested open-weight models on Track B (any plan-version), mean M3 is 0.000 ± 0.000 with two exceptions: a single 1/3 fluke for `qwen3_14b` (mean 0.17) and an isolated 1/3 for `gemma4:26b` (mean 0.33). Anthropic Opus 4.7 and Sonnet 4.6 reach 0.94 ± 0.00 on Track B; Haiku 4.5 reaches 0.67 ± 0.58 (one of three runs returns zero correct variants, instability without the recipe).

The slope plot (Figure 3) makes the inversion explicit: the open-weight curve is essentially flat at zero through Track B and Track A v1, then jumps to 1.0 at Track A v2. The Anthropic curve is flat near 1.0 across all three points. **The plan is not a productivity hack for local models; it is a capability prerequisite.**

A control (variant v0.5 = Track B + a sequence of tool names with no flags or commands) shows that telling local models *what* to do without specifying *how* is operationally identical to giving them nothing: M3 stays at 0.000 ± 0.000 across all 11 RTX 5080 models tested.

### 2.4 Mechanically extracted plans degrade weak models more than strong ones

![Figure 4. Plan-source robustness: hand-authored v1.25 vs IUC-extracted v1g.](figures/fig4_v1g_robustness.png)

If the v1.25 result reduces the v1 cliff to a single command line, a natural follow-up is whether that command line must be human-authored. Galaxy's IUC tool collection (`galaxyproject/tools-iuc`) is a community-curated registry of XML wrappers, one per tool, with `<command><![CDATA[...]]></command>` blocks that — after Cheetah templating substitution at Galaxy runtime — produce the canonical CLI invocation each tool expects.

We replaced the hand-authored lofreq snippet in v1.25 with a snippet mechanically extracted from `tools-iuc` commit `39e7456` by `scripts/galaxy_to_snippet.py` (Cheetah-strip + macro expansion + path substitution) — call this v1g. Galaxy's lofreq XML emits `--sig $value` and `--bonf $value` where `$value` is supplied at Galaxy runtime; Cheetah strips the variables, so the bare flags `--sig` and `--bonf` end up in the snippet on adjacent lines. lofreq then reads `--sig --bonf` as `--sig=--bonf`, attempts `float("--bonf")`, and exits with `ValueError`.

Figure 4 splits the test pool by self-correction capacity. Claude Opus 4.7, Claude Sonnet 4.6, qwen3.5:27b dense, qwen3.6:27b dense, and qwen3-coder:30b recognize the malformed bare flags and drop them before emitting the script — M3 stays at 1.000. Claude Haiku 4.5, qwen3.6:35b-a3b MoE, gemma4:26b, gemma4:e4b, glm-4.7-flash, and the small qwen3 models copy the snippet character-for-character and ship a script that lofreq refuses to parse — M3 drops to 0.000 ± 0.000. **Haiku regresses from M3 = 1.000 ± 0.000 on v1.25 (hand) to M3 = 0.000 ± 0.000 on v1g (IUC) with no other change.**

The Galaxy-IUC wrappers are Galaxy-runtime templates rather than self-contained CLIs. Without per-tool default-value substitution at extraction time, a deterministic transpiler cannot replace a human plan author for cheaper open-weight models. Strong models can repair noisy registry output; weak models cannot.

### 2.5 Hardware is not the bottleneck

![Figure 5. Hardware comparison.](figures/fig5_hardware.png)

Figure 5 Panel A plots mean generation time (log scale) against mean M3 for each (model × hardware) cell on v2 Track A. The Jetson AGX Orin and the RTX 5080 produce indistinguishable M3 distributions; the only sub-1.0 outliers are model-specific (`gpt-oss:20b` Harmony reasoning eating output budget, `qwen3:8b` minimum-size effects, `nemotron-3-nano`'s scripting bug). The 5080 is faster (typical 5–60 s/seed for the Group A in-VRAM models) than the Jetson (typical 15–300 s/seed), but does not score higher.

Panel B shows the Jetson sweep at 30 W power mode (original 9-of-14 perfect, green) versus MAXN (the 5 originally-failing entries retried, purple). Four of the five MAXN retries reach M3 = 1.000; only `olmo-3.1:32b` remains a real reasoning-mode timeout (excluded from the bar chart because all three of its seeds at MAXN exceed the 900 s budget and produce no run dirs). The Jetson 30 W → MAXN delta is a budget-constrained timeout effect, not a capability one. A latent harness bug — discussed in §3.3 — additionally accounts for three of the original "timeouts" being misroutes rather than real budget failures.

### Table 1 — Headline results

| model | v2 Track A M3 | v1 Track A M3 | Track B M3 | mean v2 gen (s) | M5 pass |
|---|---:|---:|---:|---:|---:|
| Claude Opus 4.7 | 1.00±0.00 (n=3) | 1.00±0.00 (n=3) | 0.94±0.00 (n=3) | 11 | 100% |
| Claude Sonnet 4.6 | 1.00±0.00 (n=3) | 1.00±0.00 (n=3) | 0.94±0.00 (n=3) | 11 | 100% |
| Claude Haiku 4.5 | 1.00±0.00 (n=3) | 1.00±0.00 (n=3) | 0.67±0.58 (n=3) | 69 | 100% |
| `qwen3.6:27b` (dense) | 1.00±0.00 (n=3) | 1.00±0.00 (n=3) | 0.33±0.58 (n=3) | 219 | 100% |
| `qwen3.6:35b` (dense) | 1.00±0.00 (n=3) | 0.00±0.00 (n=3) | — | 100 | 100% |
| `qwen3:32b` | 1.00±0.00 (n=3) | — | — | 286 | 100% |
| `qwen3:14b` | 1.00±0.00 (n=6) | 0.00±0.00 (n=6) | 0.17±0.41 (n=6) | 10 | 33% |
| `qwen3-coder:30b` | 1.00±0.00 (n=3) | 0.33±0.58 (n=3) | 0.00±0.00 (n=3) | 105 | 0% |
| `qwen3.5:27b` | 1.00±0.00 (n=3) | 0.33±0.58 (n=3) | 0.00±0.00 (n=3) | 214 | 100% |
| `qwen3.5:9b` | 1.00±0.00 (n=3) | 0.00±0.00 (n=3) | 0.00±0.00 (n=3) | 6 | 100% |
| `qwen3:8b` | 0.83±0.41 (n=6) | 0.00±0.00 (n=6) | 0.00±0.00 (n=5) | 7 | 0% |
| `qwen3.6:35b-a3b` (MoE) | 1.00±0.00 (n=3) | 0.00±0.00 (n=3) | 0.00±0.00 (n=3) | 127 | 100% |
| `gemma3:27b` | 1.00±0.00 (n=3) | — | — | 259 | 0% |
| `gemma4:26b` | 1.00±0.00 (n=3) | 0.33±0.58 (n=3) | 0.33±0.58 (n=3) | 121 | 100% |
| `gemma4:e4b` | 1.00±0.00 (n=3) | 0.00±0.00 (n=3) | 0.00±0.00 (n=3) | 9 | 100% |
| `mistral-small3.2:24b` | 1.00±0.00 (n=3) | — | — | 170 | 100% |
| `devstral-small-2:24b` | 1.00±0.00 (n=3) | — | — | 163 | 0% |
| `granite-code:34b` | 1.00±0.00 (n=3) | — | — | 203 | 33% |
| `granite4` | 1.00±0.00 (n=3) | — | — | 15 | 100% |
| `deepseek-coder-v2:16b` (MoE) | 1.00±0.00 (n=3) | — | — | 94 | 100% |
| `glm4:9b` | 1.00±0.00 (n=3) | — | — | 55 | 33% |
| `glm-4.7-flash` (MoE) | 1.00±0.00 (n=3) | 0.67±0.58 (n=3) | 0.00±0.00 (n=3) | 90 | 50% |
| `gpt-oss:20b` (MoE) | 0.67±0.58 (n=3) | 0.00±0.00 (n=2) | 0.00 (n=1) | 174 | 83% |
| `llama3.3:70b-instruct-q3_K_M` | 1.00±0.00 (n=3) | — | — | 255 | 100% |
| `nemotron-3-nano` (24 B) | 0.00±0.00 (n=3) | — | — | 22 | 0% |
| `olmo-3.1:32b` | timeout (n=3, ≥900 s) | — | — | — | — |

Mean v2 generation seconds are computed over seeds 42/43/44 of the (model × v2) cell on whichever hardware the model was tested. M5 pass is the fraction of v2 Track A runs satisfying all five script-quality flags (§4). Where a model is untested at a given (plan, track) cell the entry is "—". M3 standard deviations are computed over the n seeds in each cell; cells with `±0.58` reflect the inevitable n=3 std of a {0, 0, 1} or {1, 1, 0} pattern.

### Table 2 — Plan variants

| variant | file | bytes | derivation | hypothesis tested |
|---|---|---:|---|---|
| **Track B** | (no plan) | 0 | problem statement + tool inventory only | how much of the workflow can the model recover from internal knowledge alone? |
| v0.5 | `prompts/track_b_with_order_user.tmpl` | 1361 | Track B + a single line giving the tool order (no flags, no commands) | does sequencing without syntax help local models? |
| v1 (lean) | `plan/PLAN_v1.md` | 3118 | Opus 4.7 from `PLANNER_PROMPT.md`: numbered bullets naming tools and key flags | baseline lean plan |
| v1.25 | `plan/PLAN_v1p25.md` | 3080 | v1 + the exact `lofreq call-parallel` command line | is the v1→v2 cliff explained by a single tool's CLI? |
| v1.5 | `plan/PLAN_v1p5.md` | 1277 | v2 with every prose paragraph and "Gotchas" block deleted, code fences kept | are the prose explanations in v2 load-bearing or decorative? |
| v1g | `plan/PLAN_v1g.md` | 4187 | v1 + Galaxy-IUC-mechanical lofreq snippet (Cheetah-strip + macro expand from `tools-iuc@39e7456`) | can a tool registry replace a human plan author? |
| v2 (detailed) | `plan/PLAN.md` | 4617 | Opus 4.7 from `PLANNER_PROMPT_v2.md`: every step gives the exact command line | reference detailed plan |

### Table 3 — Failure taxonomy on v2 Track A (cells with mean M3 < 1.0)

| model | hardware | n | mean M3 | M1 pass | root cause |
|---|---|---:|---:|---:|---|
| `nemotron-3-nano` | Jetson | 3 | 0.000±0.000 | 0/3 | command-ordering bug: script calls `bgzip results/${sample}.vcf` before lofreq has produced the file; exits 255 in 8 s |
| `gpt-oss:20b` | RTX 5080 | 3 | 0.667±0.577 | varies | Harmony chain-of-thought consumes the 16384-token output budget before the script is emitted (~50% of cells) |
| `qwen3:8b` | RTX 5080 | 6 | 0.833±0.408 | varies | smallest model in the lineup; near the capability floor for v2 (1/6 cells emits an empty script) |
| `olmo-3.1:32b` | Jetson | 3 | timeout | — | reasoning-mode model; 902 s × 3 seeds at MAXN, exceeds 900 s wall budget; does not respect Ollama `/no_think` |

Three of the four failures are independent of the plan-detail axis. Two are reasoning-mode budget exhaustions (one wall-clock, one token-budget); one is a script-correctness bug in code that runs to completion. Only the script-correctness case is what most readers expect a "failure" in this benchmark to look like.

## 3. Methods

### 3.1 Dataset and ground truth

Four paired-end Illumina MiSeq samples (M117-bl, M117-ch, M117C1-bl, M117C1-ch) from [Zenodo 5119008](https://zenodo.org/records/5119008), enriched by long-range PCR for the human mitochondrial chromosome (chrM, 16,569 bp, GRCh38). Total compressed FASTQ size: ~838 KB across 8 files plus 1 reference fasta. Files are md5-verified at fetch time.

Ground truth is produced once by `ground_truth/canonical.sh`, a hand-authored bash workflow that pins to the same locked bioconda environment used by all model runs. Tools and versions are: BWA-MEM (alignment), samtools (sort/index), LoFreq (variant calling), bgzip/tabix (compression/indexing), bcftools (VCF query), SnpSift (annotation in the canonical workflow only — not required for model runs).

### 3.2 Harness

`harness/run_one.py` generates one model run, executes the resulting script in a sandboxed directory, and writes per-run JSON metadata. Scoring (`score/score_run.py`) is a separate step that compares the run's outputs against the ground truth. `harness/sweep_local.py` and `harness/matrix_5080.py` orchestrate the model and plan grids; `score/aggregate.py` produces `results.csv` from the per-run JSONs.

Generation is invoked via the local `claude` CLI for Anthropic models (no `ANTHROPIC_API_KEY` required) or via the Ollama HTTP API for local models. Sampling parameters: temperature = 0.2 for Ollama; Anthropic uses the `claude` CLI default. Ollama `num_predict` = 16384, `seed` ∈ {42, 43, 44}, `think` = `false` for all local runs except the `/think`-comparison cells in §2.5; the rationale is that thinking-mode models routinely exhaust the 900 s wall budget on this workflow's output requirements. Per-call generation budget: 900 s (15 min). Per-script execution budget: 600 s (10 min).

The user prompt is constructed from a fixed system message (`prompts/system.txt`, 1023 bytes), a per-track user template (`prompts/track_a_user.tmpl` or `prompts/track_b_user.tmpl`), the tool inventory, and the plan body. The model emits a single bash script as its raw response; the harness strips a single optional code fence around the script.

### 3.3 Scoring

Five metrics are computed per run by `score/score_run.py`:

- **M1 (Executes).** Binary. `bash run.sh` exits 0 within 600 s wall-clock.
- **M2 (Schema).** Binary. All 9 expected outputs are present (4 BAM, 4 BAI, 4 VCF.gz, 4 TBI, 1 collapsed.tsv) and `bcftools view -h` succeeds on each VCF.
- **M3 (Variant agreement).** Float ∈ [0, 1]. For each of 4 samples, the Jaccard index on the set of (chrom, pos, ref, alt) 4-tuples among PASS-or-unfiltered records, with a per-call AF tolerance of ±0.02. The four per-sample values are macro-averaged. M3 is the **primary metric**; all results sections report M3.
- **M4 (Cost and time).** Reported tuple of input/output token counts, USD cost (Anthropic only; Ollama is $0), wall-clock generation seconds, and wall-clock execution seconds.
- **M5 (Script quality).** Binary, conjunction of: `set -euo pipefail` present; no `/home/anton` paths; `shellcheck` clean (no errors); idempotent (re-run on a populated output dir exits 0).

A run with M1 = 0 returns M2 = 0 and M3 = 0 by construction. A model that uses a different but valid tool than the recipe specifies (e.g. `bcftools mpileup` instead of `lofreq`) is scored on its variant calls, not its tool selection.

### 3.4 Hardware

- **Jetson AGX Orin Developer Kit.** 64 GB unified RAM, Ampere-class GPU (sm_87), aarch64. Two power modes are exercised. The original 9-of-14 sweep ran at **30 W**; after a routing-bug fix in `harness/run_one.py` (provider dispatch was misclassifying any model tag without a `:`), the originally-failing 5 entries were retested at **MAXN**.
- **RTX 5080 desktop.** NVIDIA RTX 5080 16 GB VRAM, 125 GB system RAM, x86_64. Models ≤14 GB fit fully in VRAM (Group A); 17–23 GB models partial-offload to CPU (Group B). `/think` mode was enabled on the Group A models that fit cleanly; on Group B models /think reproducibly hit the 1800 s urlopen timeout and was disabled.

### 3.5 Statistics

Each (model × plan_version × track) cell is run with seeds 42/43/44 (occasionally additional seeds where a retry was performed). All M3 values reported in tables and figures are mean ± standard deviation across seeds, with n indicated. Per-sample Jaccards (4 samples × 3 seeds = 12 values per cell) are used only as the box-plot unit when noted; all other figures aggregate at the seed level. We report standard deviations rather than confidence intervals and avoid significance claims; n = 3 is too small for inferential statistics, but is sufficient to distinguish 0/3 from 3/3 patterns reliably.

## 4. Limitations

1. **Single workflow, single substrate.** We score one task — per-sample variant calling on a 16.6 kb mitochondrial reference — in one canonical toolchain. The recipe-detail effect documented here may not transfer to multi-step ML pipelines, statistical analyses, image processing, or workflows lacking a single canonical toolchain.
2. **Plan author is the strongest model in the test pool.** Claude Opus 4.7 wrote v1, v2, v1.25, v1.5; plan quality is therefore coupled to one model's idiom. v1g (mechanical Galaxy-IUC extraction) is a partial control demonstrating that *non*-Opus plan sources can degrade weak implementer models, but it does not control for the possibility that a different frontier author (e.g. a Gemini- or GPT-authored plan) would produce systematically different cliff structure.
3. **Q4_K_M quantization and `/no_think` are confounders for local models.** All Ollama models run quantized to ~4 bits per parameter, and all run with thinking disabled because of the 900 s wall budget. The reported M3 values therefore measure (Q4-quantized + no-think) capability, not the underlying architectures' ceilings.
4. **n = 3 seeds.** Cell-level standard deviations are informative for distinguishing stable from unstable behaviour but are not confidence intervals; we therefore report std rather than CIs, and avoid significance claims throughout.

## 5. Reproduction

```bash
git clone https://github.com/nekrut/plan-eval && cd plan-eval
bash setup/install.sh                       # miniforge + locked bioconda env (~3 GB, 2-3 min)
bash setup/fetch_data.sh                    # 9 files, ~838 KB, md5-verified
bash ground_truth/canonical.sh              # produces ground_truth/results/

# Anthropic side (uses the local claude CLI; no ANTHROPIC_API_KEY required):
python3 harness/matrix.py --tracks A
python3 score/aggregate.py

# Local-model sweep on Jetson (~3 h wall, $0):
python3 harness/sweep_local.py
python3 score/aggregate.py

# RTX 5080 matrix (11 ollama models × {v1, v2, v1.25, v1.5, v0.5, v1g} × {A, B} × 3 seeds):
python3 harness/matrix_5080.py
python3 harness/matrix_5080.py --plans v1p25,v1p5,v0p5 --include-anthropic

# Regenerate paper figures from results.csv + per-run score.json:
python3 scripts/make_figures.py --all
```

The Anthropic side authenticates through the existing `claude` CLI (Claude Code) login. The local side requires `ollama serve` with the relevant tags pulled.

### Repo layout

```
plan-eval/
├── README.md                  this file
├── results.csv                per-run flat table (n=403)
├── sweep_log.json             most recent sweep summary
├── figures/                   committed paper figures (png)
├── scripts/
│   ├── make_figures.py        regenerates figures/*.png
│   └── galaxy_to_snippet.py   Cheetah-XML → bash extractor (for v1g)
├── plan/                      6 plan variants (Table 2)
├── prompts/                   system + per-track user templates
├── harness/                   run_one, sweep_local, matrix, matrix_5080
├── score/                     score_run, aggregate
├── ground_truth/              canonical workflow + reference VCFs
├── setup/                     miniforge install + locked bioconda env
├── runs/                      Jetson v2 per-run artifacts
├── runs_v1/                   Jetson v1 per-run artifacts
└── runs_5080_v{1,2,1p25,1p5,v0p5,v1g}/   per-experiment 5080 artifacts
```

Each run dir contains `meta.json`, `usage.json`, `exec.json`, `score.json`, the verbatim `run.sh`, and `raw_response.txt`. BAMs and re-derivable VCF artifacts are gitignored.

## 6. Cost

- **Anthropic.** Total spend across the v1 and v2 matrices, the 5080 matrices (with v1.25 / v1.5 / v0.5 / v1g), and three plan-generation calls is under $3 with prompt caching enabled.
- **Local.** $0 (electricity excluded).

## 7. License

MIT (see `LICENSE`). The dataset itself (FASTQ files in Zenodo 5119008) is a separate work by A. Nekrutenko under the original Zenodo terms; this repository does not redistribute it.

## Appendix — Glossary

A reference for readers without a bioinformatics background:

- **Variant identity.** Every variant call is a 4-tuple `(chromosome, position, reference allele, alternate allele)`. Two variants are "the same" only if all four fields match.
- **PASS variants.** Most variant callers tag each call with a confidence filter; we score only calls marked `PASS` or unfiltered.
- **Allele frequency (AF).** Fraction of reads at a position that support the alternate allele. AF = 1.0 means a clean homozygous call; AF = 0.04 means ~4 % of reads carry the variant. Different callers can estimate AF slightly differently for the same position, which is why we apply a ±0.02 tolerance.
- **Jaccard index.** For two sets A and B, J = |A ∩ B| / |A ∪ B|. Bounded [0, 1]; 1.0 means the two sets are identical, 0 means disjoint.
- **chrM.** Human mitochondrial chromosome, 16,569 bp.
- **BWA-MEM.** Alignment of short Illumina reads to a reference genome.
- **LoFreq.** Variant caller emphasizing low-frequency (heteroplasmic) calls; produces a VCF.
- **MoE / A3B.** Mixture-of-Experts model with ~3 B active parameters per token (vs a "dense" model that activates all parameters every token). MoEs are faster per token than a dense model with the same total parameter count.
- **Q4_K_M.** Weight quantization to ~4 bits per parameter; standard for Ollama-hosted models, well-validated for these classes.
- **`/no_think`.** Switches off the chain-of-thought-before-answer behaviour of recent open-weight models (Qwen3 family especially) when supported. We use it for all local model runs because thinking-mode wall-clocks routinely exceed our 900 s budget on this workflow.
- **Track A vs Track B.** A = "with plan", B = "no plan" (problem statement and tool inventory only).
- **Seed.** Random-sampling seed; fixed seeds make the same prompt approximately reproducible. Three seeds per cell give a sense of variance over a single run.
