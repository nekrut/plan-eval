# plan-eval

A small benchmark comparing **Anthropic Claude models** and a **local model on a Jetson AGX Orin** at executing an Opus-authored bioinformatics plan into a runnable workflow.

The premise: *Opus 4.7 plans, a cheaper or local model implements.* This benchmark answers two questions for one concrete task class — per-sample variant calling on human mitochondrial amplicon data:

1. Which Anthropic tier is the cheapest model that can faithfully execute an Opus-authored plan?
2. Is `qwen3.6:35b` (running locally on a Jetson AGX Orin) good enough to replace the cheap-tier Anthropic model for this work?

## TL;DR

| model | track | n | M3 (variant agreement) | $/run | gen secs |
|---|---|---|---|---|---|
| Claude Haiku 4.5 | **A (with plan)** | 3 | **1.000** | $0.063 | 68 |
| Claude Sonnet 4.6 | **A (with plan)** | 3 | **1.000** | $0.049 | 24 |
| Claude Opus 4.7 | **A (with plan)** | 3 | **1.000** | $0.149 | 16 |
| Claude Haiku 4.5 | B (no plan) | 3 | 0.667 ± 0.577 | $0.065 | 64 |
| Claude Sonnet 4.6 | B (no plan) | 3 | 0.938 ± 0.000 | $0.105 | 73 |
| Claude Opus 4.7 | B (no plan) | 3 | 0.938 ± 0.000 | $0.053 | 15 |
| qwen3.6:35b (`/no_think`) | A (with plan) | 3 | **0.000 ± 0.000** | $0.000 | 121 |
| qwen3.6:35b (`/think`) | A (with plan) | — | timed out (>15 min/gen) | — | — |

**Headline:** With the Opus-authored plan, every Anthropic tier scores perfect Jaccard against the canonical variant set. Without it, Haiku is unreliable (one seed scored 0/16 variants). The local model **fails systematically** at this task — across three seeds it produced three *different* invalid `lofreq` invocations, hallucinating CLI flags freshly each time.

**Recommendation:** Use **Claude Sonnet 4.6 + a structured Opus-authored plan** for this class of work. The plan is the lever — it raises Haiku to perfect parity with Opus at one-third the cost, and prevents the high-variance failure mode Haiku exhibits without it. The local model is not a viable substitute on this Jetson for bioinformatics workflows that require obscure CLI surfaces.

## Why this dataset

[Zenodo 5119008](https://zenodo.org/records/5119008) — the *Datasets for Galaxy Collection Operations Tutorial* by A. Nekrutenko. Four paired-end Illumina MiSeq samples (~838 KB compressed), enriched by long-range PCR for human mtDNA, with a known canonical workflow documented in the [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/mitochondrial-short-variants/tutorial.html): BWA-MEM mapping → LoFreq variant calling → SnpSift annotation → collapse.

It is ideal for benchmarking because it's small, tutorial-grade, and the dataset author defined the canonical answer.

## Method

**Plan-then-implement, single-shot, two tracks.** Opus 4.7 generates one structured plan once; that plan is frozen. Each model under test receives the plan plus a problem statement and the locked tool inventory, and must emit a single self-contained `bash run.sh`. Each script runs in an identical fresh sandbox; outputs are scored against a frozen canonical run.

- **Track A** (with plan): the model gets the Opus-authored plan as authoritative.
- **Track B** (no plan): the model gets only the problem statement and tool inventory.

Why not agentic / multi-turn tool use? It conflates "can follow a plan" with "can drive tool use," tripling variance and cost without answering the user's question.

## Scoring

Five metrics per run; M3 is primary.

| metric | computation | type |
|---|---|---|
| **M1** Executes | `bash run.sh` exits 0 within 600 s | binary; gates M2–M3 |
| **M2** Schema | All expected output paths present + `bcftools` header parses | binary |
| **M3** Variant agreement | Per-sample tolerant Jaccard on `(CHROM, POS, REF, ALT)` PASS records, AF tolerance ±0.02; macro-mean across the 4 samples | continuous |
| **M4** Cost / latency | Tokens, USD, gen secs, exec secs | continuous |
| **M5** Code quality | `shellcheck` clean + `set -euo pipefail` + no hardcoded `/home/...` + idempotency check (re-run on populated dir exits 0) | binary |

Tool substitution is handled by construction: M3 compares variant tuples, not pipelines.

## Hardware

- **Jetson AGX Orin Developer Kit**: 64 GB unified RAM, Ampere-class GPU (sm_87), aarch64
- **Power mode**: 30W (MAXN unavailable without a reboot the user declined)
- **Conda env**: locked, on PATH for both canonical and model runs (see `setup/install.sh`)

## What the qwen failures look like

All three Track A `/no_think` seeds at temperature 0.2 produced different broken `lofreq` invocations:

- seed 42: `lofreq call-parallel -f $REF -r $BAM -o ...` — missing `--pp-threads`, invalid `-r`
- seed 43: `lofreq call-parallel -f $REF -i $bam -o $vcf --pp-threads ...` — invalid `-i`, lowercase `$bam`
- seed 44: `lofreq call-parallel -f $REF -d -o $VCF -r 1-16569 $BAM` — invalid `-d`, invalid `-r`

The actual invocation (per the explicit plan) is `lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o out.vcf in.bam`. The model knew the tool name and that it needs `-f` for the reference, but every other detail was hallucinated, and *differently each seed*. That pattern — same plan, three different fabrications at low temperature — indicates this is a knowledge gap, not sampling noise. Bigger budgets or slower modes don't fix knowledge gaps.

`/think` mode never finished within the 15-minute per-call budget. A previous standalone test of `qwen3.6:35b` /think on this Jetson took 6 min wall-clock to emit 2 074 reasoning tokens for a *trivial* pysam function; the larger benchmark prompt exceeded the budget repeatedly.

## Repo layout

```
plan-eval/
├── README.md                    this file
├── LICENSE                      MIT
├── results.csv                  per-run flat table
├── setup/
│   ├── install.sh              miniforge + locked bioconda env
│   ├── verify_env.sh           emits TOOL_INVENTORY string
│   └── fetch_data.sh           md5-verified Zenodo download
├── data/manifest.json           file list with md5s (data not committed)
├── ground_truth/
│   ├── canonical.sh            the answer-key workflow
│   ├── checksums.txt           content-stable VCF hashes
│   └── results/                canonical VCFs (BAMs gitignored)
├── plan/
│   ├── PLANNER_PROMPT.md       what Opus 4.7 saw to author the plan
│   └── PLAN.md                 the frozen plan injected into Track A
├── prompts/                     system + per-track user prompt templates
├── harness/
│   ├── run_one.py              generate + execute one cell
│   └── matrix.py               iterate the matrix
├── score/
│   ├── score_run.py            M1–M5 against ground truth
│   └── aggregate.py            → results.csv + summary
└── runs/<run_id>/               per-run artifacts (run.sh, score.json, exec.log, …)
```

Each `runs/<run_id>/` directory contains the exact `run.sh` the model emitted, plus per-run JSON metadata. BAMs and re-derivable artifacts are gitignored.

## Reproducing

```bash
git clone https://github.com/nekrut/plan-eval && cd plan-eval
bash setup/install.sh                 # miniforge + locked bioconda env (~3 GB, 2-3 min)
bash setup/fetch_data.sh              # 9 files, ~838 KB, md5-verified
bash ground_truth/canonical.sh        # produces ground_truth/results/
python3 harness/run_one.py --model claude-haiku-4-5 --track A --seed 42
python3 score/score_run.py runs/claude-haiku-4-5_track-A_seed-42_*/
```

To run the full matrix (Anthropic plus local):

```bash
python3 harness/matrix.py             # ~5 min Anthropic + 60+ min ollama
python3 score/aggregate.py            # → results.csv + console summary
```

The Anthropic side authenticates through the existing `claude` CLI (Claude Code) login — no `ANTHROPIC_API_KEY` needed. The local side requires `ollama serve` with `qwen3.6:35b` pulled.

## Cost

Total Anthropic spend across the whole 18-cell run plus three plan-generation calls: under **$2** with prompt caching. The `claude -p --max-budget-usd N` flag caps each call.

## Caveats

- **One task class.** This is per-sample variant calling on a 16.6 kb mitochondrial reference. Generalization to whole-genome workflows, RNA-seq, single-cell, etc., is not implied.
- **One dataset.** All four samples are amplicon-PCR mtDNA from the same study. The variants converge on three highly polymorphic positions; the dataset is intentionally easy.
- **Anthropic temperature.** Set by the `claude` CLI default (no `--temperature` flag). Anthropic's documented default is what's tuned for these models; lowering it does not reduce variance the way it does for older models.
- **Ollama temperature.** 0.2 with seeds 42, 43, 44.
- **30W on Jetson.** MAXN power mode requires a reboot on this machine; we did not enable it. Local-model wall-clock figures reflect 30W.
- **No-plan baseline noise.** The Track B `0.938` vs. `1.000` gap on Opus and Sonnet reflects a single-variant disagreement on one of the 4 samples (these models pick a slightly different caller without the plan, recovering most but not all variants at the AF tolerance).

## License

MIT. See `LICENSE`.

The dataset itself (FASTQ files in Zenodo 5119008) is a separate work by A. Nekrutenko under the original Zenodo terms; this repository does not redistribute it.
