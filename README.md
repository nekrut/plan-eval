# plan-eval

A benchmark comparing **Anthropic Claude models** and a **local model on a Jetson AGX Orin** at executing an Opus-authored bioinformatics plan into a runnable workflow — and what happens when the plan goes from *moderately detailed* (v1) to *hyper-detailed* (v2, every command and flag specified verbatim).

The premise: *Opus 4.7 plans, a cheaper or local model implements.* Two questions:

1. Which Anthropic tier is the cheapest model that can faithfully execute an Opus-authored plan?
2. Can `qwen3.6:35b` running locally on a Jetson AGX Orin replace the cheap-tier Anthropic model?

## TL;DR

The local model goes from **0/3 (v1)** to **3/3 (v2)** when the plan specifies the exact `lofreq call-parallel ...` invocation rather than just "use lofreq with --pp-threads."

| plan | model | track | n | M3 | $/run | gen secs |
|---|---|---|---|---|---|---|
| **v2 (detailed)** | qwen3.6:35b `/no_think` | A | 3 | **1.000** | $0.000 | 100 |
| **v2 (detailed)** | Claude Haiku 4.5 | A | 3 | 1.000 | $0.077 | 69 |
| **v2 (detailed)** | Claude Sonnet 4.6 | A | 3 | 1.000 | $0.064 | 11 |
| **v2 (detailed)** | Claude Opus 4.7 | A | 3 | 1.000 | $0.144 | 11 |
| v1 (lean plan) | qwen3.6:35b `/no_think` | A | 3 | **0.000** | $0.000 | 121 |
| v1 (lean plan) | Claude Haiku 4.5 | A | 3 | 1.000 | $0.063 | 68 |
| v1 (lean plan) | Claude Sonnet 4.6 | A | 3 | 1.000 | $0.049 | 24 |
| v1 (lean plan) | Claude Opus 4.7 | A | 3 | 1.000 | $0.149 | 16 |
| v1 (lean plan) | Claude Haiku 4.5 | B (no plan) | 3 | 0.667 ± 0.577 | $0.065 | 64 |
| v1 (lean plan) | Claude Sonnet 4.6 | B (no plan) | 3 | 0.938 | $0.105 | 73 |
| v1 (lean plan) | Claude Opus 4.7 | B (no plan) | 3 | 0.938 | $0.053 | 15 |

## Findings

**Plan specificity is the lever for the local model.** The v1 plan said things like "Use `lofreq call-parallel` with `--pp-threads 4` and `-f data/ref/chrM.fa`." That was enough for every Anthropic tier; qwen3.6:35b hallucinated three different invalid invocations across three seeds:

- `lofreq call-parallel -f $REF -r $BAM -o ...`         (missing `--pp-threads`, invalid `-r`)
- `lofreq call-parallel -f $REF -i $bam -o $vcf --pp-threads ...`  (invalid `-i`, lowercase `$bam`)
- `lofreq call-parallel -f $REF -d -o $VCF -r 1-16569 $BAM`  (invalid `-d`, invalid `-r`)

The v2 plan specifies the literal command: `lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/{sample}.vcf results/{sample}.bam`. With that, all three qwen seeds emitted byte-identical (modulo variable-quote style) correct invocations and scored perfect Jaccard.

**Anthropic models are insensitive to plan detail in this range.** Going from v1 to v2 didn't move the needle — every Anthropic model scored 1.000 on both. They have enough internal knowledge of lofreq's CLI surface to fill in v1's gaps; qwen3.6:35b does not.

**The plan-detail spectrum has a tradeoff.** A v1-style plan is shorter (~1 200 tokens output) and leaves room for the implementer's judgment; a v2-style plan (~2 274 output tokens) is essentially the script in prose form. With v2, the implementer's value-add shrinks toward "transliterate prose into bash." For local-model execution that's actually the point — you trade the plan-author's tokens (cheap, one-shot) for the implementer's reliability.

**No-plan baseline (Track B, v1 only).** Without any plan, Opus and Sonnet recover 0.938 (one variant difference per sample) and Haiku is unreliable (0.667 mean, σ 0.577 — one seed scored 0/16 variants). Confirming: the plan is the lever. We did not re-run Track B for v2 because Track B is plan-independent by construction.

**qwen3.6:35b `/think` mode is impractical on this Jetson.** Every `/think` attempt timed out at the 15-minute per-call budget. A standalone test of a trivial pysam function in `/think` mode took ~6 minutes wall-clock to emit ~2 000 reasoning tokens. The benchmark prompt is much larger; the resulting reasoning would exceed any reasonable budget. We did not retest `/think` for v2 — this is a hardware-class problem (Ampere GPU, 30W power), not something a more detailed plan can fix.

## Recommendation

For routine bioinformatics-plan execution where the plan author is Opus 4.7:

- **If you have an Anthropic budget**, use **Sonnet 4.6 + a moderately detailed plan** ($0.049–0.064/run, ~10 s gen). It hits 1.000 on both v1 and v2, no hyper-specification required.
- **If you want $0/run local execution on this Jetson**, use **qwen3.6:35b `/no_think` + a hyper-detailed plan** that specifies every command, flag, and filename verbatim. ~100 s gen + 10 s exec at 30W. The plan author has to do the bioinformatics knowing; the local model does the bash transliteration.

The v2 plan-authoring call cost $0.14 (Opus 4.7) and was reusable across all model runs. So the practical workflow is "spend $0.14 once on a hyper-detailed plan, then implement for free locally."

## Why this dataset

[Zenodo 5119008](https://zenodo.org/records/5119008) — *Datasets for Galaxy Collection Operations Tutorial* by A. Nekrutenko. Four paired-end Illumina MiSeq samples (~838 KB compressed), enriched by long-range PCR for human mtDNA, with a known canonical workflow documented in the [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/mitochondrial-short-variants/tutorial.html): BWA-MEM mapping → LoFreq variant calling → SnpSift annotation → collapse.

Tutorial-grade, the dataset author defined the canonical answer. Ideal for benchmarking.

## Method

**Plan-then-implement, single-shot.** Opus 4.7 generates one structured plan once; that plan is frozen. Each model under test receives the plan plus a problem statement and the locked tool inventory, and must emit a single self-contained `bash run.sh`. Each script runs in an identical fresh sandbox; outputs are scored against a frozen canonical run.

Two plan versions:

- **v1 (lean)**: a numbered list of bullet-pointed steps with named tools and key flags but no command-by-command syntax. ~1 200 output tokens. See `plan/PLAN_v1.md`.
- **v2 (detailed)**: every step gives the exact command line as a code block with all flags and arguments. ~2 274 output tokens. See `plan/PLAN.md`.

Two tracks:

- **Track A** (with plan): the model gets the plan as authoritative.
- **Track B** (no plan): the model gets only the problem statement and tool inventory. Run for v1 only since plan-version is irrelevant.

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

## Repo layout

```
plan-eval/
├── README.md                     this file
├── LICENSE                       MIT
├── results.csv                   per-run flat table; plan_version column
├── setup/
│   ├── install.sh               miniforge + locked bioconda env
│   ├── verify_env.sh            emits TOOL_INVENTORY string
│   └── fetch_data.sh            md5-verified Zenodo download
├── data/manifest.json            file list with md5s (data not committed)
├── ground_truth/
│   ├── canonical.sh             the answer-key workflow
│   ├── checksums.txt            content-stable VCF hashes
│   └── results/                 canonical VCFs (BAMs gitignored)
├── plan/
│   ├── PLANNER_PROMPT.md        v1 planner prompt
│   ├── PLANNER_PROMPT_v2.md     v2 planner prompt (hyper-detailed)
│   ├── PLAN.md                  current frozen plan (v2)
│   └── PLAN_v1.md               preserved v1 plan
├── prompts/                      system + per-track user prompt templates
├── harness/
│   ├── run_one.py               generate + execute one cell
│   └── matrix.py                iterate the matrix; --tracks A and --no-think flags
├── score/
│   ├── score_run.py             M1–M5 against ground truth
│   └── aggregate.py             scans runs/ + runs_v*/ → results.csv
├── runs/<run_id>/                v2 per-run artifacts
└── runs_v1/<run_id>/             v1 per-run artifacts (preserved)
```

Each `runs/<run_id>/` directory contains the exact `run.sh` the model emitted, plus per-run JSON metadata (`meta.json`, `usage.json`, `exec.json`, `score.json`, `raw_response.txt`, `exec.log`). BAMs and re-derivable artifacts are gitignored.

## Reproducing

```bash
git clone https://github.com/nekrut/plan-eval && cd plan-eval
bash setup/install.sh                 # miniforge + locked bioconda env (~3 GB, 2-3 min)
bash setup/fetch_data.sh              # 9 files, ~838 KB, md5-verified
bash ground_truth/canonical.sh        # produces ground_truth/results/
python3 harness/run_one.py --model claude-haiku-4-5 --track A --seed 42
python3 score/score_run.py runs/claude-haiku-4-5_track-A_seed-42_*/
```

Full matrix:

```bash
python3 harness/matrix.py --tracks A             # both Anthropic and ollama, Track A only
python3 harness/matrix.py --only-ollama --tracks A --no-think
python3 score/aggregate.py                       # → results.csv + console summary
```

The Anthropic side authenticates through the existing `claude` CLI (Claude Code) login — no `ANTHROPIC_API_KEY` needed. The local side requires `ollama serve` with `qwen3.6:35b` pulled.

## Cost

Total Anthropic spend across both v1 and v2 matrices plus three plan-generation calls: under **$3** with prompt caching. The `claude -p --max-budget-usd N` flag caps each call.

## Caveats

- **One task class.** This is per-sample variant calling on a 16.6 kb mitochondrial reference with a well-known canonical workflow. Generalization to whole-genome workflows, RNA-seq, single-cell, etc., is not implied.
- **One dataset.** All four samples are amplicon-PCR mtDNA from the same study. The variants converge on three highly polymorphic positions; the dataset is intentionally easy.
- **Anthropic temperature.** Set by the `claude` CLI default (no `--temperature` flag).
- **Ollama temperature.** 0.2 with seeds 42, 43, 44.
- **30W on Jetson.** MAXN power mode requires a reboot on this machine; we did not enable it. Local-model wall-clock figures reflect 30W.
- **No-plan baseline noise.** The v1 Track B `0.938` vs `1.000` gap on Opus and Sonnet reflects a single-variant disagreement on one of the 4 samples (these models pick a slightly different caller without the plan, recovering most but not all variants at the AF tolerance).
- **The v2 plan is close to a script in prose.** This is the *whole point* — the experimental finding is that this is what you need to drop into to get a local model working.

## License

MIT. See `LICENSE`.

The dataset itself (FASTQ files in Zenodo 5119008) is a separate work by A. Nekrutenko under the original Zenodo terms; this repository does not redistribute it.
