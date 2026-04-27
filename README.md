# plan-eval

Benchmark comparing **Anthropic Claude models** and **a wide selection of free local open-weight models** running on a Jetson AGX Orin at executing an Opus-authored bioinformatics plan.

The premise: *Opus 4.7 plans, a cheaper or local model implements.* The benchmark answers:

1. Which Anthropic tier is the cheapest model that can faithfully execute an Opus-authored plan?
2. Which **free local open-weight models** can replace the cheap-tier Anthropic model on this Jetson?
3. How does plan specificity interact with model capability?

## TL;DR

**At least 9 free local open-weight models can implement an Opus-authored hyper-detailed plan on this Jetson with perfect ground-truth Jaccard.** All-Anthropic models score perfectly with the v2 plan. With the looser v1 plan, only Anthropic models worked — all open-weight models tested produced broken `lofreq` invocations.

The local-model success depends on plan specificity: hyper-detailed plans (every command and flag verbatim) work; lean plans don't.

## Headline tables

### Free local open-weight sweep (v2 detailed plan, /no_think where applicable, n=3 seeds)

| model | M3 | mean gen secs | result | notes |
|---|---|---|---|---|
| `qwen3-coder:30b` | **1.000** | 76 | 3/3 | Qwen3 MoE A3B coder — fastest perfect |
| `glm4:9b` | **1.000** | 69 | 3/3 | small/fast Z.ai control |
| `deepseek-coder-v2:16b` | **1.000** | 108 | 3/3 | MoE A2.4B |
| `qwen3.6:35b /no_think` | **1.000** | 100 | 3/3 | (v2 baseline) |
| `devstral-small-2:24b` | **1.000** | 177 | 3/3 | Mistral coder |
| `mistral-small3.2:24b` | **1.000** | 184 | 3/3 | Mistral generalist |
| `granite-code:34b` | **1.000** | 217 | 3/3 | IBM coder |
| `gemma3:27b` | **1.000** | 273 | 3/3 | Google generalist |
| `qwen3:32b` | **1.000** | 300 | 3/3 | dense Qwen3 (vs MoE qwen3.6) |
| `gpt-oss:20b` | 1.000 | 556 | 1/3 | 2 seeds hit 900 s gen timeout |
| `glm-4.7-flash` | — | — | 0/3 | gen timeout (looks like reasoning model) |
| `olmo-3.1:32b` | — | — | 0/3 | gen timeout |
| `nemotron-3-nano` (24 B) | — | — | 0/3 | gen timeout |
| `granite4` (3 B) | — | — | 0/3 | gen timeout |
| `llama3.3:70b-instruct-q3_K_M` | — | — | — | network failure mid-pull |

### Anthropic (Track A, both plan versions)

| plan | model | M3 | $/run | gen secs |
|---|---|---|---|---|
| v2 (detailed) | Claude Opus 4.7 | 1.000 | $0.144 | 11 |
| v2 (detailed) | Claude Sonnet 4.6 | 1.000 | $0.064 | 11 |
| v2 (detailed) | Claude Haiku 4.5 | 1.000 | $0.077 | 69 |
| v1 (lean) | Claude Opus 4.7 | 1.000 | $0.149 | 16 |
| v1 (lean) | Claude Sonnet 4.6 | 1.000 | $0.049 | 24 |
| v1 (lean) | Claude Haiku 4.5 | 1.000 | $0.063 | 68 |

### No-plan baseline (v1 only — Track B is plan-independent)

| model | M3 | notes |
|---|---|---|
| Opus 4.7 | 0.938 | one-variant disagreement |
| Sonnet 4.6 | 0.938 | same |
| Haiku 4.5 | 0.667 ± 0.577 | one seed scored 0/16 — high variance without plan |
| qwen3.6:35b /no_think (v1) | 0.000 | hallucinated lofreq three different ways across seeds |

## Findings

**Plan specificity is the dominant lever.** Going from v1 (lean, ~1 200 tokens) to v2 (~2 274 tokens, every command and flag verbatim):

- Anthropic models: no change. 1.000 either way.
- Local open-weight models: **0.000 → 1.000**. The v1 plan's "use lofreq with --pp-threads 4" left enough ambiguity that qwen3.6:35b hallucinated three different invalid invocations. The v2 plan supplies the literal `lofreq call-parallel --pp-threads 4 -f data/ref/chrM.fa -o results/{sample}.vcf results/{sample}.bam`, and every working open-weight model emits the correct invocation.

**Model size is not the dividing line.** `glm4:9b` (5.5 GB on disk) and `qwen3-coder:30b` (18.6 GB) both score 1.000 with the v2 plan. The dividing line, with one exception, is whether the model finishes within the 900 s timeout — and that correlates with thinking/reasoning behavior, not parameter count.

**The four 0/3 failures are reasoning-model timeouts, not capability gaps.** `glm-4.7-flash`, `olmo-3.1:32b`, `nemotron-3-nano`, and `granite4` all crashed `urllib.request.urlopen` at our 900 s generation budget. Combined with what we already knew about `qwen3.6:35b /think` (took >>900 s on a trivial pysam prompt), this is consistent with these models emitting long internal chain-of-thought before the answer despite `think:false` in the payload — these are likely models whose reasoning isn't fully gateable from Ollama. A 1 800 s budget would probably recover them, at the cost of 30 min/cell. We left this open.

**`gpt-oss:20b` behaves like a reasoning model too.** One seed completed at 569 s; two seeds exceeded the 900 s budget. So it's borderline — a longer budget would likely recover the rest.

## Recommendation

For routine bioinformatics-plan execution where the plan author is Opus 4.7:

- **Cheapest paid (with a moderately detailed plan)**: **Sonnet 4.6 + a v1-style plan** ($0.049/run, ~10 s gen, 1.000 on both v1 and v2).
- **Cheapest free local on this Jetson (with a hyper-detailed plan)**: **`qwen3-coder:30b /no_think`** at $0/run, ~76 s gen, 1.000 across all seeds. Other comparably good free options: `glm4:9b` (faster but smaller), `deepseek-coder-v2:16b` (fast MoE), `devstral-small-2:24b`, `mistral-small3.2:24b`, `granite-code:34b`, `gemma3:27b`, `qwen3:32b`, `qwen3.6:35b`.

For free-local execution to be reliable, the plan must be **hyper-detailed**: every command line, every flag, every filename verbatim. The plan author has to do the bioinformatics knowing; the local model does the bash transliteration.

The practical workflow: **spend $0.14 once on an Opus-authored hyper-detailed plan, then implement for free locally** at ~1–5 minutes per cell on this Jetson at 30 W.

## Why this dataset

[Zenodo 5119008](https://zenodo.org/records/5119008) — *Datasets for Galaxy Collection Operations Tutorial* by A. Nekrutenko. Four paired-end Illumina MiSeq samples (~838 KB compressed), enriched by long-range PCR for human mtDNA, with a known canonical workflow documented in the [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/mitochondrial-short-variants/tutorial.html): BWA-MEM mapping → LoFreq variant calling → SnpSift annotation → collapse.

Tutorial-grade, the dataset author defined the canonical answer.

## Method

**Plan-then-implement, single-shot.** Opus 4.7 generates one structured plan once; that plan is frozen. Each model under test receives the plan plus a problem statement and the locked tool inventory, and must emit a single self-contained `bash run.sh`. Each script runs in an identical fresh sandbox; outputs are scored against a frozen canonical run.

Two plan versions:
- **v1 (lean)**: numbered list of bullet-pointed steps with named tools and key flags but no command-by-command syntax. ~1 200 output tokens. See `plan/PLAN_v1.md`.
- **v2 (detailed)**: every step gives the exact command line as a code block with all flags and arguments. ~2 274 output tokens. See `plan/PLAN.md`.

Two tracks:
- **Track A** (with plan): the model gets the plan as authoritative. The local-model sweep is Track A only.
- **Track B** (no plan): the model gets only the problem statement and tool inventory. Run for v1 Anthropic only.

Local-model sweep procedure (`harness/sweep_local.py`):
1. `ollama pull <model>`
2. For seed in {42, 43, 44}: `run_one.py --model <model> --track A --think off`. The 900 s `urlopen` budget for the generation step is hard-coded in `run_one.py`; models that exceed it appear here as "0/3 timeout."
3. Score each run.
4. `ollama rm <model>` — free disk before the next model.

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
- **Power mode**: 30 W (MAXN unavailable without a reboot the user declined)
- **Conda env**: locked, on PATH for both canonical and model runs (see `setup/install.sh`)

## Repo layout

```
plan-eval/
├── README.md                     this file
├── LICENSE                       MIT
├── results.csv                   per-run flat table; plan_version column
├── sweep_log.json                local-sweep summary (per model, per seed)
├── setup/                        miniforge + locked bioconda env, Zenodo fetch
├── data/manifest.json            file list with md5s (data not committed)
├── ground_truth/                 canonical answer-key workflow + VCFs
├── plan/
│   ├── PLANNER_PROMPT.md         v1 planner prompt
│   ├── PLANNER_PROMPT_v2.md      v2 planner prompt (hyper-detailed)
│   ├── PLAN.md                   current frozen plan (v2)
│   └── PLAN_v1.md                preserved v1 plan
├── prompts/                      system + per-track user prompt templates
├── harness/
│   ├── run_one.py                generate + execute one cell
│   ├── matrix.py                 Anthropic + ollama matrix iterator
│   └── sweep_local.py            disk-rotated local-model sweep
├── score/
│   ├── score_run.py              M1–M5 against ground truth
│   └── aggregate.py              scans runs/ + runs_v*/ → results.csv
├── runs/<run_id>/                v2 per-run artifacts
└── runs_v1/<run_id>/             v1 per-run artifacts (preserved)
```

Each run dir contains the exact `run.sh` the model emitted plus per-run JSON metadata (`meta.json`, `usage.json`, `exec.json`, `score.json`, `raw_response.txt`, `exec.log`). BAMs and re-derivable artifacts are gitignored.

## Reproducing

```bash
git clone https://github.com/nekrut/plan-eval && cd plan-eval
bash setup/install.sh                 # miniforge + locked bioconda env (~3 GB, 2-3 min)
bash setup/fetch_data.sh              # 9 files, ~838 KB, md5-verified
bash ground_truth/canonical.sh        # produces ground_truth/results/

# Anthropic side (uses claude CLI; no API key required):
python3 harness/matrix.py --tracks A
python3 score/aggregate.py

# Local-model sweep (~3 h wall, $0):
python3 harness/sweep_local.py
python3 score/aggregate.py
```

## Cost

- Total Anthropic spend across both v1 and v2 matrices plus three plan-generation calls: under **$3** with prompt caching.
- Local sweep: **$0** (electricity excluded).

## Caveats

- **One task class.** Per-sample variant calling on a 16.6 kb mitochondrial reference. Generalization to whole-genome workflows, RNA-seq, single-cell, etc., not implied.
- **One dataset.** All four samples are amplicon-PCR mtDNA from one study; intentionally easy.
- **30 W on Jetson.** MAXN power mode requires a reboot we did not perform.
- **Generation timeout 900 s** in the harness. Four open-weight models hit this on every seed; we did not retry with a higher budget. Their score of 0/3 is a "doesn't fit our 15-minute-per-cell budget" finding, not necessarily a "model can't do the task" finding.
- **The v2 plan is close to a script in prose.** This is the *whole point* — the experimental finding is that this is what you need to drop into to get a free local model working.
- **Disk pressure during the sweep.** We pull → test → `ollama rm` to keep ~17 GB free. The full 14-model sweep transferred ~210 GB sequentially.

## License

MIT. See `LICENSE`.

The dataset itself (FASTQ files in Zenodo 5119008) is a separate work by A. Nekrutenko under the original Zenodo terms; this repository does not redistribute it.
