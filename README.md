# plan-eval

Benchmark comparing **Anthropic Claude models** and **a wide selection of free local open-weight models** running on a Jetson AGX Orin at executing an Opus-authored bioinformatics plan.

## What this is, in plain English

You ask the strongest available model (Claude Opus 4.7) to write a recipe for a bioinformatics workflow. Then you hand that recipe to a cheaper or smaller model and ask it to turn the recipe into an actual runnable script. We then run the script on real data and check whether its variant calls match a known-correct answer.

The benchmark answers three questions:

1. Which Anthropic tier is the cheapest model that can faithfully execute an Opus-authored recipe?
2. Which **free, locally-runnable** models can do the same job — i.e. avoid paying for inference at all?
3. How does the *level of detail in the recipe* interact with model capability?

## TL;DR

**At least 9 free local open-weight models can implement an Opus-authored hyper-detailed recipe on this Jetson with a perfect match against the known-correct variant set.** Anthropic models match on both the lean and the detailed recipe. Free local models only work when the recipe is hyper-detailed; with a lean recipe they all fail.

## Plain-English glossary

A few terms used throughout this README and in the source files:

- **Recipe / plan** — a Markdown document describing the workflow steps. We test two versions: **v1 (lean)** = numbered bullets naming the tools and key flags, leaving the implementer some judgment; **v2 (detailed)** = every command line spelled out verbatim.
- **Implementer model** — the model under test. It receives the recipe and emits a single bash script (`run.sh`) that performs the workflow. We then *execute* the script on the data and *score* its outputs against ground truth.
- **Track A vs Track B** — A = "with the recipe", B = "no recipe, problem statement only" (a control to measure what the recipe is worth).
- **Seed** — a random number that controls the model's sampling. Same prompt + same seed ≈ same output. We run **three seeds (42, 43, 44)** for each cell so a one-off lucky or unlucky run doesn't dominate the result. *n=3* below means three runs.
- **`/think` vs `/no_think`** — some recent open-weight models (notably Qwen3) can *think out loud* before answering — emit pages of internal reasoning before the actual answer. `/no_think` turns this off. Thinking mode helps quality but is much slower; on this Jetson it routinely doesn't finish in 15 minutes per call, so we run all local models in `/no_think` (or its equivalent for non-Qwen families: setting `think: false` in the Ollama API payload).
- **MoE / A3B** — Mixture-of-Experts model with ~3B active parameters per token (vs. a "dense" model that activates all its parameters every token). MoEs are faster per token than a dense model with the same total parameter count; some of the strongest small/local options here are MoEs.
- **Q4_K_M quantization** — model weights compressed to ~4 bits per parameter to fit in memory. Some quality loss, well-validated for these model classes.
- **M3 (variant agreement)** — our primary score: across 4 samples, how many variants does the model's output share with the known-correct VCF, with a small allele-frequency tolerance. **1.000 = perfect**, 0.000 = no overlap. See "Scoring" below for the full rubric.
- **Variant calling, BWA-MEM, LoFreq, BAM/VCF, chrM** — bioinformatics: the workflow aligns short DNA reads to a reference genome (`BWA-MEM`), produces an alignment file (`BAM`), then identifies the differences from the reference (`LoFreq` → `VCF`). `chrM` = the human mitochondrial chromosome (16 569 bp).

## Headline tables

### Free local open-weight sweep
*v2 (detailed) recipe, thinking turned off, three independent runs (seeds 42/43/44) per model.*

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

### Anthropic models with the recipe (Track A, both recipe versions)

| recipe | model | M3 | $/run | gen secs |
|---|---|---|---|---|
| v2 (detailed) | Claude Opus 4.7 | 1.000 | $0.144 | 11 |
| v2 (detailed) | Claude Sonnet 4.6 | 1.000 | $0.064 | 11 |
| v2 (detailed) | Claude Haiku 4.5 | 1.000 | $0.077 | 69 |
| v1 (lean) | Claude Opus 4.7 | 1.000 | $0.149 | 16 |
| v1 (lean) | Claude Sonnet 4.6 | 1.000 | $0.049 | 24 |
| v1 (lean) | Claude Haiku 4.5 | 1.000 | $0.063 | 68 |

### What happens with NO recipe (Track B, v1 only)
*Each model is given only the problem statement and the list of available tools — no recipe.*

| model | M3 | notes |
|---|---|---|
| Opus 4.7 | 0.938 | one-variant disagreement (15 of 16 variants match) |
| Sonnet 4.6 | 0.938 | same |
| Haiku 4.5 | 0.667 ± 0.577 | one of three runs found zero correct variants — unstable without the recipe |
| qwen3.6:35b (free local) | 0.000 | invented three *different* invalid `lofreq` command lines across the three runs |

## Findings

**The level of detail in the recipe is the dominant lever.** Going from v1 (lean, ~1 200 tokens) to v2 (~2 274 tokens, every command and flag verbatim):

- Anthropic models: unchanged at 1.000. They have enough internal knowledge to fill in the v1 ambiguities themselves.
- Free local models: **0.000 → 1.000**. With the lean recipe, free local models invent broken command lines for tools they don't know well (in particular `lofreq`). With the verbatim recipe, they have nothing to invent — they just transliterate prose into bash.

**Model size is not the dividing line.** `glm4:9b` (5.5 GB on disk) and `qwen3-coder:30b` (18.6 GB) both score 1.000 with the v2 recipe. What divides "works" from "didn't finish" is whether the model can produce its answer within the 15-minute-per-call generation budget our test harness allows.

**The four 0/3 failures look like models that *want to think* but couldn't be turned off cleanly.** `glm-4.7-flash`, `olmo-3.1:32b`, `nemotron-3-nano`, and `granite4` all timed out at 15 minutes for every run. This is consistent with reasoning-style models that emit a long internal chain-of-thought before the actual answer — and don't fully respect Ollama's "no thinking" switch. A longer budget would likely recover them at ~30 min per call. We didn't retry; that's listed as a follow-up.

**`gpt-oss:20b` is borderline.** One run finished in ~9.5 minutes, two timed out at 15 — same family of issue.

## Recommendation

If Opus 4.7 is your recipe author, here are the practical answers:

- **Cheapest paid path** (recipe doesn't have to be elaborate): **Sonnet 4.6 + a lean recipe**. ~$0.05 per implementation, ~10 seconds, perfect score on every run.
- **Cheapest free path on this Jetson** (recipe must be hyper-detailed): **`qwen3-coder:30b`** in `/no_think` mode. $0 per implementation, ~75 seconds on this hardware, perfect score on every run. Alternatives that work equally well: `glm4:9b` (smaller and faster), `deepseek-coder-v2:16b` (fast Mixture-of-Experts), `devstral-small-2:24b`, `mistral-small3.2:24b`, `granite-code:34b`, `gemma3:27b`, `qwen3:32b`, `qwen3.6:35b`.

The practical workflow looks like this:

> Spend roughly **$0.14 once** on an Opus-authored hyper-detailed recipe.
> Then run that recipe through a free local model on this Jetson at $0 per execution.
> Wall-clock: 1 – 5 minutes per implementation depending on model size.

To make free-local execution reliable, the recipe has to specify *everything*: every command, every flag, every filename, character-for-character. The recipe author has to know bioinformatics; the local model just translates the prose into bash.

## Why this dataset

[Zenodo 5119008](https://zenodo.org/records/5119008) — *Datasets for Galaxy Collection Operations Tutorial* by A. Nekrutenko. Four paired-end Illumina MiSeq samples (~838 KB compressed), enriched by long-range PCR for human mtDNA, with a known canonical workflow documented in the [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/mitochondrial-short-variants/tutorial.html): BWA-MEM mapping → LoFreq variant calling → SnpSift annotation → collapse.

Tutorial-grade, the dataset author defined the canonical answer.

## Method

**Recipe-then-implement, one shot per run.** Opus 4.7 writes one recipe; that recipe is then frozen and shown to every implementer model. Each implementer must reply with a single self-contained `bash run.sh` script. The script is then executed in a fresh sandbox directory, and its outputs are compared to a known-correct reference (the "ground truth" — produced once by hand using the same locked toolchain).

Two recipe versions, written by Opus 4.7 from two different planner prompts:
- **v1 (lean)**: numbered bullets naming the tools and key flags, but no full command lines. ~1 200 output tokens. See `plan/PLAN_v1.md`.
- **v2 (detailed)**: every step gives the exact command line, flags and all, in a code fence. ~2 274 output tokens. See `plan/PLAN.md`.

Two tracks:
- **Track A** (with recipe): the model gets the recipe as authoritative. The local-model sweep is Track A only.
- **Track B** (no recipe): the model gets only the problem statement and the available-tools list — a control to measure what the recipe is worth. Run for v1 Anthropic only.

Local-model sweep procedure (`harness/sweep_local.py`):
1. `ollama pull <model>` — download.
2. Run the model three times, once per seed (42, 43, 44). Each run is a separate API call, completely independent — same prompt, different random sampling, so we get a sense of variance rather than a single lucky/unlucky data point.
3. Each call has a 15-minute generation budget. Models that don't reply within that time count as a 0/3 failure (we report those as "timed out" rather than guessing they would have been correct given more time).
4. Score each run against ground truth.
5. `ollama rm <model>` to free disk before pulling the next model. (We have ~17 GB free; the larger models are 30+ GB on disk, so we have to rotate.)

## Scoring

Five metrics computed per run. **M3 (variant agreement) is the primary one — that's what the headline tables report.**

| metric | what it measures | how it's computed |
|---|---|---|
| **M1 — Does it run?** | The script ran end-to-end without crashing within 10 minutes. | `bash run.sh` exits 0 within 600 s. Pass/fail. M1 must pass for M2 and M3 to be computed. |
| **M2 — Did it produce the expected files?** | The output structure matches the spec (right files in the right places, valid VCF headers). | Filesystem check + `bcftools view -h` succeeds on each VCF. Pass/fail. |
| **M3 — Are the variant calls right?** | For each of the 4 samples, what fraction of the variants in the model's VCF *also* appear in the known-correct VCF, allowing a small allele-frequency tolerance. Averaged across samples. **1.000 = perfect, 0.000 = no overlap.** | Per-sample Jaccard on `(chromosome, position, reference allele, alternate allele)` tuples among PASS records, with allele-frequency tolerance ±0.02; macro-mean across the 4 samples. |
| **M4 — Cost and time** | What it cost and how long it took. | Tokens, USD (Anthropic), generation seconds, execution seconds. |
| **M5 — Is the script clean?** | The script is well-formed bash and re-runnable. | `shellcheck` clean + `set -euo pipefail` present + no hardcoded user paths + re-running on a populated output dir exits 0 with no work performed. Pass/fail. |

If a model uses a *different but valid* tool than the recipe specifies (e.g. `bcftools mpileup` instead of `lofreq`), M3 still scores it correctly — the metric compares variant calls, not the pipeline that produced them.

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

- **One task class.** Per-sample variant calling on a 16.6 kb mitochondrial reference. Generalization to whole-genome workflows, RNA-seq, single-cell, etc. is not implied.
- **One dataset.** All four samples are amplicon-PCR mitochondrial DNA from one study; intentionally easy.
- **30 W power mode on the Jetson.** A higher power mode (MAXN) would require rebooting the machine, which we didn't do. So local-model wall-clock numbers are conservative — they'd be faster on a more aggressive power profile.
- **15-minute generation budget per call.** Four free-local models timed out at this limit on every run. Their 0/3 score means "did not fit our 15-minute budget on this Jetson," not "the model can't do the task in principle." A 30-minute budget would likely recover most of them.
- **The v2 recipe is close to "the script in prose."** That's the whole point — the experimental finding is that this is what you need to write to get a free local model working reliably. If you want the local model to *also* do bioinformatics reasoning, more work is needed.
- **Disk rotation during the sweep.** We pull → test → remove for each model to fit in ~17 GB free disk. The full 14-model sweep transferred ~210 GB sequentially over the network.

## License

MIT. See `LICENSE`.

The dataset itself (FASTQ files in Zenodo 5119008) is a separate work by A. Nekrutenko under the original Zenodo terms; this repository does not redistribute it.
