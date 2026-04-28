# plan-eval

Benchmark comparing **Anthropic Claude models** and **a wide selection of free local open-weight models** running on consumer hardware (Jetson AGX Orin and RTX 5080) at executing an Opus-authored bioinformatics plan.

## What this is, in plain English

You ask the strongest available model (Claude Opus 4.7) to write a recipe for a bioinformatics workflow. Then you hand that recipe to a cheaper or smaller model and ask it to turn the recipe into an actual runnable script. We then run the script on real data and check whether its variant calls match a known-correct answer.

The benchmark answers four questions:

1. Which Anthropic tier is the cheapest model that can faithfully execute an Opus-authored recipe?
2. Which **free, locally-runnable** models can do the same job — i.e. avoid paying for inference at all?
3. How does the *level of detail in the recipe* interact with model capability?
4. On a faster consumer GPU (RTX 5080, 16 GB), does the v1→v2 plan-detail effect generalize across the broader Ollama model zoo, or does the lean v1 plan suddenly work?

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
- **Variant identity** — every variant call is a 4-tuple: `(chromosome, position, reference allele, alternate allele)`. For example, `(chrM, 16519, T, C)` means "on the mitochondrial chromosome at base 16 519, the reference genome has a T but this sample has a C." Two variants are considered "the same" if and only if all four fields match.
- **PASS variants** — variant callers flag each call with a confidence filter; we only score variants the caller marked `PASS` (or unfiltered).
- **Allele frequency (AF)** — the fraction of reads at that position that support the alternate allele. AF=1.0 means a clean homozygous variant; AF=0.04 means a heteroplasmy (4% of reads carry the variant). Different callers can estimate this slightly differently for the same variant, so we allow ±0.02 wiggle room before counting two calls as different.
- **M3 (variant agreement, primary score)** — for each of the 4 samples, the **Jaccard index** on the *set* of `(chromosome, position, reference, alternate)` 4-tuples among PASS variants, allowing AF differences within ±0.02. Jaccard = (variants in BOTH the model's VCF and the known-correct VCF) / (variants in EITHER). The four per-sample Jaccards are averaged. **1.000 = every variant matches and only those variants; 0.000 = no overlap.** See "Scoring" below for the precise definition.
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

## RTX 5080 follow-up: 11 local models

Same benchmark, run on an **NVIDIA RTX 5080 (16 GB VRAM)** desktop with 125 GB system RAM. Eleven Ollama models, both plans, both tracks, three seeds each — 168 cells. `/think` was tested only on models that fit cleanly in VRAM (qwen3:8b/14b, qwen3.5:9b); on Group B partial-offload models (qwen3.5:27b, qwen3.6:27b/35b-a3b, glm-4.7-flash) the same `/think` wall-clock blowout that crippled the Jetson reproduces here, so we recorded that as a finding rather than spending ~24 hr of timeout cells.

### v1 (lean plan), Track A — only one local model breaks the wall

| model | quant size | M3 (n=3) | gen secs | notes |
|---|---|---:|---:|---|
| **qwen3.6:27b** `/no_think` | 17 GB | **1.000±0.000** | 353 | dense 27b — only model that solves the lean plan |
| glm-4.7-flash `/no_think` | 19 GB | 0.667±0.577 | 156 | 30B-class GLM, 2/3 hits |
| qwen3.5:27b `/no_think` | 17 GB | 0.333±0.577 | 276 | 1/3 |
| qwen3-coder:30b `/no_think` | 19 GB | 0.333±0.577 | 161 | 1/3 |
| gemma4:26b `/no_think` | 17 GB | 0.333±0.577 | 202 | 1/3 |
| qwen3.6:35b-a3b `/no_think` | 23 GB MoE | **0.000±0.000** | 209 | confirms Jetson — model not hardware |
| qwen3:14b /think and /no_think | 9.3 GB | 0.000±0.000 | 11–69 | reasoning didn't help |
| qwen3:8b, qwen3.5:9b, gemma4:e4b | 5–9 GB | 0.000 | 5–73 | small models lose the plan |
| gpt-oss:20b `/no_think` | 13 GB | empty_script (3 of 6 cells) | 53 | internal CoT exhausts num_predict |
| qwen3.5:9b `/think` | 6.6 GB | empty_script (6 of 6) | 74→151 | reasoning eats the entire budget, no script emitted even at 16384 num_predict |

### v2 (hyper-detailed plan), Track A — every fitted model nails it

| model | n | M3 | gen secs |
|---|---:|---:|---:|
| qwen3.6:35b-a3b `/no_think` | 3 | **1.000** | 127 |
| qwen3.6:27b `/no_think` | 3 | **1.000** | 219 |
| qwen3.5:27b `/no_think` | 3 | **1.000** | 214 |
| qwen3-coder:30b `/no_think` | 3 | **1.000** | 148 |
| glm-4.7-flash `/no_think` | 3 | **1.000** | 150 |
| gemma4:26b `/no_think` | 3 | **1.000** | 121 |
| qwen3:14b `/no_think` | 3 | **1.000** | 6 |
| qwen3:14b `/think` | 3 | **1.000** | 13 |
| qwen3.5:9b `/no_think` | 3 | **1.000** | 6 |
| gemma4:e4b `/no_think` | 3 | **1.000** | 9 |
| qwen3:8b `/think` | 3 | **1.000** | 9 |
| qwen3:8b `/no_think` | 3 | 0.667±0.577 | 6 |
| gpt-oss:20b `/no_think` | 3 | 0.667±0.577 | 53 |
| qwen3.5:9b `/think` | 3 | empty_script (all) | 151 |

### v1 + v2, Track B (no plan) — all local models 0/3

Across the 5080 model set, no Ollama model recovers the canonical mtDNA workflow without the plan. Track B sits at 0.000 for every model, sometimes with a single 1/3 fluke. This is a sharper failure than the Anthropic Track B numbers (Opus/Sonnet 0.938, Haiku 0.667 ± 0.577) — local models truly need the plan, the Anthropic models can mostly recover it from internal knowledge.

### Findings (5080)

**The v1→v2 plan-detail effect generalizes.** Every local model that *can* solve the task on the 5080 solves v2 at 1.000. The lean v1 plan continues to break almost everyone. The Jetson finding wasn't an artifact of slow hardware — it's an artifact of plan specificity, replicated across 11 different model architectures.

**Exactly one local model under 30B parameters solves the lean v1 plan: qwen3.6:27b dense.** Smaller dense Qwens (8b, 14b, 9b) don't, and neither does the 35b-a3b MoE on either hardware. This isn't a parameter-count story — qwen3.5:27b, qwen3-coder:30b, and gemma4:26b all have similar/larger budgets and only manage 1/3 each. Some training-data quirk in the 3.6-dense release puts lofreq's CLI in reach.

**`/think` mode on small models is a footgun, not an upgrade.** qwen3.5:9b `/think` failed all 12 cells: the model exhausted the 16384-token output budget on reasoning and never emitted a final script (`raw_response` was empty across all seeds; usage showed `eval_count=8192` then `eval_count=16384`). This is the same hardware-class budget exhaustion the Jetson hit on `/think` — but here it's a model-side issue, not a wall-clock one. Bigger models (qwen3:14b /think) don't fail this way but also don't outperform `/no_think`.

**`gpt-oss:20b` has the same problem without an explicit /think toggle.** Its built-in Harmony reasoning consumes the output budget and ~half its cells emit empty scripts. Effective M3 lands at 0.667 on v2 Track A despite the model being capable when it does emit a script.

**Hardware does not rescue plan specificity.** On the 5080, qwen3.6:35b-a3b on v1 still scored 0/3 — same as Jetson. Group B partial-offload speeds (100–450 s/cell) are surprisingly close to Jetson's unified-memory speeds (~100 s) because PCIe weight transfer dominates either way. The 5080's win is on Group A in-VRAM models, where 5–30 s/cell crushes anything the Jetson could do.

**Track B without a plan is a wall for local models.** This is the cleanest inversion of the Anthropic story: Opus can almost recover the workflow from problem statement alone (0.938); local 30B-class models cannot (0.000). The plan-then-implement pattern isn't a productivity hack for local models — it's a capability prerequisite.

### Recommendation (5080)

If you have a 16 GB consumer GPU and want $0/run local execution:

- **For the lean v1 style of plan**: `qwen3.6:27b /no_think` is the only reliable choice. ~5 min/run with PCIe offload, no Anthropic dependency.
- **For the hyper-detailed v2 style**: pick whichever fits cleanly in VRAM and is fastest — `gemma4:e4b` (~9 s/run, 9.6 GB), `qwen3:14b` (~6 s/run, 9.3 GB), or `qwen3.5:9b` (~6 s/run, 6.6 GB). All hit 1.000 on Track A.
- **Avoid `/think` on small models** for code-generation tasks unless you've tuned `num_predict` very high and verified the model emits a final answer.

## What part of v2 is doing the work? — three intermediates

v1 → v2 is a >2× token bump and ~1.000 score jump across the board. To isolate *which* part of v2 carries the load we ran three intermediate plan/track conditions on the same 11 local + 3 Anthropic models, three seeds each, Track A unless noted (153 cells, ~5 h on the 5080, $0).

| condition | what it adds to the previous step | tests |
|---|---|---|
| **v1.25** | v1 plus the literal `lofreq call-parallel` command (one code-fenced line) | does the lofreq positional-arg surprise *alone* explain the cliff? |
| **v1.5** | v2 stripped of every prose paragraph and "Gotchas" block — pure command lines | are the prose warnings load-bearing, or decorative? |
| **v0.5** | Track B (no plan) plus a single tool-name sequence line (`bwa → samtools → lofreq → bcftools → awk`) | does sequencing alone move the needle from a zero baseline? |

### v1.25 — one literal lofreq command unlocks 5 of 11 local models

| model | v1 | **v1.25** | v2 |
|---|---:|---:|---:|
| qwen3.6:27b dense | 1.000 | 1.000 | 1.000 |
| qwen3.5:27b dense | 0.333 | **1.000** | 1.000 |
| qwen3.6:35b-a3b MoE | **0.000** | **1.000** | 1.000 |
| gemma4:26b | 0.333 | **1.000** | 1.000 |
| qwen3-coder:30b | 0.333 | **1.000** | 1.000 |
| glm-4.7-flash | 0.667 | 0.333 | 1.000 |
| gpt-oss:20b | 0/6 (CoT) | 1/3 ok | 0.667 |
| qwen3:14b /no_think | 0.000 | 0.000 | 1.000 |
| qwen3.5:9b /no_think | 0.000 | 0.000 | 1.000 |
| qwen3:8b /no_think | 0.000 | 0.000 | 0.667 |
| gemma4:e4b | 0.000 | 0.333 | 1.000 |

Adding **one** code-fenced line — the exact `lofreq call-parallel` invocation with the BAM as a positional argument — flips four medium-large models from broken to perfect. The qwen3.6:35b-a3b MoE result is the cleanest: 0/3 on v1 → 3/3 with one extra line. The cliff for these models was a single command-syntax surprise: `lofreq` takes the BAM as a positional argument, not behind `-i`/`-b`/`-bam`, and v1's "Input: `results/{sample}.bam`" prose didn't disambiguate. Smaller models (≤14B dense) still score 0 — they need more than one fix.

### v1.5 — stripping all prose from v2 leaves M3 essentially unchanged

| model | v1.25 | **v1.5 (commands-only)** | v2 |
|---|---:|---:|---:|
| qwen3:8b /no_think | 0.000 | **1.000** | 0.667 |
| qwen3:14b /no_think | 0.000 | **1.000** | 1.000 |
| qwen3.5:9b /no_think | 0.000 | **1.000** | 1.000 |
| qwen3.5:27b | 1.000 | 1.000 | 1.000 |
| qwen3.6:27b | 1.000 | 1.000 | 1.000 |
| qwen3.6:35b-a3b | 1.000 | 1.000 | 1.000 |
| gemma4:26b | 1.000 | 1.000 | 1.000 |
| qwen3-coder:30b | 1.000 | 1.000 | 1.000 |
| glm-4.7-flash | 0.333 | **1.000** | 1.000 |
| gpt-oss:20b | 1/3 ok | 0.500 (1/2 ok, 1 err) | 0.667 |
| gemma4:e4b (4 B) | 0.333 | 0.000 | 1.000 |

v1.5 is `PLAN.md` (v2) with **every** explanatory paragraph and every "Gotchas" subsection deleted — only the numbered headings and the code-fenced commands remain. Result: 9 of 11 local models hit 1.000, including the small models that v1.25 couldn't unlock. The prose warnings ("Do NOT use `printf`/`echo -e`/`$'\t'`", "bgzip operates in place", "`%INFO/AF`, not `%AF`") that read like the load-bearing wisdom of v2 turn out to be **decorative** for this task — the verbatim commands convey the same constraints implicitly.

The two non-1.000 outliers tell their own stories: gemma4:e4b (4 B) is just too small (0/3 on both v1.5 *and* v1.25, but 1.000 on v2 — more context tokens seem to help it cohere); gpt-oss:20b's Harmony chain-of-thought continues to eat the output budget on roughly half its cells regardless of how the plan is written.

### v0.5 — tool-name ordering doesn't help (controls confirm the zero baseline)

| model | Track B (v1) | **v0.5 (B + tool order)** |
|---|---:|---:|
| 11 local models, average | 0.000 | 0.000 (one fluke each at gemma4:26b, qwen3.6:27b, gpt-oss:20b) |
| Claude Haiku 4.5 | 0.667±0.577 | **1.000** |
| Claude Sonnet 4.6 | 0.938±0.000 | 0.938±0.000 |
| Claude Opus 4.7 | 0.938±0.000 | 0.667±0.577 |

Telling local models the *order* of tools to call (without flags or commands) does nothing. The Anthropic numbers are noisy at this n=3 — they hover around the same place as plain Track B; the apparent Opus regression is a single seed-42 zero-score that lands inside the variance band. The substantive control conclusion: **sequencing alone is not what local models need**. They need every command, character-for-character. The plan's job is not to tell the model *what to do in what order* — local models can guess that. The plan's job is to literalize the syntax of every tool the model doesn't already know.

### Findings (intermediates)

**The cliff between v1 and v2 has two distinct rungs.**

1. **For ≥27 B dense local models, the cliff is a single command.** v1.25 — v1 with one extra code-fenced `lofreq call-parallel` line — hits 1.000 across all of them (qwen3.5:27b, qwen3.6:27b, qwen3.6:35b-a3b, gemma4:26b, qwen3-coder:30b). Adding any other v2 detail does nothing for them. They had already inferred BWA, samtools, bgzip, tabix and the bcftools format string from v1; the *one* tool whose CLI they couldn't reconstruct from prose was lofreq, specifically because the BAM is positional and v1 didn't say so.
2. **For ≤14 B local models, every command needs to be literalized**, but the prose around the commands does not. v1.5 (v2 minus prose) brings qwen3:8b, qwen3:14b, and qwen3.5:9b to 1.000 — same as v2 — without any of the gotchas, escape-character warnings, or guard-clause boilerplate. Smaller models have less internal CLI knowledge across the board, so every step needs the exact incantation; but they *don't* benefit from explanations of why.

**The prose in v2 is for human readers, not for local models.** This was the most surprising result: stripping every paragraph from v2, including the read-group `\t` warning that we'd considered load-bearing, leaves the score unchanged. Models read code blocks as code; the text between them is mostly ignored.

**Sequencing without syntax is worthless.** Telling a local model "call these tools in this order" without specifying how is operationally identical to giving it nothing. The information that matters is per-tool CLI specifics, not workflow shape.

### Recommendation (intermediates)

For future plans aimed at local-model implementers:

- **Always show the full command line for any tool with non-obvious CLI conventions** — positional arguments, format-string syntax, in-place behavior, escape rules. For this workflow, `lofreq call-parallel` was the single such tool.
- **Skip the prose explanations.** They cost tokens, take time to write, and don't change scores for code-emitting open-weight models.
- **Don't bother sequencing without syntax.** "Use bwa → lofreq → bcftools" is no more useful than "no plan."
- **For ≥27 B dense models, v1.25 is enough.** If you want one plan for both small and large models, write v1.5 (commands-only).

## Can a tool registry replace human plan authorship? — v1g

The v1.25 result raises a follow-up: if the cliff is one literal `lofreq` command line, does that command have to be authored by a human/Opus? **Galaxy's IUC tool collection** (`galaxyproject/tools-iuc`) is a community-curated registry of XML wrappers, one per bioinformatics tool, with `<command><![CDATA[...]]></command>` blocks that — after Cheetah templating substitution at Galaxy runtime — produce the exact CLI invocation the tool expects. lofreq's IUC wrapper has the BAM as a positional argument; samtools_sort's specifies `-@ N`; bcftools_query's records the `-f` format-string convention. If we could pull these and embed them mechanically, plan authorship would reduce to a registry lookup.

We tried this. v1g = v1 + Galaxy-IUC-derived CLI snippet for lofreq, mechanically extracted by `scripts/galaxy_to_snippet.py` from `tools-iuc` commit `39e7456` (2026-04-27).

| model | v1.25 (hand) | **v1g (IUC)** | what happened |
|---|---:|---:|---|
| Claude Opus 4.7 | 1.000 | 1.000 | self-corrected the noisy snippet |
| Claude Sonnet 4.6 | 1.000 | 1.000 | self-corrected |
| qwen3.5:27b dense | 1.000 | 1.000 | self-corrected |
| qwen3.6:27b dense | 1.000 | 1.000 | self-corrected |
| **Claude Haiku 4.5** | 1.000 | **0.000** | copied snippet literally → tool crash |
| qwen3.6:35b-a3b MoE | 1.000 | **0.000** | copied snippet literally |
| gemma4:26b | 1.000 | **0.000** | copied snippet literally |
| qwen3-coder:30b | 1.000 | 0.333 | one seed self-corrected |
| glm-4.7-flash | 0.333 | 0.000 | regressed |
| gpt-oss:20b | 1/3 ok | 0/3 (CoT issues persist) | — |
| Gemini-class small (≤14 B) | 0.000 | 0.000 | unchanged |

**The bug**: Galaxy's lofreq XML emits `--sig $value` and `--bonf $value` where `$value` is a Galaxy-runtime parameter. Cheetah strips the variables; the bare flags `--sig` and `--bonf` survive and end up in the snippet on adjacent lines. lofreq's argument parser then reads `--sig --bonf` as `--sig=--bonf`, attempts `float("--bonf")`, and dies with a Python traceback:
```
ValueError: could not convert string to float: '--bonf'
```

Sonnet, Opus, qwen3.5/3.6:27b dense, and the strongest local models recognized the bare-flag pattern as malformed and dropped both flags before emitting the script. Haiku, the MoE, and the smaller models copied the snippet character-for-character — **including the broken bare flags** — and shipped a script that lofreq refuses to parse. **Haiku regressed from 3/3 perfect on v1.25 to 0/3 on v1g** purely from snippet noise.

### Findings (Galaxy IUC as a registry)

**Galaxy IUC wrappers are not self-contained CLIs.** They are *Galaxy-runtime templates* — Cheetah scripts whose flag values are bound at execution time from XML `<param>` definitions and Galaxy's parameter-collection UI. Mechanically stripping the templating leaves bare flags that look correct at the textual level but break the underlying tool. Without Galaxy's runtime to bind values, the extracted commands are syntactically incomplete in a way that's invisible to a deterministic transpiler.

**Coverage is also worse than hoped.** Of the 8 tool steps in this workflow, only **lofreq** has a clean unconditional command core in IUC. The bwa wrappers are heavily macro-ized (most logic in `bwa_macros.xml` and `read_group_macros.xml`). samtools_faidx's command block is entirely conditional — extraction yields an empty string. bcftools_query's load-bearing `-f` format string is itself a Cheetah variable. samtools_index, bgzip, and tabix have no IUC wrappers at all (Galaxy auto-handles indexing for uploaded BAMs).

**Plan authorship via tool registry doesn't replace human authorship for this generation of models.** The strongest models can repair noisy registry output; cheaper models will faithfully emit broken commands and crash. Practically, you either need (a) a more sophisticated extractor that supplies sensible defaults for runtime-bound flags, or (b) human review of the extracted snippet before shipping it into the plan. v1.25 (hand-written) is one such "human review" pass; v1g without it loses Haiku, the 35b-a3b MoE, and gemma4:26b.

A follow-up that would change the picture: implement a per-tool default-value table in the extractor (e.g. `--sig` defaults to `0.01`, `--bonf` defaults to `dynamic`), then re-run the matrix. We've left this as a future direction — the experimental result with mechanical extraction alone is the main finding.

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
| **M3 — Are the variant calls right?** | For each of the 4 samples, the overlap between the model's variant set and the known-correct variant set. **1.000 = perfect, 0.000 = no overlap.** | **Jaccard index** on the *set* of `(chromosome, position, reference allele, alternate allele)` 4-tuples among PASS-filtered records: `J = (calls in BOTH sets) / (calls in EITHER set)`. A call counts as matched only if all four fields are identical AND the AF (allele frequency) values agree within ±0.02. Per-sample Jaccards are averaged across the 4 samples (macro-mean). |
| **M4 — Cost and time** | What it cost and how long it took. | Tokens, USD (Anthropic), generation seconds, execution seconds. |
| **M5 — Is the script clean?** | The script is well-formed bash and re-runnable. | `shellcheck` clean + `set -euo pipefail` present + no hardcoded user paths + re-running on a populated output dir exits 0 with no work performed. Pass/fail. |

If a model uses a *different but valid* tool than the recipe specifies (e.g. `bcftools mpileup` instead of `lofreq`), M3 still scores it correctly — the metric compares variant calls, not the pipeline that produced them.

## Hardware

- **Jetson AGX Orin Developer Kit**: 64 GB unified RAM, Ampere-class GPU (sm_87), aarch64; 30W power mode (MAXN unavailable without a reboot the user declined). Local-model wall-clock figures reflect 30W.
- **RTX 5080 desktop**: NVIDIA GeForce RTX 5080 16 GB VRAM, 125 GB system RAM, x86_64. Models ≤14 GB fit in VRAM (Group A); 17–23 GB models partial-offload to CPU (Group B).
- **Conda env**: locked, on PATH for both canonical and model runs (see `setup/install.sh`).

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
│   ├── PLAN_v1.md                preserved v1 plan
│   ├── PLAN_v1p25.md             v1 + literal lofreq command (intermediate)
│   ├── PLAN_v1p5.md              v2 stripped to commands only (intermediate)
│   └── PLAN_v1g.md               v1 + Galaxy-IUC-derived lofreq snippet (registry test)
├── scripts/
│   └── galaxy_to_snippet.py      Cheetah-XML → bash snippet extractor for IUC wrappers
├── prompts/                      system + per-track user prompt templates
│   └── track_b_with_order_user.tmpl   Track B + tool sequence (v0.5 control)
├── harness/
│   ├── run_one.py                generate + execute one cell
│   ├── matrix.py                 Anthropic + ollama matrix iterator (Jetson)
│   ├── sweep_local.py            disk-rotated local-model sweep (Jetson)
│   └── matrix_5080.py            iterate the 5080 matrix (11 ollama + 3 Anthropic)
├── score/
│   ├── score_run.py              M1–M5 against ground truth
│   └── aggregate.py              scans runs/ + runs_v*/ + runs_5080_v*/ → results.csv
├── runs/<run_id>/                Jetson v2 per-run artifacts
├── runs_v1/<run_id>/             Jetson v1 per-run artifacts (preserved)
├── runs_5080_v1/<run_id>/        5080 v1 per-run artifacts
├── runs_5080_v2/<run_id>/        5080 v2 per-run artifacts
├── runs_5080_v1p25/<run_id>/     5080 v1.25 (v1 + lofreq cmd) per-run artifacts
├── runs_5080_v1p5/<run_id>/      5080 v1.5 (v2 minus prose) per-run artifacts
├── runs_5080_v0p5/<run_id>/      5080 v0.5 (Track B + tool order) per-run artifacts
└── runs_5080_v1g/<run_id>/       5080 v1g (Galaxy IUC registry) per-run artifacts
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

Full matrix:

```bash
python3 harness/matrix.py --tracks A             # both Anthropic and ollama (Jetson)
python3 harness/matrix.py --only-ollama --tracks A --no-think
python3 harness/matrix_5080.py                   # 11 ollama models, v1+v2, both tracks
python3 harness/matrix_5080.py --plans v1p25,v1p5,v0p5 --include-anthropic
python3 score/aggregate.py                       # → results.csv + console summary
```

The Anthropic side authenticates through the existing `claude` CLI (Claude Code) login — no `ANTHROPIC_API_KEY` needed. The local side requires `ollama serve` with the relevant tags pulled (`qwen3.6:35b-a3b` for the Jetson story; `matrix_5080.py` will pull missing tags on first chat).


## Cost

- Total Anthropic spend across both v1 and v2 matrices plus three plan-generation calls: under **$3** with prompt caching.
- Local sweep: **$0** (electricity excluded).

## Caveats

- **One task class.** Per-sample variant calling on a 16.6 kb mitochondrial reference. Generalization to whole-genome workflows, RNA-seq, single-cell, etc. is not implied.
- **One dataset.** All four samples are amplicon-PCR mitochondrial DNA from one study; intentionally easy.
- **Anthropic temperature.** Set by the `claude` CLI default (no `--temperature` flag).
- **Ollama temperature.** 0.2 with seeds 42, 43, 44.
- **30 W power mode on the Jetson.** A higher power mode (MAXN) would require rebooting the machine, which we didn't do. So local-model wall-clock numbers are conservative — they'd be faster on a more aggressive power profile.
- **15-minute generation budget per call.** Four free-local models timed out at this limit on every run. Their 0/3 score means "did not fit our 15-minute budget on this Jetson," not "the model can't do the task in principle." A 30-minute budget would likely recover most of them.
- **Disk rotation during the Jetson sweep.** We pull → test → remove for each model to fit in ~17 GB free disk. The full 14-model sweep transferred ~210 GB sequentially over the network.
- **5080 num_predict.** Bumped from 8192 to 16384 for the 5080 matrix after qwen3.5:9b /think and gpt-oss:20b cells exhausted the original 8192-token output budget on internal reasoning. qwen3.5:9b /think kept hitting the wall even at 16384 (pure reasoning, no script emitted) — reported as a finding, not a config bug.
- **5080 /think on Group B disabled.** A single qwen3.5:27b /think cell hit the 1800 s urlopen timeout without returning anything. We disabled /think across qwen3.5:27b, qwen3.6:27b, qwen3.6:35b-a3b and glm-4.7-flash to avoid ~24 hr of timeout-doomed cells; the same Jetson `/think` failure mode reproduces here on offload-bound models.
- **No-plan baseline noise.** The v1 Track B `0.938` vs `1.000` gap on Opus and Sonnet reflects a single-variant disagreement on one of the 4 samples (these models pick a slightly different caller without the plan, recovering most but not all variants at the AF tolerance).
- **The v2 recipe is close to "the script in prose."** That's the whole point — the experimental finding is that this is what you need to write to get a free local model working reliably. If you want the local model to *also* do bioinformatics reasoning, more work is needed.

## License

MIT. See `LICENSE`.

The dataset itself (FASTQ files in Zenodo 5119008) is a separate work by A. Nekrutenko under the original Zenodo terms; this repository does not redistribute it.
