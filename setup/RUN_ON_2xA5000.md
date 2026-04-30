# Running plan-eval on a 2× NVIDIA A5000 box

This file is a self-contained recipe for an agent (or a person) starting fresh on a Linux machine with two NVIDIA A5000 GPUs (24 GB VRAM each, 48 GB total). The goal: extend the error-handling experiment (`README.md` §2.6) to the larger open-weight models that wouldn't fit on the Jetson, the RTX 5080 (16 GB VRAM), or the MacBook Air M4 (24 GB unified RAM).

The hardware unlocks four model classes that the rest of the experiment couldn't reach:

- `llama3.3:70b-instruct-q3_K_M` (34 GB) — needs both GPUs (Ollama auto-splits across them via CUDA)
- `granite-code:34b` (19 GB), `qwen3:32b` (20 GB) — fit in one A5000 each
- `gemma3:27b` (17 GB), `mistral-small3.2:24b` (15 GB) — different model families not yet in §2.6
- the existing §2.6 lineup (Opus / Sonnet / Haiku, qwen3.6:27b, qwen3.6:35b, granite4) for cross-platform confirmation

The headline question: does the **frontier-tier** of §2.6 — currently four models (Opus, Sonnet, Haiku, qwen3.6:27b) all at 21 / 14–15 / 0 / 0–1 — extend to the larger open-weight dense models, or does it stay an idiosyncrasy of qwen3.6:27b at 17 GB?

## Prerequisites (one-time, ~15 min)

```bash
sudo apt update && sudo apt install -y git curl build-essential
curl -fsSL https://ollama.com/install.sh | sh        # Linux Ollama, CUDA-backed
sudo systemctl start ollama                          # default: localhost:11434
npm install -g @anthropic-ai/claude-code             # the `claude` CLI used by the harness
```

Verify CUDA is visible to Ollama:

```bash
nvidia-smi               # both A5000s should appear
ollama --version
ollama ps                # empty is fine
```

If `ollama ps` errors with "could not connect", `sudo systemctl start ollama` and retry.

For the `claude` CLI: run `claude` once interactively to log in, then `claude --version`.

## Setup the repo (~5 min)

```bash
git clone https://github.com/nekrut/plan-eval && cd plan-eval
bash setup/install.sh                  # detects Linux + x86_64, pulls Miniforge3-Linux-x86_64.sh, creates the bench env
bash setup/fetch_data.sh               # 9 files, ~838 KB
bash ground_truth/canonical.sh         # produces ground_truth/results/
```

`install.sh` was already arch-detect-ready from the 5080 commit, plus OS-detection from the macOS work. The conda env `bench` includes bwa 0.7.18, samtools 1.21, bcftools 1.21, lofreq 2.1.5 — all native x86_64 builds from bioconda.

After ground truth: `ls ground_truth/results/` should show four `.vcf.gz` files plus their `.tbi` indexes.

## Pull ollama models (~30–40 min, ~110 GB disk)

The 48 GB of VRAM lets us run all of the §2.6 lineup plus three large dense models that previously didn't fit.

```bash
# already in §2.6 — for cross-platform confirmation:
ollama pull granite4               # 2.1 GB — defensive-scripting floor
ollama pull qwen3.6:27b            # 17 GB — current frontier-tier local
ollama pull qwen3.6:35b            # 23 GB — current mid-tier local

# new for this run — large dense models we couldn't test on smaller hardware:
ollama pull llama3.3:70b-instruct-q3_K_M    # 34 GB — needs both A5000s; the largest in the lineup
ollama pull granite-code:34b                # 19 GB — coder-tuned, IBM lineage
ollama pull qwen3:32b                       # 20 GB — tests whether qwen3.6:27b's tier is a one-off
ollama pull gemma3:27b                      # 17 GB — Google lineage, different from qwen
ollama pull mistral-small3.2:24b            # 15 GB — Mistral lineage
```

Total: ~147 GB of ollama blobs. Time: ~30–40 minutes depending on bandwidth. Free ≥ 200 GB before pulling.

## Sanity check (~2 min)

Before the long run, confirm the harness end-to-end with the smallest model:

```bash
mkdir -p runs_smoke_a5000
python3 harness/run_one.py \
    --model granite4 --track A --seed 42 \
    --plan plan/PLAN.md \
    --runs-dir runs_smoke_a5000 \
    --think off
python3 score/score_run.py runs_smoke_a5000/granite4* | head -20
```

Expect `M1=1, M3=1.0`, wall ~15–25 s.

If `M3=0.0`, look at `runs_smoke_a5000/granite4*/exec.log` — the most likely failure on a fresh machine is the conda env activation. Run `source ~/miniforge3/etc/profile.d/conda.sh && conda activate bench && which lofreq` directly to confirm.

## Run the matrix (~10–14 h overnight)

Two parallel background processes, one Anthropic, one Ollama. The Anthropic side hits the API and won't compete for compute; the ollama side saturates the GPUs. Eight ollama models is on the long side — split into two passes (smaller models first, large dense second) if you want incremental progress.

```bash
# Anthropic side — same three tiers as Jetson §2.6 for cross-platform confirmation
nohup python3 -u harness/error_matrix.py \
    --models claude-opus-4-7 claude-sonnet-4-6 claude-haiku-4-5 \
    --plans v2 v2_defensive --include-baseline \
    --log error_matrix_anthropic_a5000.jsonl \
    > /tmp/error_matrix_anthropic_a5000.log 2>&1 &
echo "anthropic pid: $!"

# Ollama side — eight models. The new ones are the science.
nohup python3 -u harness/error_matrix.py \
    --models granite4 qwen3.6:27b qwen3.6:35b qwen3:32b granite-code:34b gemma3:27b mistral-small3.2:24b llama3.3:70b-instruct-q3_K_M \
    --plans v2 v2_defensive --include-baseline \
    --log error_matrix_ollama_a5000.jsonl \
    > /tmp/error_matrix_ollama_a5000.log 2>&1 &
echo "ollama pid: $!"
```

Cell count: 234 (Anthropic) + 624 (8 ollama × 78 cells) = **858 cells**. Estimated wall: ~10–14 h depending on `llama3.3:70b` per-cell time (4× A5000 inference plus model split overhead).

If you want a shorter first pass that still answers the headline question, drop to four ollama models: `qwen3.6:27b granite-code:34b qwen3:32b llama3.3:70b-instruct-q3_K_M` — that's 312 ollama cells, ~5–7 h.

Watch progress:

```bash
wc -l error_matrix_anthropic_a5000.jsonl error_matrix_ollama_a5000.jsonl
tail -5 /tmp/error_matrix_*_a5000.log
nvidia-smi             # confirm both GPUs are loaded during ollama cells
```

Spot-check a finished cell:

```bash
python3 -c "
import json
r = [json.loads(l) for l in open('error_matrix_ollama_a5000.jsonl')]
print(f'{len(r)} cells; last 5:')
for x in r[-5:]:
    print(f'  {x[\"cell\"]:<70} {x.get(\"m_handle\",\"-\"):<10} M3={x.get(\"M3\",\"-\")}')
"
```

Anthropic spend: ~$15 across the 234 Anthropic cells (same as Jetson run).

## When it's done — push results back

```bash
git checkout -b a5000-results
git add error_matrix_anthropic_a5000.jsonl error_matrix_ollama_a5000.jsonl
# Don't add runs_inject/ — large and gitignored. The jsonl logs and the per-run score.json files have everything needed.
git status --short
git commit -m "2x A5000 error-handling matrix: $(wc -l < error_matrix_anthropic_a5000.jsonl) Anthropic + $(wc -l < error_matrix_ollama_a5000.jsonl) Ollama cells"
git push -u origin a5000-results
```

After push, the original author can pull from another machine, append the new rows to the master `error_matrix_*.jsonl` files, regenerate Figure 6 with the new rows stacked in, and update §2.6.

## Things that might trip you up

1. **`claude` not authenticated.** Run `claude` once interactively. The harness inherits the login.
2. **CUDA visibility.** `nvidia-smi` should list both A5000s. If only one shows up, check `nvidia-persistenced` is running and `lsmod | grep nvidia` shows the kernel modules. Ollama auto-uses both GPUs through CUDA — no model-side configuration needed.
3. **`llama3.3:70b` cold start.** First call to a 34 GB model takes ~30 s of model load before the first token. Subsequent calls within the same Ollama session reuse the loaded weights. The matrix harness sequences cells per-model, so the model loads once per `--models` entry.
4. **Disk pressure.** ~147 GB of ollama blobs + ~3 GB conda env + ~10 GB of run artifacts. Free ≥ 200 GB before starting.
5. **Two-GPU split for `llama3.3:70b`.** Ollama splits layers across both GPUs automatically. If `nvidia-smi` shows only one A5000 loaded during a `llama3.3:70b` cell, the model is fitting in one GPU's VRAM (less likely at 34 GB) — that's fine and faster.
6. **`gpt-oss:20b` excluded.** Its Harmony reasoning eats the output budget; on the Jetson and the 5080 it consistently emitted empty scripts on half its cells. Worth re-testing if you have headroom, but not in the headline lineup above.

## What to tell the author when results are pushed

The most useful single sentence is: "Of the new dense models tested at scale, `<list>` matched the frontier-tier (21/14/0/0); `<list>` did not." That answers the headline question of whether qwen3.6:27b's tier behavior generalizes to bigger models or stays specific to that one model.

If `llama3.3:70b` lands in the frontier tier, that's a clean ceiling result. If `qwen3:32b` matches qwen3.6:27b but `qwen3.6:35b` (already tested, doesn't) and `gemma3:27b` differ, that suggests the tier boundary is training-data-specific rather than size- or family-specific.

## What the existing repo has (for context)

- `README.md` — full paper, includes §2.6 which is what we're extending
- `harness/error_matrix.py` — the driver
- `harness/error_shims/shim.py` — the PATH shim that injects failures
- `score/score_run.py` — computes M1/M2/M3/M5 plus the error-handling triple `m_handle`/`m_recover`/`m_diagnose`
- `plan/PLAN.md` — the v2 happy-path recipe
- `plan/PLAN_v2_defensive.md` — the v2_defensive recipe (Opus-authored)
- `error_matrix_anthropic.jsonl` / `error_matrix_ollama.jsonl` — the existing Jetson results, for cross-reference
- `setup/RUN_ON_MACOS.md` — sister recipe for the M4 run, same matrix on smaller hardware

The existing tier counts to compare against:

| model | recover | partial | propagate | crash |
|---|---:|---:|---:|---:|
| Opus / Sonnet / Haiku / **qwen3.6:27b** | 21 | 14–15 | 0 | 0–1 |
| qwen3.6:35b | 7 | 29 | 0 | 0 |
| granite4 | 0 | 0 | 0 | 36 |

The 2× A5000 run will tell us where the new models fall on this spectrum.
