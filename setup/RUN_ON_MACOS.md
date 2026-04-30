# Running plan-eval on a MacBook Air M4 (24 GB)

This file is a self-contained recipe for an agent (or a person) starting fresh on a MacBook Air M4 with 24 GB unified RAM. The goal: replicate the error-handling experiment (`README.md` §2.6) on Apple silicon, so the existing Jetson + RTX 5080 result picks up a third hardware platform.

The experiment matrix from §2.6 is 5 models × 7 injection patterns × 1–2 target tools × 2 recipe variants × 3 seeds = ≈ 390 cells per model class. We'll run 3 Anthropic models + 4 fitting open-weight models, in parallel, overnight.

## Prerequisites (one-time, ~10 min)

```bash
xcode-select --install                                # if not already installed
brew install git curl ollama
npm install -g @anthropic-ai/claude-code              # the `claude` CLI used by the harness
```

After install: run `ollama serve &` (or trust that brew installed it as a launchd service — `ollama list` should return without error). Then `claude` once interactively to log in.

Verify:

```bash
git --version
curl --version | head -1
ollama --version
claude --version
```

All four should print versions. If `claude` errors with "not authenticated", run plain `claude` to start the login flow.

## Setup the repo (~5 min)

```bash
git clone https://github.com/nekrut/plan-eval && cd plan-eval
bash setup/install.sh                       # detects Darwin + arm64, pulls Miniforge3-MacOSX-arm64.sh, creates the bench env
bash setup/fetch_data.sh                    # 9 files, ~838 KB
bash ground_truth/canonical.sh              # produces ground_truth/results/
```

`install.sh` was patched to detect Darwin in commit `cb03333`. The conda env `bench` includes bwa 0.7.18, samtools 1.21, bcftools 1.21, lofreq 2.1.5, all native arm64 builds from bioconda. The fetch and ground-truth steps are platform-agnostic.

After ground truth: `ls ground_truth/results/` should show four `.vcf.gz` files plus their `.tbi` indexes.

## Pull ollama models (~10–20 min, ~52 GB disk)

24 GB unified RAM has to hold macOS, Ollama's runtime, and the model. The four models below cover the three tiers seen in §2.6 (frontier, mid-size local, MoE) plus a control:

```bash
ollama pull granite4               # 2.1 GB — fastest passer, defensive-scripting floor
ollama pull qwen3:14b              # 9 GB — clean dense fit
ollama pull qwen3-coder:30b        # 18 GB — borderline; partial offload territory
ollama pull qwen3.6:35b-a3b        # 23 GB — MoE, the data point unique to unified memory
```

If `qwen3.6:35b-a3b` evicts macOS into swap so badly that cells time out at 900 s, drop it: `ollama rm qwen3.6:35b-a3b` and exclude it from the ollama matrix below. Don't pull anything 30 GB+ dense — those won't run usefully.

## Sanity check (~2 min)

Before the long run, confirm the harness end-to-end:

```bash
mkdir -p runs_smoke_m4
python3 harness/run_one.py \
    --model granite4 --track A --seed 42 \
    --plan plan/PLAN.md \
    --runs-dir runs_smoke_m4 \
    --think off
python3 score/score_run.py runs_smoke_m4/granite4* | head -20
```

Expect `M1=1, M3=1.0`. The cell wall should be ~15–30 s.

If `M3=0.0`, check:

1. `runs_smoke_m4/granite4*/exec.log` — most likely error is conda env activation (`source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate bench`). Run that command directly and see what it says.
2. `ground_truth/results/` exists and has VCFs.
3. `runs_smoke_m4/granite4*/results/` has the model's output VCFs.

If the harness errors with a path it can't find under `/home/anton/...`, that's a regression of commit `cb03333`. Re-pull main.

## Run the matrix (~6–10 h overnight)

Two parallel background processes, one Anthropic, one Ollama. The Anthropic side hits the API and won't compete for compute; the ollama side saturates the GPU.

```bash
nohup python3 -u harness/error_matrix.py \
    --models claude-opus-4-7 claude-sonnet-4-6 claude-haiku-4-5 \
    --plans v2 v2_defensive --include-baseline \
    --log error_matrix_anthropic_m4.jsonl \
    > /tmp/error_matrix_anthropic_m4.log 2>&1 &
echo "anthropic pid: $!"

nohup python3 -u harness/error_matrix.py \
    --models granite4 qwen3:14b qwen3-coder:30b qwen3.6:35b-a3b \
    --plans v2 v2_defensive --include-baseline \
    --log error_matrix_ollama_m4.jsonl \
    > /tmp/error_matrix_ollama_m4.log 2>&1 &
echo "ollama pid: $!"
```

Cell count: 234 (Anthropic) + ~312 (4 ollama models × 78 cells) = ~546. Estimated wall: ~6–10 h.

Watch progress:

```bash
wc -l error_matrix_anthropic_m4.jsonl error_matrix_ollama_m4.jsonl
tail -5 /tmp/error_matrix_*_m4.log
```

Spot-check a cell that finished:

```bash
python3 -c "
import json
r = [json.loads(l) for l in open('error_matrix_anthropic_m4.jsonl')]
print(f'{len(r)} cells; last 5:')
for x in r[-5:]:
    print(f'  {x[\"cell\"]:<70} {x.get(\"m_handle\",\"-\"):<10} M3={x.get(\"M3\",\"-\")}')
"
```

The Anthropic spend will be roughly the same as on the Jetson run (~$15 across the 234 Anthropic cells with prompt caching).

## When it's done — push results back

```bash
git checkout -b m4-results
git add error_matrix_anthropic_m4.jsonl error_matrix_ollama_m4.jsonl
# Don't add runs_inject/ — large and gitignored. The jsonl logs have everything needed.
git status --short
git commit -m "M4 error-handling matrix: $(wc -l < error_matrix_anthropic_m4.jsonl) Anthropic + $(wc -l < error_matrix_ollama_m4.jsonl) Ollama cells"
git push -u origin m4-results
```

After push, the original author can pull the results from another machine, regenerate Figure 6 with the third platform stacked in, and update `README.md` §2.6.

## Things that might trip you up

1. **`claude` not authenticated.** Run `claude` once interactively. The harness inherits the login.
2. **macOS bash 3.2 vs the model's bash 5.** `/bin/bash` is bash 3.2 on macOS; the harness wraps with `bash -c "source ... && conda activate bench && bash run.sh"`. The outer bash 3.2 handles only `source` and `&&` (both work). The model's `run.sh` runs inside conda's bash 5 (`$HOME/miniforge3/envs/bench/bin/bash`). So 5-only features in the model's script are fine.
3. **Disk pressure.** ~52 GB ollama blobs + ~3 GB conda env + ~5 GB run artifacts. Free ≥ 70 GB before starting.
4. **Memory pressure on `qwen3.6:35b-a3b`.** macOS will swap aggressively to disk. If cells time out at 900 s, drop the model.
5. **Ollama on macOS uses Metal/MPS, not CUDA.** Same HTTP API, same model tags — the harness sees no difference.
6. **Activity Monitor.** `WindowServer` and other macOS background processes will show ~3–5 GB resident; account for that when picking models.

## What the existing repo has (for context)

- `README.md` — full paper, includes §2.6 which is the experiment we're replicating
- `harness/error_matrix.py` — the driver
- `harness/error_shims/shim.py` — the PATH shim that injects failures
- `score/score_run.py` — computes M1/M2/M3/M5 plus the error-handling triple `m_handle`/`m_recover`/`m_diagnose`
- `plan/PLAN.md` — the v2 happy-path recipe
- `plan/PLAN_v2_defensive.md` — the v2_defensive recipe (Opus-authored)
- `error_matrix_anthropic.jsonl` / `error_matrix_ollama.jsonl` — the existing Jetson results, for cross-reference

The Jetson + RTX 5080 numbers are already in the repo. The M4 run adds a third hardware column to the existing tables and Figure 6.
