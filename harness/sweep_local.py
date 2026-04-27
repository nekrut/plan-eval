#!/usr/bin/env python3
"""
Sweep across many local Ollama models with disk rotation.

For each model:
  1. ollama pull <model>
  2. for seed in {42, 43, 44}: run_one.py --track A --think off
  3. score each run
  4. ollama rm <model>   (free disk before the next model)

The full v2 plan (plan/PLAN.md) is used unchanged.

Usage:
  python3 harness/sweep_local.py
  python3 harness/sweep_local.py --models qwen3-coder:30b deepseek-coder-v2:16b
  python3 harness/sweep_local.py --keep   # don't `ollama rm` after testing
"""
from __future__ import annotations
import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
RUNS = BENCH / "runs"
RUN_ONE = BENCH / "harness" / "run_one.py"
SCORE = BENCH / "score" / "score_run.py"

SEEDS = [42, 43, 44]

# 14-model lineup. Sizes verified against ollama registry on 2026-04-26.
MODELS = [
    # Tier A — coder-tuned
    ("qwen3-coder:30b",                       18.6, "coder"),
    ("devstral-small-2:24b",                  15.2, "coder"),
    ("deepseek-coder-v2:16b",                 11.0, "coder"),
    ("granite-code:34b",                      19.2, "coder"),
    # Tier B — generalist instruct
    ("glm-4.7-flash",                         19.0, "generalist"),
    ("gpt-oss:20b",                           13.8, "generalist"),
    ("gemma3:27b",                            17.4, "generalist"),
    ("mistral-small3.2:24b",                  15.0, "generalist"),
    ("qwen3:32b",                             20.2, "generalist"),
    ("olmo-3.1:32b",                          19.5, "generalist"),
    ("nemotron-3-nano",                       24.3, "generalist"),
    ("llama3.3:70b-instruct-q3_K_M",          34.3, "generalist"),
    # Small controls
    ("glm4:9b",                                5.5, "small"),
    ("granite4",                               2.1, "small"),
]


def run(cmd, **kw):
    """Run cmd, stream stdout/stderr to our stdout, return CompletedProcess."""
    print(f"  $ {' '.join(cmd) if isinstance(cmd, list) else cmd}", flush=True)
    return subprocess.run(cmd, **kw)


def free_disk_gb() -> float:
    return shutil.disk_usage("/").free / 1e9


def model_present(tag: str) -> bool:
    p = subprocess.run(["ollama", "list"], capture_output=True, text=True)
    return tag in p.stdout


def pull(tag: str) -> dict:
    t0 = time.time()
    p = subprocess.run(["ollama", "pull", tag], capture_output=True, text=True)
    return {"ok": p.returncode == 0, "secs": time.time() - t0, "err": p.stderr[-800:] if p.returncode else ""}


def remove(tag: str) -> None:
    subprocess.run(["ollama", "rm", tag], capture_output=True, text=True)


def already_scored_count(tag: str) -> int:
    safe = tag.replace("/", "_").replace(":", "_")
    n = 0
    for d in RUNS.glob(f"{safe}_think-off_track-A_seed-*_*"):
        if (d / "score.json").exists():
            n += 1
    return n


def sweep_model(tag: str, size_gb: float, keep: bool) -> dict:
    print(f"\n[sweep] === {tag}  ({size_gb:.1f} GB)  ===", flush=True)
    print(f"[sweep] free disk before: {free_disk_gb():.1f} GB", flush=True)

    skip_existing = already_scored_count(tag)
    if skip_existing >= len(SEEDS):
        print(f"[sweep] {tag}: {skip_existing}/{len(SEEDS)} seeds already scored — skipping", flush=True)
        return {"model": tag, "skipped": True, "completed_seeds": skip_existing}

    if not model_present(tag):
        print(f"[sweep] pulling {tag}...", flush=True)
        r = pull(tag)
        if not r["ok"]:
            print(f"[sweep] pull failed for {tag} after {r['secs']:.0f}s: {r['err']}", flush=True)
            return {"model": tag, "pull_failed": True, "err": r["err"]}
        print(f"[sweep] pulled {tag} in {r['secs']:.0f}s", flush=True)
    else:
        print(f"[sweep] {tag} already present; skipping pull", flush=True)

    seeds_done = []
    for seed in SEEDS:
        print(f"[sweep] {tag}  seed={seed}", flush=True)
        cmd = [
            sys.executable, str(RUN_ONE),
            "--model", tag,
            "--track", "A",
            "--seed", str(seed),
            "--think", "off",
        ]
        t0 = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True)
        wall = time.time() - t0
        if p.returncode not in (0, 1):
            print(f"[sweep]   seed {seed} ERROR after {wall:.0f}s: {p.stderr[-400:]}", flush=True)
            continue
        run_id = p.stdout.strip().split("\n")[-1] if p.stdout.strip() else ""
        if not run_id or not (RUNS / run_id).is_dir():
            err = (p.stderr or p.stdout)[-600:]
            print(f"[sweep]   seed {seed} GENERATION_FAIL after {wall:.0f}s: {err}", flush=True)
            continue
        s = subprocess.run([sys.executable, str(SCORE), str(RUNS / run_id)],
                           capture_output=True, text=True)
        try:
            score = json.loads(s.stdout)
        except json.JSONDecodeError:
            print(f"[sweep]   seed {seed}: scoring failed: {s.stderr[:400]}", flush=True)
            continue
        seeds_done.append({"seed": seed, "run_id": run_id, "M1": score.get("M1"), "M3": score.get("M3"), "wall_s": wall})
        print(f"[sweep]   seed {seed}: M1={score.get('M1')} M3={score.get('M3')}  wall={wall:.0f}s", flush=True)

    if not keep:
        print(f"[sweep] rm {tag}", flush=True)
        remove(tag)
    print(f"[sweep] free disk after:  {free_disk_gb():.1f} GB", flush=True)
    return {"model": tag, "seeds": seeds_done}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="*", default=None,
                    help="Subset of model tags to sweep; default = all 14")
    ap.add_argument("--keep", action="store_true",
                    help="Don't `ollama rm` after testing each model")
    args = ap.parse_args()

    selected = MODELS if args.models is None else [m for m in MODELS if m[0] in args.models]
    print(f"[sweep] running {len(selected)} models, {sum(s for _, s, _ in selected):.0f} GB total downloads",
          flush=True)

    results = []
    for tag, size_gb, _kind in selected:
        results.append(sweep_model(tag, size_gb, args.keep))

    log = BENCH / "sweep_log.json"
    log.write_text(json.dumps(results, indent=2))
    print(f"\n[sweep] wrote {log}", flush=True)

    # quick summary
    print("\n[sweep] summary")
    for r in results:
        if r.get("pull_failed"):
            print(f"  PULL_FAIL  {r['model']}")
            continue
        if r.get("skipped"):
            print(f"  SKIP       {r['model']}  ({r['completed_seeds']} seeds previously scored)")
            continue
        seeds = r.get("seeds", [])
        if not seeds:
            print(f"  NO_RUNS    {r['model']}")
            continue
        m3s = [s["M3"] for s in seeds if s.get("M3") is not None]
        m1s = sum(1 for s in seeds if s.get("M1") == 1)
        m3_str = f"{sum(m3s)/len(m3s):.3f}" if m3s else "?"
        print(f"  {r['model']:<32}  M1={m1s}/{len(seeds)}  M3_mean={m3_str}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
