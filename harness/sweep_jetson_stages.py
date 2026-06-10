#!/usr/bin/env python3
"""
Run a model through the remaining Figure-1 recipe stages ON THE JETSON, so a
model can have a complete headline-heatmap row from one hardware platform
(mirrors harness/matrix_5080.py, which did this on the RTX 5080).

The v1 and v2 columns are already produced by runs_v1/ and runs/ (sweep_local).
This driver fills the rest: v1.25, v1.5, v1g (Track A), v0.5 (Track B template),
and B (Track B, no plan). Each stage writes to its own runs_jetson_<k>/ dir,
which scripts/make_figures.py PLAN_MAP attributes to the jetson platform.

Track-A vs Track-B and the v0.5 template override match matrix_5080.py exactly.

Usage:
  python3 harness/sweep_jetson_stages.py                       # gemma4:12b
  python3 harness/sweep_jetson_stages.py --model qwen3.6:27b
  python3 harness/sweep_jetson_stages.py --stages v1p25 v1p5   # subset
"""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
RUN_ONE = BENCH / "harness" / "run_one.py"
SCORE = BENCH / "score" / "score_run.py"
SEEDS = [42, 43, 44]

# (stage key, plan file, track, track-template override). Stage key -> runs_jetson_<key>/.
# Mirrors matrix_5080.py PLAN_FILES / TRACKS_BY_EXP / TRACK_TEMPLATE.
STAGES = {
    "v1p25": ("PLAN_v1p25.md", "A", None),
    "v1p5":  ("PLAN_v1p5.md",  "A", None),
    "v1g":   ("PLAN_v1g.md",   "A", None),
    "v0p5":  ("PLAN_v1.md",    "B", "track_b_with_order_user"),  # v0.5 condition
    "b":     ("PLAN.md",       "B", None),                       # no-plan "B" column
}


def already_scored(runs_dir: Path, model: str) -> int:
    safe = model.replace("/", "_").replace(":", "_")
    return sum(1 for d in runs_dir.glob(f"{safe}_*track-*_seed-*") if (d / "score.json").exists())


def run_stage(model: str, key: str) -> dict:
    plan_file, track, template = STAGES[key]
    runs_dir = BENCH / f"runs_jetson_{key}"
    runs_dir.mkdir(exist_ok=True)
    print(f"\n[stage] === {model}  {key}  (plan={plan_file}, track={track}"
          f"{', tmpl='+template if template else ''}) ===", flush=True)

    done = already_scored(runs_dir, model)
    if done >= len(SEEDS):
        print(f"[stage] {key}: {done}/{len(SEEDS)} already scored — skipping", flush=True)
        return {"stage": key, "skipped": True}

    seeds_done = []
    for seed in SEEDS:
        cmd = [
            sys.executable, str(RUN_ONE),
            "--model", model, "--track", track, "--seed", str(seed),
            "--think", "off",
            "--plan", str(BENCH / "plan" / plan_file),
            "--runs-dir", str(runs_dir),
        ]
        if template:
            cmd += ["--track-template", template]
        t0 = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True)
        wall = time.time() - t0
        rid = p.stdout.strip().split("\n")[-1] if p.stdout.strip() else ""
        if not rid or not (runs_dir / rid).is_dir():
            print(f"[stage]   seed {seed} GEN_FAIL after {wall:.0f}s: {(p.stderr or p.stdout)[-300:]}", flush=True)
            continue
        # score_run.py writes score.json itself; capture stdout only for the summary line.
        s = subprocess.run([sys.executable, str(SCORE), str(runs_dir / rid)],
                           capture_output=True, text=True)
        try:
            sc = json.loads(s.stdout)
        except json.JSONDecodeError:
            print(f"[stage]   seed {seed}: scoring failed: {s.stderr[:300]}", flush=True)
            continue
        seeds_done.append({"seed": seed, "M1": sc.get("M1"), "M3": sc.get("M3"), "wall_s": wall})
        print(f"[stage]   seed {seed}: M1={sc.get('M1')} M3={sc.get('M3')}  wall={wall:.0f}s", flush=True)
    return {"stage": key, "seeds": seeds_done}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="gemma4:12b")
    ap.add_argument("--stages", nargs="*", default=list(STAGES.keys()),
                    help=f"subset of {list(STAGES.keys())}")
    args = ap.parse_args()

    bad = [s for s in args.stages if s not in STAGES]
    if bad:
        raise SystemExit(f"unknown stage(s): {bad}; valid: {list(STAGES.keys())}")

    results = [run_stage(args.model, k) for k in args.stages]

    print("\n[stage] summary")
    for r in results:
        if r.get("skipped"):
            print(f"  SKIP  {r['stage']}")
            continue
        seeds = r.get("seeds", [])
        if not seeds:
            print(f"  NO_RUNS  {r['stage']}")
            continue
        m3s = [s["M3"] for s in seeds if s.get("M3") is not None]
        m1 = sum(1 for s in seeds if s.get("M1") == 1)
        m3_str = f"{sum(m3s)/len(m3s):.3f}" if m3s else "?"
        print(f"  {r['stage']:<6}  M1={m1}/{len(seeds)}  M3_mean={m3_str}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
