#!/usr/bin/env python3
"""
RTX 5080 matrix: 11 ollama models x {v1,v2} plans x {A,B} tracks x 3 seeds,
plus /think mode on think-capable models that fit cleanly in VRAM.

Models (group A: fits cleanly in 16 GB VRAM; group B: partial CPU offload):
    A: qwen3:8b, qwen3:14b, qwen3.5:9b, gemma4:e4b, gpt-oss:20b
    B: qwen3.5:27b, qwen3.6:27b, qwen3.6:35b-a3b, gemma4:26b,
       qwen3-coder:30b, glm-4.7-flash

/think is enabled ONLY for qwen3* family models that fit in VRAM
(qwen3:8b, qwen3:14b, qwen3.5:9b). On group B models with partial
CPU offload (qwen3.5:27b, qwen3.6:27b, qwen3.6:35b-a3b, glm-4.7-flash),
empirically /think generation does not return within a 1800 s budget —
the same hardware-class problem first observed on the Jetson AGX Orin.
gpt-oss, gemma4 and qwen3-coder do not use the /think template.

For /think cells the ollama HTTP timeout is bumped from 900 to 1800 s
to accommodate longer reasoning streams.

Outputs:
  runs_5080_v1/<run_id>/   (plan = plan/PLAN_v1.md)
  runs_5080_v2/<run_id>/   (plan = plan/PLAN.md, the v2 frozen plan)

Usage:
  python3 matrix_5080.py                # run pending cells
  python3 matrix_5080.py --dry-run      # list cells, do nothing
  python3 matrix_5080.py --plans v2     # only v2
  python3 matrix_5080.py --tracks A     # only Track A
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
TRACKS = ["A", "B"]

# (model_id, supports_think) — /think only on models that fit in VRAM.
# Group B models with partial offload time out at 1800 s for /think cells.
MODELS = [
    ("qwen3:8b",            True),
    ("qwen3:14b",           True),
    ("qwen3.5:9b",          True),
    ("qwen3.5:27b",         False),
    ("qwen3.6:27b",         False),
    ("qwen3.6:35b-a3b",     False),
    ("gemma4:e4b",          False),
    ("gemma4:26b",          False),
    ("gpt-oss:20b",         False),
    ("qwen3-coder:30b",     False),
    ("glm-4.7-flash",       False),
]

# Anthropic models (no /think axis; provider routed by claude-* prefix in run_one.py).
ANTHROPIC_MODELS = [
    "claude-haiku-4-5",
    "claude-sonnet-4-6",
    "claude-opus-4-7",
]

PLAN_FILES = {
    "v1":    BENCH / "plan" / "PLAN_v1.md",
    "v2":    BENCH / "plan" / "PLAN.md",
    "v1p25": BENCH / "plan" / "PLAN_v1p25.md",
    "v1p5":  BENCH / "plan" / "PLAN_v1p5.md",
    "v1p75": BENCH / "plan" / "PLAN_v1p75.md",
    "v1g":   BENCH / "plan" / "PLAN_v1g.md",
    # v0p5: no plan file (Track B template variant); use any path, ignored.
    "v0p5":  BENCH / "plan" / "PLAN_v1.md",
}
RUNS_DIRS = {
    k: BENCH / f"runs_5080_{k}" for k in PLAN_FILES
}

# Per-experiment track restrictions:
#   v1, v2          → both tracks (existing behaviour)
#   v1p25/p5/p75    → Track A only (the cliff lives on Track A)
#   v0p5            → Track B only (it's a Track B template variant)
TRACKS_BY_EXP = {
    "v1":    ["A", "B"],
    "v2":    ["A", "B"],
    "v1p25": ["A"],
    "v1p5":  ["A"],
    "v1p75": ["A"],
    "v1g":   ["A"],
    "v0p5":  ["B"],
}

# Template overrides (only v0p5 swaps the Track B prompt template).
TRACK_TEMPLATE = {
    "v0p5": "track_b_with_order_user",
}

GEN_TIMEOUT_NO_THINK = 900
GEN_TIMEOUT_THINK = 1800


def is_anthropic(model: str) -> bool:
    return model.startswith("claude-")


def cell_id(model: str, think: str, track: str, seed: int) -> str:
    tag = model.replace("/", "_").replace(":", "_")
    if not is_anthropic(model):
        tag += f"_think-{think}"
    return f"{tag}_track-{track}_seed-{seed}"


def already_scored(runs_dir: Path, cell: str) -> Path | None:
    if not runs_dir.exists():
        return None
    for d in runs_dir.glob(f"{cell}_*"):
        if (d / "score.json").exists():
            return d
    return None


def run_cell(model: str, think: str, track: str, seed: int, plan: str) -> dict:
    runs_dir = RUNS_DIRS[plan]
    plan_path = PLAN_FILES[plan]
    cell = cell_id(model, think, track, seed)

    existing = already_scored(runs_dir, cell)
    if existing:
        s = json.loads((existing / "score.json").read_text())
        return {"cell": cell, "plan": plan, "skipped": True,
                "run_id": existing.name, "score": s}

    runs_dir.mkdir(parents=True, exist_ok=True)
    timeout = GEN_TIMEOUT_THINK if think == "on" else GEN_TIMEOUT_NO_THINK

    cmd = [
        sys.executable, str(RUN_ONE),
        "--model", model,
        "--track", track,
        "--seed", str(seed),
        "--plan", str(plan_path),
        "--runs-dir", str(runs_dir),
        "--gen-timeout", str(timeout),
    ]
    if not is_anthropic(model):
        cmd += ["--think", think]
    if plan in TRACK_TEMPLATE:
        cmd += ["--track-template", TRACK_TEMPLATE[plan]]

    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    gen_secs = time.time() - t0

    if p.returncode not in (0, 1):
        return {"cell": cell, "plan": plan, "error": p.stderr[-2000:],
                "gen_secs": gen_secs}

    run_id = p.stdout.strip().split("\n")[-1]
    run_dir = runs_dir / run_id

    s = subprocess.run(
        [sys.executable, str(SCORE), str(run_dir)],
        capture_output=True, text=True,
    )
    score = json.loads(s.stdout) if s.returncode == 0 and s.stdout.strip() \
        else {"error": s.stderr[:500]}
    return {"cell": cell, "plan": plan, "run_id": run_id,
            "score": score, "gen_secs": gen_secs}


def parse_plans_arg(s: str) -> list[str]:
    """Accept either comma-separated ('v1,v2,v1p25') or legacy concatenation
    ('v1v2'). Unknown tokens raise."""
    if "," in s:
        toks = [t.strip() for t in s.split(",") if t.strip()]
    elif s in PLAN_FILES:
        toks = [s]
    else:
        # Legacy: 'v1', 'v2', or 'v1v2'.
        toks = []
        rest = s
        for k in ("v1", "v2"):
            if rest.startswith(k):
                toks.append(k)
                rest = rest[len(k):]
        if rest:
            toks.append(rest)  # let the validity check below catch it
    bad = [t for t in toks if t not in PLAN_FILES]
    if bad:
        raise SystemExit(f"unknown plan(s): {bad}; valid: {sorted(PLAN_FILES)}")
    return toks


def build_cells(args) -> list[tuple[str, str, str, int, str]]:
    plans = parse_plans_arg(args.plans)
    requested_tracks = set(args.tracks.upper())
    cells = []
    for plan in plans:
        allowed_tracks = [t for t in TRACKS_BY_EXP[plan] if t in requested_tracks]
        # Local ollama models
        if not args.only_anthropic:
            for model, supports_think in MODELS:
                think_modes = ["off", "on"] if (supports_think and not args.no_think) else ["off"]
                for think in think_modes:
                    for track in allowed_tracks:
                        for seed in SEEDS:
                            cells.append((model, think, track, seed, plan))
        # Anthropic models (no think axis)
        if args.include_anthropic or args.only_anthropic:
            for model in ANTHROPIC_MODELS:
                for track in allowed_tracks:
                    for seed in SEEDS:
                        cells.append((model, "off", track, seed, plan))
    return cells


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--plans", default="v1v2",
                    help="Subset, e.g. 'v2' or 'v1v2' (default both)")
    ap.add_argument("--tracks", default="AB",
                    help="Subset, e.g. 'A' or 'AB' (default both)")
    ap.add_argument("--no-think", action="store_true",
                    help="Skip /think condition entirely")
    ap.add_argument("--include-anthropic", action="store_true",
                    help="Also iterate Anthropic models (claude-haiku/sonnet/opus)")
    ap.add_argument("--only-anthropic", action="store_true",
                    help="Only iterate Anthropic models, skip local ollama")
    ap.add_argument("--models", default="",
                    help="Comma-separated subset of model ids; default = all")
    args = ap.parse_args()

    cells = build_cells(args)
    if args.models:
        wanted = set(args.models.split(","))
        cells = [c for c in cells if c[0] in wanted]

    if args.dry_run:
        print(f"[matrix_5080] {len(cells)} cells")
        for c in cells:
            print(f"  plan={c[4]}  {cell_id(c[0], c[1], c[2], c[3])}")
        return 0

    print(f"[matrix_5080] {len(cells)} cells, serialized", file=sys.stderr)
    results = []
    t_start = time.time()

    for i, (model, think, track, seed, plan) in enumerate(cells, 1):
        r = run_cell(model, think, track, seed, plan)
        results.append(r)
        elapsed = time.time() - t_start
        m3 = r.get("score", {}).get("M3", "?")
        gen = r.get("gen_secs", 0)
        tag = "skip" if r.get("skipped") else \
              "ok"   if isinstance(m3, (int, float)) else "err"
        print(f"  [{i}/{len(cells)}] [{tag}] plan={plan} "
              f"{r['cell']}: M3={m3}  gen={gen:.1f}s  total_elapsed={elapsed:.0f}s",
              file=sys.stderr)
        # persist after every cell so a crash mid-matrix doesn't lose progress
        (BENCH / "matrix_5080_log.json").write_text(json.dumps(results, indent=2))

    print(f"[matrix_5080] done in {time.time() - t_start:.0f}s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
