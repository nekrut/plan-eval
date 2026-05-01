#!/usr/bin/env python3
"""
Drive the error-handling experiment matrix:

  models × patterns × recipe × seed -> runs_inject/<run_id>/

Patterns and target tools (which patterns are valid for which tools):

  flake_first_call         bwa, lofreq
  one_sample_fails         bwa, lofreq
  silent_truncation        lofreq          (output-mutating, lofreq-only)
  stderr_warning_storm     bwa, lofreq
  slow_tool                bwa, lofreq
  wrong_format_output      lofreq          (output-mutating, lofreq-only)
  missing_lib_error        bwa, lofreq

For each cell we run run_one.py with --inject and --inject-target, then
score it. Each run produces meta.json (with inject_pattern/target),
exec.json, score.json, etc.

Usage:
  python3 harness/error_matrix.py
  python3 harness/error_matrix.py --models claude-haiku-4-5 granite4
  python3 harness/error_matrix.py --patterns silent_truncation flake_first_call
  python3 harness/error_matrix.py --plans v2 v2_defensive
  python3 harness/error_matrix.py --dry-run
"""
from __future__ import annotations
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
RUNS = BENCH / "runs_inject"
RUN_ONE = BENCH / "harness" / "run_one.py"
SCORE = BENCH / "score" / "score_run.py"

DEFAULT_MODELS = [
    "claude-opus-4-7",
    "claude-sonnet-4-6",
    "claude-haiku-4-5",
    "qwen3.6:35b",
    "granite4",
]

DEFAULT_PLANS = ["v2", "v2_defensive"]

PLAN_FILES = {
    "v2":            BENCH / "plan" / "PLAN.md",
    "v2_defensive":  BENCH / "plan" / "PLAN_v2_defensive.md",
}

# (pattern, tool) cells. Output-mutating patterns are lofreq-only.
PATTERN_TOOLS = {
    "flake_first_call":     ["bwa", "lofreq"],
    "one_sample_fails":     ["bwa", "lofreq"],
    "silent_truncation":    ["lofreq"],
    "stderr_warning_storm": ["bwa", "lofreq"],
    "slow_tool":            ["bwa", "lofreq"],
    "wrong_format_output":  ["lofreq"],
    "missing_lib_error":    ["bwa", "lofreq"],
}

DEFAULT_PATTERNS = list(PATTERN_TOOLS.keys())

SEEDS = [42, 43, 44]

# How long to wait between retries when the Claude Code CLI hits its usage cap.
# Signature: exit 1 with empty stderr → usage limit, not a model/network error.
# Retries continue indefinitely so an unattended run always completes eventually.
_RATE_LIMIT_WAIT_S = 1800  # 30 min per retry


def _is_claude_rate_limit(r: dict) -> bool:
    """Return True when a cell failed because the claude CLI exhausted its usage cap."""
    err = r.get("error", "")
    return "claude failed (exit 1):" in err and err.rstrip().endswith("exit 1):")


def run_cell(model: str, plan: str, pattern: str, target: str, seed: int, dry: bool) -> dict:
    plan_path = PLAN_FILES[plan]
    is_ollama = not model.startswith("claude-")
    cmd = [
        sys.executable, str(RUN_ONE),
        "--model", model,
        "--track", "A",
        "--seed", str(seed),
        "--plan", str(plan_path),
        "--runs-dir", str(RUNS),
        "--inject", pattern,
        "--inject-target", target,
    ]
    if is_ollama:
        cmd += ["--think", "off"]

    cell = f"{model}/{plan}/{pattern}@{target}/seed-{seed}"
    if dry:
        print(f"[DRY] {cell}: {' '.join(cmd)}")
        return {"cell": cell, "skipped": True}

    print(f"[matrix] {cell}", flush=True)
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - t0
    if p.returncode not in (0, 1):
        # generation/exec hard failure (not a non-zero from run.sh)
        return {"cell": cell, "error": p.stderr[-400:], "wall_s": wall}
    run_id = p.stdout.strip().split("\n")[-1] if p.stdout.strip() else ""
    if not run_id or not (RUNS / run_id).is_dir():
        return {"cell": cell, "error": (p.stderr or p.stdout)[-400:], "wall_s": wall}

    # Score it
    s = subprocess.run([sys.executable, str(SCORE), str(RUNS / run_id)],
                       capture_output=True, text=True)
    try:
        score = json.loads(s.stdout)
    except json.JSONDecodeError:
        return {"cell": cell, "run_id": run_id, "score_failed": s.stderr[:300], "wall_s": wall}

    eh = json.loads((RUNS / run_id / "score.json").read_text()).get("error_handling", {})
    return {
        "cell": cell,
        "run_id": run_id,
        "M1": score.get("M1"),
        "M3": score.get("M3"),
        "m_handle": eh.get("m_handle"),
        "m_recover": eh.get("m_recover"),
        "m_diagnose": eh.get("m_diagnose"),
        "n_valid": eh.get("n_valid"),
        "wall_s": wall,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="*", default=DEFAULT_MODELS)
    ap.add_argument("--patterns", nargs="*", default=DEFAULT_PATTERNS)
    ap.add_argument("--plans", nargs="*", default=DEFAULT_PLANS,
                    choices=list(PLAN_FILES.keys()))
    ap.add_argument("--seeds", nargs="*", type=int, default=SEEDS)
    ap.add_argument("--include-baseline", action="store_true",
                    help="Also run a 'none' (no injection) baseline for each model×plan")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--log", default=str(BENCH / "error_matrix_log.jsonl"))
    args = ap.parse_args()

    RUNS.mkdir(exist_ok=True)
    log_path = Path(args.log)
    cells = []

    for model in args.models:
        for plan in args.plans:
            if args.include_baseline:
                for seed in args.seeds:
                    cells.append((model, plan, "none", "lofreq", seed))
            for pat in args.patterns:
                for tool in PATTERN_TOOLS[pat]:
                    for seed in args.seeds:
                        cells.append((model, plan, pat, tool, seed))

    print(f"[matrix] {len(cells)} cells planned ({len(args.models)} models × {len(args.plans)} plans × {sum(len(PATTERN_TOOLS[p]) for p in args.patterns)} pattern-tool × {len(args.seeds)} seeds"
          f"{' + baselines' if args.include_baseline else ''})", flush=True)

    if args.dry_run:
        for m, p, pat, t, s in cells[:20]:
            print(f"  {m} {p} {pat}@{t} seed={s}")
        if len(cells) > 20:
            print(f"  ... +{len(cells)-20} more")
        return 0

    with log_path.open("a") as logf:
        for model, plan, pat, tool, seed in cells:
            is_ollama = not model.startswith("claude-")
            attempt = 0
            while True:
                r = run_cell(model, plan, pat, tool, seed, dry=False)
                if not is_ollama and _is_claude_rate_limit(r):
                    attempt += 1
                    print(
                        f"  [rate-limit] claude usage cap hit (attempt {attempt}) — "
                        f"pausing {_RATE_LIMIT_WAIT_S // 60} min...",
                        flush=True,
                    )
                    time.sleep(_RATE_LIMIT_WAIT_S)
                else:
                    break
            logf.write(json.dumps(r) + "\n")
            logf.flush()
            print(f"  -> {r.get('m_handle','-')} M3={r.get('M3','-')} valid={r.get('n_valid','-')} wall={r.get('wall_s',0):.0f}s",
                  flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
