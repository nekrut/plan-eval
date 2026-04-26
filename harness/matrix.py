#!/usr/bin/env python3
"""
Run the full benchmark matrix:
  Anthropic models: claude-opus-4-7, claude-sonnet-4-6, claude-haiku-4-5
  Local model:      qwen3.6:35b in /think and /no_think modes
  Tracks:           A (with plan), B (no plan)
  Seeds:            42, 43, 44

Anthropic cells run in parallel (3 concurrent). Ollama cells serialized.
After each cell, score it. Skip cells whose run dir already has score.json.

Usage:
  python3 matrix.py                     # run everything pending
  python3 matrix.py --only-anthropic    # skip ollama
  python3 matrix.py --only-ollama       # skip anthropic
  python3 matrix.py --dry-run           # list cells, do nothing
"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import json
import subprocess
import sys
import time
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
RUNS = BENCH / "runs"
RUN_ONE = BENCH / "harness" / "run_one.py"
SCORE = BENCH / "score" / "score_run.py"

SEEDS = [42, 43, 44]
ANTHROPIC = ["claude-opus-4-7", "claude-sonnet-4-6", "claude-haiku-4-5"]
OLLAMA = ["qwen3.6:35b"]
TRACKS = ["A", "B"]
THINK_MODES = ["on", "off"]


def cell_id(model: str, track: str, seed: int, think: str | None) -> str:
    tag = model.replace("/", "_").replace(":", "_")
    if think:
        tag += f"_think-{think}"
    return f"{tag}_track-{track}_seed-{seed}"


def already_scored(cell: str) -> Path | None:
    """Return path of an existing scored run dir matching this cell, if any."""
    for d in RUNS.glob(f"{cell}_*"):
        if (d / "score.json").exists():
            return d
    return None


def run_cell(model: str, track: str, seed: int, think: str | None) -> dict:
    cell = cell_id(model, track, seed, think)
    existing = already_scored(cell)
    if existing:
        s = json.loads((existing / "score.json").read_text())
        return {"cell": cell, "skipped": True, "run_id": existing.name, "score": s}

    cmd = [
        sys.executable, str(RUN_ONE),
        "--model", model,
        "--track", track,
        "--seed", str(seed),
    ]
    if think:
        cmd += ["--think", think]

    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    gen_secs = time.time() - t0
    if p.returncode not in (0, 1):  # 0 = ok, 1 = exec failed (still scoreable)
        return {"cell": cell, "error": p.stderr[-2000:], "gen_secs": gen_secs}

    run_id = p.stdout.strip().split("\n")[-1]
    run_dir = RUNS / run_id

    s = subprocess.run(
        [sys.executable, str(SCORE), str(run_dir)],
        capture_output=True, text=True,
    )
    score = json.loads(s.stdout) if s.returncode == 0 and s.stdout.strip() else {"error": s.stderr[:500]}
    return {"cell": cell, "run_id": run_id, "score": score, "gen_secs": gen_secs}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-anthropic", action="store_true")
    ap.add_argument("--only-ollama", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--anthropic-concurrency", type=int, default=3)
    ap.add_argument("--tracks", default="AB", help="Subset of tracks to run, e.g. 'A' or 'AB' (default)")
    ap.add_argument("--no-think", action="store_true", help="For ollama: skip the /think condition")
    args = ap.parse_args()

    tracks = [t for t in TRACKS if t in args.tracks.upper()]
    think_modes = ["off"] if args.no_think else THINK_MODES

    cells_anth = [
        (m, t, s, None)
        for m in ANTHROPIC for t in tracks for s in SEEDS
    ]
    cells_oll = [
        (m, t, s, th)
        for m in OLLAMA for t in tracks for th in think_modes for s in SEEDS
    ]

    do_anth = not args.only_ollama
    do_oll = not args.only_anthropic

    if args.dry_run:
        if do_anth:
            print(f"[anthropic] {len(cells_anth)} cells")
            for c in cells_anth:
                print("  " + cell_id(*c))
        if do_oll:
            print(f"[ollama] {len(cells_oll)} cells")
            for c in cells_oll:
                print("  " + cell_id(*c))
        return 0

    results = []

    if do_anth:
        print(f"[matrix] Anthropic: {len(cells_anth)} cells, concurrency={args.anthropic_concurrency}", file=sys.stderr)
        with cf.ThreadPoolExecutor(max_workers=args.anthropic_concurrency) as ex:
            futs = {ex.submit(run_cell, *c): c for c in cells_anth}
            for fut in cf.as_completed(futs):
                r = fut.result()
                results.append(r)
                tag = "skip" if r.get("skipped") else "ok" if "score" in r and r["score"].get("M1") is not None else "err"
                print(f"  [{tag}] {r['cell']}: {r.get('score', {}).get('M3', '?')}", file=sys.stderr)

    if do_oll:
        print(f"[matrix] Ollama: {len(cells_oll)} cells, serialized", file=sys.stderr)
        for c in cells_oll:
            r = run_cell(*c)
            results.append(r)
            tag = "skip" if r.get("skipped") else "ok" if "score" in r and r["score"].get("M1") is not None else "err"
            print(f"  [{tag}] {r['cell']}: M3={r.get('score', {}).get('M3', '?')}  gen={r.get('gen_secs', 0):.1f}s", file=sys.stderr)

    out_path = BENCH / "matrix_log.json"
    out_path.write_text(json.dumps(results, indent=2))
    print(f"[matrix] wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
