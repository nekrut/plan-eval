#!/usr/bin/env python3
"""
Aggregate all scored runs into bench/results.csv and a summary table.

Usage: python3 aggregate.py
"""
from __future__ import annotations
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev

BENCH = Path(__file__).resolve().parent.parent
RUNS = BENCH / "runs"
OUT_CSV = BENCH / "results.csv"


def parse_run_id(rid: str) -> dict:
    """qwen3.6_35b_think-on_track-A_seed-42_abcd12 -> dict"""
    m = re.match(
        r"^(?P<model>.+?)(?:_think-(?P<think>on|off))?_track-(?P<track>[AB])_seed-(?P<seed>\d+)_(?P<uuid>[0-9a-f]+)$",
        rid,
    )
    return m.groupdict() if m else {"model": rid}


def main() -> int:
    rows = []
    for d in sorted(RUNS.iterdir()):
        sj = d / "score.json"
        if not sj.exists():
            continue
        s = json.loads(sj.read_text())
        meta_path = d / "meta.json"
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        info = parse_run_id(d.name)
        rows.append({
            "run_id": d.name,
            "model": info.get("model"),
            "think": info.get("think") or "",
            "track": info.get("track"),
            "seed": int(info.get("seed", 0)),
            "M1": s.get("m1_executes"),
            "M2": s.get("m2_schema"),
            "M3": round(s.get("m3_jaccard", 0), 4),
            "M5": s.get("m5_quality"),
            "cost_usd": s.get("m4", {}).get("cost_usd"),
            "gen_seconds": round(meta.get("wall_seconds_generation") or 0, 1),
            "exec_seconds": round(s.get("m4", {}).get("wall_seconds_execution") or 0, 1),
        })

    if not rows:
        print("no scored runs found", file=sys.stderr)
        return 1

    fieldnames = list(rows[0].keys())
    with OUT_CSV.open("w") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"[aggregate] wrote {OUT_CSV} ({len(rows)} rows)", file=sys.stderr)

    # Summary by (model, think, track)
    bucket = defaultdict(list)
    for r in rows:
        key = (r["model"], r["think"], r["track"])
        bucket[key].append(r)

    fmt = "{:<28} {:<6} {:<5} {:>3} {:>6} {:>6} {:>6} {:>10} {:>8} {:>8}"
    print()
    print(fmt.format("model", "think", "track", "n", "M1", "M2", "M3", "cost_usd", "gen_s", "exec_s"))
    print("-" * 100)
    for key in sorted(bucket):
        rs = bucket[key]
        m3s = [r["M3"] for r in rs]
        m3_str = f"{mean(m3s):.3f}" + (f"±{stdev(m3s):.3f}" if len(m3s) > 1 else "")
        m1_str = f"{sum(r['M1'] for r in rs)}/{len(rs)}"
        m2_str = f"{sum(r['M2'] for r in rs)}/{len(rs)}"
        costs = [r["cost_usd"] or 0 for r in rs]
        gen_s = mean(r["gen_seconds"] for r in rs)
        exec_s = mean(r["exec_seconds"] for r in rs)
        print(fmt.format(
            (key[0] or "")[:28], key[1] or "-", key[2] or "?",
            len(rs), m1_str, m2_str, m3_str,
            f"{mean(costs):.4f}", f"{gen_s:.1f}", f"{exec_s:.1f}",
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
