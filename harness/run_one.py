#!/usr/bin/env python3
"""
Run one benchmark cell: (model_id, track, seed) -> runs/<run_id>/

Steps:
  1. Build the user message from prompts/<track>_user.tmpl, substituting
     {TOOL_INVENTORY} and {PLAN}.
  2. Generate the script via the model (claude CLI for Anthropic; ollama HTTP
     for qwen3.6).
  3. Set up sandbox: runs/<run_id>/ with data/ symlinked, results/ empty.
  4. Execute the generated script with a 600 s wall-clock budget.
  5. Persist script.sh, usage.json, meta.json, exec.json under the run dir.

Scoring is NOT done here; see score/score_run.py.

Usage:
  python3 run_one.py --model claude-haiku-4-5 --track A --seed 42
  python3 run_one.py --model qwen3.6:35b --think on --track A --seed 42
"""
from __future__ import annotations
import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
import urllib.request
import uuid
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
PROMPTS = BENCH / "prompts"
PLAN_FILE = BENCH / "plan" / "PLAN.md"
RUNS = BENCH / "runs"


def tool_inventory() -> str:
    p = subprocess.run(
        [str(BENCH / "setup" / "verify_env.sh")],
        capture_output=True, text=True, check=True,
    )
    out = p.stdout
    m = re.search(r"^TOOL_INVENTORY.*?^OK$", out, re.S | re.M)
    if not m:
        raise RuntimeError("verify_env.sh did not return expected block")
    return m.group(0).rsplit("\n", 1)[0]  # drop trailing "OK"


def build_user_message(track: str, plan_path: Path) -> str:
    tmpl = (PROMPTS / f"track_{track.lower()}_user.tmpl").read_text()
    inv = tool_inventory()
    out = tmpl.replace("{TOOL_INVENTORY}", inv)
    if "{PLAN}" in out:
        out = out.replace("{PLAN}", plan_path.read_text())
    return out


def strip_fences(s: str) -> str:
    s = s.strip()
    m = re.match(r"^```(?:bash|sh)?\s*\n(.*?)\n```\s*$", s, re.S)
    if m:
        return m.group(1).strip()
    return s


def call_claude(model: str, system_text: str, user_text: str) -> dict:
    cmd = [
        "claude", "-p",
        "--model", model,
        "--system-prompt", system_text,
        "--output-format", "json",
        "--no-session-persistence",
        "--disallowedTools", "*",
        "--max-budget-usd", "1",
    ]
    t0 = time.time()
    p = subprocess.run(cmd, input=user_text, capture_output=True, text=True)
    elapsed = time.time() - t0
    if p.returncode != 0:
        raise RuntimeError(f"claude failed (exit {p.returncode}):\n{p.stderr[:2000]}")
    try:
        d = json.loads(p.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(f"claude returned non-JSON:\n{p.stdout[:2000]}")
    return {
        "provider": "anthropic",
        "model": model,
        "script": strip_fences(d.get("result", "")),
        "usage": d.get("usage", {}),
        "cost_usd": d.get("total_cost_usd", d.get("cost_usd")),
        "duration_ms": d.get("duration_ms"),
        "wall_seconds": elapsed,
        "raw_response_preview": d.get("result", "")[:400],
    }


def call_ollama(model: str, think: bool, system_text: str, user_text: str, seed: int, gen_timeout: int = 900) -> dict:
    # /no_think is a Qwen-family control token; other model families use the
    # `think` payload field instead and treat the prefix as literal user text.
    if not think and model.lower().startswith("qwen"):
        user_text = "/no_think\n" + user_text
    payload = {
        "model": model,
        "stream": False,
        "messages": [
            {"role": "system", "content": system_text},
            {"role": "user", "content": user_text},
        ],
        "think": think,
        "options": {
            "temperature": 0.2,
            "seed": seed,
            "num_predict": 16384,
            "num_ctx": 16384,
        },
    }
    req = urllib.request.Request(
        "http://localhost:11434/api/chat",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
    )
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=gen_timeout) as r:
        d = json.loads(r.read())
    elapsed = time.time() - t0
    content = d.get("message", {}).get("content", "")
    return {
        "provider": "ollama",
        "model": model,
        "think": think,
        "script": strip_fences(content),
        "usage": {
            "prompt_eval_count": d.get("prompt_eval_count"),
            "eval_count": d.get("eval_count"),
            "prompt_eval_duration_ns": d.get("prompt_eval_duration"),
            "eval_duration_ns": d.get("eval_duration"),
            "total_duration_ns": d.get("total_duration"),
        },
        "cost_usd": 0.0,
        "duration_ms": int(elapsed * 1000),
        "wall_seconds": elapsed,
        "raw_response_preview": content[:400],
    }


def setup_sandbox(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "data").symlink_to(BENCH / "data")
    (run_dir / "results").mkdir(exist_ok=True)


def execute(run_dir: Path, script: str, budget_s: int = 600) -> dict:
    script_path = run_dir / "run.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)
    log_path = run_dir / "exec.log"
    t0 = time.time()
    activate = (
        "source $HOME/miniforge3/etc/profile.d/conda.sh && "
        "conda activate bench && exec bash run.sh"
    )
    try:
        with log_path.open("wb") as logf:
            p = subprocess.run(
                ["bash", "-c", activate],
                cwd=str(run_dir),
                stdout=logf, stderr=subprocess.STDOUT,
                timeout=budget_s,
            )
        return {
            "exit_code": p.returncode,
            "wall_seconds": time.time() - t0,
            "timed_out": False,
        }
    except subprocess.TimeoutExpired:
        return {"exit_code": -1, "wall_seconds": time.time() - t0, "timed_out": True}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--track", choices=["A", "B"], required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--think", choices=["on", "off"], default="on",
                    help="Only meaningful for ollama models")
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--plan", default=str(PLAN_FILE),
                    help="Path to PLAN.md (defaults to plan/PLAN.md)")
    ap.add_argument("--runs-dir", default=str(RUNS),
                    help="Directory to write run output (default runs/)")
    ap.add_argument("--gen-timeout", type=int, default=900,
                    help="Ollama HTTP timeout in seconds (default 900)")
    args = ap.parse_args()

    runs_root = Path(args.runs_dir)
    plan_path = Path(args.plan)

    # Anthropic model ids start with "claude-"; everything else is local (ollama).
    is_ollama = not args.model.startswith("claude-")
    cell_tag = args.model.replace("/", "_").replace(":", "_")
    if is_ollama:
        cell_tag += f"_think-{args.think}"
    run_id = args.run_id or f"{cell_tag}_track-{args.track}_seed-{args.seed}_{uuid.uuid4().hex[:6]}"
    run_dir = runs_root / run_id

    if run_dir.exists():
        print(f"[run_one] {run_dir} exists — removing for clean run", file=sys.stderr)
        shutil.rmtree(run_dir)

    system_text = (PROMPTS / "system.txt").read_text()
    user_text = build_user_message(args.track, plan_path)

    print(f"[run_one] generating script via {args.model} (track {args.track}, seed {args.seed})", file=sys.stderr)
    if is_ollama:
        gen = call_ollama(args.model, args.think == "on", system_text, user_text, args.seed, gen_timeout=args.gen_timeout)
    else:
        gen = call_claude(args.model, system_text, user_text)

    setup_sandbox(run_dir)

    meta = {
        "run_id": run_id,
        "model": args.model,
        "track": args.track,
        "seed": args.seed,
        "think": args.think if is_ollama else None,
        "provider": gen["provider"],
        "wall_seconds_generation": gen["wall_seconds"],
        "duration_ms_generation": gen["duration_ms"],
    }
    (run_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    (run_dir / "usage.json").write_text(json.dumps({
        "usage": gen["usage"],
        "cost_usd": gen["cost_usd"],
    }, indent=2))
    (run_dir / "raw_response.txt").write_text(gen["raw_response_preview"])

    if not gen["script"]:
        print("[run_one] EMPTY script returned; aborting before exec", file=sys.stderr)
        (run_dir / "exec.json").write_text(json.dumps(
            {"exit_code": None, "skipped": "empty_script"}, indent=2))
        return 2

    print(f"[run_one] executing run.sh in {run_dir} (budget 600s)", file=sys.stderr)
    exec_res = execute(run_dir, gen["script"], budget_s=600)
    (run_dir / "exec.json").write_text(json.dumps(exec_res, indent=2))

    print(f"[run_one] done: {run_dir}  exit={exec_res['exit_code']}  exec_wall={exec_res['wall_seconds']:.1f}s", file=sys.stderr)
    print(run_id)  # stdout = run_id for matrix.py
    return 0 if exec_res["exit_code"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
