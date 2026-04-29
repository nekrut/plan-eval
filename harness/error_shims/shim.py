#!/usr/bin/env python3
"""
Error-injection shim. Stands in front of a real tool on PATH and either
delegates cleanly or simulates a configured failure mode.

Configured by env vars set by the harness before running run.sh:
  EVAL_INJECT_PATTERN  pattern name (one of the keys below)
  EVAL_INJECT_TARGET   tool the pattern applies to ('bwa' or 'lofreq')
  EVAL_INJECT_STATE    state dir for cross-call counters (one per run)
  EVAL_REAL_BIN_DIR    where to find the real binaries (typically the
                       conda env's bin/)

Invoked as:
  shim.py <tool_name> <real-args...>

Patterns (tool applicability noted; '-' means pattern is a no-op for
that tool and we just delegate):

  none                  no injection (baseline)
  flake_first_call      first invocation exits 1; subsequent succeed
  one_sample_fails      exits 1 when sample 'M117C1-ch' appears in args
  silent_truncation     (lofreq only) run real tool, then truncate
                         output VCF to 0 bytes
  stderr_warning_storm  emit 200 'WARNING' lines on stderr, then delegate
  slow_tool             sleep 30s, then delegate
  wrong_format_output   (lofreq only) run real tool, then strip every
                         variant line from output VCF (header-only file)
  missing_lib_error     exit 127 with 'error while loading shared libs'
"""
from __future__ import annotations
import os
import re
import sys
import time
import subprocess
from pathlib import Path

INJECT_TARGET_SAMPLE = "M117C1-ch"
SLOW_SECONDS = int(os.environ.get("EVAL_SLOW_SECONDS", "30"))
NOISE_LINES = 200
LIB_ERROR = "error while loading shared libraries: libhts.so.3: cannot open shared object file: No such file or directory"


def real_bin(tool: str) -> str:
    real_dir = os.environ.get("EVAL_REAL_BIN_DIR", "/home/anton/miniforge3/envs/bench/bin")
    return str(Path(real_dir) / tool)


def state_path(tool: str, key: str) -> Path:
    state = Path(os.environ.get("EVAL_INJECT_STATE", "/tmp/eval-shim-state"))
    state.mkdir(parents=True, exist_ok=True)
    return state / f"{tool}.{key}"


def inc_counter(tool: str, key: str) -> int:
    """Atomically bump a counter file and return the new value."""
    p = state_path(tool, key)
    n = 0
    if p.exists():
        try:
            n = int(p.read_text().strip())
        except ValueError:
            n = 0
    n += 1
    p.write_text(str(n))
    return n


def find_o_value(args: list[str]) -> str | None:
    for i, a in enumerate(args):
        if a == "-o" and i + 1 < len(args):
            return args[i + 1]
        if a.startswith("--out=") or a.startswith("-o="):
            return a.split("=", 1)[1]
        if a == "--out" and i + 1 < len(args):
            return args[i + 1]
    return None


def truncate_to_header(path: Path) -> None:
    """Keep only header (`##` and `#CHROM`) lines; drop all data lines."""
    if not path.exists():
        return
    body = path.read_text()
    kept = [ln for ln in body.splitlines() if ln.startswith("#")]
    path.write_text("\n".join(kept) + ("\n" if kept else ""))


def truncate_to_zero(path: Path) -> None:
    if path.exists():
        path.write_text("")


def run_real(tool: str, args: list[str]) -> int:
    try:
        return subprocess.run([real_bin(tool)] + args, check=False).returncode
    except FileNotFoundError as e:
        print(f"shim.py: real tool not found at {real_bin(tool)}: {e}", file=sys.stderr)
        return 127


def main() -> int:
    if len(sys.argv) < 2:
        print("shim.py: expected <tool_name> as argv[1]", file=sys.stderr)
        return 2
    tool = sys.argv[1]
    args = sys.argv[2:]

    pattern = os.environ.get("EVAL_INJECT_PATTERN", "none")
    target = os.environ.get("EVAL_INJECT_TARGET", "")

    # Pattern only applies to its target tool. Others delegate cleanly.
    if pattern != "none" and target and target != tool:
        return run_real(tool, args)

    if pattern == "none":
        return run_real(tool, args)

    if pattern == "flake_first_call":
        n = inc_counter(tool, "calls")
        if n == 1:
            print(f"{tool}: transient I/O error (eval-injected, call #1)", file=sys.stderr)
            return 1
        return run_real(tool, args)

    if pattern == "one_sample_fails":
        joined = " ".join(args)
        if INJECT_TARGET_SAMPLE in joined:
            print(f"{tool}: error processing {INJECT_TARGET_SAMPLE} (eval-injected)", file=sys.stderr)
            return 1
        return run_real(tool, args)

    if pattern == "stderr_warning_storm":
        for i in range(NOISE_LINES):
            print(f"{tool}: WARNING [{i:03d}] eval-injected stderr noise — safe to ignore", file=sys.stderr)
        return run_real(tool, args)

    if pattern == "slow_tool":
        time.sleep(SLOW_SECONDS)
        return run_real(tool, args)

    if pattern == "missing_lib_error":
        print(f"{real_bin(tool)}: {LIB_ERROR}", file=sys.stderr)
        return 127

    if pattern == "silent_truncation":
        if tool != "lofreq":
            return run_real(tool, args)
        rc = run_real(tool, args)
        out = find_o_value(args)
        if out:
            truncate_to_zero(Path(out))
        return rc

    if pattern == "wrong_format_output":
        if tool != "lofreq":
            return run_real(tool, args)
        rc = run_real(tool, args)
        out = find_o_value(args)
        if out:
            truncate_to_header(Path(out))
        return rc

    print(f"shim.py: unknown EVAL_INJECT_PATTERN={pattern!r}", file=sys.stderr)
    return run_real(tool, args)


if __name__ == "__main__":
    sys.exit(main())
