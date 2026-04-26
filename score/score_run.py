#!/usr/bin/env python3
"""
Score one benchmark run against bench/ground_truth/results/.

Metrics:
  M1 executes:        bash run.sh exited 0 within 600 s        (binary)
  M2 schema:          all expected output paths present + bcftools header parses (binary)
  M3 jaccard:         per-sample tolerant Jaccard (CHROM,POS,REF,ALT) with AF tol ±0.02,
                      PASS only, macro-mean across 4 samples   (continuous, primary)
  M4 cost/latency:    wall_seconds (gen + exec), tokens, USD   (continuous, reported)
  M5 code quality:    shellcheck pass, set -euo pipefail, no /home/anton, idempotent (binary checks)

Usage: python3 score_run.py runs/<run_id>
"""
from __future__ import annotations
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parent.parent
GT = BENCH / "ground_truth" / "results"
SAMPLES = ["M117-bl", "M117-ch", "M117C1-bl", "M117C1-ch"]
EXPECTED = (
    [f"{s}.bam" for s in SAMPLES]
    + [f"{s}.bam.bai" for s in SAMPLES]
    + [f"{s}.vcf.gz" for s in SAMPLES]
    + [f"{s}.vcf.gz.tbi" for s in SAMPLES]
    + ["collapsed.tsv"]
)
AF_TOL = 0.02


def conda_run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    """Run a command inside the bench conda env."""
    full = f"source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate bench && {' '.join(cmd)}"
    return subprocess.run(["bash", "-lc", full], capture_output=True, text=True, **kw)


def m1_executes(run_dir: Path) -> tuple[int, dict]:
    e = json.loads((run_dir / "exec.json").read_text())
    ok = int(e.get("exit_code") == 0 and not e.get("timed_out"))
    return ok, {"exit_code": e.get("exit_code"), "timed_out": e.get("timed_out")}


def m2_schema(run_dir: Path) -> tuple[int, dict]:
    res = run_dir / "results"
    missing = [f for f in EXPECTED if not (res / f).exists()]
    if missing:
        return 0, {"missing": missing}
    # bcftools must accept every VCF
    for s in SAMPLES:
        p = conda_run(["bcftools", "view", "-h", str(res / f"{s}.vcf.gz")])
        if p.returncode != 0:
            return 0, {"bcftools_header_failed": s, "stderr": p.stderr[:300]}
    return 1, {}


def parse_variants(vcf: Path) -> list[tuple]:
    """Returns list of (chrom, pos, ref, alt, af) for PASS records."""
    p = conda_run([
        "bcftools", "query", "-f",
        r"'%FILTER\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n'",
        str(vcf),
    ])
    if p.returncode != 0:
        return []
    out = []
    for line in p.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        flt, chrom, pos, ref, alt, af = parts[:6]
        if flt not in ("PASS", "."):  # "." = no filter applied = treat as pass
            continue
        try:
            af_f = float(af) if af not in (".", "") else 0.0
        except ValueError:
            af_f = 0.0
        out.append((chrom, int(pos), ref, alt, af_f))
    return out


def tolerant_jaccard(model_vars: list[tuple], canon_vars: list[tuple]) -> float:
    """Match on (chrom,pos,ref,alt); count as match if AF within ±AF_TOL."""
    canon_by_key = {(c, p, r, a): af for (c, p, r, a, af) in canon_vars}
    model_by_key = {(c, p, r, a): af for (c, p, r, a, af) in model_vars}
    keys = set(canon_by_key) | set(model_by_key)
    if not keys:
        return 1.0
    matches = 0
    for k in keys:
        if k in canon_by_key and k in model_by_key:
            if abs(canon_by_key[k] - model_by_key[k]) <= AF_TOL:
                matches += 1
    return matches / len(keys)


def m3_jaccard(run_dir: Path) -> tuple[float, dict]:
    per_sample = {}
    for s in SAMPLES:
        model_vcf = run_dir / "results" / f"{s}.vcf.gz"
        canon_vcf = GT / f"{s}.vcf.gz"
        if not model_vcf.exists():
            per_sample[s] = 0.0
            continue
        m = parse_variants(model_vcf)
        c = parse_variants(canon_vcf)
        per_sample[s] = tolerant_jaccard(m, c)
    macro = sum(per_sample.values()) / len(per_sample) if per_sample else 0.0
    return macro, {"per_sample": per_sample}


def m4_costlatency(run_dir: Path) -> dict:
    meta = json.loads((run_dir / "meta.json").read_text())
    usage = json.loads((run_dir / "usage.json").read_text())
    exec_ = json.loads((run_dir / "exec.json").read_text())
    return {
        "wall_seconds_generation": meta.get("wall_seconds_generation"),
        "wall_seconds_execution": exec_.get("wall_seconds"),
        "cost_usd": usage.get("cost_usd"),
        "usage": usage.get("usage"),
    }


def m5_quality(run_dir: Path) -> tuple[int, dict]:
    script = (run_dir / "run.sh").read_text() if (run_dir / "run.sh").exists() else ""
    if not script:
        return 0, {"reason": "no script"}
    checks = {}
    checks["set_euo_pipefail"] = bool(re.search(r"^\s*set\s+-euo\s+pipefail", script, re.M))
    checks["no_home_anton"] = "/home/anton" not in script
    sc = conda_run(["shellcheck", "-S", "warning", str(run_dir / "run.sh")])
    checks["shellcheck_clean"] = sc.returncode == 0
    checks["shellcheck_output"] = sc.stdout[:600] if sc.stdout else ""
    # Idempotency: only test if M1 already passed
    e = json.loads((run_dir / "exec.json").read_text())
    if e.get("exit_code") == 0:
        activate = (
            "source $HOME/miniforge3/etc/profile.d/conda.sh && "
            "conda activate bench && exec bash run.sh"
        )
        rerun = subprocess.run(
            ["bash", "-c", activate], cwd=str(run_dir),
            capture_output=True, text=True, timeout=300,
        )
        checks["idempotent"] = rerun.returncode == 0
        checks["idempotent_stderr"] = rerun.stderr[:300] if rerun.returncode != 0 else ""
    else:
        checks["idempotent"] = None
    score = int(
        checks["set_euo_pipefail"] and checks["no_home_anton"] and checks["shellcheck_clean"]
        and (checks["idempotent"] is True)
    )
    return score, checks


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    args = ap.parse_args()
    run_dir = Path(args.run_dir).resolve()

    out = {"run_id": run_dir.name}
    out["m1_executes"], out["m1_detail"] = m1_executes(run_dir)
    if out["m1_executes"]:
        out["m2_schema"], out["m2_detail"] = m2_schema(run_dir)
    else:
        out["m2_schema"], out["m2_detail"] = 0, {"reason": "M1 failed"}
    if out["m2_schema"]:
        out["m3_jaccard"], out["m3_detail"] = m3_jaccard(run_dir)
    else:
        out["m3_jaccard"], out["m3_detail"] = 0.0, {"reason": "M2 failed"}
    out["m4"] = m4_costlatency(run_dir)
    out["m5_quality"], out["m5_detail"] = m5_quality(run_dir)

    score_path = run_dir / "score.json"
    score_path.write_text(json.dumps(out, indent=2))

    print(json.dumps({
        "run_id": out["run_id"],
        "M1": out["m1_executes"],
        "M2": out["m2_schema"],
        "M3": round(out["m3_jaccard"], 3),
        "M5": out["m5_quality"],
        "cost_usd": out["m4"]["cost_usd"],
        "exec_s": round(out["m4"]["wall_seconds_execution"] or 0, 1),
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
