#!/usr/bin/env bash
#
# Reproduce the gemma4:12b capability result (README §2.1) on a MacBook Air M4.
#
# Expected outcome (matches the Jetson run committed in 3d64220):
#   gemma4:12b  M1=3/3  M3=1.000  on the v2 plan, track A, think off, 3 seeds.
#
# This is the *capability* sweep (sweep_local.py), not the §2.6 error matrix.
# It is non-destructive: it writes per-run artifacts under runs/ and prints the
# score; it does NOT touch the committed results.csv.
#
# Prereqs (see the rest of RUN_ON_MACOS.md for the one-time setup):
#   - conda env `bench` with bwa/samtools/bcftools/lofreq  (bash setup/install.sh)
#   - data fetched                                          (bash setup/fetch_data.sh)
#   - ground truth built                                    (bash ground_truth/canonical.sh)
#   - `ollama serve` running, recent version (>= 0.30)
#   - `claude` CLI not required for this open-weight-only run
#
# Usage:
#   bash setup/reproduce_gemma4_m4.sh                 # gemma4:12b (default)
#   bash setup/reproduce_gemma4_m4.sh gemma4:12b lfm2.5:8b   # also reproduce the lfm2.5 think-leak floor
#
set -euo pipefail

cd "$(dirname "$0")/.."          # repo root, regardless of where invoked from
MODELS="${*:-gemma4:12b}"

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
fail() { printf '\033[31mFAIL: %s\033[0m\n' "$*" >&2; exit 1; }

say "Preflight"

command -v ollama >/dev/null || fail "ollama not found. brew install ollama"
ollama list >/dev/null 2>&1   || fail "ollama server not reachable. Run: ollama serve &"

# gemma4:12b was published ~2026-06-04 and 412s ('requires a newer version of
# Ollama') on builds older than ~0.30. Warn early instead of after a long wait.
ver="$(ollama --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
if [ -n "$ver" ]; then
  major_minor="$(printf '%s' "$ver" | cut -d. -f1-2)"
  case "$major_minor" in
    0.2*|0.1*|0.0*) printf '\033[33mWARNING: ollama %s is old; gemma4:12b may 412. Run: brew upgrade ollama\033[0m\n' "$ver" ;;
    *) echo "ollama $ver — ok" ;;
  esac
fi

[ -f plan/PLAN.md ] || fail "plan/PLAN.md missing — wrong directory or incomplete clone."
ls ground_truth/results/*.vcf.gz >/dev/null 2>&1 \
  || fail "ground truth missing. Run: bash ground_truth/canonical.sh"

say "Pulling models: $MODELS"
for m in $MODELS; do
  if ollama list | grep -q "^${m%%:*}"; then
    echo "$m already present"
  else
    ollama pull "$m" || fail "pull failed for $m (if this is a 412, run: brew upgrade ollama)"
  fi
done

say "Running sweep (3 seeds x track A x v2 plan, think off)"
# --keep: don't `ollama rm` afterwards, so a re-run is instant.
python3 harness/sweep_local.py --models $MODELS --keep

say "Result (parsed straight from per-run score.json)"
python3 - "$MODELS" <<'PY'
import sys, json, glob, os, statistics
models = sys.argv[1].split()
expected = {"gemma4:12b": 1.000, "lfm2.5:8b": 0.000}
rc = 0
for m in models:
    safe = m.replace("/", "_").replace(":", "_")
    scores = []
    for d in sorted(glob.glob(f"runs/{safe}_think-off_track-A_seed-*")):
        sj = os.path.join(d, "score.json")
        if not os.path.exists(sj):
            continue
        j = json.load(open(sj))
        scores.append((j.get("m1_executes"), j.get("m3_jaccard")))
    if not scores:
        print(f"  {m:<14} no scored runs found"); rc = 1; continue
    m1 = sum(1 for x, _ in scores if x == 1)
    m3 = statistics.mean(y for _, y in scores)
    line = f"  {m:<14} M1={m1}/{len(scores)}  M3_mean={m3:.3f}"
    if m in expected:
        ok = abs(m3 - expected[m]) < 1e-6
        line += f"   (Jetson: {expected[m]:.3f})  {'MATCH' if ok else 'DIFFERS'}"
        if not ok: rc = 1
    print(line)
sys.exit(rc)
PY

say "Done. Per-run artifacts are under runs/  (results.csv was not modified)."
