#!/usr/bin/env bash
set -euo pipefail

# Asserts every benchmark tool is on PATH at the pinned version.
# Stdout is the TOOL_INVENTORY string injected verbatim into model prompts.
# Exits non-zero if any tool is missing or version mismatches.

# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate bench

fail=0
report() {
  local name="$1" expected="$2" actual="$3"
  if [[ "$actual" == *"$expected"* ]]; then
    printf "  %-12s %s  [ok]\n" "$name" "$actual"
  else
    printf "  %-12s %s  [EXPECTED %s]\n" "$name" "$actual" "$expected"
    fail=1
  fi
}

echo "TOOL_INVENTORY (linux-aarch64, conda env 'bench')"
echo "================================================="

# bwa prints version on stderr to a no-arg invocation; capture line 3
report bwa        "0.7.18"   "$(bwa 2>&1 | grep -oE 'Version: [0-9.-]+' | head -1)"
report samtools   "1.21"     "$(samtools --version | head -1)"
report bcftools   "1.21"     "$(bcftools --version | head -1)"
# htslib version comes via samtools/bcftools build line; check tabix as proxy
report tabix      "1.21"     "$(tabix --version 2>&1 | head -1 || true)"
report lofreq     "2.1.5"    "$(lofreq version 2>&1 | head -1)"
report SnpSift    "5.2"      "$(SnpSift 2>&1 | grep -oE 'SnpSift version [0-9a-z.]+' | head -1 || true)"
report snpEff     "5.2"      "$(snpEff -version 2>&1 | head -1 || true)"
report fastqc     "0.12.1"   "$(fastqc --version 2>&1 | head -1)"
report seqkit     "2.8"      "$(seqkit version 2>&1 | head -1)"
report snakemake  "8.20"     "$(snakemake --version 2>&1 | head -1)"
report shellcheck "0.10"     "$(shellcheck --version | grep -oE 'version: [0-9.]+' | head -1)"
report java       "21"       "$(java -version 2>&1 | head -1)"

echo "================================================="
if [[ $fail -ne 0 ]]; then
  echo "FAIL: one or more tools missing or wrong version" >&2
  exit 1
fi
echo "OK"
