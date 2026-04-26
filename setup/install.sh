#!/usr/bin/env bash
set -euo pipefail

# Idempotent installer: miniforge3-aarch64 + locked bioconda env named `bench`.
# Re-running skips already-completed steps.

MINIFORGE_DIR="$HOME/miniforge3"
ENV_NAME="bench"
INSTALLER="/tmp/Miniforge3-Linux-aarch64.sh"
URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh"

if [[ ! -d "$MINIFORGE_DIR" ]]; then
  echo "[install] downloading miniforge..."
  curl -fsSL -o "$INSTALLER" "$URL"
  echo "[install] installing miniforge to $MINIFORGE_DIR (no rc-file modification)..."
  bash "$INSTALLER" -b -p "$MINIFORGE_DIR"
  rm -f "$INSTALLER"
else
  echo "[install] miniforge already present at $MINIFORGE_DIR — skipping."
fi

# shellcheck disable=SC1091
source "$MINIFORGE_DIR/etc/profile.d/conda.sh"

conda config --set channel_priority strict --file "$MINIFORGE_DIR/.condarc" || true

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[install] env '$ENV_NAME' already exists — skipping create."
else
  echo "[install] creating env '$ENV_NAME' (this may take 2-3 min)..."
  conda create -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda \
    bwa=0.7.18 \
    samtools=1.21 \
    bcftools=1.21 \
    htslib=1.21 \
    lofreq=2.1.5 \
    snpsift=5.2 \
    snpeff=5.2 \
    fastqc=0.12.1 \
    seqkit=2.8 \
    snakemake-minimal=8.20 \
    shellcheck=0.10 \
    openjdk=21 \
    python=3.12
fi

echo "[install] done. activate with: source $MINIFORGE_DIR/etc/profile.d/conda.sh && conda activate $ENV_NAME"
