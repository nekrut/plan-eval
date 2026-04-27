#!/usr/bin/env bash
set -euo pipefail

# Define paths relative to current working directory
RAW_DIR="data/raw"
REF_DIR="data/ref"
RES_DIR="results"
REF_GENOME="${REF_DIR}/chrM.fa"

# Create results directory if it doesn't exist
mkdir -p "${RES_DIR}"

# Function to check if a file exists and is not empty
check_file() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "Error: Required file not found: $file" >&2
        exit 1
    fi
}

# Function to run a command only if output file doesn't exist (idempotency)
run_if_missing() {
    local output_file="$1"
    local cmd="$2"
    if [[ ! -f "$output_file" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")
    local cmd="$1"
    shift
    local outputs=("$@")
    if [[ ! -f "${outputs[0]}" ]]; then
        eval "$cmd"
    fi
}

# Function to run a command only if output file doesn't exist (idempotency) for multiple outputs
run_if_missing_multi() {
    local output_files=("$@")