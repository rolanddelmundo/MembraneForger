#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

mf_cpu_threads() {
  if command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.logicalcpu 2>/dev/null && return
  fi
  getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1
}

: "${MEMBRANEFORGER_LOCAL_THREADS:=$(mf_cpu_threads)}"
export OMP_NUM_THREADS="$MEMBRANEFORGER_LOCAL_THREADS"
export OPENMM_CPU_THREADS="$MEMBRANEFORGER_LOCAL_THREADS"
export VECLIB_MAXIMUM_THREADS="$MEMBRANEFORGER_LOCAL_THREADS"
export OPENBLAS_NUM_THREADS="$MEMBRANEFORGER_LOCAL_THREADS"
export MKL_NUM_THREADS="$MEMBRANEFORGER_LOCAL_THREADS"
export NUMEXPR_NUM_THREADS="$MEMBRANEFORGER_LOCAL_THREADS"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  if command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
  else
    echo "ERROR: Python executable not found; set PYTHON_BIN" >&2
    exit 127
  fi
fi
export PYTHONPATH="$ROOT_DIR${PYTHONPATH:+:$PYTHONPATH}"
exec "$PYTHON_BIN" -m membraneforger.cli stage stage4 "$@"
