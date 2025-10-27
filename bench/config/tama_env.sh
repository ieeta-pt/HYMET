#!/usr/bin/env bash
# Source this file to configure TAMA integration for the HYMET benchmark harness.
# Edit TAMA_ROOT to match your local installation of the TAMA repository.

_tama_env_this="${BASH_SOURCE[0]}"
if [[ "${_tama_env_this}" != */* ]]; then
  _tama_env_this="./${_tama_env_this}"
fi
_tama_env_dir="$(cd "$(dirname "${_tama_env_this}")" && pwd)"
_bench_dir="$(cd "${_tama_env_dir}/.." && pwd)"

export BENCH_ROOT="${BENCH_ROOT:-${_bench_dir}}"
export HYMET_ROOT="${HYMET_ROOT:-$(cd "${BENCH_ROOT}/.." && pwd)}"

export PATH="/root/.local/share/mamba/bin:${PATH}"

if [[ -z "${TAMA_ROOT:-}" ]]; then
  export TAMA_ROOT="/data/tools/TAMA"
fi

export TAMA_PARAM_DIR="${BENCH_ROOT}/config/tama_params"

# Default to per-sample parameter files generated under $TAMA_PARAM_DIR; the runner
# will append the sample name automatically if TAMA_PARAM_FILE is unset.
unset TAMA_PARAM_FILE

# Keep intermediates only when explicitly requested to avoid large disk usage.
export KEEP_TAMA_WORK="${KEEP_TAMA_WORK:-0}"
