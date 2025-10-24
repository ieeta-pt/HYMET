#!/usr/bin/env bash
# Reproducible CAMI suite runner.
#
# This wrapper executes bin/hymet bench with user-selected tool panels and
# stores *every* artefact under results/<scenario>/<suite>/run_<timestamp>/.
# Each run folder contains:
#   raw/      – full bench outputs grouped by mode (contigs, reads, ...)
#   tables/   – summary TSVs copied from each mode
#   figures/  – regenerated figures (per mode)
#   metadata.json – manifest, tools, commands, commit hash

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYMET_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_ROOT="${HYMET_ROOT}/results"
CONFIG_FILE="${SCRIPT_DIR}/config/cami_suite.cfg"

MANIFEST="${HYMET_ROOT}/bench/cami_manifest.tsv"
SCENARIO="cami"
SUITE_NAME="custom"
SCENARIO_SET=0
SUITE_SET=0
MODES="contigs,reads"
THREADS="${THREADS:-16}"
CACHE_ROOT=""
DRY_RUN=0
BENCH_EXTRA=""
READ_CHUNK_SIZE="${READ_CHUNK_SIZE:-250}"
READ_MIN_CHUNK="${READ_MIN_CHUNK:-125}"
CONTIG_TOOLS=""
READ_TOOLS=""
SUITE_PATH_ARG=""

usage(){
  cat <<'USAGE'
Usage: run_cami_suite.sh [options]

Options:
  --manifest PATH        Manifest TSV (default: bench/cami_manifest.tsv)
  --scenario NAME        Results namespace (default: cami)
  --suite NAME           Suite name within the scenario (default: custom)
  --modes LIST           Comma-separated list (default: contigs,reads)
  --threads N            Thread count (default: env THREADS or 16)
  --cache-root PATH      Override CACHE_ROOT for HYMET runs
  --contig-tools LIST    Tool list for contig mode (override default)
  --read-tools LIST      Tool list for read mode (override default)
  --bench-extra "ARGS"   Extra args forwarded to bin/hymet bench
  --suite-path REL/PATH  Custom results/<path> target (instead of scenario/suite)
  --dry-run              Record metadata without executing the benchmark
  --config PATH          Alternate config file
  -h, --help             Show this message
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --scenario) SCENARIO="$2"; SCENARIO_SET=1; shift 2;;
    --suite) SUITE_NAME="$2"; SUITE_SET=1; shift 2;;
    --modes) MODES="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --cache-root) CACHE_ROOT="$2"; shift 2;;
    --contig-tools) CONTIG_TOOLS="$2"; shift 2;;
    --read-tools) READ_TOOLS="$2"; shift 2;;
    --bench-extra) BENCH_EXTRA="$2"; shift 2;;
    --suite-path) SUITE_PATH_ARG="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    --config) CONFIG_FILE="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

if [[ -f "${CONFIG_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${CONFIG_FILE}"
  if [[ ${SUITE_SET} -eq 0 && -n "${SUITE_NAME_DEFAULT:-}" ]]; then
    SUITE_NAME="${SUITE_NAME_DEFAULT}"
  fi
  if [[ ${SCENARIO_SET} -eq 0 && -n "${SUITE_SCENARIO_DEFAULT:-}" ]]; then
    SCENARIO="${SUITE_SCENARIO_DEFAULT}"
  fi
fi

MANIFEST="$(cd "${HYMET_ROOT}" && python3 -c 'import pathlib, sys; p=pathlib.Path(sys.argv[1]).resolve(); print(p, end="")' "${MANIFEST}")"
[[ -s "${MANIFEST}" ]] || { echo "[suite] manifest not found: ${MANIFEST}" >&2; exit 1; }

CONTIG_TOOLS="${CONTIG_TOOLS:-${CONTIG_TOOLS_DEFAULT:-hymet,kraken2}}"
READ_TOOLS="${READ_TOOLS:-${READ_TOOLS_DEFAULT:-hymet_reads}}"
IFS=',' read -r -a MODE_LIST <<< "${MODES}"

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
if [[ -n "${SUITE_PATH_ARG}" ]]; then
  if [[ "${SUITE_PATH_ARG}" = /* ]]; then
    RUN_DIR="${SUITE_PATH_ARG}/run_${RUN_STAMP}"
  else
    RUN_DIR="${RESULTS_ROOT}/${SUITE_PATH_ARG}/run_${RUN_STAMP}"
  fi
else
  RUN_DIR="${RESULTS_ROOT}/${SCENARIO}/${SUITE_NAME}/run_${RUN_STAMP}"
fi
RAW_ROOT="${RUN_DIR}/raw"
TABLES_ROOT="${RUN_DIR}/tables"
FIGURES_ROOT="${RUN_DIR}/figures"
mkdir -p "${RAW_ROOT}" "${TABLES_ROOT}" "${FIGURES_ROOT}"

commands_log=()
if [[ "${DRY_RUN}" -eq 0 ]]; then
  for mode in "${MODE_LIST[@]}"; do
    mode_trim="${mode// /}"
    case "${mode_trim}" in
      contigs)
        tools="${CONTIG_TOOLS}"
        env_prefix=(env THREADS="${THREADS}")
        ;;
      reads)
        tools="${READ_TOOLS}"
        env_prefix=(env THREADS="${THREADS}" READ_CHUNK_SIZE="${READ_CHUNK_SIZE}" READ_MIN_CHUNK="${READ_MIN_CHUNK}")
        ;;
      *)
        echo "[suite] unsupported mode: ${mode_trim}" >&2
        continue
        ;;
    esac
    if [[ -n "${CACHE_ROOT}" ]]; then
      env_prefix+=(CACHE_ROOT="${CACHE_ROOT}")
    fi
    MODE_RAW="${RAW_ROOT}/${mode_trim}"
    mkdir -p "${MODE_RAW}"
    cmd=("${HYMET_ROOT}/bin/hymet" bench --manifest "${MANIFEST}" --tools "${tools}")
    if [[ -n "${BENCH_EXTRA}" ]]; then
      # shellcheck disable=SC2206
      extra=( ${BENCH_EXTRA} )
      cmd+=("${extra[@]}")
    fi
    commands_log+=("mode=${mode_trim} tools=${tools} command=${env_prefix[*]} BENCH_OUT_ROOT=${MODE_RAW} ${cmd[*]}")
    env BENCH_OUT_ROOT="${MODE_RAW}" "${env_prefix[@]}" "${cmd[@]}"

    mode_tables="${TABLES_ROOT}/${mode_trim}"
    mode_figs="${FIGURES_ROOT}/${mode_trim}"
    mkdir -p "${mode_tables}" "${mode_figs}"
    for table in summary_per_tool_per_sample.tsv leaderboard_by_rank.tsv runtime_memory.tsv contig_accuracy_per_tool.tsv manifest.snapshot.tsv; do
      src="${MODE_RAW}/${table}"
      [[ -f "${src}" ]] && cp "${src}" "${mode_tables}/"
    done
    cp "${MODE_RAW}/"fig_*.png "${mode_figs}/" 2>/dev/null || true
  done
else
  commands_log+=("dry_run modes=${MODES} contig_tools=${CONTIG_TOOLS} read_tools=${READ_TOOLS}")
fi

if ((${#commands_log[@]})); then
  COMMANDS_TEXT=$(printf "%s\n" "${commands_log[@]}")
else
  COMMANDS_TEXT=""
fi

RUN_DIR="${RUN_DIR}" \
RUN_STAMP="${RUN_STAMP}" \
SCENARIO="${SCENARIO}" \
SUITE_NAME="${SUITE_NAME}" \
MANIFEST="${MANIFEST}" \
DRY_RUN="${DRY_RUN}" \
THREADS="${THREADS}" \
CONTIG_TOOLS="${CONTIG_TOOLS}" \
READ_TOOLS="${READ_TOOLS}" \
MODES="${MODES}" \
COMMANDS_LOG="${COMMANDS_TEXT}" \
python3 - <<'PY'
import json, pathlib, os
from datetime import datetime
root = pathlib.Path(os.environ['RUN_DIR'])
meta = {
    "scenario": os.environ['SCENARIO'],
    "suite": os.environ['SUITE_NAME'],
    "run_id": root.name,
    "timestamp": os.environ['RUN_STAMP'],
    "manifest": os.environ['MANIFEST'],
    "dry_run": bool(int(os.environ['DRY_RUN'])),
    "threads": int(os.environ['THREADS']),
    "contig_tools": os.environ['CONTIG_TOOLS'],
    "read_tools": os.environ['READ_TOOLS'],
    "modes": os.environ['MODES'].split(','),
    "commands": os.environ['COMMANDS_LOG'].split('\n') if os.environ['COMMANDS_LOG'] else [],
    "git_commit": os.popen('git rev-parse HEAD').read().strip()
}
(root / 'metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
PY
