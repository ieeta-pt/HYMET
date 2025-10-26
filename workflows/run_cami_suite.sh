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
PRIMARY_MODE=""
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
  --primary-mode NAME    Mode whose outputs populate top-level figures/tables (default: first mode)
  --bench-extra "ARGS"   Extra args forwarded to bin/hymet bench
  --suite-path REL/PATH  Custom results/<path> target (instead of scenario/suite)
  --dry-run              Record metadata without executing the benchmark
  --config PATH          Alternate config file
  -h, --help             Show this message
USAGE
  exit 1
}

copy_runtime_with_mode(){
  local src="$1"
  local dst="$2"
  local mode="$3"
  [[ -f "${src}" ]] || return 0
  python3 - "$src" "$dst" "$mode" <<'PY'
import csv, sys, pathlib
src_path, dst_path, mode = sys.argv[1:4]
src = pathlib.Path(src_path)
dst = pathlib.Path(dst_path)
if not src.is_file():
    raise SystemExit(0)
with src.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames or []
if not rows:
    raise SystemExit(0)
if "mode" not in fieldnames:
    idx = 2 if len(fieldnames) >= 2 else len(fieldnames)
    fieldnames = fieldnames[:idx] + ["mode"] + fieldnames[idx:]
    for row in rows:
        row["mode"] = mode
else:
    for row in rows:
        row["mode"] = row.get("mode") or mode
dst.parent.mkdir(parents=True, exist_ok=True)
with dst.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "") for key in fieldnames})
PY
}

append_runtime_aggregate(){
  local src="$1"
  local dst="$2"
  [[ -f "${src}" ]] || return 0
  mkdir -p "$(dirname "${dst}")"
  if [[ ! -s "${dst}" ]]; then
    cp "${src}" "${dst}"
  else
    tail -n +2 "${src}" >> "${dst}"
  fi
}

append_table_with_mode(){
  local src="$1"
  local dst="$2"
  local mode="$3"
  [[ -f "${src}" ]] || return 0
  python3 - "$src" "$dst" "$mode" <<'PY'
import csv, sys, pathlib
src_path, dst_path, mode = sys.argv[1:4]
src = pathlib.Path(src_path)
dst = pathlib.Path(dst_path)
if not src.is_file():
    raise SystemExit(0)
with src.open(newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows = list(reader)
    fieldnames = reader.fieldnames or []
if not rows:
    raise SystemExit(0)
if "mode" not in fieldnames:
    if "tool" in fieldnames:
        idx = fieldnames.index("tool") + 1
    else:
        idx = len(fieldnames)
    fieldnames = fieldnames[:idx] + ["mode"] + fieldnames[idx:]
    for row in rows:
        row["mode"] = mode
else:
    for row in rows:
        row["mode"] = row.get("mode") or mode
dst.parent.mkdir(parents=True, exist_ok=True)
write_header = not dst.exists()
with dst.open("a", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    if write_header:
        writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "") for key in fieldnames})
PY
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
    --primary-mode) PRIMARY_MODE="$2"; shift 2;;
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
if [[ -z "${PRIMARY_MODE}" && ${#MODE_LIST[@]} -gt 0 ]]; then
  PRIMARY_MODE="${MODE_LIST[0]// /}"
fi
PRIMARY_MODE="${PRIMARY_MODE// /}"

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
primary_mirrored=0
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
    env_prefix+=(HYMET_BENCH_MODE="${mode_trim}")
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
    for table in summary_per_tool_per_sample.tsv leaderboard_by_rank.tsv contig_accuracy_per_tool.tsv manifest.snapshot.tsv; do
      src="${MODE_RAW}/${table}"
      if [[ -f "${src}" ]]; then
        cp "${src}" "${mode_tables}/"
        case "${table}" in
          summary_per_tool_per_sample.tsv|leaderboard_by_rank.tsv|contig_accuracy_per_tool.tsv)
            append_table_with_mode "${src}" "${TABLES_ROOT}/combined/${table}" "${mode_trim}"
            ;;
        esac
      fi
    done
    cp "${MODE_RAW}/"fig_*.png "${mode_figs}/" 2>/dev/null || true

    runtime_src="${MODE_RAW}/runtime_memory.tsv"
    if [[ -f "${runtime_src}" ]]; then
      copy_runtime_with_mode "${runtime_src}" "${mode_tables}/runtime_memory.tsv" "${mode_trim}"
      append_runtime_aggregate "${mode_tables}/runtime_memory.tsv" "${TABLES_ROOT}/runtime_memory.tsv"
    fi

    if [[ ${primary_mirrored} -eq 0 && "${mode_trim}" == "${PRIMARY_MODE}" ]]; then
      for table in summary_per_tool_per_sample.tsv leaderboard_by_rank.tsv contig_accuracy_per_tool.tsv manifest.snapshot.tsv; do
        src="${MODE_RAW}/${table}"
        [[ -f "${src}" ]] && cp "${src}" "${TABLES_ROOT}/"
      done
      cp "${MODE_RAW}/"fig_*.png "${FIGURES_ROOT}/" 2>/dev/null || true
      primary_mirrored=1
    fi
  done
else
  commands_log+=("dry_run modes=${MODES} contig_tools=${CONTIG_TOOLS} read_tools=${READ_TOOLS}")
fi

if [[ ${primary_mirrored} -eq 0 && "${DRY_RUN}" -eq 0 ]]; then
  echo "[suite] WARNING: primary mode '${PRIMARY_MODE}' not found; top-level figures/tables were not populated" >&2
fi

if [[ "${DRY_RUN}" -eq 0 ]]; then
  combined_tables="${TABLES_ROOT}/combined"
  if [[ -d "${combined_tables}" ]]; then
    mkdir -p "${FIGURES_ROOT}/combined"
    python3 "${SCRIPT_DIR}/plot/make_combined_figures.py" \
      --tables "${TABLES_ROOT}" \
      --outdir "${FIGURES_ROOT}/combined" \
      --suite "${SUITE_NAME}" \
      --run "${RUN_DIR}" \
      || echo "[suite] WARNING: combined figure generation failed" >&2
  fi
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
PRIMARY_MODE="${PRIMARY_MODE}" \
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
    "primary_mode": os.environ['PRIMARY_MODE'],
    "commands": os.environ['COMMANDS_LOG'].split('\n') if os.environ['COMMANDS_LOG'] else [],
    "git_commit": os.popen('git rev-parse HEAD').read().strip()
}
(root / 'metadata.json').write_text(json.dumps(meta, indent=2) + '\n')
PY
