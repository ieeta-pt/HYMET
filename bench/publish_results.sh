#!/usr/bin/env bash
# Snapshot bench/out into results/<scenario>/<suite>/run_<timestamp>/.

set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${REPO_ROOT}/bench"
OUT_DIR="${BENCH_DIR}/out"

SCENARIO="${PUBLISH_SCENARIO:-cami}"
SUITE="${PUBLISH_SUITE:-canonical}"
RUN_STAMP="${PUBLISH_RUN_STAMP:-$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="${PUBLISH_RUN_DIR:-${REPO_ROOT}/results/${SCENARIO}/${SUITE}/run_${RUN_STAMP}}"
MANIFEST="${PUBLISH_MANIFEST:-${BENCH_DIR}/cami_manifest.tsv}"
THREADS_META="${PUBLISH_THREADS:-${THREADS:-8}}"
TOOLS_META="${PUBLISH_TOOLS:-all}"
MODES_META="${PUBLISH_MODES:-contigs}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-dir) RUN_DIR="$2"; shift 2;;
    --scenario) SCENARIO="$2"; shift 2;;
    --suite) SUITE="$2"; shift 2;;
    --run-stamp) RUN_STAMP="$2"; shift 2;;
    --manifest) MANIFEST="$2"; shift 2;;
    --threads) THREADS_META="$2"; shift 2;;
    --tools) TOOLS_META="$2"; shift 2;;
    --modes) MODES_META="$2"; shift 2;;
    -h|--help)
      cat <<'USAGE'
Usage: publish_results.sh [options]

Options:
  --run-dir DIR       Target run directory (default: results/<scenario>/<suite>/run_<timestamp>)
  --scenario NAME     Scenario namespace (default: cami)
  --suite NAME        Suite name (default: canonical)
  --run-stamp TS      Timestamp suffix (default: current UTC)
  --manifest PATH     Manifest recorded in metadata
  --threads N         Thread count for metadata
  --tools LIST        Tool list for metadata
  --modes LIST        Modes captured in metadata
USAGE
      exit 0
      ;;
    *) echo "[publish] Unknown option: $1" >&2; exit 1;;
  esac
done

export REPO_ROOT RUN_DIR SCENARIO SUITE RUN_STAMP MANIFEST THREADS_META TOOLS_META MODES_META

if [[ "${SKIP_RECOMPUTE:-0}" != "1" ]]; then
  python3 "${BENCH_DIR}/aggregate_metrics.py" --bench-root "${BENCH_DIR}" --outdir out
  python3 "${BENCH_DIR}/plot/make_figures.py" --bench-root "${BENCH_DIR}" --outdir out
fi

[[ -d "${OUT_DIR}" ]] || { echo "[publish] bench/out/ missing; nothing to package" >&2; exit 1; }

RAW_DIR="${RUN_DIR}/raw"
TABLES_DIR="${RUN_DIR}/tables"
FIGURES_DIR="${RUN_DIR}/figures"
mkdir -p "${RAW_DIR}" "${TABLES_DIR}" "${FIGURES_DIR}"

echo "[publish] Snapshotting bench/out into ${RAW_DIR}"
if command -v rsync >/dev/null 2>&1; then
  rsync -a --delete "${OUT_DIR}/" "${RAW_DIR}/"
else
  echo "[publish] rsync unavailable; falling back to cp -a" >&2
  rm -rf "${RAW_DIR}"
  mkdir -p "${RAW_DIR}"
  cp -a "${OUT_DIR}/." "${RAW_DIR}/"
fi

declare -a TABLE_FILES=(
  "summary_per_tool_per_sample.tsv"
  "leaderboard_by_rank.tsv"
  "contig_accuracy_per_tool.tsv"
  "runtime_memory.tsv"
  "manifest.snapshot.tsv"
)
declare -a FIGURE_FILES=(
  "fig_f1_by_rank_lines.png"
  "fig_l1_braycurtis_lines.png"
  "fig_accuracy_by_rank_lines.png"
  "fig_per_sample_f1_stack.png"
  "fig_cpu_time_by_tool.png"
  "fig_peak_memory_by_tool.png"
  "fig_runtime_cpu_mem.png"
)

copy_if_exists() {
  local src="$1"
  local dst="$2"
  if [[ -f "${src}" ]]; then
    install -m 0644 -D "${src}" "${dst}"
  fi
}

for rel in "${TABLE_FILES[@]}"; do
  copy_if_exists "${OUT_DIR}/${rel}" "${TABLES_DIR}/${rel}"
done
for rel in "${FIGURE_FILES[@]}"; do
  copy_if_exists "${OUT_DIR}/${rel}" "${FIGURES_DIR}/${rel}"
done

python3 - <<'PY'
import json, os, pathlib, subprocess
repo_root = pathlib.Path(os.environ["REPO_ROOT"])
run_dir = pathlib.Path(os.environ["RUN_DIR"])
meta = {
    "scenario": os.environ["SCENARIO"],
    "suite": os.environ["SUITE"],
    "run_id": run_dir.name,
    "timestamp": os.environ["RUN_STAMP"],
    "manifest": str(pathlib.Path(os.environ["MANIFEST"]).resolve()),
    "threads": int(os.environ["THREADS_META"]),
    "tools": [tool for tool in os.environ["TOOLS_META"].split(",") if tool],
    "modes": [mode for mode in os.environ["MODES_META"].split(",") if mode],
    "source": "bench/run_all_cami.sh",
    "git_commit": subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    ).stdout.strip() or "unknown",
}
(run_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
PY

echo "[publish] Results packaged in ${RUN_DIR}"
