#!/usr/bin/env bash
# Regenerate benchmark aggregates and publish snapshots into results/bench.

set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${REPO_ROOT}/bench"
OUT_DIR="${BENCH_DIR}/out"
PUBLISH_DIR="${REPO_ROOT}/results/bench"

if [[ "${SKIP_RECOMPUTE:-0}" != "1" ]]; then
  python3 "${BENCH_DIR}/aggregate_metrics.py" --bench-root "${BENCH_DIR}" --outdir out
  python3 "${BENCH_DIR}/plot/make_figures.py" --bench-root "${BENCH_DIR}" --outdir out
fi

mkdir -p "${PUBLISH_DIR}"

publish_file() {
  local src="$1"
  local dst="$2"
  if [[ -f "${src}" ]]; then
    install -m 0644 -D "${src}" "${dst}"
  fi
}

declare -a PUBLISH_FILES=(
  "summary_per_tool_per_sample.tsv"
  "leaderboard_by_rank.tsv"
  "contig_accuracy_per_tool.tsv"
  "runtime_memory.tsv"
  "fig_f1_by_rank.png"
  "fig_f1_by_rank_lines.png"
  "fig_l1_braycurtis.png"
  "fig_l1_braycurtis_lines.png"
  "fig_accuracy_by_rank.png"
  "fig_accuracy_by_rank_lines.png"
  "fig_per_sample_f1_stack.png"
  "fig_cpu_time_by_tool.png"
  "fig_peak_memory_by_tool.png"
  "fig_runtime_cpu_mem.png"
)

for rel in "${PUBLISH_FILES[@]}"; do
  publish_file "${OUT_DIR}/${rel}" "${PUBLISH_DIR}/${rel}"
done

echo "[publish] Synced aggregates into ${PUBLISH_DIR}"
