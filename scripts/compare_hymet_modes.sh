#!/usr/bin/env bash
# Automate HYMET contig-vs-read benchmarking on the CAMI suite and generate comparison figures.

set -Eeuo pipefail

log() {
  printf '[%(%F %T)T] %s\n' -1 "$*"
}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_ROOT="${ROOT}/bench"
RUN_SCRIPT="${BENCH_ROOT}/run_all_cami.sh"

[[ -x "${RUN_SCRIPT}" ]] || { log "Cannot find bench runner at ${RUN_SCRIPT}"; exit 1; }

THREADS="${THREADS:-8}"
SCENARIO="${SCENARIO:-cami}"
CONTIG_SUITE="${CONTIG_SUITE:-mode_contigs}"
READ_SUITE="${READ_SUITE:-mode_reads}"
COMPARE_SUITE="${COMPARE_SUITE:-reads_vs_contigs}"

latest_run_dir() {
  local suite="$1"
  local base="${ROOT}/results/${SCENARIO}/${suite}"
  [[ -d "${base}" ]] || return 1
  ls -dt "${base}"/run_* 2>/dev/null | head -n1
}

log "Running HYMET (contig mode)"
PUBLISH_SCENARIO="${SCENARIO}" \
PUBLISH_SUITE="${CONTIG_SUITE}" \
"${RUN_SCRIPT}" --tools hymet --threads "${THREADS}" "$@"
CONTIG_RUN="$(latest_run_dir "${CONTIG_SUITE}")"
[[ -n "${CONTIG_RUN}" && -d "${CONTIG_RUN}" ]] || { log "Contig run directory not found"; exit 1; }

log "Running HYMET (read mode)"
PUBLISH_SCENARIO="${SCENARIO}" \
PUBLISH_SUITE="${READ_SUITE}" \
"${RUN_SCRIPT}" --tools hymet_reads --threads "${THREADS}" "$@"
READ_RUN="$(latest_run_dir "${READ_SUITE}")"
[[ -n "${READ_RUN}" && -d "${READ_RUN}" ]] || { log "Read run directory not found"; exit 1; }

log "Contig run: ${CONTIG_RUN}"
log "Read run:   ${READ_RUN}"

timestamp="$(date -u +%Y%m%dT%H%M%SZ)"
OUT_BASE="${ROOT}/results/${SCENARIO}/${COMPARE_SUITE}/run_${timestamp}"
TABLES_DIR="${OUT_BASE}/tables"
FIGURES_DIR="${OUT_BASE}/figures"
mkdir -p "${TABLES_DIR}" "${FIGURES_DIR}"

contig_tables="${CONTIG_RUN}/tables"
read_tables="${READ_RUN}/tables"

merge_table() {
  local name="$1"
  local dest="${TABLES_DIR}/${name}"
  local header_written=0

  rm -f "${dest}"

  for entry in "contigs:hymet_contigs:${contig_tables}" "reads:hymet_reads:${read_tables}"; do
    IFS=':' read -r mode label source_dir <<<"${entry}"
    local src="${source_dir}/${name}"
    [[ -s "${src}" ]] || continue

    python3 - "$src" "$dest" "$label" "$mode" "$header_written" <<'PY'
import csv, sys
src, dest, label, mode, header_written = sys.argv[1:6]
header_written = int(header_written)

with open(src, newline='', encoding='utf-8') as fh:
    reader = csv.reader(fh, delimiter='\t')
    rows = list(reader)

if not rows:
    sys.exit(0)

header = rows[0]
try:
    tool_idx = header.index('tool')
except ValueError:
    tool_idx = None

try:
    mode_idx = header.index('mode')
except ValueError:
    header.append('mode')
    mode_idx = len(header) - 1
    for row in rows[1:]:
        row.append('')

with open(dest, 'a', newline='', encoding='utf-8') as out:
    writer = csv.writer(out, delimiter='\t')
    if not header_written:
        writer.writerow(header)
    for row in rows[1:]:
        if tool_idx is not None:
            row[tool_idx] = label
        row[mode_idx] = mode
        writer.writerow(row)
PY

    header_written=1
  done
}

merge_table summary_per_tool_per_sample.tsv
merge_table contig_accuracy_per_tool.tsv
merge_table runtime_memory.tsv
merge_table leaderboard_by_rank.tsv

if [[ -s "${contig_tables}/manifest.snapshot.tsv" ]]; then
  cp -f "${contig_tables}/manifest.snapshot.tsv" "${TABLES_DIR}/manifest.snapshot.tsv"
fi

log "Generating comparison figures in ${FIGURES_DIR}"
python3 "${BENCH_ROOT}/plot/make_figures.py" \
  --bench-root "${BENCH_ROOT}" \
  --tables "${TABLES_DIR}" \
  --outdir "${FIGURES_DIR}"

cat >"${OUT_BASE}/metadata.json" <<EOF
{
  "contig_run": "${CONTIG_RUN}",
  "read_run": "${READ_RUN}",
  "comparison_suite": "${SCENARIO}/${COMPARE_SUITE}",
  "generated": "${timestamp}",
  "threads": ${THREADS}
}
EOF

log "Comparison complete → ${OUT_BASE}"
