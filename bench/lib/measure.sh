#!/usr/bin/env bash
# Wrap a command with /usr/bin/time -v and append structured metrics.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

TOOL=""
SAMPLE=""
STAGE="overall"
OUT_FILE="${BENCH_ROOT}/out/runtime_memory.tsv"

usage(){
  cat <<'EOF'
Usage: measure.sh --tool TOOL --sample SAMPLE [--stage STAGE] [--out FILE] -- COMMAND...
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tool) TOOL="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --stage) STAGE="$2"; shift 2;;
    --out) OUT_FILE="$2"; shift 2;;
    --) shift; break;;
    *) usage;;
  esac
done

if [[ -z "${TOOL}" || -z "${SAMPLE}" ]]; then
  usage
fi

if [[ $# -lt 1 ]]; then
  usage
fi

append_runtime_header "${OUT_FILE}"
TMP_LOG="$(mktemp)"
trap 'rm -f "${TMP_LOG}"' EXIT

log "Running (${TOOL}/${SAMPLE}/${STAGE}) → ${*}"
set +e
{ /usr/bin/time -v "$@" ; } 2> >(tee "${TMP_LOG}" >&2)
STATUS=$?
set -e

python3 - "${TMP_LOG}" "${OUT_FILE}" "${SAMPLE}" "${TOOL}" "${STAGE}" <<'PY'
import csv, sys, math
time_log, out_path, sample, tool, stage = sys.argv[1:6]
metrics = {
    "User time (seconds)": 0.0,
    "System time (seconds)": 0.0,
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "0:00.00",
    "Maximum resident set size (kbytes)": 0.0,
    "File system inputs": 0.0,
    "File system outputs": 0.0,
}

def parse_wall(value: str) -> float:
    value = value.strip()
    if not value:
        return 0.0
    if value.isdigit():
        return float(value)
    parts = value.split(":")
    parts = [float(p.replace(",", ".")) for p in parts]
    if len(parts) == 3:
        h, m, s = parts
    elif len(parts) == 2:
        h = 0.0
        m, s = parts
    else:
        return float(parts[0])
    return h * 3600.0 + m * 60.0 + s

with open(time_log) as fh:
    for raw in fh:
        if ": " not in raw:
            continue
        key, val = raw.split(": ", 1)
        key = key.strip()
        val = val.strip()
        if key in metrics:
            metrics[key] = val

wall = parse_wall(str(metrics["Elapsed (wall clock) time (h:mm:ss or m:ss)"]))
user = float(str(metrics["User time (seconds)"]).replace(",", ".") or 0.0)
sys_time = float(str(metrics["System time (seconds)"]).replace(",", ".") or 0.0)
rss_gb = float(str(metrics["Maximum resident set size (kbytes)"]).replace(",", ".") or 0.0) / (1024.0 * 1024.0)
io_in = float(str(metrics["File system inputs"]).replace(",", ".") or 0.0) / (1024.0 * 1024.0)
io_out = float(str(metrics["File system outputs"]).replace(",", ".") or 0.0) / (1024.0 * 1024.0)

with open(out_path, "a", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow([
        sample, tool, stage,
        f"{wall:.3f}", f"{user:.3f}", f"{sys_time:.3f}",
        f"{rss_gb:.3f}", f"{io_in:.3f}", f"{io_out:.3f}",
    ])
PY

if [[ "${STATUS}" -ne 0 ]]; then
  log "Command for ${TOOL}/${SAMPLE}/${STAGE} exited with status ${STATUS}"
fi

exit "${STATUS}"
