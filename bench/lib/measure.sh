#!/usr/bin/env bash
# Wrap a command with /usr/bin/time -v and append structured metrics.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

TOOL=""
SAMPLE=""
STAGE="overall"
OUT_FILE="${BENCH_ROOT}/out/runtime_memory.tsv"
LOCAL_FILE=""

usage(){
  cat <<'EOF'
Usage: measure.sh --tool TOOL --sample SAMPLE [--stage STAGE] [--out FILE] [--local FILE] -- COMMAND...
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tool) TOOL="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --stage) STAGE="$2"; shift 2;;
    --out) OUT_FILE="$2"; shift 2;;
    --local) LOCAL_FILE="$2"; shift 2;;
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
if [[ -n "${LOCAL_FILE}" ]]; then
  append_runtime_header "${LOCAL_FILE}"
fi
TMP_LOG="$(mktemp)"
trap 'rm -f "${TMP_LOG}"' EXIT

log "Running (${TOOL}/${SAMPLE}/${STAGE}) → ${*}"
MODE_TAG="${HYMET_BENCH_MODE:-${BENCH_MODE:-}}"
START_TS="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
CMD_STR="$*"
set +e
{ /usr/bin/time -v "$@" ; } 2> >(tee "${TMP_LOG}" >&2)
STATUS=$?
set -e
END_TS="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

python3 - "${TMP_LOG}" "${OUT_FILE}" "${LOCAL_FILE}" "${SAMPLE}" "${TOOL}" "${MODE_TAG}" "${STAGE}" "${START_TS}" "${END_TS}" "${STATUS}" "${CMD_STR}" <<'PY'
import csv, sys
time_log, global_out, local_out, sample, tool, mode, stage, started, finished, status, command = sys.argv[1:12]
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

row = [
    sample,
    tool,
    mode or "",
    stage,
    started,
    finished,
    f"{wall:.3f}",
    f"{user:.3f}",
    f"{sys_time:.3f}",
    f"{rss_gb:.3f}",
    f"{io_in:.3f}",
    f"{io_out:.3f}",
    command,
    status,
]

def append_row(path):
    if not path:
        return
    with open(path, "a", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(row)

append_row(global_out)
append_row(local_out)
PY

if [[ "${STATUS}" -ne 0 ]]; then
  log "Command for ${TOOL}/${SAMPLE}/${STAGE} exited with status ${STATUS}"
fi

exit "${STATUS}"
