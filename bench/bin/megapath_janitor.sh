#!/usr/bin/env bash
# Periodically prune abandoned MegaPath-Nano split/temp files in per-sample run dirs.
# Safe defaults:
#  - Only touches files older than AGE_MIN minutes (default: 45)
#  - Only under HYMET/bench/out/*/megapath_nano/run/tmp (and select tmp.*.tmp in run/)
#  - Logs actions to HYMET/bench/out/logs/megapath_janitor.log

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RUN_ROOT="${RUN_ROOT:-${BENCH_ROOT}/out}"
AGE_MIN="${AGE_MIN:-45}"
INTERVAL_SEC="${INTERVAL_SEC:-300}"
LOG_DIR="${BENCH_ROOT}/out/logs"
LOG_FILE="${LOG_DIR}/megapath_janitor.log"
LOCK_FILE="${LOG_DIR}/megapath_janitor.lock"

mkdir -p "${LOG_DIR}"

log(){ printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOG_FILE}"; }

# Ensure single instance
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  log "Another janitor instance is running; exiting"
  exit 0
fi

log "Starting janitor: RUN_ROOT=${RUN_ROOT}, AGE_MIN=${AGE_MIN} min, INTERVAL=${INTERVAL_SEC}s"

while true; do
  start_ts=$(date +%s)
  # Targets within each sample run dir
  # 1) minimap2 split files (we set --split-prefix to mm2_split under run/tmp)
  # 2) generic tmp.*.tmp files in run/ (in case of legacy runs)
  removed=0
  bytes=0

  # mm2 split files under run/tmp
  while IFS= read -r -d '' f; do
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    rm -f -- "$f" 2>/dev/null || true
    removed=$((removed+1))
    bytes=$((bytes+sz))
  done < <(find "${RUN_ROOT}" -type f \
            -path '*/megapath_nano/run/tmp/*' \
            -name 'mm2_split*' -mmin +"${AGE_MIN}" -print0 2>/dev/null)

  # Legacy tmp files directly under run/
  while IFS= read -r -d '' f; do
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    rm -f -- "$f" 2>/dev/null || true
    removed=$((removed+1))
    bytes=$((bytes+sz))
  done < <(find "${RUN_ROOT}" -maxdepth 4 -type f \
            -path '*/megapath_nano/run/*' \
            -name 'tmp.*.tmp' -mmin +"${AGE_MIN}" -print0 2>/dev/null)

  # Also trim empty dirs left behind under run/tmp
  find "${RUN_ROOT}" -type d -path '*/megapath_nano/run/tmp/*' -empty -delete 2>/dev/null || true

  # Human-readable bytes
  hr_bytes="${bytes}B"
  python3 - <<'PY' "${bytes}" || true
import sys
b=int(sys.argv[1])
u=['B','KB','MB','GB','TB']
i=0
v=float(b)
while v>=1024 and i<4:
    v/=1024;i+=1
print(f"{v:.2f} {u[i]}")
PY
  hr=$(python3 - <<'PY' "${bytes}" 2>/dev/null || true)
import sys
b=int(sys.argv[1]);u=['B','KB','MB','GB','TB'];i=0;v=float(b)
while v>=1024 and i<4:
    v/=1024;i+=1
print(f"{v:.2f} {u[i]}")
PY
  hr=${hr:-${hr_bytes}}

  end_ts=$(date +%s)
  took=$((end_ts-start_ts))
  log "Sweep: removed ${removed} files (~${hr}) in ${took}s"
  sleep "${INTERVAL_SEC}"
done

