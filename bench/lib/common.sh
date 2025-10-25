#!/usr/bin/env bash
# Shared helpers for bench scripts.

if [[ -n "${BENCH_COMMON_SOURCED:-}" ]]; then
  return 0
fi
BENCH_COMMON_SOURCED=1

: "${BENCH_CALLER_PWD:=$(pwd -P)}"

set -o errexit
set -o nounset
set -o pipefail

_bench__this="${BASH_SOURCE[0]}"
if [[ "${_bench__this}" != */* ]]; then
  _bench__this="./${_bench__this}"
fi
_bench__dir="$(cd "$(dirname "${_bench__this}")" && pwd)"
export BENCH_ROOT="$(cd "${_bench__dir}/.." && pwd)"
export HYMET_ROOT="$(cd "${BENCH_ROOT}/.." && pwd)"

log(){ printf '[%(%F %T)T] %s\n' -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

ensure_dir(){
  local path="$1"
  mkdir -p "${path}"
}

resolve_path(){
  local input="${1:-}"
  local base="${2:-${BENCH_CALLER_PWD}}"
  python3 - "$input" "$base" <<'PY'
import os, pathlib, sys
value, base = sys.argv[1], sys.argv[2]
if not value:
    print("", end="")
    raise SystemExit
value = os.path.expanduser(value)
path = pathlib.Path(value)
if path.is_absolute():
    print(str(path.resolve()), end="")
    raise SystemExit
base_path = pathlib.Path(base or ".").resolve()
print(str((base_path / path).resolve()), end="")
PY
}

normalize_metadata_json(){
  local json_path="$1"
  local base_dir="${2:-$(dirname "$1")}"
  if [[ ! -f "${json_path}" ]]; then
    return 0
  fi
  python3 - "$json_path" "$base_dir" <<'PY'
import json, os, sys

json_path, base_dir = sys.argv[1:3]
base_dir = os.path.abspath(base_dir)

def convert(value):
    if isinstance(value, dict):
        return {key: convert(val) for key, val in value.items()}
    if isinstance(value, list):
        return [convert(item) for item in value]
    if isinstance(value, str) and os.path.isabs(value):
        try:
            return os.path.relpath(value, base_dir)
        except ValueError:
            return value
    return value

with open(json_path, "r", encoding="utf-8") as handle:
    data = json.load(handle)

data = convert(data)

with open(json_path, "w", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
}

need_cmd(){
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    log "ERROR: missing required command '${cmd}'"
    return 1
  fi
}

append_runtime_header(){
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    ensure_dir "$(dirname "${path}")"
    cat <<'EOF' >"${path}"
sample	tool	mode	stage	started_at	finished_at	wall_seconds	user_seconds	sys_seconds	max_rss_gb	io_input_mb	io_output_mb	command	exit_code
EOF
  fi
}

manifest_split_line(){
  # Usage: manifest_split_line "<line>"
  local line="$1"
  python3 - <<'PY' -- "${line}"
import csv, sys
line = sys.argv[1]
row = next(csv.reader([line], delimiter="\t"))
sys.stdout.write("\0".join(row))
PY
}
