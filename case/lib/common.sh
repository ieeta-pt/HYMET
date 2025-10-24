#!/usr/bin/env bash
# Shared helpers for case-study scripts.

if [[ -n "${CASE_COMMON_SOURCED:-}" ]]; then
  return 0
fi
CASE_COMMON_SOURCED=1

set -o errexit
set -o nounset
set -o pipefail

_case__this="${BASH_SOURCE[0]}"
if [[ "${_case__this}" != */* ]]; then
  _case__this="./${_case__this}"
fi
_case__dir="$(cd "$(dirname "${_case__this}")" && pwd)"
export CASE_ROOT="$(cd "${_case__dir}/.." && pwd)"
export HYMET_ROOT="$(cd "${CASE_ROOT}/.." && pwd)"
export BENCH_ROOT="$(cd "${CASE_ROOT}/../bench" && pwd)"

log(){ printf '[%(%F %T)T] %s\n' -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }

ensure_dir(){
  local path="$1"
  mkdir -p "${path}"
}

resolve_path(){
  local input="$1"
  python3 - "$input" "$CASE_ROOT" <<'PY'
import os, sys
value, case_root = sys.argv[1], sys.argv[2]
if not value:
    print("", end="")
elif os.path.isabs(value):
    print(os.path.normpath(value), end="")
else:
    print(os.path.normpath(os.path.join(case_root, value)), end="")
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

manifest_split_line(){
  local line="$1"
  python3 - <<'PY' -- "${line}"
import csv, sys
SEP = "\x1f"
if len(sys.argv) < 3:
    raise SystemExit("manifest_split_line requires a non-empty line")
line = sys.argv[2]
row = next(csv.reader([line], delimiter='\t'))
sys.stdout.write(SEP.join(row))
PY
}
