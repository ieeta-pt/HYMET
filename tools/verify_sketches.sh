#!/usr/bin/env bash
# Compare local HYMET sketch files with the published Zenodo record.
#
# Usage:
#   tools/verify_sketches.sh [--record 17428354] [--base-url https://zenodo.org] \
#                            [--local-dir data]
#
# Options:
#   --record ID       Zenodo record/DOI numeric ID (default: 17428354)
#   --base-url URL    Base URL for the API (default: https://zenodo.org)
#   --local-dir DIR   Directory containing sketch*.msh (default: HYMET/data)
#   -h, --help        Show this message

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYMET_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

RECORD_ID=17428354
BASE_URL="https://zenodo.org"
LOCAL_DIR="${HYMET_ROOT}/data"

usage(){
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --record) RECORD_ID="$2"; shift 2;;
    --base-url) BASE_URL="$2"; shift 2;;
    --local-dir) LOCAL_DIR="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1" >&2; usage;;
  }
done

[[ -d "$LOCAL_DIR" ]] || { echo "[verify-sketches] Local directory not found: $LOCAL_DIR" >&2; exit 2; }

TMP_META=$(mktemp)
TMP_REMOTE_SUM=$(mktemp)
TMP_LOCAL_SUM=$(mktemp)
trap 'rm -f "$TMP_META" "$TMP_REMOTE_SUM" "$TMP_LOCAL_SUM"' EXIT

curl -fsSL "$BASE_URL/api/records/$RECORD_ID" -o "$TMP_META"
curl -fsSL "$BASE_URL/records/$RECORD_ID/files/sketch_sha256.txt?download=1" -o "$TMP_REMOTE_SUM"

sha256sum "$LOCAL_DIR"/sketch*.msh > "$TMP_LOCAL_SUM"

python3 - "$TMP_META" "$TMP_REMOTE_SUM" "$TMP_LOCAL_SUM" "$LOCAL_DIR" <<'PY'
import json
import sys
from pathlib import Path

meta_path, remote_sum_path, local_sum_path, local_dir = sys.argv[1:5]

with open(meta_path, 'r', encoding='utf-8') as fh:
    meta = json.load(fh)

remote_files = {f['key']: f for f in meta.get('files', [])}

remote_hashes = {}
with open(remote_sum_path, 'r', encoding='utf-8') as fh:
    for line in fh:
        parts = line.strip().split()
        if len(parts) >= 2:
            rel = parts[-1]
            name = Path(rel).name
            remote_hashes[name] = parts[0].lower()

local_hashes = {}
with open(local_sum_path, 'r', encoding='utf-8') as fh:
    for line in fh:
        parts = line.strip().split()
        if len(parts) >= 2:
            local_hashes[Path(parts[-1]).name] = parts[0].lower()

expected = [name for name in remote_hashes if name.startswith('sketch')]

all_ok = True
rows = []
for name in sorted(expected):
    remote_hash = remote_hashes.get(name)
    local_hash = local_hashes.get(name)
    local_path = Path(local_dir) / name
    local_exists = local_path.exists()
    remote_size = remote_files.get(name, {}).get('size')
    local_size = local_path.stat().st_size if local_exists else None
    hash_ok = (remote_hash is not None and local_hash == remote_hash)
    size_ok = (remote_size is not None and local_size == remote_size)
    if not (local_exists and hash_ok and size_ok):
        all_ok = False
    rows.append((name, remote_size, local_size, hash_ok, size_ok, local_exists))

links = meta.get('links', {})
landing = links.get('html') or links.get('self_html') or links.get('self')

print(f"Zenodo record: {meta.get('metadata', {}).get('title', '<unknown>')} (ID {meta.get('id')})")
print(f"Landing page: {landing}")
print()
print(f"{'File':<15} {'Remote Size':>12} {'Local Size':>12}  HashOK  SizeOK  Exists")
print('-' * 70)
for name, rsize, lsize, hash_ok, size_ok, exists in rows:
    rsize_str = f"{rsize:,}" if rsize is not None else 'n/a'
    lsize_str = f"{lsize:,}" if lsize is not None else 'missing'
    print(f"{name:<15} {rsize_str:>12} {lsize_str:>12}  {str(hash_ok):<6} {str(size_ok):<6} {str(exists):<6}")

print('\nRemote SHA256:')
for name in sorted(expected):
    print(f"  {name}: {remote_hashes.get(name, 'n/a')}")

print('\nLocal SHA256:')
for name in sorted(expected):
    print(f"  {name}: {local_hashes.get(name, 'missing')}")

print()
if all_ok:
    print('All comparisons PASS ✅')
else:
    print('Comparisons FAIL ❌ (see above)')
    sys.exit(1)
PY
