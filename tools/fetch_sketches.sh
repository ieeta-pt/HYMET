#!/usr/bin/env bash
# Fetch HYMET Mash sketch databases with optional checksum verification.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYMET_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DEST_DIR="${HYMET_ROOT}/data"
CHECKSUMS_TSV="${HYMET_ROOT}/data/sketch_checksums.tsv"
# Published archive hosting the Mash sketches (override via HYMET_SKETCH_BASE_URL or --base-url)
DEFAULT_BASE_URL="https://zenodo.org/records/17428354/files"
BASE_URL="${HYMET_SKETCH_BASE_URL:-$DEFAULT_BASE_URL}"
SKETCHES=(sketch1.msh sketch2.msh sketch3.msh)
SKIP_VERIFY=0

usage(){
  cat <<'USAGE'
Usage: tools/fetch_sketches.sh [--dest DIR] [--base-url URL] [--checksums FILE] [--skip-verify]

Notes:
  - If --base-url is omitted, uses HYMET_SKETCH_BASE_URL or the default Zenodo record.
  - Expected files: sketch1.msh, sketch2.msh, sketch3.msh
  - If a checksum for a file is present in the TSV, it will be verified.
    Missing checksum entries result in a warning (not an error) unless verification is enforced.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest) DEST_DIR="$2"; shift 2;;
    --base-url) BASE_URL="$2"; shift 2;;
    --checksums) CHECKSUMS_TSV="$2"; shift 2;;
    --skip-verify) SKIP_VERIFY=1; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

mkdir -p "${DEST_DIR}"

if [[ -z "${BASE_URL}" ]]; then
  echo "[fetch-sketches] ERROR: --base-url not set and HYMET_SKETCH_BASE_URL is empty." >&2
  echo "Set HYMET_SKETCH_BASE_URL or pass --base-url (e.g., a Zenodo/GitHub release URL)." >&2
  exit 2
fi

have_cmd(){ command -v "$1" >/dev/null 2>&1; }

download(){
  local url="$1"; local out="$2"
  if have_cmd curl; then
    curl -L --fail --retry 5 --retry-delay 2 -o "${out}" "${url}"
  elif have_cmd aria2c; then
    aria2c -x 8 -s 8 -o "${out}" "${url}"
  elif have_cmd wget; then
    wget -O "${out}" "${url}"
  else
    echo "[fetch-sketches] ERROR: need aria2c, curl or wget." >&2
    return 1
  fi
}

read_checksum(){
  local file="$1"
  local tsv="$2"
  [[ -s "${tsv}" ]] || return 1
  awk -v f="${file}" 'BEGIN{FS="\t"} NR>1{ if($1==f){print $2; exit} }' "${tsv}" || true
}

verify_file(){
  local path="$1"; local expected="$2"
  if [[ -z "${expected}" ]]; then
    echo "[fetch-sketches] WARNING: no checksum for $(basename "${path}") — skipping verification" >&2
    return 0
  fi
  local actual
  actual=$(sha256sum "${path}" | awk '{print $1}')
  if [[ "${actual}" != "${expected}" ]]; then
    echo "[fetch-sketches] ERROR: checksum mismatch for $(basename "${path}")" >&2
    echo "  expected: ${expected}" >&2
    echo "  actual:   ${actual}" >&2
    return 1
  fi
  return 0
}

for name in "${SKETCHES[@]}"; do
  dst="${DEST_DIR}/${name}"
  url="${BASE_URL%/}/${name}"
  echo "[fetch-sketches] Fetching ${name} → ${dst}"

  tmp="${dst}.tmp.$$"
  download "${url}" "${tmp}"

  if [[ "${SKIP_VERIFY}" -eq 0 ]]; then
    sum=$(read_checksum "${name}" "${CHECKSUMS_TSV}" || true)
    if ! verify_file "${tmp}" "${sum}"; then
      rm -f "${tmp}"
      echo "[fetch-sketches] Aborting due to checksum error." >&2
      exit 3
    fi
  fi

  mv -f "${tmp}" "${dst}"
  echo "[fetch-sketches] OK: ${dst}"
done

echo "[fetch-sketches] All sketches downloaded."
