#!/usr/bin/env bash
# Idempotent helper to fetch CAMI assets referenced in the manifest.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${SCRIPT_DIR}/cami_manifest.tsv"
DRY_RUN=0
MAX_SAMPLES=0

usage(){
  cat <<'USAGE'
Usage: fetch_cami.sh [--manifest TSV] [--dry-run] [--max-samples N]

The manifest may optionally include columns with URLs named:
  contigs_url, truth_contigs_url, truth_profile_url
Only missing files with matching URLs are downloaded. Existing files are left untouched.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    --max-samples) MAX_SAMPLES="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"
MANIFEST_DIR="$(dirname "${MANIFEST}")"
[[ -s "${MANIFEST}" ]] || die "Manifest not found: ${MANIFEST}"

python3 - "$MANIFEST" "$MAX_SAMPLES" <<'PY' | while IFS=$'	' read -r path url; do
import csv, sys
manifest = sys.argv[1]
max_samples = int(sys.argv[2]) if sys.argv[2] and sys.argv[2] != "0" else None
rows = []
with open(manifest, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for idx, row in enumerate(reader):
        if max_samples is not None and idx >= max_samples:
            break
        paths = {
            "contigs_fa": row.get("contigs_fa", ""),
            "truth_contigs_tsv": row.get("truth_contigs_tsv", ""),
            "truth_profile_tsv": row.get("truth_profile_tsv", ""),
        }
        urls = {
            "contigs_fa": row.get("contigs_url", ""),
            "truth_contigs_tsv": row.get("truth_contigs_url", ""),
            "truth_profile_tsv": row.get("truth_profile_url", ""),
        }
        for key, path in paths.items():
            if not path:
                continue
            url = urls.get(key) or ""
            print(f"{path}\t{url}")
PY
  [[ -z "${path}" ]] && continue
  abs_path="$(resolve_path "${path}" "${MANIFEST_DIR}")"
  url="${url}"
  if [[ -s "${abs_path}" ]]; then
    log "exists: ${abs_path}"
    continue
  fi
  if [[ -z "${url}" ]]; then
    log "missing and no URL: ${abs_path}"
    continue
  fi
  ensure_dir "$(dirname "${abs_path}")"
  if [[ ${DRY_RUN} -eq 1 ]]; then
    log "would download ${url} -> ${abs_path}"
    continue
  fi
  log "downloading ${url}"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 4 -s 4 -o "${abs_path}" "${url}" || die "aria2c failed for ${url}"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 -o "${abs_path}" "${url}" || die "curl failed for ${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${abs_path}" "${url}" || die "wget failed for ${url}"
  else
    die "No downloader (aria2c/curl/wget) available for ${url}"
  fi
done
