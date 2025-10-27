#!/usr/bin/env bash
# Download the contig FASTA assets referenced in case/manifest.tsv.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_ROOT="${SCRIPT_DIR}"
DEST_ROOT="/data/case"

usage(){
  cat <<'USAGE'
Usage: fetch_case_data.sh [--dest DIR] [sample ...]

Samples:
  zymo_mc   ZymoBIOMICS mock community (Loman Lab assembly)
  gut_case  MGnify human stool metagenome (Processed contigs)

Specify one or more sample IDs or omit to fetch all.
USAGE
  exit 1
}

ensure_dir(){
  mkdir -p "$1"
}

download(){
  local url="$1"
  local dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "[fetch] Exists: ${dest}" >&2
    return 0
  fi
  echo "[fetch] Downloading ${url} → ${dest}" >&2
  ensure_dir "$(dirname "${dest}")"
  curl -L --fail --retry 3 --retry-delay 2 -o "${dest}" "${url}"
}

SHA256(){
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    echo "(sha256 tool not available)"
  fi
}

DEST="${DEST_ROOT}"
ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest)
      [[ $# -ge 2 ]] || usage
      DEST="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    *)
      ARGS+=("$1")
      shift
      ;;
  esac
done

ensure_dir "${DEST}"

fetch_zymo(){
  local target="${DEST}/zymo_mc_contigs.fna"
  local url="http://nanopore.s3.climb.ac.uk/mockcommunity/v3/7cd60d3b-eafb-48d1-9aab-c8701232f2f8.ctg.cns.fa"
  download "${url}" "${target}"
  if [[ -s "${target}" ]]; then
    echo "[fetch] zymo_mc sha256: $(SHA256 "${target}")"
  fi
}

fetch_gut(){
  local gz="${DEST}/gut_case_contigs.fna.gz"
  local dest="${DEST}/gut_case_contigs.fna"
  local url="https://www.ebi.ac.uk/metagenomics/api/v1/analyses/MGYA00794604/file/ERZ24911249_FASTA.fasta.gz"
  if [[ ! -s "${dest}" ]]; then
    download "${url}" "${gz}"
    echo "[fetch] Decompressing ${gz}" >&2
    gunzip -f "${gz}"
  else
    echo "[fetch] Exists: ${dest}" >&2
  fi
  if [[ -s "${dest}" ]]; then
    echo "[fetch] gut_case sha256: $(SHA256 "${dest}")"
  fi
}

samples=()
if [[ ${#ARGS[@]} -eq 0 ]]; then
  samples=(zymo_mc gut_case)
else
  samples=(${ARGS[@]})
fi

for sample in "${samples[@]}"; do
  case "${sample}" in
    zymo_mc)
      fetch_zymo
      ;;
    gut_case)
      fetch_gut
      ;;
    *)
      echo "[fetch] Unknown sample: ${sample}" >&2
      usage
      ;;
  esac
done

echo "[fetch] All done. Files staged under ${DEST}"
