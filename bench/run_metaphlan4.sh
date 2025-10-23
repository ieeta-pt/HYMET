#!/usr/bin/env bash
# Run MetaPhlAn4 on contigs and convert output to CAMI format.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
THREADS="${METAPHLAN_THREADS:-${THREADS}}"
EXTRA_OPTS="${METAPHLAN_OPTS:-}"
TAXDIR="${TAXDIR:-${HYMET_ROOT}/taxonomy_files}"
METAPHLAN_CMD="${METAPHLAN_CMD:-metaphlan}"
DB_DIR="${METAPHLAN_DB_DIR:-${BENCH_ROOT}/db/metaphlan}"
INDEX_NAME="${METAPHLAN_INDEX:-mpa_vJun23_CHOCOPhlAnSGB_202307}"

usage(){
  cat <<'USAGE'
Usage: run_metaphlan4.sh --sample ID --contigs FASTA [--out DIR] [--threads N] [--opts "..."]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --opts) EXTRA_OPTS="$2"; shift 2;;
    --taxdir) TAXDIR="$2"; shift 2;;
    --cmd) METAPHLAN_CMD="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/metaphlan4}"
ensure_dir "${OUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "MetaPhlAn input FASTA missing (${CONTIGS_ABS})"
fi

if ! command -v "${METAPHLAN_CMD}" >/dev/null 2>&1; then
  die "${METAPHLAN_CMD} executable not found. Please install MetaPhlAn4 and set METAPHLAN_CMD if needed."
fi

if [[ ! -d "${DB_DIR}" ]]; then
  die "MetaPhlAn database directory missing (${DB_DIR})."
fi

if [[ ! -s "${DB_DIR}/${INDEX_NAME}.pkl" ]]; then
  die "MetaPhlAn index ${INDEX_NAME} not found under ${DB_DIR}."
fi

RAW_PROFILE="${OUT_DIR}/metaphlan_profile.tsv"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"
IFS=' ' read -r -a EXTRA_ARGS <<< "${EXTRA_OPTS:-}"
if [[ ${#EXTRA_ARGS[@]} -eq 0 || -z "${EXTRA_ARGS[0]:-}" ]]; then
  EXTRA_ARGS=()
fi
if [[ ${#EXTRA_ARGS[@]} -eq 0 ]]; then
  EXTRA_ARGS=(--split_reads)
fi

run_metaphlan_once(){
  local nproc="$1"; shift
  local label="$1"; shift
  local -a extra=("$@")
  local -a cmd=(
    "${METAPHLAN_CMD}"
    "${CONTIGS_ABS}"
    "--input_type" "fasta"
    "--nproc" "${nproc}"
    "--db_dir" "${DB_DIR}"
    "-x" "${INDEX_NAME}"
  )
  if [[ ${#extra[@]} -gt 0 ]]; then
    cmd+=("${extra[@]}")
  fi
  cmd+=("-o" "${RAW_PROFILE}")
  local extra_pretty="${extra[*]:-∅}"
  log "MetaPhlAn4 ${label} run: threads=${nproc}, extra_opts=${extra_pretty}"
  "${cmd[@]}"
}

log "Running MetaPhlAn4 for ${SAMPLE}"
PRIMARY_THREADS="${THREADS}"
PRIMARY_EXTRA=("${EXTRA_ARGS[@]}")
USED_THREADS="${PRIMARY_THREADS}"
USED_EXTRA_ARGS=("${PRIMARY_EXTRA[@]}")

if ! run_metaphlan_once "${PRIMARY_THREADS}" "primary" "${PRIMARY_EXTRA[@]}"; then
  rm -f "${RAW_PROFILE}"
  FALLBACK_THREADS="${PRIMARY_THREADS}"
  if [[ "${FALLBACK_THREADS}" -gt 4 ]]; then
    FALLBACK_THREADS=4
  fi
  if [[ "${FALLBACK_THREADS}" -lt 1 ]]; then
    FALLBACK_THREADS=1
  fi
  FALLBACK_EXTRA=("${PRIMARY_EXTRA[@]}")
  if [[ " ${FALLBACK_EXTRA[*]} " != *" --split_reads "* ]]; then
    FALLBACK_EXTRA+=(--split_reads)
  fi
  primary_str="${PRIMARY_EXTRA[*]}"
  fallback_str="${FALLBACK_EXTRA[*]}"
  if [[ "${FALLBACK_THREADS}" -eq "${PRIMARY_THREADS}" && "${fallback_str}" == "${primary_str}" ]]; then
    die "MetaPhlAn failed. Try setting METAPHLAN_OPTS=\"--split_reads\" and reducing METAPHLAN_THREADS."
  fi
  log "MetaPhlAn primary run failed; retrying with threads=${FALLBACK_THREADS}, extra_opts=${fallback_str}"
  if ! run_metaphlan_once "${FALLBACK_THREADS}" "fallback" "${FALLBACK_EXTRA[@]}"; then
    die "MetaPhlAn failed after fallback attempt"
  fi
  USED_THREADS="${FALLBACK_THREADS}"
  USED_EXTRA_ARGS=("${FALLBACK_EXTRA[@]}")
fi
USED_EXTRA_OPTS="${USED_EXTRA_ARGS[*]}"

python3 "${SCRIPT_DIR}/convert/metaphlan4_to_cami.py" \
  --input "${RAW_PROFILE}" \
  --out "${PROFILE_CAMI}" \
  --sample-id "${SAMPLE}" \
  --tool "metaphlan4" \
  --taxdb "${TAXDIR}"

METAPHLAN_VER=$(${METAPHLAN_CMD} --version 2>/dev/null | head -n1 | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g' || echo "unknown")
cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "metaphlan4", "profile": "${PROFILE_CAMI}", "raw_profile": "${RAW_PROFILE}", "metaphlan_cmd": "${METAPHLAN_CMD}", "extra_opts": "${USED_EXTRA_OPTS}", "threads": "${USED_THREADS}", "db_dir": "${DB_DIR}", "index": "${INDEX_NAME}", "metaphlan_version": "${METAPHLAN_VER}"}
EOF
