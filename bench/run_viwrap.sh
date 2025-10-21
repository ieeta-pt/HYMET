#!/usr/bin/env bash
# Wrapper to execute geNomad (ViWrap) on assembled contigs within HYMET.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
VIWRAP_ENV_PREFIX="${VIWRAP_ENV_PREFIX:-/opt/envs/genomad}"
VIWRAP_DB_DIR="${VIWRAP_DB_DIR:-}"
VIWRAP_CMD="${VIWRAP_CMD:-genomad}"
VIWRAP_EXTRA_OPTS="${VIWRAP_EXTRA_OPTS:-}"
VIWRAP_SCORE_CUTOFF="${VIWRAP_SCORE_CUTOFF:-0.5}"

usage(){
  cat <<'USAGE'
Usage: run_viwrap.sh --sample ID --contigs FASTA --db DIR [--out DIR] [--threads N]

Environment variables:
  VIWRAP_ENV_PREFIX    Path to geNomad conda/mamba environment (default: /opt/envs/genomad).
  VIWRAP_DB_DIR        geNomad database directory (required).
  VIWRAP_CMD           geNomad executable (default: genomad).
  VIWRAP_EXTRA_OPTS    Extra options appended to `genomad end-to-end`.
  VIWRAP_SCORE_CUTOFF  Minimum virus_score retained for CAMI export (default: 0.5).
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --db) VIWRAP_DB_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/viwrap}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "Input contigs FASTA missing (${CONTIGS_ABS})"

VIWRAP_DB_DIR="$(resolve_path "${VIWRAP_DB_DIR}")"
[[ -d "${VIWRAP_DB_DIR}" ]] || die "geNomad database directory not found (${VIWRAP_DB_DIR})"

[[ -d "${VIWRAP_ENV_PREFIX}" ]] || die "geNomad environment prefix not found (${VIWRAP_ENV_PREFIX})"
need_cmd micromamba

run_in_env(){
  micromamba run -p "${VIWRAP_ENV_PREFIX}" "$@"
}

LIB_NAME="${SAMPLE//[^A-Za-z0-9_]/_}"
[[ -n "${LIB_NAME}" ]] || LIB_NAME="sample"
INPUT_FASTA="${RUN_DIR}/${LIB_NAME}.fna"

if [[ ! -e "${INPUT_FASTA}" ]]; then
  ln -sf "${CONTIGS_ABS}" "${INPUT_FASTA}"
fi

SUMMARY_DIR="${RUN_DIR}/${LIB_NAME}_summary"
SUMMARY_FILE="${SUMMARY_DIR}/${LIB_NAME}_virus_summary.tsv"

if [[ ! -s "${SUMMARY_FILE}" ]]; then
  log "Running geNomad for ${SAMPLE}"
  GENOMAD_ARGS=(
    "${VIWRAP_CMD}"
    end-to-end
    "${INPUT_FASTA}"
    "${RUN_DIR}"
    "${VIWRAP_DB_DIR}"
    --threads "${THREADS}"
  )
  if [[ -n "${VIWRAP_EXTRA_OPTS}" ]]; then
    # shellcheck disable=SC2206
    EXTRA=( ${VIWRAP_EXTRA_OPTS} )
    GENOMAD_ARGS+=("${EXTRA[@]}")
  fi
  run_in_env "${GENOMAD_ARGS[@]}"
else
  log "Reusing existing geNomad outputs for ${SAMPLE}"
fi

[[ -s "${SUMMARY_FILE}" ]] || die "geNomad virus summary not found (${SUMMARY_FILE})"

if [[ -z "${TAXONKIT_DB:-}" && -d "${HYMET_ROOT}/taxonomy_files" ]]; then
  export TAXONKIT_DB="${HYMET_ROOT}/taxonomy_files"
fi

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"

python3 "${SCRIPT_DIR}/convert/viwrap_to_cami.py" \
  --input "${SUMMARY_FILE}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool viwrap \
  --taxdb "${TAXONKIT_DB:-}" \
  --score-cutoff "${VIWRAP_SCORE_CUTOFF}" \
  --classified-out "${CLASSIFIED_DST}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "viwrap", "profile": "${PROFILE_DST}", "contigs": "${CONTIGS_ABS}", "db_dir": "${VIWRAP_DB_DIR}", "threads": "${THREADS}", "score_cutoff": "${VIWRAP_SCORE_CUTOFF}"}
EOF
