#!/usr/bin/env bash
# Wrapper to execute CAMITAX within the HYMET benchmark harness.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
CAMITAX_DB="${CAMITAX_DB:-}"
CAMITAX_PIPELINE="${CAMITAX_PIPELINE:-CAMI-challenge/CAMITAX}"
CAMITAX_PROFILE="${CAMITAX_PROFILE:-docker}"
CAMITAX_INPUT_EXT="${CAMITAX_INPUT_EXT:-fna}"
CAMITAX_EXTRA_OPTS="${CAMITAX_EXTRA_OPTS:-}"
CAMITAX_REPORT_NAME="${CAMITAX_REPORT_NAME:-camitax.tsv}"
NEXTFLOW_CMD="${NEXTFLOW_CMD:-nextflow}"
KEEP_WORK="${KEEP_CAMITAX_WORK:-0}"

usage(){
  cat <<'USAGE'
Usage: run_camitax.sh --sample ID --contigs FASTA [--out DIR] [--threads N] [--db PATH]

Environment variables:
  CAMITAX_DB           Path to the CAMITAX database directory (required).
  CAMITAX_PIPELINE     Nextflow pipeline identifier (default: CAMI-challenge/CAMITAX).
  CAMITAX_PROFILE      Nextflow profile to use (default: docker).
  CAMITAX_INPUT_EXT    File suffix for input genomes (default: fna).
  CAMITAX_EXTRA_OPTS   Additional options appended to the nextflow command.
  CAMITAX_REPORT_NAME  Expected TSV file emitted by CAMITAX (default: camitax.tsv).
  NEXTFLOW_CMD         Nextflow executable (default: nextflow).
  KEEP_CAMITAX_WORK    If set to 1, retain intermediate work directories.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --db) CAMITAX_DB="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/camitax}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"
INPUT_DIR="${RUN_DIR}/input"
ensure_dir "${INPUT_DIR}"
NF_WORK_DIR="${RUN_DIR}/work"
CAMITAX_DATA_DIR="${RUN_DIR}/data"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "Input contigs FASTA missing (${CONTIGS_ABS})"

CAMITAX_DB="$(resolve_path "${CAMITAX_DB}")"
[[ -d "${CAMITAX_DB}" ]] || die "CAMITAX database directory not found (${CAMITAX_DB}). Set --db or CAMITAX_DB."

need_cmd "${NEXTFLOW_CMD}"

rm -rf "${CAMITAX_DATA_DIR}" "${NF_WORK_DIR}"

INPUT_FASTA="${INPUT_DIR}/${SAMPLE}.${CAMITAX_INPUT_EXT}"
ln -sf "${CONTIGS_ABS}" "${INPUT_FASTA}"

log "Running CAMITAX for ${SAMPLE}"
log "[debug] pipeline=${CAMITAX_PIPELINE}, profile=${CAMITAX_PROFILE}, threads=${THREADS}"

export NXF_DEFAULT_CACHEDIR="${NF_WORK_DIR}"

cmd=(
  "${NEXTFLOW_CMD}" run "${CAMITAX_PIPELINE}"
  -work-dir "${NF_WORK_DIR}"
  -profile "${CAMITAX_PROFILE}"
  --db "${CAMITAX_DB}"
  --i "${INPUT_DIR}"
  --x "${CAMITAX_INPUT_EXT}"
)
if [[ -n "${CAMITAX_EXTRA_OPTS}" ]]; then
  # shellcheck disable=SC2206
  extra=( ${CAMITAX_EXTRA_OPTS} )
  cmd+=("${extra[@]}")
fi

set +e
pushd "${RUN_DIR}" >/dev/null
"${cmd[@]}"
status=$?
popd >/dev/null
set -e
[[ ${status} -eq 0 ]] || die "CAMITAX execution failed with exit code ${status}"

REPORT_CANDIDATE="${CAMITAX_DATA_DIR}/${CAMITAX_REPORT_NAME}"
if [[ ! -s "${REPORT_CANDIDATE}" ]]; then
  REPORT_CANDIDATE="$(find "${CAMITAX_DATA_DIR}" -maxdepth 2 -type f -name "*.tsv" | head -n 1 || true)"
fi
[[ -n "${REPORT_CANDIDATE}" && -s "${REPORT_CANDIDATE}" ]] || die "Failed to locate CAMITAX report under ${CAMITAX_DATA_DIR}"

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"

python3 "${SCRIPT_DIR}/convert/camitax_to_cami.py" \
  --input "${REPORT_CANDIDATE}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool camitax

if [[ "${KEEP_WORK}" -ne 1 ]]; then
  rm -rf "${NF_WORK_DIR}" || true
  rm -rf "${CAMITAX_DATA_DIR}" || true
fi

cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "camitax",
  "profile": "${PROFILE_DST}",
  "input_fasta": "${CONTIGS_ABS}",
  "threads": "${THREADS}"
}
EOF

normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
