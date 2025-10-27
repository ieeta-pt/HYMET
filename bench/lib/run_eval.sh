#!/usr/bin/env bash
# Helper to run HYMET/tools/eval_cami.py with consistent paths.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common.sh"

SAMPLE=""
TOOL=""
PRED_PROFILE=""
PRED_CONTIGS=""
PRED_FASTA=""
TRUTH_PROFILE=""
TRUTH_CONTIGS=""
TRUTH_FASTA=""
PAF=""
OUTDIR=""
THREADS="${THREADS:-8}"

usage(){
  cat <<'EOF'
run_eval.sh --sample ID --tool TOOL --pred-profile FILE --truth-profile FILE [options]
Options:
  --pred-contigs FILE
  --pred-fasta FILE
  --truth-contigs FILE
  --truth-fasta FILE
  --paf FILE
  --outdir DIR            (defaults to bench/out/<sample>/<tool>/eval)
  --threads N
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --tool) TOOL="$2"; shift 2;;
    --pred-profile) PRED_PROFILE="$2"; shift 2;;
    --pred-contigs) PRED_CONTIGS="$2"; shift 2;;
    --pred-fasta) PRED_FASTA="$2"; shift 2;;
    --truth-profile) TRUTH_PROFILE="$2"; shift 2;;
    --truth-contigs) TRUTH_CONTIGS="$2"; shift 2;;
    --truth-fasta) TRUTH_FASTA="$2"; shift 2;;
    --paf) PAF="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) usage;;
  esac
done

if [[ -z "${SAMPLE}" || -z "${TOOL}" || -z "${PRED_PROFILE}" || -z "${TRUTH_PROFILE}" ]]; then
  usage
fi

if [[ -z "${OUTDIR}" ]]; then
  if [[ -n "${BENCH_OUT_ROOT:-}" ]]; then
    OUTDIR="$(resolve_path "${BENCH_OUT_ROOT}")/${SAMPLE}/${TOOL}/eval"
  else
    OUTDIR="${BENCH_ROOT}/out/${SAMPLE}/${TOOL}/eval"
  fi
fi

ensure_dir "${OUTDIR}"

PRED_PROFILE="$(resolve_path "${PRED_PROFILE}")"
TRUTH_PROFILE="$(resolve_path "${TRUTH_PROFILE}")"
PRED_CONTIGS="$(resolve_path "${PRED_CONTIGS}")"
TRUTH_CONTIGS="$(resolve_path "${TRUTH_CONTIGS}")"
PRED_FASTA="$(resolve_path "${PRED_FASTA}")"
TRUTH_FASTA="$(resolve_path "${TRUTH_FASTA}")"
PAF="$(resolve_path "${PAF}")"

if [[ ! -s "${PRED_PROFILE}" ]]; then
  log "WARNING: ${TOOL}/${SAMPLE} missing predicted profile (${PRED_PROFILE}); skipping evaluation."
  exit 0
fi

if [[ ! -s "${TRUTH_PROFILE}" ]]; then
  log "WARNING: ${TOOL}/${SAMPLE} missing truth profile (${TRUTH_PROFILE}); skipping evaluation."
  exit 0
fi

# Align predicted superkingdom taxids with the truth profile before scoring.
FIX_SUPERKINGDOM="${BENCH_ROOT}/tools/fix_superkingdom_taxids.py"
if [[ -f "${FIX_SUPERKINGDOM}" ]]; then
  python3 "${FIX_SUPERKINGDOM}" --profile "${PRED_PROFILE}" --truth-profile "${TRUTH_PROFILE}" || log "WARNING: superkingdom normalization failed for ${TOOL}/${SAMPLE}"
fi

log "Evaluating ${TOOL} for ${SAMPLE}"

python3 "${HYMET_ROOT}/tools/eval_cami.py" \
  --pred-profile "${PRED_PROFILE}" \
  --truth-profile "${TRUTH_PROFILE}" \
  --pred-contigs "${PRED_CONTIGS}" \
  --truth-contigs "${TRUTH_CONTIGS}" \
  --pred-fasta "${PRED_FASTA}" \
  --truth-fasta "${TRUTH_FASTA}" \
  --paf "${PAF}" \
  --outdir "${OUTDIR}" \
  --threads "${THREADS}"
