#!/usr/bin/env bash
# Run CAMI benchmark across samples and tools with optional database builds.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

MANIFEST="${SCRIPT_DIR}/cami_manifest.tsv"
TOOLS_REQUEST="all"
BUILD_DBS=1
THREADS="${THREADS:-8}"
RESUME=0
MAX_SAMPLES=0

usage(){
  cat <<'USAGE'
Usage: run_all_cami.sh [--manifest TSV] [--tools list] [--no-build] [--threads N] [--resume] [--max-samples N]

Options:
  --manifest PATH     Manifest TSV (default: bench/cami_manifest.tsv)
  --tools LIST        Comma-separated tools (default: all)
  --no-build          Skip database build step
  --threads N         Thread count passed to runners (default: env THREADS or 8)
  --resume            Keep existing runtime log and outputs (default: overwrite runtime log)
  --max-samples N     Limit number of samples processed from manifest (0 = all)
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --tools) TOOLS_REQUEST="$2"; shift 2;;
    --no-build) BUILD_DBS=0; shift;;
    --threads) THREADS="$2"; shift 2;;
    --resume) RESUME=1; shift;;
    --max-samples) MAX_SAMPLES="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"
[[ -s "${MANIFEST}" ]] || die "Manifest not found: ${MANIFEST}"

DEFAULT_TOOLS=(hymet kraken2 centrifuge ganon2 sourmash_gather metaphlan4 camitax phabox phyloflash viwrap)
declare -A TOOL_SCRIPTS=(
  [hymet]="${SCRIPT_DIR}/run_hymet.sh"
  [kraken2]="${SCRIPT_DIR}/run_kraken2.sh"
  [centrifuge]="${SCRIPT_DIR}/run_centrifuge.sh"
  [ganon2]="${SCRIPT_DIR}/run_ganon2.sh"
  [sourmash_gather]="${SCRIPT_DIR}/run_sourmash_gather.sh"
  [metaphlan4]="${SCRIPT_DIR}/run_metaphlan4.sh"
  [camitax]="${SCRIPT_DIR}/run_camitax.sh"
  [basta]="${SCRIPT_DIR}/run_basta.sh"
  [tama]="${SCRIPT_DIR}/run_tama.sh"
  [phabox]="${SCRIPT_DIR}/run_phabox.sh"
  [phyloflash]="${SCRIPT_DIR}/run_phyloflash.sh"
  [viwrap]="${SCRIPT_DIR}/run_viwrap.sh"
)
declare -A TOOL_BUILDERS=(
  [kraken2]="${SCRIPT_DIR}/db/build_kraken2.sh"
  [centrifuge]="${SCRIPT_DIR}/db/build_centrifuge.sh"
  [ganon2]="${SCRIPT_DIR}/db/build_ganon2.sh"
  [sourmash_gather]="${SCRIPT_DIR}/db/build_sourmash.sh"
)

IFS=',' read -r -a TOOLS <<< "${TOOLS_REQUEST}"
if [[ ${#TOOLS[@]} -eq 1 && "${TOOLS[0]}" == "all" ]]; then
  TOOLS=(${DEFAULT_TOOLS[@]})
fi

MEASURE="${SCRIPT_DIR}/lib/measure.sh"
[[ -x "${MEASURE}" ]] || die "measure.sh not executable: ${MEASURE}"

OUT_ROOT="${SCRIPT_DIR}/out"
ensure_dir "${OUT_ROOT}"
if [[ ${RESUME} -eq 0 ]]; then
  rm -f "${OUT_ROOT}/runtime_memory.tsv"
fi

if [[ ${BUILD_DBS} -eq 1 ]]; then
  for tool in "${TOOLS[@]}"; do
    builder="${TOOL_BUILDERS[$tool]:-}"
    if [[ -n "${builder}" ]]; then
      log "Ensuring database for ${tool}"
      if [[ ! -x "${builder}" ]]; then
        die "Builder script missing for ${tool}: ${builder}"
      fi
      bash "${builder}" || die "Database build failed for ${tool}"
    fi
  done
fi

processed=0
while IFS=$'\t' read -r sample_id contigs truth_contigs truth_profile rest; do
  [[ -z "${sample_id}" || "${sample_id}" =~ ^# ]] && continue
  if [[ "${sample_id}" == "sample_id" ]]; then
    continue
  fi
  if [[ ${MAX_SAMPLES} -gt 0 && ${processed} -ge ${MAX_SAMPLES} ]]; then
    break
  fi
  processed=$((processed+1))

  SAMPLE_DIR="${OUT_ROOT}/${sample_id}"
  ensure_dir "${SAMPLE_DIR}"

  contigs_abs="$(resolve_path "${contigs}")"
  truth_profile_abs="$(resolve_path "${truth_profile}")"
  truth_contigs_abs="$(resolve_path "${truth_contigs}")"

  if [[ ! -s "${contigs_abs}" ]]; then
    log "WARNING: sample ${sample_id} missing contigs at ${contigs_abs}; skipping"
    continue
  fi

  for tool in "${TOOLS[@]}"; do
    script="${TOOL_SCRIPTS[$tool]:-}"
    if [[ -z "${script}" || ! -x "${script}" ]]; then
      log "WARNING: tool ${tool} not configured; skipping"
      continue
    fi

    log "[sample=${sample_id}] Running tool ${tool}"
    run_cmd=("${script}" --sample "${sample_id}" --contigs "${contigs_abs}" --threads "${THREADS}")
    "${MEASURE}" --tool "${tool}" --sample "${sample_id}" --stage run -- "${run_cmd[@]}" || log "WARNING: ${tool} failed for ${sample_id}"

    TOOL_DIR="${OUT_ROOT}/${sample_id}/${tool}"
    pred_profile="${TOOL_DIR}/profile.cami.tsv"
    pred_contigs=""
    pred_paf=""
    case "${tool}" in
      hymet)
        pred_contigs="${TOOL_DIR}/classified_sequences.tsv"
        pred_paf="${TOOL_DIR}/resultados.paf"
        ;;
      kraken2|centrifuge|ganon2)
        pred_contigs="${TOOL_DIR}/classified_sequences.tsv"
        ;;
    esac

    if [[ ! -s "${pred_profile}" ]]; then
      log "WARNING: ${tool} produced no profile for ${sample_id}; skipping evaluation"
      continue
    fi

    if [[ ! -s "${truth_profile_abs}" ]]; then
      log "WARNING: truth profile missing for ${sample_id}; skipping evaluation"
      continue
    fi

    eval_cmd=(
      "${SCRIPT_DIR}/lib/run_eval.sh"
      --sample "${sample_id}"
      --tool "${tool}"
      --pred-profile "${pred_profile}"
      --truth-profile "${truth_profile_abs}"
      --pred-contigs "${pred_contigs}"
      --truth-contigs "${truth_contigs_abs}"
      --pred-fasta "${contigs_abs}"
      --paf "${pred_paf}"
      --threads "${THREADS}"
    )
    "${MEASURE}" --tool "${tool}" --sample "${sample_id}" --stage eval -- "${eval_cmd[@]}" || log "WARNING: evaluation failed for ${tool}/${sample_id}"
  done
done < "${MANIFEST}"

if [[ -s "${OUT_ROOT}/runtime_memory.tsv" ]]; then
  log "Aggregating metrics"
  python3 "${SCRIPT_DIR}/aggregate_metrics.py" --bench-root "${SCRIPT_DIR}" --outdir "out"
  python3 "${SCRIPT_DIR}/plot/make_figures.py" --bench-root "${SCRIPT_DIR}" --outdir "out" || log "WARNING: plotting step failed"
fi

log "Benchmark completed. Outputs under ${OUT_ROOT}"
