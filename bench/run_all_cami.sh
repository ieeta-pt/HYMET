#!/usr/bin/env bash
# Run CAMI benchmark across samples and tools with optional database builds.

# Re-exec with Bash if the script was invoked by /bin/sh or another shell.
if [ -z "${BASH_VERSION:-}" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

set -Eeuo pipefail

if (( BASH_VERSINFO[0] < 4 )); then
  echo "ERROR: bench/run_all_cami.sh requires Bash >= 4" >&2
  exit 1
fi

# Ensure deterministic Python behavior for any helper scripts invoked
export PYTHONHASHSEED="0"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
export PATH="${SCRIPT_DIR}/bin:${PATH}"

MANIFEST="${SCRIPT_DIR}/cami_manifest.tsv"
TOOLS_REQUEST="all"
BUILD_DBS=1
THREADS="${THREADS:-8}"
RESUME=0
MAX_SAMPLES=0
SUITE_SCENARIO="${PUBLISH_SCENARIO:-cami}"
SUITE_NAME="${PUBLISH_SUITE:-canonical}"
SUITE_NAME_SOURCE="default"
if [[ -n "${PUBLISH_SUITE:-}" ]]; then
  SUITE_NAME_SOURCE="env"
fi
SUITE_PATH=""
PUBLISH_RESULTS=1

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
  --scenario NAME     Results namespace (default: cami)
  --suite NAME        Suite name within the scenario (default: canonical)
  --suite-path REL    Explicit path under results/ (overrides scenario/suite)
  --no-publish        Skip publishing a run directory under results/
USAGE
  exit 1
}

normalize_tool_list(){
  if [[ $# -eq 0 ]]; then
    echo ""
    return
  fi
  printf '%s\n' "$@" | LC_ALL=C sort -u | paste -sd',' -
}

infer_suite_from_tools(){
  local requested="$1"
  shift || true
  local actual_norm
  actual_norm="$(normalize_tool_list "$@")"
  local contig_norm
  contig_norm="$(normalize_tool_list "${DEFAULT_CONTIG_TOOLS[@]}")"
  local canonical_norm
  canonical_norm="$(normalize_tool_list "${DEFAULT_FULL_PANEL[@]}")"
  case "${requested}" in
    contigs)
      echo "contig_full"
      return
      ;;
    ""|all)
      echo "canonical"
      return
      ;;
  esac
  if [[ -n "${actual_norm}" && "${actual_norm}" == "${contig_norm}" ]]; then
    echo "contig_full"
    return
  fi
  if [[ -n "${actual_norm}" && "${actual_norm}" == "${canonical_norm}" ]]; then
    echo "canonical"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2;;
    --tools) TOOLS_REQUEST="$2"; shift 2;;
    --no-build) BUILD_DBS=0; shift;;
    --threads) THREADS="$2"; shift 2;;
    --resume) RESUME=1; shift;;
    --max-samples) MAX_SAMPLES="$2"; shift 2;;
    --scenario) SUITE_SCENARIO="$2"; shift 2;;
    --suite) SUITE_NAME="$2"; SUITE_NAME_SOURCE="cli"; shift 2;;
    --suite-path) SUITE_PATH="$2"; shift 2;;
    --no-publish) PUBLISH_RESULTS=0; shift;;
    -h|--help) usage;;
    *) usage;;
  esac
done

MANIFEST="$(resolve_path "${MANIFEST}")"
MANIFEST_DIR="$(dirname "${MANIFEST}")"
[[ -s "${MANIFEST}" ]] || die "Manifest not found: ${MANIFEST}"

DEFAULT_CONTIG_TOOLS=(hymet kraken2 centrifuge ganon2 viwrap tama squeezemeta megapath_nano)
DEFAULT_FULL_PANEL=(hymet kraken2 centrifuge ganon2 sourmash_gather metaphlan4 camitax phabox phyloflash viwrap squeezemeta megapath_nano snakemags)
declare -A TOOL_SCRIPTS=(
  [hymet]="${SCRIPT_DIR}/run_hymet.sh"
  [hymet_reads]="${SCRIPT_DIR}/run_hymet.sh"
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
  [squeezemeta]="${SCRIPT_DIR}/run_squeezemeta.sh"
  [megapath_nano]="${SCRIPT_DIR}/run_megapath_nano.sh"
  [snakemags]="${SCRIPT_DIR}/run_snakemags.sh"
)
declare -A TOOL_BUILDERS=(
  [kraken2]="${SCRIPT_DIR}/db/build_kraken2.sh"
  [centrifuge]="${SCRIPT_DIR}/db/build_centrifuge.sh"
  [ganon2]="${SCRIPT_DIR}/db/build_ganon2.sh"
  [sourmash_gather]="${SCRIPT_DIR}/db/build_sourmash.sh"
)

IFS=',' read -r -a TOOLS <<< "${TOOLS_REQUEST}"
if [[ ${#TOOLS[@]} -eq 1 ]]; then
  case "${TOOLS[0]}" in
    ""|all)
      TOOLS=("${DEFAULT_FULL_PANEL[@]}")
      ;;
    contigs)
      TOOLS=("${DEFAULT_CONTIG_TOOLS[@]}")
      ;;
    reads)
      TOOLS=("hymet_reads")
      ;;
  esac
fi

if [[ "${SUITE_NAME_SOURCE}" == "default" && -z "${SUITE_PATH}" ]]; then
  suite_guess="$(infer_suite_from_tools "${TOOLS_REQUEST}" "${TOOLS[@]}")"
  if [[ -n "${suite_guess}" && "${suite_guess}" != "${SUITE_NAME}" ]]; then
    tools_label="${TOOLS_REQUEST:-$(IFS=','; echo "${TOOLS[*]}")}"
    log "Auto-selecting suite ${suite_guess} for tools: ${tools_label}"
    SUITE_NAME="${suite_guess}"
  fi
fi

MEASURE="${SCRIPT_DIR}/lib/measure.sh"
[[ -x "${MEASURE}" ]] || die "measure.sh not executable: ${MEASURE}"

OUT_ROOT="${BENCH_OUT_ROOT:-${SCRIPT_DIR}/out}"
RUNTIME_TSV="${OUT_ROOT}/runtime_memory.tsv"
ensure_dir "${OUT_ROOT}"
if [[ ${RESUME} -eq 0 ]]; then
  rm -f "${RUNTIME_TSV}"
fi

# Persist a copy of the manifest used for this run
cp -f "${MANIFEST}" "${OUT_ROOT}/manifest.snapshot.tsv" || true

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

  contigs_abs="$(resolve_path "${contigs}" "${MANIFEST_DIR}")"
  truth_profile_abs="$(resolve_path "${truth_profile}" "${MANIFEST_DIR}")"
  truth_contigs_abs="$(resolve_path "${truth_contigs}" "${MANIFEST_DIR}")"

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
    TOOL_DIR="${SAMPLE_DIR}/${tool}"
    ensure_dir "${TOOL_DIR}"
    TOOL_RUNTIME="${TOOL_DIR}/runtime_memory.tsv"
    run_cmd=("${script}" --sample "${sample_id}" --contigs "${contigs_abs}" --threads "${THREADS}")
    case "${tool}" in
      hymet)
        run_cmd+=("--out" "${TOOL_DIR}")
        ;;
      hymet_reads)
        run_cmd+=("--out" "${TOOL_DIR}" "--mode" "reads")
        ;;
    esac
    "${MEASURE}" \
      --tool "${tool}" \
      --sample "${sample_id}" \
      --stage run \
      --out "${RUNTIME_TSV}" \
      --local "${TOOL_RUNTIME}" \
      -- "${run_cmd[@]}" || log "WARNING: ${tool} failed for ${sample_id}"

    pred_profile="${TOOL_DIR}/profile.cami.tsv"
    pred_contigs=""
    pred_paf=""
    case "${tool}" in
      hymet|hymet_reads)
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
    "${MEASURE}" \
      --tool "${tool}" \
      --sample "${sample_id}" \
      --stage eval \
      --out "${RUNTIME_TSV}" \
      --local "${TOOL_RUNTIME}" \
      -- "${eval_cmd[@]}" || log "WARNING: evaluation failed for ${tool}/${sample_id}"
  done
done < "${MANIFEST}"

if [[ -s "${RUNTIME_TSV}" ]]; then
  log "Aggregating metrics"
  python3 "${SCRIPT_DIR}/aggregate_metrics.py" --bench-root "${SCRIPT_DIR}" --outdir "out"
  python3 "${SCRIPT_DIR}/plot/make_figures.py" --bench-root "${SCRIPT_DIR}" --outdir "out" || log "WARNING: plotting step failed"

  if [[ ${PUBLISH_RESULTS} -eq 1 ]]; then
    if [[ -n "${SUITE_PATH}" ]]; then
      if [[ "${SUITE_PATH}" = /* ]]; then
        run_dir_base="${SUITE_PATH}"
      else
        run_dir_base="${REPO_ROOT}/results/${SUITE_PATH}"
      fi
    else
      run_dir_base="${REPO_ROOT}/results/${SUITE_SCENARIO}/${SUITE_NAME}"
    fi
    RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
    export PUBLISH_RUN_DIR="${run_dir_base}/run_${RUN_STAMP}"
    export PUBLISH_RUN_STAMP="${RUN_STAMP}"
    export PUBLISH_SCENARIO="${SUITE_SCENARIO}"
    export PUBLISH_SUITE="${SUITE_NAME}"
    export PUBLISH_MANIFEST="${MANIFEST}"
    export PUBLISH_THREADS="${THREADS}"
    export PUBLISH_TOOLS="$(IFS=','; echo "${TOOLS[*]}")"
    export PUBLISH_MODES="contigs"
    SKIP_RECOMPUTE=1 "${SCRIPT_DIR}/publish_results.sh" || log "WARNING: failed to publish results snapshot"
  else
    log "Skipping publish step (user disabled)"
  fi
fi

log "Benchmark completed. Outputs under ${OUT_ROOT}"
