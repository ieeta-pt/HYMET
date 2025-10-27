#!/usr/bin/env bash
# Run all published HYMET case-study suites end-to-end (HYMET + figures + runtime charts).

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HYMET_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CASE_ROOT="${HYMET_ROOT}/case"
BENCH_ROOT="${HYMET_ROOT}/bench"
MAKE_FIGURES="${BENCH_ROOT}/plot/make_figures.py"
CASE_PLOTTER="${CASE_ROOT}/plot_case.py"

THREADS="${THREADS:-16}"
SCENARIO="cases"
SKIP_FIGURES=0
DRY_RUN=0
declare -a REQUESTED_SUITES=()

usage() {
  cat <<'USAGE'
Usage: run_cases_full.sh [--threads N] [--scenario NAME] [--suite NAME ...]
                         [--skip-figures] [--dry-run]

Runs the HYMET case-study workflow for each configured suite (canonical, gut, zymo by default),
publishes the results under results/<scenario>/<suite>/run_<timestamp>/, and regenerates both the
case-specific figures and the shared runtime/memory charts from the produced tables.

Options:
  --threads N        Number of threads for HYMET (default: THREADS env or 16)
  --scenario NAME    Results namespace under results/<scenario>/... (default: cases)
  --suite NAME       Restrict execution to the given suite; repeat for multiple suites.
                     (default: run all known suites - canonical, gut, zymo)
  --skip-figures     Run HYMET but skip post-processing/figure regeneration.
  --dry-run          Print the actions that would be taken without executing.
  -h, --help         Show this help message.
USAGE
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="$2"; shift 2;;
    --scenario) SCENARIO="$2"; shift 2;;
    --suite) REQUESTED_SUITES+=("$2"); shift 2;;
    --suites) IFS=',' read -ra suites <<<"$2"; REQUESTED_SUITES+=("${suites[@]}"); shift 2;;
    --skip-figures) SKIP_FIGURES=1; shift;;
    --dry-run) DRY_RUN=1; shift;;
    -h|--help) usage;;
    *) echo "Unknown option: $1" >&2; usage;;
  esac
done

declare -A SUITE_MANIFESTS=(
  [canonical]="${CASE_ROOT}/manifest.tsv"
  [gut]="${CASE_ROOT}/manifest_gut.tsv"
  [zymo]="${CASE_ROOT}/manifest_zymo.tsv"
)

if [[ ${#REQUESTED_SUITES[@]} -eq 0 ]]; then
  REQUESTED_SUITES=(canonical gut zymo)
fi

log() { printf '[cases] %s\n' "$*" >&2; }

ensure_file() {
  local path="$1"
  [[ -s "$path" ]] || { log "ERROR: required file missing: $path"; return 1; }
}

run_suite() {
  local suite="$1"
  local manifest="${SUITE_MANIFESTS[$suite]}"
  if [[ -z "$manifest" ]]; then
    log "Skipping unknown suite '${suite}'. Known suites: ${!SUITE_MANIFESTS[*]}"
    return 0
  fi
  ensure_file "$manifest" || return 1

  log "Running HYMET case suite '${suite}' (manifest: ${manifest})"
  local run_cmd=(env CASE_SUITE="$suite" CASE_SCENARIO="$SCENARIO" THREADS="$THREADS" "${CASE_ROOT}/run_case.sh" --manifest "$manifest")
  if [[ $DRY_RUN -eq 1 ]]; then
    printf '[dry-run]'
    for arg in "${run_cmd[@]}"; do printf ' %q' "$arg"; done
    printf '\n'
    return 0
  fi
  "${run_cmd[@]}"

  local suite_root="${HYMET_ROOT}/results/${SCENARIO}/${suite}"
  if [[ ! -d "$suite_root" ]]; then
    log "WARNING: Expected results directory not found: $suite_root"
    return 0
  fi
  local latest_run
  latest_run="$(ls -1dt "${suite_root}"/run_* 2>/dev/null | head -n 1 || true)"
  if [[ -z "$latest_run" ]]; then
    log "WARNING: No run directories detected under ${suite_root}"
    return 0
  fi
  log "Latest run directory: ${latest_run}"

  if [[ $SKIP_FIGURES -eq 1 ]]; then
    log "Skipping figure regeneration for ${suite}"
    return 0
  fi

  local tables_dir="${latest_run}/tables"
  local figures_dir="${latest_run}/figures"
  local raw_dir="${latest_run}/raw"

  if [[ -d "$tables_dir" ]]; then
    log "Refreshing runtime/memory figures from ${tables_dir}"
    python3 "$MAKE_FIGURES" \
      --bench-root "$BENCH_ROOT" \
      --tables "$tables_dir" \
      --outdir "$figures_dir"
  else
    log "WARNING: Tables directory missing (${tables_dir}); skipping runtime plots."
  fi

  if [[ -d "$raw_dir" ]]; then
    log "Refreshing case-study figures from ${raw_dir}"
    python3 "$CASE_PLOTTER" \
      --case-root "$raw_dir" \
      --figures-dir "$figures_dir"
  else
    log "WARNING: Raw directory missing (${raw_dir}); skipping case figures."
  fi
}

for suite in "${REQUESTED_SUITES[@]}"; do
  run_suite "$suite"
done
