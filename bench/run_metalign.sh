#!/usr/bin/env bash
# Run Metalign on synthetic reads derived from CAMI contigs and produce CAMI-formatted profile.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
EXTRA_OPTS="${METALIGN_OPTS:-}"
DB_DIR_DEFAULT="${BENCH_ROOT}/db/metalign/data"
DB_DIR="${METALIGN_DB_DIR:-${DB_DIR_DEFAULT}}"
ENV_PREFIX="${METALIGN_ENV_PREFIX:-/opt/envs/metalign}"
PRESET="${METALIGN_PRESET:-}"  # e.g., --precise or --sensitive

usage(){
  cat <<'USAGE'
Usage: run_metalign.sh --sample ID --contigs FASTA [--out DIR] [--threads N]

Environment toggles:
  METALIGN_DB_DIR     Path to Metalign data directory (default: bench/db/metalign/data)
  METALIGN_OPTS       Extra options passed through to Metalign (e.g., "--length_normalize")
  METALIGN_PRESET     One of "--precise" or "--sensitive" (empty for default)
  METALIGN_ENV_PREFIX Micromamba/Conda prefix for installing/running metalign (default: /opt/conda)
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/metalign}"
RUN_DIR="${OUT_DIR}/run"
ensure_dir "${OUT_DIR}"; ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "Metalign input FASTA missing (${CONTIGS_ABS})"
fi

if [[ ! -d "${DB_DIR}" ]]; then
  die "Metalign database directory not found (${DB_DIR}). Run bench/db/build_metalign.sh or set METALIGN_DB_DIR."
fi

USE_MM=0
if command -v micromamba >/dev/null 2>&1 && [[ -d "${ENV_PREFIX}" ]]; then
  USE_MM=1
fi

# Synthesize single-end reads from contigs for the read-based Metalign pipeline.
SYN_FASTQ="${RUN_DIR}/synthetic_reads.fastq"
python3 "${SCRIPT_DIR}/tools/contigs_to_reads.py" \
  --contigs "${CONTIGS_ABS}" \
  --out "${SYN_FASTQ}" \
  --chunk-size "${METALIGN_CHUNK_SIZE:-250}" \
  --min-chunk "${METALIGN_MIN_CHUNK:-100}"

RAW_OUT="${OUT_DIR}/metalign_abundances.tsv"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"

IFS=' ' read -r -a EXTRA_ARGS <<< "${EXTRA_OPTS:-}"
if [[ ${#EXTRA_ARGS[@]} -eq 0 || -z "${EXTRA_ARGS[0]:-}" ]]; then
  EXTRA_ARGS=()
fi
if [[ -n "${PRESET}" ]]; then
  EXTRA_ARGS=("${EXTRA_ARGS[@]}" "${PRESET}")
fi

if [[ ${USE_MM} -eq 1 ]]; then
  run_cmd=("micromamba" "run" "-p" "${ENV_PREFIX}" "metalign" "${SYN_FASTQ}" "${DB_DIR}" "--threads" "${THREADS}" "--output" "${RAW_OUT}")
else
  run_cmd=("metalign" "${SYN_FASTQ}" "${DB_DIR}" "--threads" "${THREADS}" "--output" "${RAW_OUT}")
fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  run_cmd+=("${EXTRA_ARGS[@]}")
fi

log "Running Metalign for ${SAMPLE} (threads=${THREADS})"
set +e
"${run_cmd[@]}"
status=$?
set -e
if [[ ${status} -ne 0 ]]; then
  # Fallback to Python module invocation if console_script not present
  log "Metalign CLI failed (status=${status}); trying Python fallback"
  if command -v micromamba >/dev/null 2>&1 && [[ -d "${ENV_PREFIX}" ]]; then
    micromamba run -p "${ENV_PREFIX}" python - <<'PY' "${SYN_FASTQ}" "${DB_DIR}" "${THREADS}" "${RAW_OUT}" "${EXTRA_OPTS}" "${PRESET}"
import os, shlex, subprocess, sys
fq, db, threads, outp, extra_opts, preset = sys.argv[1:]
args = [sys.executable, '-m', 'metalign.metalign', fq, db, '--threads', threads, '--output', outp]
if preset:
    args.append(preset)
if extra_opts:
    args.extend(shlex.split(extra_opts))
subprocess.check_call(args)
PY
  else
    python3 - <<'PY' "${SYN_FASTQ}" "${DB_DIR}" "${THREADS}" "${RAW_OUT}" "${EXTRA_OPTS}" "${PRESET}"
import os, shlex, subprocess, sys
fq, db, threads, outp, extra_opts, preset = sys.argv[1:]
script = 'metalign.py'
if not os.path.exists(script):
    # Try to import module as a last resort
    args = [sys.executable, '-m', 'metalign.metalign', fq, db, '--threads', threads, '--output', outp]
else:
    args = [sys.executable, script, fq, db, '--threads', threads, '--output', outp]
if preset:
    args.append(preset)
if extra_opts:
    args.extend(shlex.split(extra_opts))
subprocess.check_call(args)
PY
  fi
fi

if [[ ! -s "${RAW_OUT}" ]]; then
  die "Metalign did not produce an abundance file (${RAW_OUT})"
fi

python3 "${SCRIPT_DIR}/convert/metalign_to_cami.py" \
  --input "${RAW_OUT}" \
  --out "${PROFILE_CAMI}" \
  --sample-id "${SAMPLE}" \
  --tool "metalign"

# Metalign is read-based; no per-contig classifications available. Emit header-only TSV.
CLASSIFIED_TSV="${OUT_DIR}/classified_sequences.tsv"
printf "Query\tTaxID\n" >"${CLASSIFIED_TSV}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "metalign", "profile": "${PROFILE_CAMI}", "abundances_raw": "${RAW_OUT}", "db_dir": "${DB_DIR}", "synthetic_fastq": "${SYN_FASTQ}", "threads": "${THREADS}", "opts": "${EXTRA_OPTS}", "preset": "${PRESET}"}
EOF
