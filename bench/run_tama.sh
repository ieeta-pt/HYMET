#!/usr/bin/env bash
# Wrapper to execute TAMA consensus taxonomy assignments within the HYMET benchmark.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"

TAMA_ROOT="${TAMA_ROOT:-}"
TAMA_PARAM_FILE="${TAMA_PARAM_FILE:-}"
TAMA_DBNAME="${TAMA_DBNAME:-tama}"
TAMA_RANK="${TAMA_RANK:-species}"
TAMA_TOOLS="${TAMA_TOOLS:-centrifuge,kraken}"
TAMA_META_THRESHOLD="${TAMA_META_THRESHOLD:-0.34}"
TAMA_WEIGHT_CLARK="${TAMA_WEIGHT_CLARK:-0.9374}"
TAMA_WEIGHT_CENTRIFUGE="${TAMA_WEIGHT_CENTRIFUGE:-0.9600}"
TAMA_WEIGHT_KRAKEN="${TAMA_WEIGHT_KRAKEN:-0.9362}"
TAMA_SINGLE_READS="${TAMA_SINGLE_READS:-}"
TAMA_READS1="${TAMA_READS1:-}"
TAMA_READS2="${TAMA_READS2:-}"
TAMA_ENV_FILE="${TAMA_ENV_FILE:-}"
TAMA_CMD="${TAMA_CMD:-perl TAMA.pl}"
KEEP_WORK="${KEEP_TAMA_WORK:-0}"
TAMA_PARAM_DIR="$(resolve_path "${TAMA_PARAM_DIR:-${SCRIPT_DIR}/config/tama_params}")"

usage(){
  cat <<'USAGE'
Usage: run_tama.sh --sample ID [--contigs FASTA] [--out DIR] [--threads N] [--tama-root DIR] [--param FILE]

Environment variables / options:
  TAMA_ROOT               Path to the TAMA repository (required).
  TAMA_PARAM_FILE         Preconfigured parameter file; if unset, one is generated.
  TAMA_DBNAME             Database name registered within TAMA (default: tama).
  TAMA_TOOLS              Comma-separated list of component tools (default: centrifuge,kraken).
  TAMA_RANK               Target taxonomic rank (default: species).
  TAMA_META_THRESHOLD     Meta-analysis threshold (default: 0.34).
  TAMA_WEIGHT_CLARK       Weight assigned to CLARK (default: 0.9374; ignored when CLARK is disabled).
  TAMA_WEIGHT_CENTRIFUGE  Weight assigned to Centrifuge (default: 0.9600).
  TAMA_WEIGHT_KRAKEN      Weight assigned to Kraken (default: 0.9362).
  TAMA_SINGLE_READS       Single-end read file(s) to analyse (comma separated).
  TAMA_READS1             Paired-end forward reads (comma separated).
  TAMA_READS2             Paired-end reverse reads (comma separated).
  TAMA_ENV_FILE           Optional env.sh file to source before invoking TAMA.
  TAMA_CMD                Command used to launch TAMA (default: "perl TAMA.pl").
  KEEP_TAMA_WORK          Set to 1 to retain intermediate work directories.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --tama-root) TAMA_ROOT="$2"; shift 2;;
    --param) TAMA_PARAM_FILE="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/tama}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"
WORK_DIR="${RUN_DIR}/work"
ensure_dir "${WORK_DIR}"
PARAM_DST="${OUT_DIR}/params.txt"

if [[ -n "${CONTIGS}" ]]; then
  CONTIGS="$(resolve_path "${CONTIGS}")"
  [[ -e "${CONTIGS}" ]] || die "Provided contigs file not found (${CONTIGS})"
fi

TAMA_ROOT="$(resolve_path "${TAMA_ROOT}")"
[[ -d "${TAMA_ROOT}" ]] || die "TAMA_ROOT must point to the TAMA repository (got: ${TAMA_ROOT})"

if [[ -z "${TAMA_PARAM_FILE}" && -n "${TAMA_PARAM_DIR:-}" ]]; then
  TAMA_PARAM_FILE="${TAMA_PARAM_DIR%/}/${SAMPLE}.params.txt"
fi

if [[ -z "${TAMA_PARAM_FILE}" ]]; then
  if [[ -z "${TAMA_SINGLE_READS}" && ( -z "${TAMA_READS1}" || -z "${TAMA_READS2}" ) && -z "${CONTIGS}" ]]; then
    die "Provide a TAMA parameter file (--param/TAMA_PARAM_FILE) or input reads via TAMA_SINGLE_READS / TAMA_READS1+TAMA_READS2 (or supply --contigs to fall back on)."
  fi
fi

PARAM_PATH="${WORK_DIR}/params.txt"

if [[ -n "${TAMA_PARAM_FILE}" ]]; then
  TAMA_PARAM_FILE="$(resolve_path "${TAMA_PARAM_FILE}")"
  [[ -f "${TAMA_PARAM_FILE}" ]] || die "TAMA parameter file not found (${TAMA_PARAM_FILE})"
  cp "${TAMA_PARAM_FILE}" "${PARAM_PATH}"
else
  ensure_dir "$(dirname "${PARAM_PATH}")"
  {
    printf "[Project]\n"
    printf "\$PROJECTNAME=%s\n\n" "${SAMPLE}"
    printf "[Basic_options]\n"
    printf "\$TOOL=%s\n" "${TAMA_TOOLS}"
    printf "\$RANK=%s\n" "${TAMA_RANK}"
    printf "\$META-THRESHOLD=%s\n" "${TAMA_META_THRESHOLD}"
    printf "\$WEIGHT-CLARK=%s\n" "${TAMA_WEIGHT_CLARK}"
    printf "\$WEIGHT-centrifuge=%s\n" "${TAMA_WEIGHT_CENTRIFUGE}"
    printf "\$WEIGHT-kraken=%s\n\n" "${TAMA_WEIGHT_KRAKEN}"
    printf "[Database]\n"
    printf "\$DBNAME=%s\n\n" "${TAMA_DBNAME}"
    printf "[Input]\n"
    printf ">%s\n" "${SAMPLE}"
    if [[ -n "${TAMA_READS1}" && -n "${TAMA_READS2}" ]]; then
      printf "\$PAIRED1=%s\n" "${TAMA_READS1}"
      printf "\$PAIRED2=%s\n" "${TAMA_READS2}"
    fi
    if [[ -n "${TAMA_SINGLE_READS}" ]]; then
      printf "\$SINGLE=%s\n" "${TAMA_SINGLE_READS}"
    elif [[ -z "${TAMA_READS1}" || -z "${TAMA_READS2}" ]]; then
      # Fall back to treating contigs as single-end reads if supplied.
      if [[ -n "${CONTIGS}" ]]; then
        printf "\$SINGLE=%s\n" "${CONTIGS}"
      else
        die "Unable to determine input reads for TAMA"
      fi
    fi
    printf "\n[Preprocessing]\n"
    printf "\$TRIMMOMATIC-RUN=false\n"
    printf "\$BAYESHAMMER-RUN=false\n"
  } > "${PARAM_PATH}"
fi
cp "${PARAM_PATH}" "${PARAM_DST}"

if [[ -n "${TAMA_ENV_FILE}" ]]; then
  if [[ "${TAMA_ENV_FILE}" == /* ]]; then
    TAMA_ENV_FILE="$(resolve_path "${TAMA_ENV_FILE}")"
  else
    TAMA_ENV_FILE="${TAMA_ROOT}/${TAMA_ENV_FILE}"
  fi
elif [[ -f "${TAMA_ROOT}/src/env.sh" ]]; then
  TAMA_ENV_FILE="${TAMA_ROOT}/src/env.sh"
fi
if [[ -n "${TAMA_ENV_FILE}" && ! -f "${TAMA_ENV_FILE}" ]]; then
  die "TAMA_ENV_FILE points to a missing file (${TAMA_ENV_FILE})"
fi

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"
METADATA_PATH="${OUT_DIR}/metadata.json"

log "Running TAMA for ${SAMPLE}"
pushd "${TAMA_ROOT}" >/dev/null
if [[ -n "${TAMA_ENV_FILE}" ]]; then
  # shellcheck disable=SC1090
  source "${TAMA_ENV_FILE}"
fi

export PATH="/root/.local/share/mamba/bin:${PATH}"
export PERL5LIB="/root/.local/share/mamba/lib/perl5/5.32/vendor_perl${PERL5LIB:+:${PERL5LIB}}"

read -r -a tama_cmd <<<"${TAMA_CMD}"
tama_cmd+=(-p "${THREADS}" -o "${WORK_DIR}" --param "${PARAM_PATH}")

set +e
"${tama_cmd[@]}"
status=$?
set -e
popd >/dev/null
[[ ${status} -eq 0 ]] || die "TAMA execution failed with exit code ${status}"

mapfile -t ABUNDANCE_FILES < <(find "${WORK_DIR}" -maxdepth 6 -type f -name "abundance_profile*.out" | sort)
[[ ${#ABUNDANCE_FILES[@]} -gt 0 ]] || die "Unable to locate TAMA abundance_profile*.out under ${WORK_DIR}"
ABUNDANCE_SRC="${ABUNDANCE_FILES[0]}"

mapfile -t READ_CLASSI_FILES < <(find "${WORK_DIR}" -maxdepth 6 -type f -name "read_classi*.out" | sort)
READ_CLASSI_SRC=""
if [[ ${#READ_CLASSI_FILES[@]} -gt 0 ]]; then
  READ_CLASSI_SRC="${READ_CLASSI_FILES[0]}"
fi

rm -f "${CLASSIFIED_DST}"

python3 "${SCRIPT_DIR}/convert/tama_to_cami.py" \
  --profile "${ABUNDANCE_SRC}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool tama \
  --rank "${TAMA_RANK}" \
  --taxdb "${TAXONKIT_DB:-}" \
  ${READ_CLASSI_SRC:+--read-classi "${READ_CLASSI_SRC}"} \
  --classified-out "${CLASSIFIED_DST}"

CLASSIFIED_REF="${CLASSIFIED_DST}"
if [[ ! -f "${CLASSIFIED_DST}" ]]; then
  CLASSIFIED_REF=""
fi

cat > "${METADATA_PATH}" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "tama",
  "profile": "${PROFILE_DST}",
  "classified_sequences": "${CLASSIFIED_REF}",
  "param_file": "${PARAM_DST}",
  "tama_root": "${TAMA_ROOT}",
  "threads": "${THREADS}"
}
EOF

normalize_metadata_json "${METADATA_PATH}" "${OUT_DIR}"

if [[ "${KEEP_WORK}" -ne 1 ]]; then
  rm -rf "${WORK_DIR}" || true
fi
