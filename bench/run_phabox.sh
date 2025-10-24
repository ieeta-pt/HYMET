#!/usr/bin/env bash
# Wrapper to execute PhaBOX predictions within the HYMET benchmark harness.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
PHABOX_CMD="${PHABOX_CMD:-phabox2}"
PHABOX_DB_DIR="${PHABOX_DB_DIR:-}"
PHABOX_TASK="${PHABOX_TASK:-phagcn}"
PHABOX_EXTRA_OPTS="${PHABOX_EXTRA_OPTS:-}"
PHABOX_WORKDIR="${PHABOX_WORKDIR:-}"

ensure_prodigal_gv(){
  if command -v prodigal-gv >/dev/null 2>&1; then
    return 0
  fi

  if command -v micromamba >/dev/null 2>&1; then
    log "Installing prodigal-gv via micromamba (bioconda channel)"
    micromamba install -y -p /opt/conda -c conda-forge -c bioconda prodigal-gv || die "Failed to install prodigal-gv with micromamba"
  else
    die "prodigal-gv not found on PATH and micromamba unavailable. Install it manually via 'conda install -c bioconda prodigal-gv'."
  fi

  command -v prodigal-gv >/dev/null 2>&1 || die "prodigal-gv still missing after attempted installation"
}

ensure_taxonkit(){
  if command -v taxonkit >/dev/null 2>&1; then
    return 0
  fi

  if command -v micromamba >/dev/null 2>&1; then
    log "Installing taxonkit via micromamba (bioconda channel)"
    micromamba install -y -p /opt/conda -c conda-forge -c bioconda taxonkit || die "Failed to install taxonkit with micromamba"
  else
    die "taxonkit not found and micromamba unavailable. Install it manually via 'conda install -c bioconda taxonkit'."
  fi

  command -v taxonkit >/dev/null 2>&1 || die "taxonkit still missing after attempted installation"
}

usage(){
  cat <<'USAGE'
Usage: run_phabox.sh --sample ID --contigs FASTA --db DIR [--out DIR] [--threads N]

Environment variables:
  PHABOX_CMD         Command used to invoke PhaBOX (default: phabox2).
  PHABOX_DB_DIR      Directory with the PhaBOX database (required).
  PHABOX_TASK        PhaBOX task to run (default: phagcn).
  PHABOX_EXTRA_OPTS  Extra arguments appended to the command.
  PHABOX_WORKDIR     Optional working directory (e.g., PhaBOX installation root).
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --db) PHABOX_DB_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" && -n "${PHABOX_DB_DIR}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/phabox}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
PHABOX_DB_DIR="$(resolve_path "${PHABOX_DB_DIR}")"

[[ -s "${CONTIGS_ABS}" ]] || die "PhaBOX input FASTA missing (${CONTIGS_ABS})"
[[ -d "${PHABOX_DB_DIR}" ]] || die "PhaBOX database directory not found (${PHABOX_DB_DIR})"

ensure_prodigal_gv
ensure_taxonkit

CONVERTED_FASTA="${RUN_DIR}/contigs_for_phabox.fasta"
ID_MAP="${RUN_DIR}/id_mapping.tsv"

python3 - <<'PY' "${CONTIGS_ABS}" "${CONVERTED_FASTA}" "${ID_MAP}"
import sys

src_path, dst_path, map_path = sys.argv[1:4]

def write_record(handle, idx, seq):
    handle.write(f">{idx}\n")
    for i in range(0, len(seq), 80):
        handle.write(seq[i:i+80] + "\n")

count = 0

with open(src_path, "r", encoding="utf-8", errors="ignore") as src, \
     open(dst_path, "w", encoding="utf-8") as dst, \
     open(map_path, "w", encoding="utf-8") as mapping:
    mapping.write("New_ID\tOriginal_ID\n")
    current_id = None
    seq_lines = []
    for line in src:
        if line.startswith(">"):
            if current_id is not None:
                count += 1
                write_record(dst, count, "".join(seq_lines))
                mapping.write(f"{count}\t{current_id}\n")
            current_id = line[1:].strip().split()[0]
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if current_id is not None:
        count += 1
        write_record(dst, count, "".join(seq_lines))
        mapping.write(f"{count}\t{current_id}\n")
PY

PHABOX_OUT="${RUN_DIR}/phabox_out"
ensure_dir "${PHABOX_OUT}"

IFS=' ' read -r -a PHABOX_CMD_ARR <<<"${PHABOX_CMD}"
[[ ${#PHABOX_CMD_ARR[@]} -gt 0 ]] || die "PHABOX_CMD is empty"

run_command(){
  if [[ -n "${PHABOX_WORKDIR}" ]]; then
    (cd "${PHABOX_WORKDIR}" && "$@")
  else
    "$@"
  fi
}

log "Running PhaBOX for ${SAMPLE}"
if ! command -v "${PHABOX_CMD_ARR[0]}" >/dev/null 2>&1; then
  die "PhaBOX command '${PHABOX_CMD_ARR[0]}' not found on PATH. Adjust PHABOX_CMD or add it to PATH."
fi
PHABOX_ARGS=(
  "${PHABOX_CMD_ARR[@]}"
  --task "${PHABOX_TASK}"
  --contigs "${CONVERTED_FASTA}"
  --outpth "${PHABOX_OUT}"
  --threads "${THREADS}"
  --dbdir "${PHABOX_DB_DIR}"
)

if [[ -n "${PHABOX_EXTRA_OPTS}" ]]; then
  # shellcheck disable=SC2206
  EXTRA=( ${PHABOX_EXTRA_OPTS} )
  PHABOX_ARGS+=("${EXTRA[@]}")
fi

if ! run_command "${PHABOX_ARGS[@]}"; then
  die "PhaBOX execution failed for sample ${SAMPLE}"
fi

PREDICTION_FILE="$(find "${PHABOX_OUT}" -type f -name '*phagcn_prediction*.tsv' -o -name 'phagcn_prediction.tsv' | head -n 1 || true)"
[[ -n "${PREDICTION_FILE}" && -s "${PREDICTION_FILE}" ]] || die "Unable to locate phagcn prediction file under ${PHABOX_OUT}"

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"

if [[ -z "${TAXONKIT_DB:-}" && -d "${HYMET_ROOT}/taxonomy_files" ]]; then
  export TAXONKIT_DB="${HYMET_ROOT}/taxonomy_files"
fi

CONVERTER="${SCRIPT_DIR}/convert/phabox_to_cami.py"
if [[ ! -f "${CONVERTER}" ]]; then
  ALT="$(python3 - <<'PY' "${HYMET_ROOT}"
import os, sys
root = sys.argv[1]
candidate = os.path.normpath(os.path.join(root, "..", "convert", "phabox_to_cami.py"))
print(candidate if os.path.isfile(candidate) else "")
PY
)"
  if [[ -n "${ALT}" ]]; then
    CONVERTER="${ALT}"
  fi
fi
[[ -f "${CONVERTER}" ]] || die "Unable to locate phabox_to_cami.py converter script"

python3 "${CONVERTER}" \
  --input "${PREDICTION_FILE}" \
  --id-map "${ID_MAP}" \
  --out "${PROFILE_DST}" \
  --classified-out "${CLASSIFIED_DST}" \
  --sample-id "${SAMPLE}" \
  --tool phabox \
  --taxdb "${TAXONKIT_DB:-}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "phabox",
  "profile": "${PROFILE_DST}",
  "contigs": "${CLASSIFIED_DST}",
  "prediction_file": "${PREDICTION_FILE}",
  "db_dir": "${PHABOX_DB_DIR}",
  "threads": "${THREADS}"
}
EOF

normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
