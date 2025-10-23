#!/usr/bin/env bash
# Run Centrifuge on a CAMI sample and convert to CAMI profile.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
DB_DIR="${BENCH_ROOT}/db/centrifuge"
THREADS="${THREADS:-8}"

usage(){
  cat <<'USAGE'
Usage: run_centrifuge.sh --sample ID --contigs FASTA [--db DIR] [--out DIR] [--threads N]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --db) DB_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/centrifuge}"
ensure_dir "${OUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "Centrifuge input FASTA missing (${CONTIGS_ABS})"
fi

if ! command -v centrifuge >/dev/null 2>&1; then
  if command -v micromamba >/dev/null 2>&1; then
    micromamba install -y -p /opt/conda -c conda-forge -c bioconda -c defaults centrifuge
  else
    die "centrifuge executable not found and micromamba unavailable."
  fi
fi

INDEX_PREFIX="${DB_DIR}/centrifuge"
if [[ ! -s "${INDEX_PREFIX}.1.cf" ]]; then
  die "Centrifuge index missing (${INDEX_PREFIX}.1.cf)."
fi

OUTPUT="${OUT_DIR}/centrifuge.output.tsv"
REPORT="${OUT_DIR}/centrifuge.report.tsv"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_TSV="${OUT_DIR}/classified_sequences.tsv"
OUTPUT_META="${OUTPUT}"

log "Running Centrifuge for ${SAMPLE}"
centrifuge \
  -x "${INDEX_PREFIX}" \
  -U "${CONTIGS_ABS}" \
  -f \
  -S "${OUTPUT}" \
  --report-file "${REPORT}" \
  -p "${THREADS}"

python3 "${SCRIPT_DIR}/convert/centrifuge_to_cami.py" \
  --report "${REPORT}" \
  --out "${PROFILE_CAMI}" \
  --sample-id "${SAMPLE}" \
  --tool "centrifuge"

python3 - "${OUTPUT}" "${CLASSIFIED_TSV}" <<'PY'
import csv, sys
inp, outp = sys.argv[1:3]
with open(outp, "w", newline="") as out:
    wr = csv.writer(out, delimiter="\t")
    wr.writerow(["Query", "TaxID", "Score"])
    with open(inp) as fin:
        next(fin, None)
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            query = parts[0]
            taxid = parts[2]
            score = parts[3]
            wr.writerow([query, taxid, score])
PY

if [[ "${KEEP_CENTRIFUGE_RAW:-0}" -eq 0 ]]; then
  rm -f "${OUTPUT}"
  OUTPUT_META=""
fi

CENTRIFUGE_VER=$(centrifuge --version 2>/dev/null | head -n1 | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g' || echo "unknown")
cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "centrifuge", "profile": "${PROFILE_CAMI}", "classified": "${CLASSIFIED_TSV}", "report": "${REPORT}", "db_prefix": "${INDEX_PREFIX}", "threads": "${THREADS}", "output_raw": "${OUTPUT_META}", "centrifuge_version": "${CENTRIFUGE_VER}"}
EOF
