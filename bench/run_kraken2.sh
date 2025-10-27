#!/usr/bin/env bash
# Run Kraken2 + Bracken on a CAMI sample and produce CAMI-formatted profile.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
DB_DIR="${BENCH_ROOT}/db/kraken2"
THREADS="${THREADS:-8}"
CONFIDENCE="${KRAKEN2_CONFIDENCE:-0.0}"
READ_LEN="${BRACKEN_READ_LEN:-150}"

usage(){
  cat <<'USAGE'
Usage: run_kraken2.sh --sample ID --contigs FASTA [--db DIR] [--out DIR] [--threads N]
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
    --confidence) CONFIDENCE="$2"; shift 2;;
    --read-len) READ_LEN="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/kraken2}"
ensure_dir "${OUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "Kraken2 input FASTA missing (${CONTIGS_ABS})"
fi

if [[ ! -d "${DB_DIR}" ]]; then
  die "Kraken2 database directory not found (${DB_DIR})"
fi

if ! command -v kraken2 >/dev/null 2>&1; then
  if command -v micromamba >/dev/null 2>&1; then
    micromamba install -y -p /opt/conda -c conda-forge -c bioconda -c defaults kraken2
  else
    die "kraken2 executable not found and micromamba unavailable."
  fi
fi

REPORT="${OUT_DIR}/kraken.report.tsv"
OUTPUT="${OUT_DIR}/kraken.output.tsv"
BRACKEN_TABLE="${OUT_DIR}/bracken_species.tsv"
BRACKEN_REPORT="${OUT_DIR}/bracken_species.kreport"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_TSV="${OUT_DIR}/classified_sequences.tsv"
OUTPUT_META="${OUTPUT}"

log "Running Kraken2 for ${SAMPLE}"
kraken2 \
  --db "${DB_DIR}" \
  --threads "${THREADS}" \
  --use-names \
  --confidence "${CONFIDENCE}" \
  --output "${OUTPUT}" \
  --report "${REPORT}" \
  "${CONTIGS_ABS}"

log "Running Bracken refinement"
if command -v bracken >/dev/null 2>&1; then
  if ! bracken \
    -d "${DB_DIR}" \
    -i "${REPORT}" \
    -o "${BRACKEN_TABLE}" \
    -w "${BRACKEN_REPORT}" \
    -r "${READ_LEN}" \
    -l S; then
    log "WARNING: bracken failed; falling back to raw Kraken2 report."
    BRACKEN_REPORT="${REPORT}"
  fi
else
  log "WARNING: bracken executable not available; falling back to raw Kraken2 report."
  BRACKEN_REPORT="${REPORT}"
fi

if [[ ! -s "${BRACKEN_REPORT}" ]]; then
  if [[ -s "${REPORT}" ]]; then
    log "WARNING: expected Bracken output missing; using raw Kraken2 report."
    BRACKEN_REPORT="${REPORT}"
  else
    log "WARNING: no Bracken or Kraken report available; aborting."
    exit 1
  fi
fi

python3 "${SCRIPT_DIR}/convert/kraken2_to_cami.py" \
  --report "${BRACKEN_REPORT}" \
  --out "${PROFILE_CAMI}" \
  --sample-id "${SAMPLE}" \
  --tool "kraken2+bracken"

python3 - "${OUTPUT}" "${CLASSIFIED_TSV}" <<'PY'
import csv, re, sys

def extract_taxid(value: str) -> str:
    if not value:
        return ""
    m = re.search(r"taxid\s*(\d+)", value, re.IGNORECASE)
    if m:
        return m.group(1)
    m = re.findall(r"\d+", value)
    return m[-1] if m else ""

inp, outp = sys.argv[1:3]
with open(outp, "w", newline="", encoding="utf-8") as out:
    wr = csv.writer(out, delimiter="\t")
    wr.writerow(["Query", "TaxID"])
    with open(inp, encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != "C":
                continue
            query = parts[1].strip()
            taxid = extract_taxid(parts[2])
            if not taxid:
                continue
            wr.writerow([query, taxid])
PY

if [[ "${KEEP_KRAKEN2_RAW:-0}" -eq 0 ]]; then
  rm -f "${OUTPUT}" "${BRACKEN_TABLE}"
  OUTPUT_META=""
fi

# Write per-run provenance metadata (include versions when available)
KRAKEN_VER=$(kraken2 --version 2>/dev/null | tr '\n' ' ' | sed 's/\s\+/ /g' || echo "unknown")
BRACKEN_VER=$(bracken -v 2>/dev/null | tr '\n' ' ' | sed 's/\s\+/ /g' || echo "unknown")
cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "kraken2+bracken",
  "profile": "${PROFILE_CAMI}",
  "classified": "${CLASSIFIED_TSV}",
  "report": "${BRACKEN_REPORT}",
  "db_dir": "${DB_DIR}",
  "threads": "${THREADS}",
  "output_raw": "${OUTPUT_META}",
  "kraken2_version": "${KRAKEN_VER}",
  "bracken_version": "${BRACKEN_VER}"
}
EOF

normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
