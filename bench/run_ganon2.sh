#!/usr/bin/env bash
# Run ganon2 classification and convert to CAMI profile.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
DB_DIR="${BENCH_ROOT}/db/ganon2"
THREADS="${THREADS:-8}"
REL_CUTOFF="${GANON_REL_CUTOFF:-0}"
REL_FILTER="${GANON_REL_FILTER:-1}"

usage(){
  cat <<'USAGE'
Usage: run_ganon2.sh --sample ID --contigs FASTA [--db DIR] [--out DIR] [--threads N]
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

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/ganon2}"
ensure_dir "${OUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "ganon2 input FASTA missing (${CONTIGS_ABS})"
fi

if ! command -v ganon >/dev/null 2>&1; then
  if command -v micromamba >/dev/null 2>&1; then
    micromamba install -y -p /opt/conda -c conda-forge -c bioconda -c defaults ganon
  else
    die "ganon executable not found and micromamba unavailable."
  fi
fi

INDEX_PREFIX="${DB_DIR}/ganon"
if [[ ! -s "${INDEX_PREFIX}.ibf" && ! -s "${INDEX_PREFIX}.hibf" ]]; then
  die "ganon database index missing (${INDEX_PREFIX}.ibf/.hibf)."
fi

ASSIGN_PREFIX="${OUT_DIR}/ganon"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_TSV="${OUT_DIR}/classified_sequences.tsv"
REPORT_FILE=""
ASSIGN_FILE="${ASSIGN_PREFIX}.one"

log "Running ganon classify for ${SAMPLE}"
ganon classify \
  --db-prefix "${INDEX_PREFIX}" \
  --single-reads "${CONTIGS_ABS}" \
  --threads "${THREADS}" \
  --rel-cutoff "${REL_CUTOFF}" \
  --rel-filter "${REL_FILTER}" \
  --output-prefix "${ASSIGN_PREFIX}" \
  --output-one \
  --skip-report \
  --multiple-matches lca

[[ -s "${ASSIGN_FILE}" ]] || die "ganon classify did not produce ${ASSIGN_FILE}"

log "Generating ganon report"
ganon report \
  --input "${ASSIGN_PREFIX}.rep" \
  --db-prefix "${INDEX_PREFIX}" \
  --output-prefix "${OUT_DIR}/ganon" \
  --output-format bioboxes \
  --report-type reads \
  --min-count 0 \
  --max-count 0 \
  --split-hierarchy

GANON_CAMI_FILE=$(ls "${OUT_DIR}"/ganon*.tre 2>/dev/null | head -n 1 || true)
[[ -n "${GANON_CAMI_FILE}" && -s "${GANON_CAMI_FILE}" ]] || die "ganon report failed to produce CAMI profile"
cp -f "${GANON_CAMI_FILE}" "${PROFILE_CAMI}"
REPORT_FILE="${GANON_CAMI_FILE}"

python3 - "${ASSIGN_FILE}" "${CLASSIFIED_TSV}" <<'PY'
import csv, re, sys

inp, outp = sys.argv[1:3]

def extract_taxid(token: str) -> str:
    if not token:
        return "0"
    token = token.strip()
    parts = token.split("|")
    if parts and parts[0].isdigit():
        return parts[0]
    matches = re.findall(r"\d+", token)
    return matches[0] if matches else "0"

with open(outp, "w", newline="", encoding="utf-8") as out:
    wr = csv.writer(out, delimiter="\t")
    wr.writerow(["Query", "TaxID"])
    with open(inp, encoding="utf-8", errors="ignore") as fin:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            query = parts[0]
            taxid = extract_taxid(parts[1])
            wr.writerow([query, taxid])
PY

ASSIGN_META="${ASSIGN_FILE}"
if [[ "${KEEP_GANON_RAW:-0}" -eq 0 ]]; then
  rm -f "${ASSIGN_FILE}" "${ASSIGN_PREFIX}.rep"
  ASSIGN_META=""
fi

GANON_VER=$(ganon --version 2>/dev/null | head -n1 | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g' || echo "unknown")
cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "ganon2", "profile": "${PROFILE_CAMI}", "classified": "${CLASSIFIED_TSV}", "assignments": "${ASSIGN_META}", "report": "${REPORT_FILE}", "db_prefix": "${INDEX_PREFIX}", "threads": "${THREADS}", "rel_cutoff": "${REL_CUTOFF}", "rel_filter": "${REL_FILTER}", "ganon_version": "${GANON_VER}"}
EOF
