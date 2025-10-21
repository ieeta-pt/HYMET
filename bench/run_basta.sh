#!/usr/bin/env bash
# Wrapper to execute BASTA taxonomy assignments within the HYMET benchmark.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
BLAST_DB="${BASTA_BLAST_DB:-}"
BLASTX_CMD="${BLASTX_CMD:-blastx}"
BLASTX_EXTRA_OPTS="${BLASTX_EXTRA_OPTS:-}"
BASTA_CMD="${BASTA_CMD:-basta}"
BASTA_TAXON_MODE="${BASTA_TAXON_MODE:-uni}"
BASTA_EXTRA_OPTS="${BASTA_EXTRA_OPTS:-}"

usage(){
  cat <<'USAGE'
Usage: run_basta.sh --sample ID --contigs FASTA --blast-db DBPATH [--out DIR] [--threads N]

Environment variables:
  BLASTX_CMD         blastx executable (default: blastx)
  BLASTX_EXTRA_OPTS  Additional options appended to blastx command.
  BASTA_CMD          BASTA command (default: basta)
  BASTA_TAXON_MODE   BASTA taxonomic mode (default: uni)
  BASTA_EXTRA_OPTS   Extra arguments passed to BASTA.
  BASTA_BLAST_DB     Path to BLAST database (can also be set via --blast-db).
  TAXONKIT_DB        Optional path for TaxonKit database.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --blast-db) BLAST_DB="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/basta}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "BASTA input FASTA missing (${CONTIGS_ABS})"

BLAST_DB="$(resolve_path "${BLAST_DB}")"
[[ -n "${BLAST_DB}" ]] || die "Protein database path must be provided via --blast-db or BASTA_BLAST_DB."

DIAMOND_AVAILABLE=0
if [[ "${BASTA_USE_DIAMOND:-0}" -eq 1 || -n "${BASTA_DIAMOND_DB:-}" ]]; then
  if command -v diamond >/dev/null 2>&1; then
    DIAMOND_AVAILABLE=1
  else
    log "WARNING: diamond not found on PATH; falling back to ${BLASTX_CMD}"
  fi
elif command -v diamond >/dev/null 2>&1; then
  DIAMOND_AVAILABLE=1
fi

if [[ ${DIAMOND_AVAILABLE} -eq 1 ]]; then
  DIAMOND_DB="${BASTA_DIAMOND_DB:-${BLAST_DB}}"
  [[ -f "${DIAMOND_DB}" || -f "${DIAMOND_DB}.dmnd" ]] || die "DIAMOND database not found at ${DIAMOND_DB}"
else
  DIAMOND_DB=""
fi

if [[ ${DIAMOND_AVAILABLE} -eq 0 ]]; then
  if [[ ! -e "${BLAST_DB}" && ! -e "${BLAST_DB}.pin" ]]; then
    die "BLAST database not found at ${BLAST_DB}"
  fi
  need_cmd "${BLASTX_CMD%% *}"
fi

need_cmd "${BASTA_CMD%% *}"

BLAST_OUT="${RUN_DIR}/blast.tsv"
BASTA_PREFIX="${RUN_DIR}/basta"

if [[ ${DIAMOND_AVAILABLE} -eq 1 ]]; then
  log "Running DIAMOND blastx for ${SAMPLE}"
  if [[ -n "${DIAMOND_EXTRA_OPTS:-}" ]]; then
    read -r -a diamond_extra <<<"${DIAMOND_EXTRA_OPTS}"
  else
    diamond_extra=()
  fi
  diamond blastx \
    --query "${CONTIGS_ABS}" \
    --db "${DIAMOND_DB}" \
    --out "${BLAST_OUT}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --max-target-seqs "${BLAST_MAX_TARGETS:-50}" \
    --threads "${THREADS}" \
    "${diamond_extra[@]}"
else
  log "Running BLASTX for ${SAMPLE}"
  if [[ -n "${BLASTX_EXTRA_OPTS:-}" ]]; then
    read -r -a blast_extra <<<"${BLASTX_EXTRA_OPTS}"
  else
    blast_extra=()
  fi
  "${BLASTX_CMD}" \
    -query "${CONTIGS_ABS}" \
    -db "${BLAST_DB}" \
    -out "${BLAST_OUT}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -max_target_seqs "${BLAST_MAX_TARGETS:-50}" \
    -num_threads "${THREADS}" \
    "${blast_extra[@]}"
fi

[[ -s "${BLAST_OUT}" ]] || die "BLASTX did not produce any results for ${SAMPLE}"

log "Running BASTA for ${SAMPLE}"
if [[ -n "${BASTA_EXTRA_OPTS:-}" ]]; then
  read -r -a basta_extra <<<"${BASTA_EXTRA_OPTS}"
else
  basta_extra=()
fi
"${BASTA_CMD}" sequence "${BLAST_OUT}" "${BASTA_PREFIX}" "${BASTA_TAXON_MODE}" "${basta_extra[@]}"

declare -a taxonomy_candidates=(
  "${BASTA_PREFIX}.taxonomy"
  "${BASTA_PREFIX}.taxonomy.txt"
  "${BASTA_PREFIX}.taxonomy.tsv"
  "${BASTA_PREFIX}"
)
BASTA_TAXONOMY=""
for candidate in "${taxonomy_candidates[@]}"; do
  if [[ -s "${candidate}" ]]; then
    BASTA_TAXONOMY="${candidate}"
    break
  fi
done
[[ -n "${BASTA_TAXONOMY}" ]] || die "BASTA taxonomy output not found. Checked: ${taxonomy_candidates[*]}"

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"

python3 "${SCRIPT_DIR}/convert/basta_to_cami.py" \
  --input "${BASTA_TAXONOMY}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool basta \
  --taxdb "${TAXONKIT_DB:-}" \
  --classified-out "${CLASSIFIED_DST}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "basta", "profile": "${PROFILE_DST}", "contigs": "${CLASSIFIED_DST}", "blast_output": "${BLAST_OUT}", "blast_db": "${BLAST_DB}", "basta_taxonomy": "${BASTA_TAXONOMY}"}
EOF
