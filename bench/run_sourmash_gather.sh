#!/usr/bin/env bash
# Run sourmash gather against the shared signature index and produce a CAMI profile.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
DB_DIR="${BENCH_ROOT}/db/sourmash"
THREADS="${THREADS:-8}"
KSIZE="${SOURMASH_KSIZE:-31}"
SCALED="${SOURMASH_SCALED:-1000}"
TOP_HITS="${SOURMASH_TOP_HITS:-500}"
SEQMAP=""
TAXDIR="${TAXDIR:-${HYMET_ROOT}/taxonomy_files}"

usage(){
  cat <<'USAGE'
Usage: run_sourmash_gather.sh --sample ID --contigs FASTA [--db DIR] [--out DIR] [--seqmap FILE] [--threads N]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --db) DB_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --seqmap) SEQMAP="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --ksize) KSIZE="$2"; shift 2;;
    --scaled) SCALED="$2"; shift 2;;
    --top-hits) TOP_HITS="$2"; shift 2;;
    --taxdir) TAXDIR="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/sourmash_gather}"
ensure_dir "${OUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
DB_DIR="$(resolve_path "${DB_DIR}")"

if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "sourmash gather input FASTA missing (${CONTIGS_ABS})"
fi

if ! command -v sourmash >/dev/null 2>&1; then
  if command -v micromamba >/dev/null 2>&1; then
    micromamba install -y -p /opt/conda -c conda-forge -c bioconda -c defaults sourmash
  else
    die "sourmash executable not found and micromamba unavailable."
  fi
fi

SIGFILE="${DB_DIR}/reference.k${KSIZE}.scaled${SCALED}.sig"
INDEX_FILE="${DB_DIR}/reference_k${KSIZE}.sbt.zip"
if [[ ! -s "${SIGFILE}" ]]; then
  die "sourmash signature file missing (${SIGFILE}). Run bench/db/build_sourmash.sh first."
fi
if [[ ! -s "${INDEX_FILE}" ]]; then
  die "sourmash index file missing (${INDEX_FILE}). Run bench/db/build_sourmash.sh first."
fi

REF_FASTA="${HYMET_ROOT}/data/downloaded_genomes/combined_genomes.fasta"
if [[ ! -s "${REF_FASTA}" ]]; then
  die "Reference FASTA missing (${REF_FASTA}). Run HYMET pipeline once or build databases first."
fi

if [[ -z "${SEQMAP}" ]]; then
  for candidate in \
    "${BENCH_ROOT}/db/centrifuge/seqid2taxid.map" \
    "${BENCH_ROOT}/db/ganon2/seqid2taxid.map" \
    "${BENCH_ROOT}/db/shared/seqid2taxid.map"; do
    if [[ -s "${candidate}" ]]; then
      SEQMAP="${candidate}"
      break
    fi
  done
fi

if [[ -z "${SEQMAP}" ]]; then
  SHARED_DIR="${BENCH_ROOT}/db/shared"
  ensure_dir "${SHARED_DIR}"
  IDMAP="${SHARED_DIR}/id2taxid.tsv"
  if [[ ! -s "${IDMAP}" ]]; then
    python3 "${HYMET_ROOT}/tools/build_id_map.py" "${HYMET_ROOT}/data/detailed_taxonomy.tsv" "${IDMAP}"
  fi
  SEQMAP="${SHARED_DIR}/seqid2taxid.map"
  if [[ ! -s "${SEQMAP}" || "${IDMAP}" -nt "${SEQMAP}" || "${HYMET_ROOT}/data/downloaded_genomes/combined_genomes.fasta" -nt "${SEQMAP}" ]]; then
    python3 "${BENCH_ROOT}/lib/make_seqid_map.py" \
      --fasta "${REF_FASTA}" \
      --taxonomy-map "${IDMAP}" \
      --out "${SEQMAP}"
  fi
fi

SEQMAP="$(resolve_path "${SEQMAP}")"
if [[ ! -s "${SEQMAP}" ]]; then
  die "seqid2taxid map not found (${SEQMAP})."
fi

GATHER_CSV="${OUT_DIR}/sourmash.gather.csv"
PROFILE_CAMI="${OUT_DIR}/profile.cami.tsv"
QUERY_SIG="${OUT_DIR}/query.k${KSIZE}.scaled${SCALED}.sig"

log "Sketching query signature for ${SAMPLE}"
sourmash sketch dna "${CONTIGS_ABS}" \
  -p "k=${KSIZE},scaled=${SCALED}" \
  --name-from-first \
  --output "${QUERY_SIG}" \
  --force

log "Running sourmash gather for ${SAMPLE}"
gather_cmd=(sourmash gather "${QUERY_SIG}" "${INDEX_FILE}" --ksize "${KSIZE}" --threshold-bp 0 --output "${GATHER_CSV}")
if [[ "${TOP_HITS}" != "all" ]]; then
  if sourmash gather --help 2>&1 | grep -q -- "--num-results"; then
    gather_cmd+=(--num-results "${TOP_HITS}")
  fi
fi
"${gather_cmd[@]}"

python3 "${SCRIPT_DIR}/convert/sourmash_gather_to_cami.py" \
  --gather "${GATHER_CSV}" \
  --seqmap "${SEQMAP}" \
  --taxdb "${TAXDIR}" \
  --sample-id "${SAMPLE}" \
  --tool "sourmash_gather" \
  --out "${PROFILE_CAMI}"

SOURMASH_VER=$(sourmash --version 2>/dev/null | head -n1 | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g' || echo "unknown")
cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "sourmash_gather",
  "profile": "${PROFILE_CAMI}",
  "gather_csv": "${GATHER_CSV}",
  "seqmap": "${SEQMAP}",
  "index": "${INDEX_FILE}",
  "query_sig": "${QUERY_SIG}",
  "ksize": "${KSIZE}",
  "scaled": "${SCALED}",
  "top_hits": "${TOP_HITS}",
  "threads": "${THREADS}",
  "sourmash_version": "${SOURMASH_VER}"
}
EOF

normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
