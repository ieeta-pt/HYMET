#!/usr/bin/env bash
# Wrapper to execute phyloFlash within the HYMET benchmark harness.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
PHYLOFLASH_DB_DIR="${PHYLOFLASH_DB_DIR:-}"
PHYLOFLASH_ENV_PREFIX="${PHYLOFLASH_ENV_PREFIX:-/opt/envs/phyloflash}"
PHYLOFLASH_CMD="${PHYLOFLASH_CMD:-${PHYLOFLASH_ENV_PREFIX}/bin/phyloFlash.pl}"
PHYLOFLASH_EXTRA_OPTS="${PHYLOFLASH_EXTRA_OPTS:-}"
MIN_FRAGMENT_LEN="${PHYLOFLASH_MIN_FRAGMENT_LEN:-60}"
FRAGMENT_LEN="${PHYLOFLASH_FRAGMENT_LEN:-250}"

usage(){
  cat <<'USAGE'
Usage: run_phyloflash.sh --sample ID --contigs FASTA --db DIR [--out DIR] [--threads N]

Environment variables:
  PHYLOFLASH_DB_DIR        Path to the phyloFlash database directory (required).
  PHYLOFLASH_ENV_PREFIX    Conda/mamba environment prefix containing phyloFlash (default: /opt/envs/phyloflash).
  PHYLOFLASH_CMD           phyloFlash executable (default: phyloFlash.pl within the environment).
  PHYLOFLASH_EXTRA_OPTS    Extra options passed to phyloFlash.pl.
  PHYLOFLASH_FRAGMENT_LEN  Fragment length (bp) for synthetic reads (default: 250).
  PHYLOFLASH_MIN_FRAGMENT_LEN  Minimum fragment length to retain (default: 60).
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --db) PHYLOFLASH_DB_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/phyloflash}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "Input contigs FASTA missing (${CONTIGS_ABS})"

PHYLOFLASH_DB_DIR="$(resolve_path "${PHYLOFLASH_DB_DIR}")"
[[ -d "${PHYLOFLASH_DB_DIR}" ]] || die "phyloFlash database directory not found (${PHYLOFLASH_DB_DIR}). Set --db or PHYLOFLASH_DB_DIR."

[[ -d "${PHYLOFLASH_ENV_PREFIX}" ]] || die "phyloFlash environment prefix not found (${PHYLOFLASH_ENV_PREFIX})"
need_cmd micromamba

run_in_env(){
  micromamba run -p "${PHYLOFLASH_ENV_PREFIX}" "$@"
}

LIB_NAME="${SAMPLE//[^A-Za-z0-9_]/_}"
[[ -n "${LIB_NAME}" ]] || LIB_NAME="sample"

GFF_PATH="${RUN_DIR}/rrna.gff"
RRNA_FASTA="${RUN_DIR}/rrna_sequences.fna"
READS_FASTQ="${RUN_DIR}/rrna_reads.fastq"

log "Extracting rRNA annotations with barrnap"
if ! run_in_env barrnap --threads "${THREADS}" "${CONTIGS_ABS}" > "${GFF_PATH}"; then
  die "barrnap failed on ${CONTIGS_ABS}"
fi

if [[ ! -s "${GFF_PATH}" ]]; then
  log "No rRNA features detected; emitting empty phyloFlash profile"
  PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
  CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"
  EMPTY_CSV="${RUN_DIR}/empty_phyloflash.csv"
  : > "${EMPTY_CSV}"
  python3 "${SCRIPT_DIR}/convert/phyloflash_to_cami.py" \
    --input "${EMPTY_CSV}" \
    --out "${PROFILE_DST}" \
    --sample-id "${SAMPLE}" \
    --tool phyloflash \
    --taxdb "${TAXONKIT_DB:-}" \
    --classified-out "${CLASSIFIED_DST}" || true
  cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "phyloflash",
  "profile": "${PROFILE_DST}",
  "contigs": "${CONTIGS_ABS}",
  "db_dir": "${PHYLOFLASH_DB_DIR}",
  "threads": "${THREADS}",
  "note": "No rRNA features detected"
}
EOF
  normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
  exit 0
fi

log "Extracting rRNA sequences with bedtools"
run_in_env bedtools getfasta -fi "${CONTIGS_ABS}" -bed "${GFF_PATH}" -fo "${RRNA_FASTA}"
[[ -s "${RRNA_FASTA}" ]] || die "Failed to extract rRNA sequences for ${SAMPLE}"

log "Creating synthetic reads (${FRAGMENT_LEN} bp) from rRNA sequences"
python3 - "${RRNA_FASTA}" "${READS_FASTQ}" "${FRAGMENT_LEN}" "${MIN_FRAGMENT_LEN}" <<'PY'
import sys
from pathlib import Path

src, dst, frag_len_str, min_len_str = sys.argv[1:5]
fragment_len = int(frag_len_str)
min_len = int(min_len_str)
source = Path(src)
target = Path(dst)

if not source.exists() or source.stat().st_size == 0:
    sys.exit("No source sequences to fragment")

def flush_record(identifier, seq, handle):
    if not seq:
        return
    seq = seq.upper()
    idx = 0
    frag_idx = 0
    length = len(seq)
    while idx < length:
        frag = seq[idx : idx + fragment_len]
        if len(frag) < min_len:
            break
        frag_idx += 1
        handle.write(f"@{identifier}_{frag_idx}\n{frag}\n+\n{'I' * len(frag)}\n")
        idx += fragment_len

with source.open("r") as fin, target.open("w") as fout:
    current_id = None
    chunks = []
    for line in fin:
        if line.startswith(">"):
            if current_id is not None:
                flush_record(current_id, "".join(chunks), fout)
            current_id = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if current_id is not None:
        flush_record(current_id, "".join(chunks), fout)
PY

[[ -s "${READS_FASTQ}" ]] || die "Synthetic read file is empty for ${SAMPLE}"

PHYLOFLASH_ARGS=(
  "${PHYLOFLASH_CMD}"
  -lib "${LIB_NAME}"
  -read1 "$(basename "${READS_FASTQ}")"
  -CPUs "${THREADS}"
  -dbhome "${PHYLOFLASH_DB_DIR}"
)

if [[ -n "${PHYLOFLASH_EXTRA_OPTS}" ]]; then
  # shellcheck disable=SC2206
  extra=( ${PHYLOFLASH_EXTRA_OPTS} )
  PHYLOFLASH_ARGS+=("${extra[@]}")
fi

log "Running phyloFlash for ${SAMPLE}"
pushd "${RUN_DIR}" >/dev/null
if ! run_in_env "${PHYLOFLASH_ARGS[@]}"; then
  popd >/dev/null
  die "phyloFlash execution failed for ${SAMPLE}"
fi
popd >/dev/null

NTU_FILE="${RUN_DIR}/${LIB_NAME}.phyloFlash.NTUfull_abundance.csv"
[[ -s "${NTU_FILE}" ]] || die "phyloFlash output missing (${NTU_FILE})"

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"

if [[ -z "${TAXONKIT_DB:-}" && -d "${HYMET_ROOT}/taxonomy_files" ]]; then
  export TAXONKIT_DB="${HYMET_ROOT}/taxonomy_files"
fi

python3 "${SCRIPT_DIR}/convert/phyloflash_to_cami.py" \
  --input "${NTU_FILE}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool phyloflash \
  --taxdb "${TAXONKIT_DB:-}" \
  --classified-out "${CLASSIFIED_DST}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "phyloflash",
  "profile": "${PROFILE_DST}",
  "contigs": "${CONTIGS_ABS}",
  "reads": "${READS_FASTQ}",
  "db_dir": "${PHYLOFLASH_DB_DIR}",
  "threads": "${THREADS}"
}
EOF

normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
