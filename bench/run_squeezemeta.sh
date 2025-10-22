#!/usr/bin/env bash
# Wrapper to execute SqueezeMeta sequential workflow on CAMI contigs.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"
SQUEEZEMETA_ENV_PREFIX="${SQUEEZEMETA_ENV_PREFIX:-/opt/envs/squeezemeta}"
SQUEEZEMETA_DB_DIR="${SQUEEZEMETA_DB_DIR:-}"
SQUEEZEMETA_EXTRA_OPTS="${SQUEEZEMETA_EXTRA_OPTS:-}"
SQUEEZEMETA_SYNTH_FRAG_LEN="${SQUEEZEMETA_SYNTH_FRAG_LEN:-150}"

usage(){
  cat <<'USAGE'
Usage: run_squeezemeta.sh --sample ID --contigs FASTA --db DIR [--out DIR] [--threads N]

Environment variables:
  SQUEEZEMETA_ENV_PREFIX   Path to SqueezeMeta conda environment (default: /opt/envs/squeezemeta).
  SQUEEZEMETA_DB_DIR       Directory with SqueezeMeta databases (required).
  SQUEEZEMETA_EXTRA_OPTS   Extra options appended to SqueezeMeta.pl invocation.
  SQUEEZEMETA_SYNTH_FRAG_LEN  Fragment length used for synthetic reads (default: 150).
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --db) SQUEEZEMETA_DB_DIR="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/squeezemeta}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "Input contigs FASTA missing (${CONTIGS_ABS})"

SQUEEZEMETA_DB_DIR="$(resolve_path "${SQUEEZEMETA_DB_DIR}")"
[[ -d "${SQUEEZEMETA_DB_DIR}" ]] || die "SqueezeMeta database directory not found (${SQUEEZEMETA_DB_DIR})"

[[ -d "${SQUEEZEMETA_ENV_PREFIX}" ]] || die "SqueezeMeta environment prefix not found (${SQUEEZEMETA_ENV_PREFIX})"
need_cmd micromamba

run_in_env(){
  micromamba run -p "${SQUEEZEMETA_ENV_PREFIX}" "$@"
}

LIB_NAME="${SAMPLE//[^A-Za-z0-9_]/_}"
[[ -n "${LIB_NAME}" ]] || LIB_NAME="sample"

READS_DIR="${RUN_DIR}/reads"
ensure_dir "${READS_DIR}"
SAMPLE_FASTA="${READS_DIR}/${LIB_NAME}.fasta"

if [[ ! -s "${SAMPLE_FASTA}" ]]; then
  python3 - "${CONTIGS_ABS}" "${SAMPLE_FASTA}" "${SQUEEZEMETA_SYNTH_FRAG_LEN}" <<'PY'
import sys
from pathlib import Path

src_path, dst_path, frag_len_str = sys.argv[1:4]
fragment_len = int(frag_len_str)
src = Path(src_path)
dst = Path(dst_path)

def emit_fragments(sequence):
    seq = sequence.upper()
    for i in range(0, len(seq), fragment_len):
        frag = seq[i : i + fragment_len]
        if len(frag) < fragment_len:
            break
        yield frag

with src.open("r") as fin, dst.open("w") as fout:
    current_id = None
    seq_chunks = []
    for line in fin:
        if line.startswith(">"):
            if current_id is not None:
                sequence = "".join(seq_chunks)
                for idx, frag in enumerate(emit_fragments(sequence), start=1):
                    fout.write(f">{current_id}_{idx}\n{frag}\n")
            current_id = line[1:].strip().split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if current_id is not None:
        sequence = "".join(seq_chunks)
        for idx, frag in enumerate(emit_fragments(sequence), start=1):
            fout.write(f">{current_id}_{idx}\n{frag}\n")
PY
fi

SAMPLE_R1="${READS_DIR}/${LIB_NAME}_R1.fastq"
SAMPLE_R2="${READS_DIR}/${LIB_NAME}_R2.fastq"

if [[ ! -s "${SAMPLE_R1}" ]]; then
  python3 - "${SAMPLE_FASTA}" "${SAMPLE_R1}" "${SAMPLE_R2}" <<'PY'
import sys
from pathlib import Path

src_path, r1_path, r2_path = map(Path, sys.argv[1:4])

def to_fastq(fasta_path, fastq_path):
    with fasta_path.open("r") as fin, fastq_path.open("w") as fout:
        current_id = None
        seq = []
        for line in fin:
            if line.startswith(">"):
                if current_id is not None:
                    sequence = "".join(seq)
                    fout.write(f"@{current_id}\n{sequence}\n+\n{'I'*len(sequence)}\n")
                current_id = line[1:].strip().split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if current_id is not None:
            sequence = "".join(seq)
            fout.write(f"@{current_id}\n{sequence}\n+\n{'I'*len(sequence)}\n")

to_fastq(src_path, r1_path)
to_fastq(src_path, r2_path)
PY
fi

SAMPLES_TXT="${RUN_DIR}/samples.txt"
cat > "${SAMPLES_TXT}" <<EOF
${LIB_NAME}	${SAMPLE_R1}	${SAMPLE_R2}
EOF

PROJECT_NAME="${LIB_NAME}"
PROJECT_DIR="${RUN_DIR}/${PROJECT_NAME}"

if [[ ! -d "${PROJECT_DIR}" ]]; then
  ensure_dir "${PROJECT_DIR}"
fi

SQUEEZEMETA_ARGS=(
  "SqueezeMeta.pl"
  -m sequential
  -p "${PROJECT_NAME}"
  -s "${SAMPLES_TXT}"
  -f "${READS_DIR}"
  -t "${THREADS}"
  -extassembly "${CONTIGS_ABS}"
  -dbdir "${SQUEEZEMETA_DB_DIR}"
)

if [[ -n "${SQUEEZEMETA_EXTRA_OPTS}" ]]; then
  # shellcheck disable=SC2206
  EXTRA=( ${SQUEEZEMETA_EXTRA_OPTS} )
  SQUEEZEMETA_ARGS+=("${EXTRA[@]}")
fi

PROJECT_RESULTS="${PROJECT_DIR}/results"

if [[ ! -d "${PROJECT_RESULTS}" ]]; then
  log "Running SqueezeMeta for ${SAMPLE}"
  run_in_env "${SQUEEZEMETA_ARGS[@]}"
else
  log "Reusing existing SqueezeMeta project at ${PROJECT_DIR}"
fi

SUMMARY_FILE="$(find "${PROJECT_RESULTS}" -maxdepth 2 -type f -name 'contig_taxonomy.summary' | head -n 1 || true)"
[[ -n "${SUMMARY_FILE}" && -s "${SUMMARY_FILE}" ]] || die "Unable to locate contig_taxonomy.summary under ${PROJECT_RESULTS}"

if [[ -z "${TAXONKIT_DB:-}" && -d "${HYMET_ROOT}/taxonomy_files" ]]; then
  export TAXONKIT_DB="${HYMET_ROOT}/taxonomy_files"
fi

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"

python3 "${SCRIPT_DIR}/convert/squeezemeta_to_cami.py" \
  --input "${SUMMARY_FILE}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool squeezemeta \
  --taxdb "${TAXONKIT_DB:-}" \
  --classified-out "${CLASSIFIED_DST}"

cat > "${OUT_DIR}/metadata.json" <<EOF
{"sample_id": "${SAMPLE}", "tool": "squeezemeta", "profile": "${PROFILE_DST}", "contigs": "${CONTIGS_ABS}", "db_dir": "${SQUEEZEMETA_DB_DIR}", "threads": "${THREADS}", "project_dir": "${PROJECT_DIR}"}
EOF
