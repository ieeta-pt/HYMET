#!/usr/bin/env bash
# Lightweight SnakeMAGs-style classification harness for CAMI benchmarking.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"

SNAKEMAGS_SNAKEFILE="${SNAKEMAGS_SNAKEFILE:-${SCRIPT_DIR}/snakes/snakemags_classify.smk}"
SNAKEMAGS_SNKMK_CMD="${SNAKEMAGS_SNKMK_CMD:-snakemake}"
SNAKEMAGS_GTDB="${SNAKEMAGS_GTDB:-}"
SNAKEMAGS_MIN_CONTIG="${SNAKEMAGS_MIN_CONTIG:-5000}"
SNAKEMAGS_KEEP_WORK="${SNAKEMAGS_KEEP_WORK:-0}"
SNAKEMAGS_EXTRA_ARGS="${SNAKEMAGS_EXTRA_ARGS:-}"

usage(){
  cat <<'USAGE'
Usage: run_snakemags.sh --sample ID --contigs FASTA [--out DIR] [--threads N]

Environment variables:
  SNAKEMAGS_GTDB          Path to GTDB-Tk data directory (required).
  SNAKEMAGS_SNAKEFILE     Snakemake workflow file (default: snakes/snakemags_classify.smk).
  SNAKEMAGS_SNKMK_CMD     Command used to invoke Snakemake (default: snakemake).
  SNAKEMAGS_MIN_CONTIG    Minimum contig length (bp) to treat as a MAG (default: 5000).
  SNAKEMAGS_KEEP_WORK     Set to 1 to retain the intermediate work directory.
  SNAKEMAGS_EXTRA_ARGS    Extra arguments appended to the Snakemake command.
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) usage;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage
[[ -n "${SNAKEMAGS_GTDB}" ]] || die "SNAKEMAGS_GTDB must point to a GTDB-Tk data directory"

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/snakemags}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
MAG_DIR="${RUN_DIR}/mags"
WORK_DIR="${RUN_DIR}/work"
OUTPUT_DIR="${RUN_DIR}/output"
CONFIG_PATH="${RUN_DIR}/config.yaml"
MAPPING_PATH="${RUN_DIR}/mag_mapping.tsv"
SUMMARY_PATH="${OUTPUT_DIR}/gtdbtk.summary.tsv"

ensure_dir "${RUN_DIR}"
ensure_dir "${MAG_DIR}"
ensure_dir "${WORK_DIR}"
ensure_dir "${OUTPUT_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "SnakeMAGs input FASTA missing (${CONTIGS_ABS})"

SNAKEMAGS_SNAKEFILE="$(resolve_path "${SNAKEMAGS_SNAKEFILE}")"
[[ -f "${SNAKEMAGS_SNAKEFILE}" ]] || die "SnakeMAGs snakefile not found (${SNAKEMAGS_SNAKEFILE})"

SNAKEMAGS_GTDB="$(resolve_path "${SNAKEMAGS_GTDB}")"
[[ -d "${SNAKEMAGS_GTDB}" ]] || die "GTDB directory not found (${SNAKEMAGS_GTDB})"

log "Preparing MAG surrogates for SnakeMAGs (${SAMPLE})"
if ! MAG_COUNT=$(python3 - <<'PY' "${CONTIGS_ABS}" "${MAG_DIR}" "${MAPPING_PATH}" "${SNAKEMAGS_MIN_CONTIG}"
import re
import sys
from pathlib import Path

contigs_path, mag_dir, map_path, min_len_raw = sys.argv[1:5]
min_len = max(int(min_len_raw), 1)

def iter_fasta(path: Path):
    name = None
    seq_parts = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper()
                name = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if name is not None:
        yield name, "".join(seq_parts).upper()

contigs = Path(contigs_path)
mag_dir_path = Path(mag_dir)
mag_dir_path.mkdir(parents=True, exist_ok=True)

mapping = Path(map_path)
mapping.parent.mkdir(parents=True, exist_ok=True)

def sanitise(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    if not cleaned:
        cleaned = "mag"
    return cleaned[:120]

seen = set()
written = 0

with mapping.open("w", encoding="utf-8") as map_handle:
    map_handle.write("mag_id\tcontig_id\tlength\n")
    for contig_id, seq in iter_fasta(contigs):
        length = len(seq)
        if length < min_len:
            continue
        base = sanitise(contig_id)
        candidate = base
        suffix = 1
        while candidate in seen:
            suffix += 1
            candidate = f"{base}_{suffix}"
        seen.add(candidate)
        mag_path = mag_dir_path / f"{candidate}.fa"
        with mag_path.open("w", encoding="utf-8") as mag_handle:
            mag_handle.write(f">{contig_id}\n")
            for idx in range(0, length, 80):
                mag_handle.write(seq[idx:idx+80] + "\n")
        map_handle.write(f"{candidate}\t{contig_id}\t{length}\n")
        written += 1

print(written)
PY
); then
  die "Failed to construct MAG surrogates for SnakeMAGs"
fi
MAG_COUNT="${MAG_COUNT//$'\n'/}"

if [[ -z "${MAG_COUNT}" || "${MAG_COUNT}" == "0" ]]; then
  log "No contigs met the minimum length (${SNAKEMAGS_MIN_CONTIG} bp); producing empty profile."
  : > "${SUMMARY_PATH}"
else
  PPLACER_THREADS="${THREADS}"
  if ! [[ "${PPLACER_THREADS}" =~ ^[0-9]+$ ]]; then
    PPLACER_THREADS=4
  fi
  if [[ "${PPLACER_THREADS}" -gt 4 ]]; then
    PPLACER_THREADS=4
  fi

  cat > "${CONFIG_PATH}" <<EOF
sample: "${SAMPLE}"
mag_dir: "${MAG_DIR}"
work_dir: "${WORK_DIR}"
out_dir: "${OUTPUT_DIR}"
gtdb: "${SNAKEMAGS_GTDB}"
threads: ${THREADS}
pplacer_threads: ${PPLACER_THREADS}
EOF

  IFS=' ' read -r -a SNKMK_CMD_ARR <<<"${SNAKEMAGS_SNKMK_CMD}"
  [[ ${#SNKMK_CMD_ARR[@]} -gt 0 ]] || die "SNAKEMAGS_SNKMK_CMD is empty"
  read -r -a SNKMK_EXTRA_ARR <<<"${SNAKEMAGS_EXTRA_ARGS}"

  log "Running SnakeMAGs Snakemake workflow for ${SAMPLE}"
  pushd "${RUN_DIR}" >/dev/null
  set +e
  # ensure GTDBTK_DATA_PATH is visible to all spawned jobs
  export GTDBTK_DATA_PATH="${SNAKEMAGS_GTDB}"
  "${SNKMK_CMD_ARR[@]}" \
    --cores "${THREADS}" \
    --snakefile "${SNAKEMAGS_SNAKEFILE}" \
    --configfile "${CONFIG_PATH}" \
    --directory "${RUN_DIR}" \
    --rerun-incomplete \
    --nolock \
    --quiet \
    "${SNKMK_EXTRA_ARR[@]}"
  status=$?
  set -e
  popd >/dev/null
  [[ ${status} -eq 0 ]] || die "SnakeMAGs Snakemake execution failed with exit code ${status}"

  [[ -f "${SUMMARY_PATH}" ]] || : > "${SUMMARY_PATH}"
fi

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"
CONVERTER="${SCRIPT_DIR}/convert/snakemags_to_cami.py"

python3 "${CONVERTER}" \
  --summary "${SUMMARY_PATH}" \
  --mapping "${MAPPING_PATH}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool snakemags \
  --taxdb "${TAXONKIT_DB:-}" \
  --classified-out "${CLASSIFIED_DST}"

if [[ "${SNAKEMAGS_KEEP_WORK}" != "1" ]]; then
  rm -rf "${RUN_DIR:?}"
fi

log "SnakeMAGs processing completed for ${SAMPLE}"
