#!/usr/bin/env bash
# Wrapper to execute MegaPath-Nano taxonomy profiling within the HYMET benchmark harness.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
OUT_DIR=""
THREADS="${THREADS:-8}"

MPN_ROOT="${MPN_ROOT:-}"
MPN_CMD="${MPN_CMD:-python3}"
MPN_SCRIPT="${MPN_SCRIPT:-megapath_nano.py}"
MPN_CHUNK_SIZE="${MPN_CHUNK_SIZE:-50000}"
MPN_MIN_CHUNK="${MPN_MIN_CHUNK:-5000}"
MPN_KEEP_WORK="${MPN_KEEP_WORK:-0}"
MPN_EXTRA_ARGS="${MPN_EXTRA_ARGS:-}"

usage(){
  cat <<'USAGE'
Usage: run_megapath_nano.sh --sample ID --contigs FASTA [--out DIR] [--threads N]

Environment variables:
  MPN_ROOT         Path to the MegaPath-Nano repository (required).
  MPN_CMD          Command prefix used to invoke Python (default: python3).
  MPN_SCRIPT       Entry script name (default: megapath_nano.py).
  MPN_CHUNK_SIZE   Chunk size (bp) when slicing contigs into synthetic reads (default: 10000).
  MPN_MIN_CHUNK    Minimum chunk size threshold (default: 2500).
  MPN_KEEP_WORK    Set to 1 to retain the raw MegaPath output directory.
  MPN_EXTRA_ARGS   Additional arguments passed verbatim to MegaPath-Nano.
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
[[ -n "${MPN_ROOT}" ]] || die "MPN_ROOT must point to the MegaPath-Nano repository"

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/megapath_nano}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"

RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
[[ -s "${CONTIGS_ABS}" ]] || die "MegaPath-Nano input FASTA missing (${CONTIGS_ABS})"

MPN_ROOT="$(resolve_path "${MPN_ROOT}")"
[[ -d "${MPN_ROOT}" ]] || die "MPN_ROOT directory not found (${MPN_ROOT})"

if [[ ! -f "${MPN_ROOT}/bin/${MPN_SCRIPT}" ]]; then
  die "MegaPath-Nano script not found at ${MPN_ROOT}/bin/${MPN_SCRIPT}"
fi

log "Preparing synthetic reads for MegaPath-Nano (${SAMPLE})"
QUERY_FASTQ="${RUN_DIR}/input_reads.fastq"
ID_MAP="${RUN_DIR}/read_to_contig.tsv"

if ! READ_COUNT=$(python3 - <<'PY' "${CONTIGS_ABS}" "${QUERY_FASTQ}" "${ID_MAP}" "${MPN_CHUNK_SIZE}" "${MPN_MIN_CHUNK}"
import sys
from pathlib import Path

contigs_path, fastq_path, map_path, chunk_raw, min_chunk_raw = sys.argv[1:6]
chunk_size = max(int(chunk_raw), 1)
min_chunk = max(int(min_chunk_raw), 1)

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
if not contigs.is_file():
    raise SystemExit(f"Input FASTA not found: {contigs}")

fastq_file = Path(fastq_path)
fastq_file.parent.mkdir(parents=True, exist_ok=True)
map_file = Path(map_path)

total_reads = 0
with fastq_file.open("w", encoding="utf-8") as fastq, map_file.open("w", encoding="utf-8") as mapping:
    mapping.write("read_id\tcontig_id\n")
    for contig_id, seq in iter_fasta(contigs):
        if not seq:
            continue
        seq = seq.replace("N", "A")
        length = len(seq)
        if length <= chunk_size:
            read_id = f"{contig_id}|0|{length}"
            fastq.write(f"@{read_id}\n{seq}\n+\n{'I' * length}\n")
            mapping.write(f"{read_id}\t{contig_id}\n")
            total_reads += 1
            continue
        start = 0
        while start < length:
            end = min(start + chunk_size, length)
            chunk = seq[start:end]
            if len(chunk) < min_chunk and start != 0:
                break
            read_id = f"{contig_id}|{start}|{end}"
            fastq.write(f"@{read_id}\n{chunk}\n+\n{'I' * len(chunk)}\n")
            mapping.write(f"{read_id}\t{contig_id}\n")
            total_reads += 1
            start += chunk_size

if total_reads == 0:
    raise SystemExit("No reads generated from input contigs.")
print(total_reads)
PY
); then
  die "Failed to synthesise reads for MegaPath-Nano"
fi

READ_COUNT="${READ_COUNT//$'\n'/}"

if [[ -z "${READ_COUNT}" || "${READ_COUNT}" == "0" ]]; then
  die "MegaPath-Nano synthetic read generation produced no reads"
fi

log "MegaPath-Nano synthetic read count: ${READ_COUNT}"

WORK_TMP="${RUN_DIR}/tmp"
WORK_RAM="${RUN_DIR}/ram"
MPN_OUTPUT="${RUN_DIR}/mpn_output"
ensure_dir "${WORK_TMP}"
ensure_dir "${WORK_RAM}"
ensure_dir "${MPN_OUTPUT}"

export OMP_NUM_THREADS="${THREADS}"
export PATH="${MPN_ROOT}/bin:${PATH}"

IFS=' ' read -r -a MPN_CMD_ARR <<<"${MPN_CMD}"
[[ ${#MPN_CMD_ARR[@]} -gt 0 ]] || die "MPN_CMD is empty"

read -r -a MPN_EXTRA_ARR <<<"${MPN_EXTRA_ARGS}"

MPN_SCRIPT_ABS="${MPN_ROOT}/bin/${MPN_SCRIPT}"
ALN_THREADS="${THREADS}"
if ! [[ "${ALN_THREADS}" =~ ^[0-9]+$ ]]; then
  ALN_THREADS=8
fi
if [[ "${ALN_THREADS}" -gt 64 ]]; then
  ALN_THREADS=64
fi

log "Running MegaPath-Nano for sample ${SAMPLE}"
# Run from the per-sample run directory so any stray temp files land under RUN_DIR
# instead of polluting the tool's bin directory.
export TMPDIR="${WORK_TMP}"
pushd "${RUN_DIR}" >/dev/null
set +e
"${MPN_CMD_ARR[@]}" "${MPN_SCRIPT_ABS}" \
  --query "${QUERY_FASTQ}" \
  --output_prefix "${SAMPLE}" \
  --output_folder "${MPN_OUTPUT}" \
  --temp_folder "${WORK_TMP}" \
  --RAM_folder "${WORK_RAM}" \
  --max_aligner_thread "${ALN_THREADS}" \
  --taxon_module_only \
  --no-adaptor_trimming \
  --no-read_trimming \
  --no-read_filter \
  --no-human_filter \
  --no-decoy_filter \
  --no-variable_region_adjustment \
  --no-spike_filter \
  --no-closing_spike_filter \
  --no-human_repetitive_region_filter \
  --no-microbe_repetitive_region_filter \
  --no-short_alignment_filter \
  --no-noise_projection \
  --no-similar_species_marker \
  --no-output_PAF \
  --output_noise_stat \
  --no-output_separate_noise_bed \
  --no-output_human_stat \
  --no-output_decoy_stat \
  --no-output_id_signal \
  --no-output_raw_signal \
  --no-output_per_read_data \
  --no-output_quality_score_histogram \
  --no-output_read_length_histogram \
  --no-output_genome_set \
  "${MPN_EXTRA_ARR[@]}"
status=$?
set -e
popd >/dev/null
[[ ${status} -eq 0 ]] || die "MegaPath-Nano execution failed with exit code ${status}"

ASSEMBLY_STAT="${MPN_OUTPUT}/${SAMPLE}.microbe_stat"
SEQ_STAT="${MPN_OUTPUT}/${SAMPLE}.microbe_stat_by_sequence_id_assembly_info"
if [[ ! -s "${SEQ_STAT}" ]]; then
  mapfile -t candidates < <(find "${MPN_OUTPUT}" -maxdepth 1 -type f -name '*.microbe_stat_by_sequence_id_assembly_info' | sort)
  if [[ ${#candidates[@]} -gt 0 ]]; then
    SEQ_STAT="${candidates[0]}"
  fi
fi
[[ -s "${SEQ_STAT}" ]] || die "MegaPath-Nano did not produce microbe_stat_by_sequence_id_assembly_info"

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"
CONVERTER="${SCRIPT_DIR}/convert/megapath_nano_to_cami.py"

python3 "${CONVERTER}" \
  --input "${SEQ_STAT}" \
  --assembly-stat "${ASSEMBLY_STAT}" \
  --out "${PROFILE_DST}" \
  --sample-id "${SAMPLE}" \
  --tool megapath_nano \
  --taxdb "${TAXONKIT_DB:-}" \
  --classified-out "${CLASSIFIED_DST}" \
  --id-map "${ID_MAP}"

if [[ "${MPN_KEEP_WORK}" != "1" ]]; then
  rm -rf "${MPN_OUTPUT:?}" "${WORK_TMP:?}" "${WORK_RAM:?}"
  rm -f "${QUERY_FASTQ}" "${ID_MAP}"
fi

log "MegaPath-Nano completed for ${SAMPLE}"
