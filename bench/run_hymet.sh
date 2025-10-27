#!/usr/bin/env bash
# Run HYMET classifier on a sample (wrapper for run_hymet_cami.sh).

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"

SAMPLE=""
CONTIGS=""
READS=""
MODE="contigs"
OUT_DIR=""
THREADS="${THREADS:-8}"

usage(){
  cat <<'USAGE'
Usage: run_hymet.sh --sample ID --contigs FASTA [--reads FASTQ] [--mode contigs|reads]
                    [--out DIR] [--threads N]
USAGE
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) SAMPLE="$2"; shift 2;;
    --contigs) CONTIGS="$2"; shift 2;;
    --reads) READS="$2"; shift 2;;
    --mode) MODE="$2"; shift 2;;
    --out) OUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "${SAMPLE}" && -n "${CONTIGS}" ]] || usage

case "$MODE" in
  contigs|reads) ;;
  *) die "Unsupported mode '${MODE}' (expected contigs or reads)";;
esac

if [[ -n "${READS}" && "${MODE}" != "reads" ]]; then
  MODE="reads"
fi

OUT_DIR="${OUT_DIR:-${BENCH_ROOT}/out/${SAMPLE}/hymet}"
OUT_DIR="$(resolve_path "${OUT_DIR}")"
ensure_dir "${OUT_DIR}"
RUN_DIR="${OUT_DIR}/run"
ensure_dir "${RUN_DIR}"

CONTIGS_ABS="$(resolve_path "${CONTIGS}")"
if [[ ! -s "${CONTIGS_ABS}" ]]; then
  die "Input contigs FASTA missing (${CONTIGS_ABS})"
fi

READS_ABS=""
GENERATED_READS=""
if [[ -n "${READS}" ]]; then
  READS_ABS="$(resolve_path "${READS}")"
  [[ -s "${READS_ABS}" ]] || die "Reads file missing (${READS_ABS})"
elif [[ "${MODE}" == "reads" ]]; then
  READS_ABS="${RUN_DIR}/input_reads.fastq"
  python3 "${BENCH_ROOT}/tools/contigs_to_reads.py" \
    --contigs "${CONTIGS_ABS}" \
    --out "${READS_ABS}" \
    --chunk-size "${READ_CHUNK_SIZE:-250}" \
    --min-chunk "${READ_MIN_CHUNK:-100}"
  GENERATED_READS="${READS_ABS}"
fi

CACHE_ROOT_EFFECTIVE="${CACHE_ROOT:-${HYMET_ROOT}/data/downloaded_genomes/cache_bench}"
CAND_MAX_EFFECTIVE="${CAND_MAX:-1500}"
SPECIES_DEDUP_EFFECTIVE="${SPECIES_DEDUP:-1}"
ASSEMBLY_SUMMARY_DIR_EFFECTIVE="${ASSEMBLY_SUMMARY_DIR:-${HYMET_ROOT}/data/downloaded_genomes/assembly_summaries}"
ensure_dir "${OUT_DIR}/logs"
CAND_LIMIT_LOG_PATH="${OUT_DIR}/logs/candidate_limit.log"

log "Running HYMET classifier for ${SAMPLE} (mode=${MODE})"
if [[ "${MODE}" == "reads" ]]; then
  env \
    INPUT_MODE="reads" \
    INPUT_READS="${READS_ABS}" \
    OUTDIR="${RUN_DIR}" \
    THREADS="${THREADS}" \
    ROOT="${HYMET_ROOT}" \
    CACHE_ROOT="${CACHE_ROOT_EFFECTIVE}" \
    CAND_MAX="${CAND_MAX_EFFECTIVE}" \
    SPECIES_DEDUP="${SPECIES_DEDUP_EFFECTIVE}" \
    ASSEMBLY_SUMMARY_DIR="${ASSEMBLY_SUMMARY_DIR_EFFECTIVE}" \
    CAND_LIMIT_LOG="${CAND_LIMIT_LOG_PATH}" \
    bash "${HYMET_ROOT}/run_hymet_cami.sh"
else
  env \
    INPUT_MODE="contigs" \
    INPUT_FASTA="${CONTIGS_ABS}" \
    OUTDIR="${RUN_DIR}" \
    THREADS="${THREADS}" \
    ROOT="${HYMET_ROOT}" \
    CACHE_ROOT="${CACHE_ROOT_EFFECTIVE}" \
    CAND_MAX="${CAND_MAX_EFFECTIVE}" \
    SPECIES_DEDUP="${SPECIES_DEDUP_EFFECTIVE}" \
    ASSEMBLY_SUMMARY_DIR="${ASSEMBLY_SUMMARY_DIR_EFFECTIVE}" \
    CAND_LIMIT_LOG="${CAND_LIMIT_LOG_PATH}" \
    bash "${HYMET_ROOT}/run_hymet_cami.sh"
fi

PROFILE_SRC="${RUN_DIR}/hymet.sample_0.cami.tsv"
CLASSIFIED_SRC="${RUN_DIR}/classified_sequences.tsv"
PAF_SRC="${RUN_DIR}/work/resultados.paf"
log "[debug] contents of ${RUN_DIR}:"
ls -l "${RUN_DIR}" || true

PROFILE_DST="${OUT_DIR}/profile.cami.tsv"
CLASSIFIED_DST="${OUT_DIR}/classified_sequences.tsv"
PAF_DST="${OUT_DIR}/resultados.paf"

if [[ -s "${PROFILE_SRC}" ]]; then
  cp -f "${PROFILE_SRC}" "${PROFILE_DST}"
else
  die "Expected HYMET CAMI profile not found at ${PROFILE_SRC}"
fi

CLASSIFIED_META=""
if [[ -s "${CLASSIFIED_SRC}" ]]; then
  cp -f "${CLASSIFIED_SRC}" "${CLASSIFIED_DST}"
  CLASSIFIED_META="${CLASSIFIED_DST}"
fi

PAF_META=""
if [[ -s "${PAF_SRC}" ]]; then
  cp -f "${PAF_SRC}" "${PAF_DST}"
  PAF_META="${PAF_DST}"
fi

if [[ "${KEEP_HYMET_WORK:-0}" -eq 0 ]]; then
  rm -f "${RUN_DIR}/work/reference.mmi"
  rm -f "${CLASSIFIED_SRC}" "${PAF_SRC}" "${PROFILE_SRC}"
  if [[ -n "${GENERATED_READS}" ]]; then
    rm -f "${GENERATED_READS}"
  fi
fi

git_commit="$(git -C "${HYMET_ROOT}" rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
git_dirty="$(git -C "${HYMET_ROOT}" status --porcelain 2>/dev/null || echo '')"
if [[ -n "${git_dirty}" && "${git_commit}" != "unknown" ]]; then git_commit="${git_commit}-dirty"; fi

tool_label="hymet"
if [[ "${MODE}" == "reads" ]]; then
  tool_label="hymet_reads"
fi

cat > "${OUT_DIR}/metadata.json" <<EOF
{
  "sample_id": "${SAMPLE}",
  "tool": "${tool_label}",
  "hymet_commit": "${git_commit}",
  "threads": ${THREADS},
  "cache_root": "${CACHE_ROOT_EFFECTIVE}",
  "cand_max": ${CAND_MAX_EFFECTIVE},
  "species_dedup": ${SPECIES_DEDUP_EFFECTIVE},
  "assembly_summary_dir": "${ASSEMBLY_SUMMARY_DIR_EFFECTIVE}",
  "profile": "${PROFILE_DST}",
  "contigs": "${CLASSIFIED_META}",
  "paf": "${PAF_META}",
  "input_mode": "${MODE}",
  "input_fasta": "${CONTIGS_ABS}",
  "input_reads": "${READS_ABS}",
  "run_dir": "${RUN_DIR}"
}
EOF

normalize_metadata_json "${OUT_DIR}/metadata.json" "${OUT_DIR}"
