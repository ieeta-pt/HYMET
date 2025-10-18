#!/usr/bin/env bash
# run_hymet_cami.sh — HYMET end‑to‑end: Mash → Download → Minimap2 → Classify → CAMI export
#
# Usage examples:
#   bash run_hymet_cami.sh
#   INPUT_FASTA=/data/cami/sample_0.fna OUTDIR=/data/hymet_out/sample_0 THREADS=16 bash run_hymet_cami.sh
#
# Notes:
# - Works with the fixed classifier in scripts/classification_cami.py. If that fails (empty result),
#   a robust fallback (first-hit mapping via Identifiers) kicks in automatically.
# - Avoids copying a file onto itself and includes clear diagnostics.
#
# Reviewer comment to retain in repo docs:
# * Benchmarking is currently based solely on the authors' own simulated data. This is insufficient to validate
#   performance claims. Standard datasets such as those from the CAMI challenge should be incorporated.

set -Eeuo pipefail

# ----------------- params -----------------
SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-${SCRIPT_ROOT}}"

INPUT_FASTA="${INPUT_FASTA:-/data/cami/sample_0.fna}"
OUTDIR="${OUTDIR:-/data/hymet_out/sample_0}"
THREADS="${THREADS:-8}"
CAND_MAX="${CAND_MAX:-5000}"      # max candidates after Mash
SPECIES_DEDUP="${SPECIES_DEDUP:-0}"
ASSEMBLY_SUMMARY_DIR="${ASSEMBLY_SUMMARY_DIR:-${ROOT}/data/downloaded_genomes/assembly_summaries}"
CAND_LIMIT_LOG="${CAND_LIMIT_LOG:-}"
SPLIT_IDX="${SPLIT_IDX:-2g}"      # minimap2 -I chunk size for lower RAM
MASH_THRESH="${MASH_THRESH:-0.9}" # Mash screen threshold
FORCE_DOWNLOAD="${FORCE_DOWNLOAD:-0}"
CACHE_ROOT="${CACHE_ROOT:-data/downloaded_genomes/cache}"
export PATH="/opt/conda/bin:$PATH"
export TMPDIR="${TMPDIR:-/data/tmp}"
export OMP_NUM_THREADS="$THREADS"
export TAXONKIT_DB="$ROOT/taxonomy_files"
umask 0022
# -----------------------------------------

log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }
trap 'die "failed at line $LINENO"' ERR

# ensure required dirs exist
mkdir -p "$TMPDIR"

cd "$ROOT" || die "cannot cd to $ROOT"

# ---- checks ----
[ -s "$INPUT_FASTA" ] || die "missing FASTA $INPUT_FASTA"
for f in data/sketch1.msh data/sketch2.msh data/sketch3.msh data/detailed_taxonomy.tsv data/taxonomy_hierarchy.tsv \
         scripts/mash.sh scripts/minimap2.sh scripts/classification_cami.py; do
  [ -s "$f" ] || die "missing $f"
done
[ -d taxonomy_files ] || { log "taxonomy_files missing → running ./config.pl"; ./config.pl; }

# Check for required external tools
for tool in "${ROOT}/tools/hymet2cami.py" "${ROOT}/tools/build_id_map.py" "${ROOT}/tools/mini_classify.py"; do
  [ -f "$tool" ] || die "missing required tool: $tool"
done

# ---- IO layout ----
WORKDIR="$OUTDIR/work"
rm -rf output 2>/dev/null || true
mkdir -p "$OUTDIR" "$WORKDIR" "$OUTDIR/logs" input cache data/downloaded_genomes
ln -s "$WORKDIR" output 2>/dev/null || true
cp -f "$INPUT_FASTA" input/sample_0.fna

# ---- deps ----
command -v minimap2 >/dev/null || micromamba install -y -p /opt/conda -c bioconda minimap2
command -v mash     >/dev/null || micromamba install -y -p /opt/conda -c bioconda mash
command -v taxonkit >/dev/null || micromamba install -y -p /opt/conda -c bioconda taxonkit

# ---- make minimap2 memory-friendly (idempotent) ----
# ensure all invocations inside scripts/minimap2.sh carry -I$SPLIT_IDX
if ! grep -Eq 'minimap2 .* -I[0-9]+' scripts/minimap2.sh; then
  sed -i "s/minimap2 -d /minimap2 -I${SPLIT_IDX} -d /g" scripts/minimap2.sh || true
  sed -i "s/minimap2 -t /minimap2 -I${SPLIT_IDX} -t /g" scripts/minimap2.sh || true
fi

# 1) Mash → selected_genomes
if [ ! -s output/selected_genomes.txt ]; then
  log "Mash screen (sketch1)"
  bash scripts/mash.sh input data/sketch1.msh \
    output/screen.tab output/filtered_screen.tab output/sorted_screen.tab \
    output/top_hits.tab output/selected_genomes.txt "$MASH_THRESH"
  log "Mash screen (sketch2)"
  bash scripts/mash.sh input data/sketch2.msh \
    output/gtdb_screen.tab output/gtdb_filtered.tab output/gtdb_sorted.tab \
    output/gtdb_top_hits.tab output/gtdb_selected_genomes.txt "$MASH_THRESH" || true
  cat output/gtdb_selected_genomes.txt >> output/selected_genomes.txt 2>/dev/null || true
  log "Mash screen (sketch3)"
  bash scripts/mash.sh input data/sketch3.msh \
    output/custom_screen.tab output/custom_filtered.tab output/custom_sorted.tab \
    output/custom_top_hits.tab output/custom_selected_genomes.txt "$MASH_THRESH" || true
  cat output/custom_selected_genomes.txt >> output/selected_genomes.txt 2>/dev/null || true
  sort -u -o output/selected_genomes.txt output/selected_genomes.txt
fi

# 2) deduplicate/limit candidates if requested
LIMIT_TMP="output/selected_genomes.limited.txt"
LIMIT_ARGS=(
  --selected output/selected_genomes.txt
  --output "${LIMIT_TMP}"
  --max "${CAND_MAX}"
  --assembly-dir "${ASSEMBLY_SUMMARY_DIR}"
  --score-file output/sorted_screen.tab
)
if [ -s output/gtdb_sorted.tab ]; then
  LIMIT_ARGS+=(--score-file output/gtdb_sorted.tab)
fi
if [ -s output/custom_sorted.tab ]; then
  LIMIT_ARGS+=(--score-file output/custom_sorted.tab)
fi
if [ "${SPECIES_DEDUP}" -eq 1 ]; then
  LIMIT_ARGS+=(--dedupe)
fi
if [ -n "${CAND_LIMIT_LOG}" ]; then
  LIMIT_ARGS+=(--log "${CAND_LIMIT_LOG}")
fi
python3 scripts/limit_candidates.py "${LIMIT_ARGS[@]}"
mv "${LIMIT_TMP}" output/selected_genomes.txt

LINES=$(wc -l < output/selected_genomes.txt || echo 0)
[ "$LINES" -gt 0 ] || die "candidate list empty after applying limit"

#  cache key derived from selected genomes
ensure_cache_dirs(){
  local dir="$1"
  mkdir -p "${dir}"
}

ensure_cache_dirs "${CACHE_ROOT}"
CACHE_KEY="$(sha1sum output/selected_genomes.txt | awk '{print $1}')"
[ -n "${CACHE_KEY}" ] || die "Failed to compute cache key from selected genomes."
CACHE_DIR="${CACHE_ROOT}/${CACHE_KEY}"
CACHE_FASTA="${CACHE_DIR}/combined_genomes.fasta"
CACHE_TAX="${CACHE_DIR}/detailed_taxonomy.tsv"
CACHE_DL="${CACHE_DIR}/download_cache"
CACHE_MMI="${CACHE_DIR}/reference.mmi"
ensure_cache_dirs "${CACHE_DIR}"
ensure_cache_dirs "${CACHE_DL}"
log "cache key ${CACHE_KEY} → ${CACHE_DIR}"

# 3) download refs + concat
if [ "${FORCE_DOWNLOAD}" -eq 1 ]; then
  log "FORCE_DOWNLOAD=1 → clearing cached reference"
  rm -f "${CACHE_FASTA}" "${CACHE_TAX}" "${CACHE_MMI}"
fi
if [ ! -s "${CACHE_FASTA}" ]; then
  log "downloadDB.py (cache ${CACHE_KEY})"
  python3 scripts/downloadDB.py \
    output/selected_genomes.txt \
    "${CACHE_DIR}" \
    "${CACHE_TAX}" \
    "${CACHE_DL}"
else
  log "cache hit for ${CACHE_KEY}; reusing ${CACHE_FASTA}"
fi
mkdir -p "${ROOT}/data/downloaded_genomes"
mkdir -p "${ROOT}/data"
ln -sf "$(readlink -f "${CACHE_FASTA}")" "${ROOT}/data/downloaded_genomes/combined_genomes.fasta"
ln -sf "$(readlink -f "${CACHE_TAX}")" "${ROOT}/data/detailed_taxonomy.tsv"

# 4) minimap2 index+map
if [ ! -s output/resultados.paf ]; then
  log "minimap2 index+map"
  bash scripts/minimap2.sh input "${CACHE_FASTA}" "${CACHE_MMI}" output/resultados.paf
  ln -sf "$(readlink -f "${CACHE_MMI}")" output/reference.mmi
fi

# 5) classify (primary path)
log "classification_cami.py on raw PAF"
python3 scripts/classification_cami.py \
  --paf output/resultados.paf \
  --taxonomy data/detailed_taxonomy.tsv \
  --hierarchy data/taxonomy_hierarchy.tsv \
  --output output/classified_sequences.tsv \
  --processes "$THREADS" || true

ROWS=$(wc -l < output/classified_sequences.tsv 2>/dev/null || echo 0)
if [ "$ROWS" -lt 2 ]; then
  log "primary classification empty → running robust fallback (first-hit via Identifiers)"

  # Use existing fallback tools
  BUILD_ID_MAP="${ROOT}/tools/build_id_map.py"
  MINI_CLASSIFY="${ROOT}/tools/mini_classify.py"
  
  [ -f "$BUILD_ID_MAP" ] || die "missing build_id_map.py: $BUILD_ID_MAP"
  [ -f "$MINI_CLASSIFY" ] || die "missing mini_classify.py: $MINI_CLASSIFY"

  python3 "$BUILD_ID_MAP" data/detailed_taxonomy.tsv "$WORKDIR/id_to_taxid.tsv" | tee "$OUTDIR/logs/step_fallback_build_map.log"
  python3 "$MINI_CLASSIFY" output/resultados.paf "$WORKDIR/id_to_taxid.tsv" "$WORKDIR/fallback_classified.tsv" | tee "$OUTDIR/logs/step_fallback_classify.log"

  # Convert mini_classify output to HYMET format if needed
  if [ -s "$WORKDIR/fallback_classified.tsv" ]; then
    python3 - "$WORKDIR/fallback_classified.tsv" data/taxonomy_hierarchy.tsv output/classified_sequences.tsv <<'PY'
import csv, sys

src, hier, dest = sys.argv[1:4]
rank_order = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

def determine_level(lineage: str) -> str:
    if not lineage or lineage.lower() == 'unknown':
        return 'root'
    level = 'root'
    for part in lineage.split(';'):
        part = part.strip()
        if ':' not in part:
            continue
        rank, _ = part.split(':', 1)
        rank = rank.strip().lower()
        if rank in rank_order:
            level = rank
    return level

hierarchy = {}
with open(hier, newline='') as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        tid = (row.get('TaxID') or '').strip()
        lineage = (row.get('Lineage') or '').strip() or 'Unknown'
        if tid:
            hierarchy[tid] = lineage

with open(src, newline='') as inp, open(dest, 'w', newline='') as outp:
    reader = csv.DictReader(inp, delimiter='\t')
    writer = csv.writer(outp, delimiter='\t')
    writer.writerow(['Query', 'Lineage', 'Taxonomic Level', 'TaxID', 'Confidence'])
    for row in reader:
        query = row.get('qname') or row.get('Query') or ''
        taxid = (row.get('taxid') or '').strip() or 'Unknown'
        lineage = hierarchy.get(taxid, 'Unknown')
        level = determine_level(lineage)
        writer.writerow([query, lineage, level, taxid, '1.0000'])
PY
  fi

  ROWS=$(wc -l < output/classified_sequences.tsv 2>/dev/null || echo 0)
  [ "$ROWS" -ge 2 ] || die "classification still empty after fallback"
fi

# 6) persist result (avoid copy-on-self)
SRC="$(readlink -f output/classified_sequences.tsv)"
DST="$(readlink -f "$OUTDIR/classified_sequences.tsv")"
if [ "$SRC" != "$DST" ]; then cp -f "$SRC" "$DST"; fi
[ -s "$OUTDIR/classified_sequences.tsv" ] || die "empty classification"

# 7) CAMI export (robust to 'rank:name' or plain 'name' lineage)
HYMET2CAMI_SCRIPT="${ROOT}/tools/hymet2cami.py"
[ -f "$HYMET2CAMI_SCRIPT" ] || die "missing hymet2cami script: $HYMET2CAMI_SCRIPT"

python3 "$HYMET2CAMI_SCRIPT" "$OUTDIR/classified_sequences.tsv" > "$OUTDIR/hymet.sample_0.cami.tsv"

log "OK: $OUTDIR/classified_sequences.tsv"
log "OK: $OUTDIR/hymet.sample_0.cami.tsv"

# small preview
echo "[preview] classified_sequences.tsv:"
head -n 5 "$OUTDIR/classified_sequences.tsv" || true
echo "[preview] CAMI profile:"
head -n 5 "$OUTDIR/hymet.sample_0.cami.tsv" || true
