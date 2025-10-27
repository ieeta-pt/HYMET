#!/usr/bin/env bash
# benchmark_cami.sh — run HYMET as-is, then evaluate vs CAMI (no classifier changes)
set -Eeuo pipefail

# ---- user-tunable (defaults match your layout) ----
ROOT="${ROOT:-/data/HYMET}"
INPUT_FASTA="${INPUT_FASTA:-/data/cami/sample_0.fna}"
OUTDIR="${OUTDIR:-/data/hymet_out/sample_0}"
THREADS="${THREADS:-16}"

# CAMI ground truth (both variants accepted; first one found is used)
TRUTH_CONTIGS_1="${TRUTH_CONTIGS_1:-/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping.tsv}"
TRUTH_CONTIGS_2="${TRUTH_CONTIGS_2:-/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/gsa_mapping_new.tsv}"
TRUTH_PROFILE="${TRUTH_PROFILE:-/data/cami/sample_0/taxonomic_profile_0.txt}"

# ---------------------------------------------------
log(){ printf "[%(%F %T)T] %s\n" -1 "$*"; }
die(){ log "ERROR: $*"; exit 1; }
trap 'die "failed at line $LINENO"' ERR

[ -s "$INPUT_FASTA" ] || die "missing INPUT_FASTA: $INPUT_FASTA"
[ -d "$ROOT" ] || die "missing HYMET root: $ROOT"

cd "$ROOT"

# 1) Run HYMET end-to-end (your existing script; classifier unchanged)
log "Running HYMET (unchanged)…"
INPUT_FASTA="$INPUT_FASTA" OUTDIR="$OUTDIR" THREADS="$THREADS" bash "$ROOT/run_hymet_cami.sh"

# 2) Choose which CAMI ground-truth contig mapping to use
TRUTH_CONTIGS=""
if [ -s "$TRUTH_CONTIGS_1" ]; then TRUTH_CONTIGS="$TRUTH_CONTIGS_1"; fi
if [ -z "$TRUTH_CONTIGS" ] && [ -s "$TRUTH_CONTIGS_2" ]; then TRUTH_CONTIGS="$TRUTH_CONTIGS_2"; fi
[ -n "$TRUTH_CONTIGS" ] || die "no CAMI contig mapping found"

# 3) Evaluate (pure post-processing; does NOT alter HYMET outputs)
log "Evaluating vs CAMI ground truth…"
python3 "$ROOT/tools/eval_cami.py" \
  --threads "$THREADS" \
  --pred-profile "$OUTDIR/hymet.sample_0.cami.tsv" \
  --truth-profile "$TRUTH_PROFILE" \
  --pred-contigs "$OUTDIR/work/classified_sequences.tsv" \
  --truth-contigs "$TRUTH_CONTIGS" \
  --pred-fasta "/data/cami/sample_0.fna" \
  --truth-fasta "/data/cami/sample_0/2017.12.29_11.37.26_sample_0/contigs/anonymous_gsa.fasta" \
  --taxdb "$ROOT/taxonomy_files" \
  --taxmap "$ROOT/data/detailed_taxonomy.tsv" \
  --paf "$OUTDIR/work/resultados.paf" \
  --outdir "$OUTDIR/eval"

log "Done. Key files:"
echo " - $OUTDIR/eval/summary.txt"
echo " - $OUTDIR/eval/profile_summary.tsv"
echo " - $OUTDIR/eval/contigs_per_rank.tsv"
