# HYMET CAMI Benchmark Notes

This file tracks the current status of the CAMI benchmark harness and how to reproduce the runs after the per-run reference cache changes.

## 1. Environment

- HYMET root: `$(pwd)/HYMET`
- Cache root (configurable): `data/downloaded_genomes/cache_bench` (set via `CACHE_ROOT`)
- Threads: default 16, override with `THREADS`
- `REF_FASTA` must point to a shared reference FASTA (we use `bench/refsets/combined_subset.fasta`).
- Candidate selection: harness defaults to species-level deduplication with `CAND_MAX=1500` (override with `SPECIES_DEDUP=0` or a different `CAND_MAX`).
- Cache hygiene: run `python bench/tools/prune_cache.py --max-age-days 30 --max-size-gb 150` to trim stale `data/downloaded_genomes/cache_bench` entries.

`run_hymet_cami.sh` hashes the Mash-selected accession list and stores downloads in `data/downloaded_genomes/cache/<sha1>/`, so subsequent runs reuse the same minimap index and FASTA. For day-to-day usage, invoke the unified CLI instead of the raw scripts:

```bash
# single run
bin/hymet run --contigs /path/to/contigs.fna --out /path/to/output --threads 16

# CAMI benchmark
bin/hymet bench --manifest bench/cami_manifest.tsv --tools hymet,kraken2,centrifuge
```

## 2. Reproducible run recipe

```bash
cd HYMET/bench

# one-time helper: derive a taxonomy table for the shared FASTA
python bench/tools/make_refset_taxonomy.py \
  --fasta bench/refsets/combined_subset.fasta \
  --taxonkit-db taxonomy_files \
  --output data/detailed_taxonomy.tsv

THREADS=16 \
CACHE_ROOT=data/downloaded_genomes/cache_bench \
REF_FASTA=$(pwd)/refsets/combined_subset.fasta \
./run_all_cami.sh
```

- Outputs land in `bench/out/<sample>/<tool>/`. Aggregated metrics are written to:
  - `bench/out/summary_per_tool_per_sample.tsv`
  - `bench/out/leaderboard_by_rank.tsv`
  - `bench/out/contig_accuracy_per_tool.tsv`
  - `bench/out/runtime_memory.tsv`
  - Figures: `bench/out/fig_accuracy_by_rank.png`, `bench/out/fig_f1_by_rank.png`, `bench/out/fig_l1_braycurtis.png`, `bench/out/fig_per_sample_f1_stack.png`, `bench/out/fig_cpu_time_by_tool.png`, `bench/out/fig_peak_memory_by_tool.png` (mirrored under `results/bench/` alongside `fig_contig_accuracy_heatmap.png` for repo-level access)
- Cache keys are logged for each HYMET invocation; omit `FORCE_DOWNLOAD` to reuse them. Remove old entries in `data/downloaded_genomes/cache_bench/` when disk space gets tight.
- MetaPhlAn 4 retries automatically with `--split_reads` and ≤4 threads if the primary run fails, which eliminates the previous Bowtie2 broken pipe. Use `METAPHLAN_OPTS`/`METAPHLAN_THREADS` to override as needed.

`run_all_cami.sh` triggers `aggregate_metrics.py` and `plot/make_figures.py` automatically at the end of a successful run, so no extra commands are required to refresh the TSVs and figures listed above.

## 3. Latest metrics snapshot (current settings)

| Sample              | Tool   | Rank    | F1 (%) | Notes |
|---------------------|--------|---------|-------:|-------|
| `cami_i_hc`         | HYMET  | species | 52.94  | 9 TP / 11 FP / 5 FN — tougher filters trim long-tail false positives but some bleed-through remains. |
| `cami_i_lc`         | HYMET  | species | 60.00  | 6 TP / 6 FP / 2 FN — compact reference keeps both precision and recall balanced. |
| `cami_i_mc`         | HYMET  | species | 55.56  | 5 TP / 1 FP / 7 FN — higher precision, but recall is limited by contigs without confident hits. |
| `cami_ii_marine`    | HYMET  | species | 78.26  | 9 TP / 3 FP / 2 FN — high-complexity marine panel benefits most from the filtered candidate list. |
| `cami_ii_mousegut`  | HYMET  | species | 60.87  | 7 TP / 2 FP / 7 FN — balanced hit list after deduplicating candidates. |
| `cami_ii_strainmadness` | HYMET | species | 50.00  | 6 TP / 7 FP / 5 FN — strain crowding still causes near-neighbour swaps at species rank. |
| `cami_sample_0`     | HYMET  | species | 63.64  | 14 TP / 9 FP / 7 FN — hard coverage cutoffs suppress most spurious matches from the mega-mix. |

Refer to `bench/out/summary_per_tool_per_sample.tsv` for per-rank detail and `bench/out/leaderboard_by_rank.tsv` for rank-wise means (7 CAMI samples).

This configuration used the following HYMET parameters:

```
CAND_MAX=200 SPECIES_DEDUP=1 HYMET_REL_COV_THRESHOLD=0.2 HYMET_ABS_COV_THRESHOLD=0.02 \
HYMET_TAXID_MIN_SUPPORT=1 HYMET_TAXID_MIN_WEIGHT=0
```

Candidate logs (`out/<sample>/hymet/logs/candidate_limit.log`) confirm the pruning: `cami_sample_0` keeps 200 of 37,556 Mash hits, while smaller panels such as `cami_i_lc` retain their full 147 deduplicated candidates. Run metadata and resource usage live in `bench/out/runtime_memory.tsv`.

### Figure interpretations
- **fig_f1_by_rank.png** – HYMET averages ~93% F1 through order, ~88% at family, and settles near 60% at species across the seven CAMI samples.
- **fig_accuracy_by_rank.png / fig_contig_accuracy_heatmap.png** – Contig-level accuracy stays above 90% down to genus and remains in the mid-80% range at species despite the stringent filters.
- **fig_l1_braycurtis.png** – Mean abundance error trends downward with rank (~28 pct-pts at class/order, ~52 pct-pts at species) reflecting the precision/recall trade-off.
- **fig_per_sample_f1_stack.png** – Highlights which CAMI subsets gain the most from the tightened filtering (marine and mock communities show the biggest lift).
- **fig_cpu_time_by_tool.png / fig_peak_memory_by_tool.png** – HYMET completes in roughly 1–3.5 wall minutes per sample; peak RSS ranges from ~1.1 GB on CAMI I to ~17.4 GB on CAMI II marine/strainmadness.

### Tool-specific notes
- **Kraken2/Bracken** – The rebuilt Bracken database (`database150mers.kmer_distrib`) now feeds the evaluation, lifting mean species-level F1 to ~55% (precision 69%, recall 47%).
- **MetaPhlAn4** – Lineage conversion now consumes MetaPhlAn’s taxid hierarchy directly, producing populated CAMI profiles and ~75% mean species F1 across the seven CAMI samples.
- **sourmash_gather** – Profiles are now rolled up across the taxonomy so intermediate ranks appear in the tables/plots. Several CAMI samples still have zero F1 below phylum simply because gather reports no deeper hits; the zeros now reflect the underlying predictions rather than missing rows.
- **Centrifuge & Ganon2** – Both tools complete successfully, but high abundance error remains without additional filtering; consult `summary_per_tool_per_sample.tsv` for per-rank deltas.
- **BASTA (DIAMOND backend)** – Now executed against a UniProt Swiss-Prot subset converted to DIAMOND; run times range from ~7 s on the CAMI I panels to ~5.5 min for `cami_sample_0`. The converter emits CAMI-compliant profiles and contig assignments, delivering high precision at upper ranks (superkingdom/phylum ≥100% on most samples) while keeping the overall benchmark turnaround on par with the other profilers. Species-level recall remains bounded by protein coverage but is captured in the updated summary tables.
- **PhaBOX** – Integrated via the bench runner to re-label contigs, execute `phabox2`, and parse `phagcn_prediction.tsv`. Outputs feed CAMI evaluation (profile + contig assignments). Runtime depends on the PhaBOX database size, but using the CLI natively keeps wall-clock comparable to other profilers when DIAMOND acceleration is available.
- **Viral-only tools** – PhaBOX (and any other viral classifiers) still score near-zero on the bacterial CAMI datasets; with the roll-up in place, those zeros are genuine false positives/negatives rather than a converter gap.
