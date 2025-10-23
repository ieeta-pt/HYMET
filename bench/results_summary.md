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

## 3. Latest results (aggregated across CAMI samples)

This section summarises the current benchmark after refreshing all aggregates and figures from `bench/out/`. See the referenced figures and TSVs for complete details.

Species-rank F1 (mean across samples):
- TAMA ≈ 79.5%
- MetaPhlAn4 ≈ 75.5%
- Kraken2 ≈ 69.4%
- HYMET ≈ 62.6%
- MegaPath‑Nano ≈ 45.6%
- sourmash gather ≈ 42.9%
- SnakeMAGs ≈ 35.7%
- phyloFlash ≈ 27.6%
- BASTA ≈ 22.4%, Centrifuge ≈ 22.6%, Ganon2 ≈ 23.8%
- CAMITAX, PhaBOX, SqueezeMeta, ViWrap: ≈ 0.0% at species on these bacterial CAMI panels (no/little species‑level signal under the current converters/run‑paths).

Higher ranks (mean F1):
- Genus: HYMET (~89.5%) leads, with MegaPath‑Nano (~86%) and TAMA (~83.5%) close; Kraken2 (~78.7%) and MetaPhlAn4 (~74.7%) follow.
- Family: HYMET (~98%) tops, MetaPhlAn4 (~92%), MegaPath‑Nano (~92.8%), TAMA (~90%), Kraken2 (~87%).
- Superkingdom: Perfect 100% appears for several tools (e.g. BASTA, Ganon2, phyloFlash) as expected at this coarse level.

Contig‑level accuracy (species):
- HYMET ≈ 86.1% average; Kraken2 ≈ 73.9%.
- Other tools either lack contig assignments in our converters or score near zero at species for these datasets.

Abundance error trends (L1 total variation and Bray–Curtis):
- Error increases toward the species rank for all profilers; line plots show a consistent monotonic trend across ranks.
- The top profilers at genus/family keep substantially lower error; the relative ordering broadly mirrors F1.

Runtime and peak memory (means across “run” stages):
- Fast/light: CAMITAX (~0.0 min), ganon2 (~0.13 min, ~0.18 GB), sourmash gather (~0.17 min, ~0.81 GB).
- Mid‑pack: Kraken2 (~0.42 min, ~11.15 GB), HYMET (~1.96 min, ~17.36 GB), phyloFlash (~0.95 min, ~4.34 GB), Centrifuge (~1.13 min, ~0.32 GB), MegaPath‑Nano (~0.79 min, ~11.49 GB), BASTA (~4.36 min, ~2.28 GB), MetaPhlAn4 (~4.63 min, ~18.76 GB).
- Heavier: SnakeMAGs (~17.4 min, ~28.9 GB), ViWrap (~87.9 min, ~18.4 GB).

### F1 by rank

![F1 by rank (lines)](../results/bench/fig_f1_by_rank_lines.png)

- End labels show per‑tool trajectories across ranks; the top‑3 are highlighted.
- Relative ordering remains stable across ranks; HYMET’s advantage increases toward family/genus.

### Abundance error (L1 & Bray–Curtis)

![Abundance error (lines)](../results/bench/fig_l1_braycurtis_lines.png)

- Lines emphasise the monotonic increase in error from superkingdom → species.
- Ordering mirrors the F1 plots, underscoring precision/recall trade‑offs at deep ranks.

### Contig accuracy by rank

![Contig accuracy (lines)](../results/bench/fig_accuracy_by_rank_lines.png)

- Line view of contig accuracy reinforces HYMET’s gap at genus/species.

### Per‑sample stacked F1 (species)

![Per‑sample F1 stack](../results/bench/fig_per_sample_f1_stack.png)

- Shows dataset sensitivity: marine/mock communities benefit most from richer references and consensus schemes.

### CPU time by tool

![CPU time by tool](../results/bench/fig_cpu_time_by_tool.png)

- Fast/light: CAMITAX, ganon2, sourmash gather.
- Mid‑range: Kraken2, HYMET, phyloFlash, Centrifuge, MegaPath‑Nano, BASTA, MetaPhlAn4.
- Heavy: SnakeMAGs, ViWrap.

### Peak memory by tool

![Peak memory by tool](../results/bench/fig_peak_memory_by_tool.png)

- Peak RSS spans sub‑GB (ganon2) to ~29 GB (SnakeMAGs); HYMET ~17 GB, Kraken2/MetaPhlAn4 ~11–19 GB.

Tables (CSV/TSV):
- Per‑sample, per‑rank metrics: `bench/out/summary_per_tool_per_sample.tsv`
- Rank‑wise leaderboard (means): `bench/out/leaderboard_by_rank.tsv`
- Contig accuracy per rank/tool: `bench/out/contig_accuracy_per_tool.tsv`
- Runtime/memory per stage: `bench/out/runtime_memory.tsv`

This configuration used the following HYMET parameters:

```
CAND_MAX=200 SPECIES_DEDUP=1 HYMET_REL_COV_THRESHOLD=0.2 HYMET_ABS_COV_THRESHOLD=0.02 \
HYMET_TAXID_MIN_SUPPORT=1 HYMET_TAXID_MIN_WEIGHT=0
```

Candidate logs (`out/<sample>/hymet/logs/candidate_limit.log`) confirm the pruning: `cami_sample_0` keeps 200 of 37,556 Mash hits, while smaller panels such as `cami_i_lc` retain their full 147 deduplicated candidates. Run metadata and resource usage live in `bench/out/runtime_memory.tsv`.

### Figure interpretations
See the discussion sections following each figure above.

### Tool-specific notes
- **Kraken2/Bracken** – The rebuilt Bracken database (`database150mers.kmer_distrib`) now feeds the evaluation, lifting mean species-level F1 to ~55% (precision 69%, recall 47%).
- **MetaPhlAn4** – Lineage conversion now consumes MetaPhlAn’s taxid hierarchy directly, producing populated CAMI profiles and ~75% mean species F1 across the seven CAMI samples.
- **sourmash_gather** – Profiles are now rolled up across the taxonomy so intermediate ranks appear in the tables/plots. Several CAMI samples still have zero F1 below phylum simply because gather reports no deeper hits; the zeros now reflect the underlying predictions rather than missing rows.
- **Centrifuge & Ganon2** – Both tools complete successfully, but high abundance error remains without additional filtering; consult `summary_per_tool_per_sample.tsv` for per-rank deltas.
- **BASTA (DIAMOND backend)** – Now executed against a UniProt Swiss-Prot subset converted to DIAMOND; run times range from ~7 s on the CAMI I panels to ~5.5 min for `cami_sample_0`. The converter emits CAMI-compliant profiles and contig assignments, delivering high precision at upper ranks (superkingdom/phylum ≥100% on most samples) while keeping the overall benchmark turnaround on par with the other profilers. Species-level recall remains bounded by protein coverage but is captured in the updated summary tables.
- **PhaBOX** – Integrated via the bench runner to re-label contigs, execute `phabox2`, and parse `phagcn_prediction.tsv`. Outputs feed CAMI evaluation (profile + contig assignments). Runtime depends on the PhaBOX database size, but using the CLI natively keeps wall-clock comparable to other profilers when DIAMOND acceleration is available.
- **Viral-only tools** – PhaBOX (and any other viral classifiers) still score near-zero on the bacterial CAMI datasets; with the roll-up in place, those zeros are genuine false positives/negatives rather than a converter gap.
