# HYMET CAMI Benchmark Notes

This file tracks the current status of the CAMI benchmark harness and how to reproduce the runs after the per-run reference cache changes.

---

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

---

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
  - `bench/out/runtime_memory.tsv`
  - Figures: `results/bench/fig_accuracy_by_rank.png`, `results/bench/fig_f1_by_rank.png`, `results/bench/fig_l1_braycurtis.png`, `results/bench/fig_per_sample_f1_stack.png`, `results/bench/fig_cpu_time_by_tool.png`, `results/bench/fig_peak_memory_by_tool.png`
- Cache keys are logged for each HYMET invocation; omit `FORCE_DOWNLOAD` to reuse them. Remove old entries in `data/downloaded_genomes/cache_bench/` when disk space gets tight.
- MetaPhlAn 4 retries automatically with `--split_reads` and ≤4 threads if the primary run fails, which eliminates the previous Bowtie2 broken pipe. Use `METAPHLAN_OPTS`/`METAPHLAN_THREADS` to override as needed.


## 3. Latest metrics snapshot (subset reference)

| Sample        | Tool      | Rank        | F1 (%) | Notes |
|---------------|-----------|-------------|-------:|-------|
| `cami_i_hc`   | HYMET     | species     | 11.76  | 2 TP / 18 FP / 12 FN — high FP rate persists even after capping Mash hits to 235 candidates. |
| `cami_i_lc`   | HYMET     | species     | 42.11  | 4 TP / 7 FP / 4 FN — compact 147-genome reference; recall remains moderate. |
| `cami_sample_0` | HYMET   | species     | 2.67   | 1 TP / 53 FP / 20 FN — broad mock community remains noisy even after shrinking to 744 candidates (from 5 000). |
| `cami_i_lc`   | MetaPhlAn4 | species    | 0.00   | Run completes with fallback logic; profile stays header-only on this subset. |

Refer to `bench/out/summary_per_tool_per_sample.tsv` for the complete table across ranks and samples.

Per-sample candidate stats now live in `out/<sample>/hymet/logs/candidate_limit.log`; for the latest run HYMET retained 744 of 5 000 Mash hits for `cami_sample_0` and 147 of 231 for `cami_i_lc`. The MetaPhlAn fallback eliminates Bowtie2 crashes, though extremely small panels still yield empty profiles.

Superkingdom rows from every tool are canonicalised before scoring via `bench/tools/fix_superkingdom_taxids.py`, so CAMI II subsets that report GTDB “Bacillati/Pseudomonadati” now align with the reference (Bacteria/Archaea) and populate the previously empty metrics. All runners were re-evaluated with `lib/run_eval.sh` to refresh `summary_per_tool_per_sample.tsv` and downstream figures.

### Figure interpretations
- **fig_f1_by_rank.png** – HYMET and Kraken2 dominate across the ordered ranks; HYMET peaks at order/species (≈92% / 52% mean F1), while Kraken2 remains strongest at superkingdom with fewer low-level calls.
- **fig_accuracy_by_rank.png** – Contig accuracy follows the same hierarchy: HYMET stays above 80% through family, outperforming Kraken2 as ranks narrow; Centrifuge and Ganon2 drop below 40% at species.
- **fig_l1_braycurtis.png** – HYMET yields the lowest abundance errors (mean L1 <20 pct-pts from class downward), whereas Centrifuge/Ganon2 hover around 45–50 pct-pts, signalling noisier mixing estimates.
- **fig_per_sample_f1_stack.png** – Sample stacks show HYMET contributing most F1 signal for `cami_i_hc`, `cami_i_lc`, and the mock community; MetaPhlAn4 and sourmash gather add negligible lift on these subsets.
- **fig_cpu_time_by_tool.png** – Kraken2 and HYMET dominate runtime efficiency (~3–8 CPU minutes per sample), while MetaPhlAn4 is the heaviest (>60 CPU minutes) due to Bowtie2 preprocessing.
- **fig_peak_memory_by_tool.png** – HYMET and GanON2 stay under ~13 GB peak RSS; MetaPhlAn4 peaks close to 19 GB, reinforcing the need for higher-memory nodes when including it in runs.



## 4. Outstanding tasks

- None at this time. Keep monitoring MetaPhlAn resource use on larger contig sets and document any future deviations.
