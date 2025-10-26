# HYMET Case Study & Database Ablation Toolkit

This directory documents the tooling used to run HYMET on real-world contig datasets and to measure graceful degradation when reference sequences are progressively removed.

## Layout

```
case/
├── manifest.tsv          # Case-study samples (sample_id, contigs_fa, truth_contigs_tsv, truth_profile_tsv, expected_taxa, citation)
├── fetch_case_data.sh    # Helper to download the referenced contig FASTA files
├── run_case.sh           # Execute HYMET (and optional MetaPhlAn sanity checks)
├── run_ablation.sh       # Run HYMET across reference ablation levels
├── ablate_db.py          # Produce ablated FASTA references
├── lib/
│   ├── common.sh         # Shared helpers (path resolution, logging)
│   └── measure.sh        # Runtime/memory tracker (wrapper over /usr/bin/time -v)
└── .gitignore            # Ignores heavy outputs (out/, tmp/, ablation/)
```

## Quick Start

1. **Populate the manifest** (`manifest.tsv`)

   Provide the contig FASTA and optional metadata. Supply ground-truth paths when available to enable F1/misassignment metrics during ablations:
   ```tsv
   sample_id	contigs_fa	truth_contigs_tsv	truth_profile_tsv	expected_taxa	citation
   gut_case	/data/case/gut_assembly.fna			"Bacteroides fragilis;Escherichia coli"	"Doe et al., Microbiome (2024)"
   ```

2. **Fetch data (optional helper)**

   ```bash
   cd HYMET/case
   ./fetch_case_data.sh              # downloads both samples to /data/case
   ./fetch_case_data.sh zymo_mc      # download only the mock community
   ./fetch_case_data.sh --dest /tmp/case gut_case
   ```

3. **Run the case study**

   ```bash
   cd HYMET/case
   THREADS=16 ./run_case.sh --sanity-metaphlan
   # or explicitly pin the suite:
   THREADS=16 ./run_case.sh --suite reviewer_panel --scenario cases
   # run every bundled suite (canonical + gut + zymo) end-to-end:
   ./run_cases_full.sh --threads 16
   ```

   By default every invocation publishes into `results/cases/<suite>/run_<timestamp>/`:
   - `raw/<sample>/hymet/` – `profile.cami.tsv`, `classified_sequences.tsv`, logs, and per-sample metadata.
   - `raw/<sample>/top_taxa.tsv` – Top‑N taxa summary.
   - `raw/<sample>/metaphlan/` – Optional sanity profiles + comparison tables.
   - `raw/runtime_memory.tsv` – Wall/user/sys time plus RSS for each stage.
   - `tables/top_taxa_summary.tsv`, `tables/runtime_memory.tsv`, `tables/metaphlan_metrics.tsv` (when applicable).
   - `figures/` – Runtime/memory bar charts and taxa panels regenerated via `case/plot_case.py`.
   Use `--out /custom/path` (or `--no-publish`) if you need the legacy `case/out` workspace for ad‑hoc experiments.

4. **Ablation experiment**

   Remove increasing fractions of dominant taxa from the shared FASTA and re-run HYMET:
   ```bash
   ./run_ablation.sh \
     --taxa 1239,976 \       # TaxIDs to ablate (e.g., Bacillota, Bacteroidota)
     --levels 0,0.25,0.5,0.75,1.0 \
     --threads 16
   ```

   Each run now lands in `results/ablation/<suite>/run_<timestamp>/`. Within `raw/` you will find the ablated FASTA sets, HYMET outputs per `level_xxx`, per-level evaluation reports, and a unified `runtime_memory.tsv`. The `tables/` folder holds copies of `ablation_summary.tsv`, `ablation_eval_summary.tsv`, and runtime stats, while `figures/` contains the fallback curves produced by `case/plot_ablation.py`. Use `--out` to bypass publishing and reuse the legacy `case/ablation/` scratch space when needed.

## Outputs

- `results/cases/<suite>/run_<timestamp>/raw/<sample>/hymet/` – HYMET predictions per case-study sample.
- `results/cases/<suite>/run_<timestamp>/raw/runtime_memory.tsv` – Wall/user/sys time plus memory and I/O (including MetaPhlAn stages when enabled).
- `results/cases/<suite>/run_<timestamp>/tables/top_taxa_summary.tsv` – Aggregated top‑taxa snapshot for all samples.
- `results/cases/<suite>/run_<timestamp>/figures/` – Runtime/memory plot, taxa panels, and overlap heatmap.
- `results/ablation/<suite>/run_<timestamp>/raw/refsets/*.fasta` – Ablated references (percentage encoded in the filename).
- `results/ablation/<suite>/run_<timestamp>/raw/<sample>/level_<label>/hymet/` – HYMET outputs and optional `eval/` metrics per ablation level.
- `results/ablation/<suite>/run_<timestamp>/tables/ablation_summary.tsv` (plus `ablation_eval_summary.tsv` when truth is provided).
- `results/ablation/<suite>/run_<timestamp>/figures/` – Rank fallback curve, stacked assignment chart, and optional F1-by-rank plot.

## Notes

- The scripts reuse the CAMI harness (`bench/run_hymet.sh`) to guarantee identical classifier behaviour.
- MetaPhlAn sanity checks are optional; set `--sanity-metaphlan` and ensure its database is installed.
- Ablation temporarily replaces `HYMET/data/downloaded_genomes/combined_genomes.fasta`; the `run_ablation.sh` script backs up and restores the original file automatically.
- All heavy artefacts are ignored by git (`out/`, `ablation/`, `tmp/`).
- `run_ablation.sh` forwards the `--seed` parameter to `ablate_db.py` (default 1337) so sequence removal is reproducible.
- `run_cases_full.sh` can execute the canonical, gut, and zymo suites (or a subset via `--suite`) and automatically re-run both `case/plot_case.py` and `bench/plot/make_figures.py --tables results/.../tables` for each published run.

These outputs support documenting real-data performance and robustness under incomplete reference databases.
