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
   ```

   Outputs land in `case/out/<sample>/`:
   - `hymet/profile.cami.tsv`, `hymet/classified_sequences.tsv`
   - `top_taxa.tsv` (Top N taxa by abundance)
   - Optional `metaphlan/profile.tsv`, `metaphlan/comparison.tsv`, and `metaphlan/metrics.tsv` (symmetric KL divergence + Spearman rank correlation)
   - `metadata.json` summarising inputs
   - `out/runtime_memory.tsv` capturing wall time, CPU time, RSS, and I/O

4. **Ablation experiment**

   Remove increasing fractions of dominant taxa from the shared FASTA and re-run HYMET:
   ```bash
   ./run_ablation.sh \
     --taxa 1239,976 \       # TaxIDs to ablate (e.g., Bacillota, Bacteroidota)
     --levels 0,0.25,0.5,0.75,1.0 \
     --threads 16
   ```

   For each level, `ablate_db.py` creates FASTA files under `case/ablation/refsets/`, HYMET runs, and `case/ablation_summary.tsv` records how assignments fall back from species → higher ranks. When truth files are provided, `bench/lib/run_eval.sh` is invoked automatically and per-rank F1/precision/recall plus contig misassignment percentages are appended to `case/ablation_eval_summary.tsv`. Runtime stats are appended to `case/out/runtime_memory.tsv` with stage labels `ablation_<level>` and `ablation_eval_<level>`.

## Outputs

- `case/out/runtime_memory.tsv` – Absolute wall/user/sys time plus memory and I/O for HYMET (and optional MetaPhlAn) case runs.
- `case/out/<sample>/top_taxa.tsv` – Top-N ranked taxa summary.
- `case/out/<sample>/metaphlan/metrics.tsv` – Symmetric KL divergence + Spearman rank correlation between HYMET and MetaPhlAn (when the sanity check is enabled).
- `case/ablation/refsets/*.fasta` – Ablated references (percentage encoded in the filename).
- `case/ablation/<sample>/level_<label>/hymet/` – HYMET outputs for each ablation level (with optional `eval/` when truth paths are supplied).
- `case/ablation_summary.tsv` – Per-level fallback stats (% assignments retained at species/genus/family, etc.).
- `case/ablation_eval_summary.tsv` – Per-rank F1/precision/recall and contig misassignment percentages (if truth data available).
- `case/ablation/figures/` – Rank fallback curve, stacked assignment chart, and optional F1-by-rank plot.

## Notes

- The scripts reuse the CAMI harness (`bench/run_hymet.sh`) to guarantee identical classifier behaviour.
- MetaPhlAn sanity checks are optional; set `--sanity-metaphlan` and ensure its database is installed.
- Ablation temporarily replaces `HYMET/data/downloaded_genomes/combined_genomes.fasta`; the `run_ablation.sh` script backs up and restores the original file automatically.
- All heavy artefacts are ignored by git (`out/`, `ablation/`, `tmp/`).
- `run_ablation.sh` forwards the `--seed` parameter to `ablate_db.py` (default 1337) so sequence removal is reproducible.

These outputs support documenting real-data performance and robustness under incomplete reference databases.
