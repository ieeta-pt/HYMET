# HYMET (Hybrid Metagenomic Tool)

[![Conda Version](https://anaconda.org/bioconda/hymet/badges/version.svg)](https://anaconda.org/bioconda/hymet)
[![Downloads](https://anaconda.org/bioconda/hymet/badges/downloads.svg)](https://anaconda.org/bioconda/hymet)
[![Platforms](https://anaconda.org/bioconda/hymet/badges/platforms.svg)](https://anaconda.org/bioconda/hymet)
[![License](https://anaconda.org/bioconda/hymet/badges/license.svg)](https://anaconda.org/bioconda/hymet)
[![Latest Release Date](https://anaconda.org/bioconda/hymet/badges/latest_release_date.svg)](https://anaconda.org/bioconda/hymet)

HYMET performs contig-level metagenomic classification by combining Mash-based candidate selection, minimap2 alignment, and a weighted-LCA resolver. The repository includes the classifier, the CAMI benchmark harness, real-data case-study tooling, and auxiliary scripts.

## Feature Snapshot

- **Candidate filtering** – Mash containment scores cap the number of references passed to minimap2 (200 by default).
- **CLI workflows** – `bin/hymet` provides `run`, `bench`, `case`, `ablation`, `truth build-zymo`, and `legacy` subcommands with consistent metadata outputs.
- **Benchmark automation** – The CAMI harness produces evaluation tables, runtime logs, and figures from a single driver script.
- **Case-study tooling** – Dedicated scripts execute MGnify and Zymo contig workflows and perform reference ablation experiments.
- **Deployment options** – Install via Bioconda, Docker/Singularity images, or a source checkout with the supplied environment file.

<p align="center">
  <img src="results/cami/canonical/run_20251024T171353Z/raw/fig_f1_by_rank_lines.png" alt="HYMET F1 by taxonomic rank" width="49%"> 
  <img src="results/cases/canonical/run_20251018T220513Z/figures/fig_case_top_taxa_panels.png" alt="Case study top taxa" width="44%">
</p>

<!-- Detailed benchmark figures and discussion live in bench/results_summary.md -->

## What’s Included

| Directory | Purpose |
|-----------|---------|
| `bench/` | CAMI benchmark harness, database builders, evaluation, and plotting scripts. |
| `case/` | Real-data case study runner plus reference ablation tooling. |
| `workflows/` | High-level runners (e.g., CAMI suite) that stage artefacts under `results/<scenario>/<suite>/`. |
| `bin/` | Python CLI entry points (preferred interface for new workflows). |
| `scripts/` | Legacy Perl/Bash helpers retained for reproducibility (`main.pl`, `config.pl`, etc.). |
| `testdataset/` | Utilities to assemble small synthetic evaluation sets. |
| `data/`, `taxonomy_files/` | Expected locations for downloaded references and taxonomy dumps. |

## Quick Start Commands

```bash
# Single-sample classification
your/env/bin/hymet run   --contigs /path/to/sample.fna   --out out/sample   --threads 16

# CAMI benchmark (HYMET + baselines)
your/env/bin/hymet bench   --manifest bench/cami_manifest.tsv   --tools hymet,kraken2,centrifuge,ganon2,sourmash_gather,metaphlan4   --threads 16

# Case-study bundle (MGnify gut + Zymo mock community)
your/env/bin/hymet case   --manifest case/manifest_zymo.tsv   --threads 8

# Reference ablation experiment
your/env/bin/hymet ablation   --sample zymo_mc   --taxa 1423,562   --levels 0,0.5,1.0   --threads 4

# Refresh supplementary tables and figures
your/env/bin/hymet artifacts
```

`bin/hymet` auto-detects `HYMET_ROOT`. Export it explicitly (`export HYMET_ROOT=/path/to/HYMET`) if you prefer running from arbitrary directories. The legacy Perl entry point remains available as `bin/hymet legacy -- …`.

## Using HYMET Beyond Benchmarks

HYMET’s `run` subcommand is the supported way to classify your own assemblies or read sets; the CAMI harness is just a bundled example. A typical ad-hoc run looks like this:

```bash
export HYMET_ROOT=/path/to/HYMET
${HYMET_ENV}/bin/hymet run \
  --contigs my_assembly.fna \
  --out results/my_assembly \
  --threads 32 \
  --cand-max 500 \
  --species-dedup
```

1. **Prepare inputs** – Provide a contig FASTA via `--contigs` or a read FASTA/FASTQ via `--reads`. Place the pre-built Mash sketches under `HYMET/data/` (see *Preparing Data* below) and point `taxonomy_files/` at a fresh NCBI dump. HYMET will download any missing reference genomes into `CACHE_ROOT` on demand.
2. **Launch classification** – `hymet run` stages the sample under `OUTDIR`, copies your FASTA/FASTQ into `OUTDIR/input/`, screens candidates with Mash, downloads the needed references (or reuses a cache hit), aligns with minimap2, and resolves calls with the weighted LCA resolver. Tweak `--cand-max`, `--species-dedup`, `--threads`, `--cache-root`, or `--assembly-summary-dir` to fit your hardware and naming conventions. Passing `--keep-work` retains intermediates under `OUTDIR/work/` for debugging.
3. **Consume outputs** – Every run writes:
   - `OUTDIR/classified_sequences.tsv` – one row per contig/read with lineage, rank, TaxID, and confidence (compatible with downstream Krona/plots).
   - `OUTDIR/hymet.sample_0.cami.tsv` – CAMI-format profile built from the classification table (rename or symlink as desired for multi-sample batches).
   - `OUTDIR/metadata.json` – reproducibility snapshot (HYMET commit, sketch checksums, cache key, tool versions, tunables).
   - `OUTDIR/logs/` plus `OUTDIR/work/` – diagnostics and reusable intermediates. Set `KEEP_HYMET_WORK=1` if you want to inspect minimap2 PAFs or the selected genome list later.

For throughput runs, iterate over samples in a simple shell loop or a workflow manager, changing only `--contigs/--reads` and `--out` per sample; the cache key (hash of `selected_genomes.txt`) lets multiple runs share downloaded references automatically.

### Reproducing CAMI suites

- Follow the detailed playbook in `docs/reproducibility.md` for the original manuscript run. The published artefacts now live under `results/cami/canonical/run_<timestamp>/` (raw outputs, tables, figures, metadata).
- Use `workflows/run_cami_suite.sh` to stage new CAMI experiments. Every invocation creates `results/<scenario>/<suite>/run_<timestamp>/` and fills it with `raw/`, `tables/`, `figures/`, and `metadata.json`. For example:

  ```bash
  THREADS=8 CACHE_ROOT=data/downloaded_genomes/cache_bench \
  workflows/run_cami_suite.sh \
    --scenario cami \
    --suite contig_full \
    --contig-tools hymet,kraken2,centrifuge,ganon2,viwrap,tama,squeezemeta,megapath_nano \
    --read-tools hymet_reads
  ```

  All artefacts for that run appear in `results/cami/contig_full/run_<timestamp>/`; nothing under `bench/out/` or `results/cami/canonical/` is touched.

## Installation Options

| Method | Command | Notes |
|--------|---------|-------|
| **Bioconda / mamba** | `mamba install -c bioconda hymet` | Installs the CLI and dependencies into the active environment. |
| **Docker** | `docker build -t hymet .`<br>`docker run --rm -it hymet hymet --help` | Image bundles the benchmark harness; bind data/cache directories as needed. |
| **Singularity / Apptainer** | `apptainer build hymet.sif Singularity.def`<br>`apptainer exec hymet.sif hymet --help` | Mirrors the Docker build for HPC clusters. |
| **Source checkout** | `git clone https://github.com/ieeta-pt/HYMET.git`<br>`cd HYMET`<br>`mamba env create -f environment.yml` | Recommended for development; activate the environment before using `bin/hymet`. For exact pins, use `environment.lock.yml`. |

## Preparing Data

1. **References & taxonomy** – Download the Mash sketches from the Zenodo archive (`10.5281/zenodo.17428354`) with:
   ```bash
   tools/fetch_sketches.sh    # defaults to the Zenodo record + checksum verification
   tools/verify_sketches.sh   # optional: confirm local files match the archive
   ```
   Place NCBI taxonomy dumps under `HYMET/taxonomy_files/`. Builders in `bench/db/` derive tool-specific indices on demand.
2. **CAMI subsets** – Use `bench/fetch_cami.sh` (supports `--dry-run`) to download the contigs listed in `bench/cami_manifest.tsv`.
3. **Case-study contigs** – `case/fetch_case_data.sh` retrieves the MGnify gut assembly and Zymo mock community.
4. **Truth tables** – CAMI truth lives under `bench/data/`; case-study truth files (including curated Zymo panels) live under `case/truth/`.

## Outputs at a Glance

- Canonical CAMI artefacts: `results/cami/canonical/run_<timestamp>/` (raw benchmark outputs, summary tables, figures, metadata).
- Reviewer suites: `results/cami/<suite>/run_<timestamp>/` (one folder per run). Raw per-tool outputs are grouped by mode; derived tables/figures live alongside metadata for immediate inspection.
- Case studies and ablations follow the same pattern under `results/cases/…/run_<timestamp>/` and `results/ablation/…/run_<timestamp>/`.
- Bench runs no longer write to `bench/out/` (all artefacts land in `results/<scenario>/<suite>/run_<timestamp>/…`). Case workflows follow the same convention: `case/run_case.sh` (and `bin/hymet case`) publish into `results/cases/<suite>/run_<timestamp>/…` unless you override the destination with `--out`.
- Use `python bench/plot/make_figures.py --bench-root <raw_dir> --outdir <target>` to regenerate figures for any run.

## Documentation & Reporting

- CAMI harness details: `bench/README.md`, latest metrics: `bench/results_summary.md`.
- Case-study workflows: `case/README.md`, results recap: `case/results_summary.md`.
- Reproducibility playbook: `docs/reproducibility.md`.

## Repository Layout

```
HYMET/
├── bin/                 # CLI entry points (Python)
├── bench/               # CAMI harness (runners, builders, plots)
├── case/                # Real-data case study + ablation toolkit
├── docs/                # Additional guides
├── results/             # Canonical artefacts (cami, cases, ablation, …)
├── workflows/           # Repro runners that populate results/<scenario>/<suite>/
├── scripts/             # Legacy helpers (Perl/Bash)
├── testdataset/         # Synthetic dataset utilities
└── data/, taxonomy_files/, …  # Downloaded references and taxonomy dumps
```

The maintained workflow is through the Python CLI. Legacy scripts (`config.pl`, `main.pl`, `scripts/*.sh`) are retained for historical pipelines but no longer required for fresh runs.

## Support & Citation

- Open issues and feature requests on the GitHub tracker.
- Cite HYMET using `CITATION.cff` in this repository.
