## Systematic Review Tool Integration

This file tracks the porting effort of tools from `/data/Syst_Review` into the HYMET benchmark harness.

### Completed
- **CAMITAX** – Runner added at `run_camitax.sh`; converter lives in `convert/camitax_to_cami.py`; registered inside `run_all_cami.sh`.
- **BASTA** – Runner available at `run_basta.sh`; converter located in `convert/basta_to_cami.py`; supports DIAMOND-backed searches (`BASTA_USE_DIAMOND=1`, `BASTA_DIAMOND_DB=/path/to/db.dmnd`) to accelerate CAMI benchmarking while producing CAMI-format profiles and contig tables.
- **PhaBOX** – Integrated via `run_phabox.sh`; converter `convert/phabox_to_cami.py` parses `phagcn_prediction.tsv`, restores original contig IDs, and emits CAMI profiles/contig tables. Supports DIAMOND by default (`PHABOX_CMD`, `PHABOX_DB_DIR`, `PHABOX_EXTRA_OPTS`).
- **TAMA** – Wrapper at `run_tama.sh` prepares/executes TAMA, converts `abundance_profile*.out` via `convert/tama_to_cami.py`, emits both CAMI profile and `classified_sequences.tsv`, and relies on reproducible params stored under `config/tama_params/`. Synthetic single-end reads are generated from CAMI contigs (`tools/contigs_to_reads.py → data/tama_reads/`) so Kraken/Centrifuge receive FASTQ input. Current params run the Centrifuge+Kraken consensus (CLARK skipped to avoid the 150 GB RAM requirement).
- **PhyloFlash** – Runner at `run_phyloflash.sh` synthesizes rRNA-derived reads from contigs, invokes `phyloFlash.pl`, and converts `*.phyloFlash.NTUfull_abundance.csv` with `convert/phyloflash_to_cami.py`. Requires database path via `PHYLOFLASH_DB_DIR` and optionally reuses `/opt/envs/phyloflash` environment.
- **ViWrap / geNomad** – `run_viwrap.sh` links contigs, executes `genomad end-to-end`, and converts `*_virus_summary.tsv` via `convert/viwrap_to_cami.py`. Expects the geNomad conda env at `/opt/envs/genomad` (override with `VIWRAP_ENV_PREFIX`) and the downloaded database under `/data/ref/viwrap/genomad_db`.

### Pending
- **MegaPath-Nano (MPN)** – Needs GPU-aware environment notes, runner stub, and conversion of abundance outputs.
- **SnakeMAGs** – Provide Snakemake invocation through a stable entrypoint and convert resulting MAG classifications.
- **SqueezeMeta** – Wrap the existing pipeline, translate its `sum_tsv` profile into CAMI format, and surface contig annotations.
- **ViWrap** – Add Viral workflow runner and adapt taxonomic outputs for CAMI metrics.

### Shared Implementation Notes
- Follow `run_camitax.sh` for argument parsing and metadata handling.
- Place converters in `bench/convert/` and register new runners in `run_all_cami.sh`.
- Each runner should emit `profile.cami.tsv` and, when available, `classified_sequences.tsv`.
- Use environment variables to point to pre-built databases rather than hard-coded paths.
