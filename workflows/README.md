# HYMET Workflows

The `workflows/` directory contains optional orchestration scripts that build on the core HYMET utilities (`bin/hymet`, `bench/run_all_cami.sh`, etc.) without modifying their layouts. They provide reproducible entry points for manuscript-driven experiments while keeping automation separate from the benchmarking harness.

## Contents

| Path | Description |
|------|-------------|
| `run_cami_suite.sh` | Multi-mode CAMI benchmark runner (contigs + synthetic reads). |
| `config/cami_suite.cfg` | Default scenario/suite and tool lists. |
| `docs/cami_suite_protocol.md` | Step-by-step reproduction protocol with provenance checklist. |

## Usage

```bash
cd HYMET
THREADS=16 CACHE_ROOT=data/downloaded_genomes/cache_bench \
workflows/run_cami_suite.sh \
  --scenario cami \
  --suite contig_full
```
Each invocation creates `results/<scenario>/<suite>/run_<timestamp>/` with subfolders for raw outputs, tables, figures, and metadata. The runner sets `BENCH_OUT_ROOT` automatically, so the core harness never writes into `bench/out/`. See `docs/cami_suite_protocol.md` for the full reproducibility checklist and packaging instructions.
