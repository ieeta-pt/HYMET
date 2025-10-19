import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "bin" / "hymet"
SIM_MUT = ROOT / "testdataset" / "simulate_mutations.py"


def run_cli(*args):
    cmd = [str(CLI), *args]
    subprocess.run(cmd, check=True, cwd=ROOT)


def test_run_dry_run():
    contigs = ROOT / "bench" / "data" / "cami_i_lc" / "contigs.fna"
    outdir = ROOT / "out" / "ci"  # won't be used due to --dry-run
    run_cli(
        "run",
        "--contigs",
        str(contigs),
        "--out",
        str(outdir),
        "--threads",
        "1",
        "--dry-run",
    )


def test_bench_dry_run():
    manifest = ROOT / "bench" / "cami_manifest.tsv"
    run_cli(
        "bench",
        "--manifest",
        str(manifest),
        "--tools",
        "hymet",
        "--max-samples",
        "1",
        "--dry-run",
    )


def test_case_dry_run():
    manifest = ROOT / "case" / "manifest_zymo.tsv"
    run_cli(
        "case",
        "--manifest",
        str(manifest),
        "--dry-run",
    )


def test_ablation_dry_run():
    run_cli(
        "ablation",
        "--sample",
        "zymo_mc",
        "--taxa",
        "1423,562",
        "--levels",
        "0,1",
        "--fasta",
        str(ROOT / "case" / "truth" / "zymo_refs" / "zymo_refs.fna.gz"),
        "--seqmap",
        str(ROOT / "case" / "truth" / "zymo_refs" / "seqid2taxid.tsv"),
        "--dry-run",
    )


def test_truth_build_zymo_dry_run():
    run_cli(
        "truth",
        "build-zymo",
        "--contigs",
        str(ROOT / "bench" / "data" / "cami_i_lc" / "contigs.fna"),
        "--paf",
        str(ROOT / "case" / "truth" / "zymo_mc" / "zymo_mc_vs_refs.paf"),
        "--out-contigs",
        str(ROOT / "out" / "dummy_truth.tsv"),
        "--out-profile",
        str(ROOT / "out" / "dummy_truth.cami.tsv"),
        "--dry-run",
    )


def test_legacy_dry_run():
    run_cli(
        "legacy",
        "--dry-run",
        "--",
        "--input_dir",
        "/tmp/input",
    )


def test_simulate_mutations_deterministic(tmp_path):
    input_fasta = ROOT / "tests" / "data" / "mutation_input.fna"
    expected_fasta = ROOT / "tests" / "data" / "mutation_expected.fna"
    output_fasta = tmp_path / "mutated.fna"
    cmd = [
        sys.executable,
        str(SIM_MUT),
        "--fasta",
        str(input_fasta),
        "--output",
        str(output_fasta),
        "--sub-rate",
        "0.3",
        "--indel-rate",
        "0.05",
        "--max-indel-length",
        "2",
        "--seed",
        "1337",
    ]
    subprocess.run(cmd, check=True, cwd=ROOT)
    assert output_fasta.read_text() == expected_fasta.read_text()
