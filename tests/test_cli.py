import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "bin" / "hymet"
SIM_MUT = ROOT / "testdataset" / "simulate_mutations.py"
PRELOAD_CACHE = ROOT / "case" / "tools" / "preload_cache_from_fasta.py"
CONTIGS_TO_READS = ROOT / "bench" / "tools" / "contigs_to_reads.py"


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


def test_preload_cache_from_fasta(tmp_path):
    cache_dir = tmp_path / "cache"
    fasta = tmp_path / "refs.fna"
    fasta.write_text(
        ">seq1\nACGTACGT\n>seq2\nGGNNA\n",
        encoding="utf-8",
    )
    seqmap = tmp_path / "seqmap.tsv"
    seqmap.write_text("seq1\t123\nseq2\t456\n", encoding="utf-8")

    cache_dir.mkdir()
    (cache_dir / "reference.mmi").write_text("placeholder", encoding="utf-8")

    cmd = [
        sys.executable,
        str(PRELOAD_CACHE),
        "--cache-dir",
        str(cache_dir),
        "--fasta",
        str(fasta),
        "--seqmap",
        str(seqmap),
        "--taxid-prefix",
        "Test",
    ]
    subprocess.run(cmd, check=True, cwd=ROOT)

    combined = (cache_dir / "combined_genomes.fasta").read_text(encoding="utf-8")
    assert ">seq1" in combined and "ACGTACGT" in combined
    assert ">seq2" in combined and "GGNNA" in combined

    taxonomy_lines = (cache_dir / "detailed_taxonomy.tsv").read_text(encoding="utf-8").strip().splitlines()
    assert taxonomy_lines[0] == "GCF\tTaxID\tIdentifiers"
    assert len(taxonomy_lines) == 3
    assert "Test_123" in taxonomy_lines[1]
    assert not (cache_dir / "reference.mmi").exists()


def test_contigs_to_reads(tmp_path):
    contigs = tmp_path / "contigs.fna"
    contigs.write_text(
        ">contig1\nACGTACGT\n>contig2\nNNAC\n",
        encoding="utf-8",
    )
    out_fastq = tmp_path / "reads.fastq"
    cmd = [
        sys.executable,
        str(CONTIGS_TO_READS),
        "--contigs",
        str(contigs),
        "--out",
        str(out_fastq),
        "--chunk-size",
        "3",
        "--min-chunk",
        "2",
    ]
    subprocess.run(cmd, check=True, cwd=ROOT)

    lines = out_fastq.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) % 4 == 0
    sequences = lines[1::4]
    qualities = lines[3::4]
    assert sequences, "Expected at least one synthetic read"
    for seq, qual in zip(sequences, qualities):
        assert len(seq) == len(qual)
        assert set(qual) == {"I"}
        assert all(base in {"A", "C", "G", "T"} for base in seq)
