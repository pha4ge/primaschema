import subprocess

from pathlib import Path

import pytest

import primaschema.lib as lib


data_dir = Path("test/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def test_cli_version():
    run("primaschema --version")


def test_hash_ref():
    assert (
        lib.hash_ref(
            "test/data/primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0/reference.fasta"
        )
        == "primaschema:b1acd7163146bf17"
    )


def test_cli_hash_ref():
    run_cmd = run(
        "primaschema hash-ref primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0/reference.fasta"
    )
    assert "primaschema:b1acd7163146bf17" in run_cmd.stdout


def test_cli_hash_primer_bed():
    run_cmd = run(
        "primaschema hash-bed primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    assert "primaschema:3ef3e7bb23008684" in run_cmd.stdout


def test_cli_scheme_bed():
    run_cmd = run(
        "primaschema hash-bed primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/scheme.bed"
    )
    assert "primaschema:3ef3e7bb23008684" in run_cmd.stdout


def test_artic_v41_scheme_hash_matches_primer_hash():
    scheme_bed_hash = lib.hash_scheme_bed(
        "test/data/primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/scheme.bed",
        "test/data/primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/reference.fasta",
    )
    primer_bed_hash = lib.hash_primer_bed(
        "test/data/primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    assert scheme_bed_hash == primer_bed_hash


def test_valid_eden_v1():
    lib.validate(
        data_dir / "primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0",
    )
    lib.validate(
        data_dir / "primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0",
        full=True,
    )


def test_valid_artic_v41():
    lib.validate(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0",
    )


def test_checksum_case_normalisation():
    assert lib.hash_bed(
        data_dir / "primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0/primer.bed"
    ) == lib.hash_bed(data_dir / "different-case/eden.modified.primer.bed")


def test_cli_valid_recursive():
    run("primaschema validate --recursive primer-schemes")


def test_valid_rebuild():
    lib.validate(
        data_dir / "primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0",
        rebuild=True,
    )


def test_hash_bed():
    lib.hash_bed(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    lib.hash_bed(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/scheme.bed"
    )


def test_build():
    run("primaschema build primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0")
    run("rm -rf artic-v4.1.0")


def test_build_recursive():
    lib.build(data_dir / "primer-schemes", recursive=True)
    run("rm -rf built")


def test_build_manifest():
    lib.build_manifest(root_dir=data_dir / "primer-schemes")
    run("rm -rf built index.json", cwd="./")


def test_primer_bed_to_scheme_bed():
    scheme_bed_path = (
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/scheme.bed"
    )
    primer_bed_path = (
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    bed_str = lib.convert_primer_bed_to_scheme_bed(bed_path=primer_bed_path)
    with open(scheme_bed_path) as fh:
        expected_bed_str = fh.read()
    assert bed_str == expected_bed_str


def test_scheme_bed_to_primer_bed():
    scheme_bed_path = (
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/scheme.bed"
    )
    primer_bed_path = (
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    reference_path = (
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/reference.fasta"
    )
    bed_str = lib.convert_scheme_bed_to_primer_bed(
        bed_path=scheme_bed_path, fasta_path=reference_path
    )
    with open(primer_bed_path) as fh:
        expected_bed_str = fh.read()
    assert bed_str == expected_bed_str


def test_diff():
    run_cmd = run(
        "primaschema diff primer-schemes/schemes/sars-cov-2/midnight/1200/v1.0.0/primer.bed primer-schemes/schemes/sars-cov-2/midnight/1200/v2.0.0/primer.bed"
    )
    assert (
        """SARS-CoV-2_28_LEFT_2""" in run_cmd.stdout.strip()
        and len(run_cmd.stdout.strip().split("\n")) == 2
    )


def test_calculate_intervals():
    all_intervals = lib.amplicon_intervals(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    assert "MN908947.3" in all_intervals
    intervals = all_intervals["MN908947.3"]
    assert "SARS-CoV-2_99" in intervals
    assert intervals["SARS-CoV-2_99"] == (29452, 29854)


def test_print_intervals():
    run_cmd = run(
        "primaschema show-intervals primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    assert """MN908947.3\t29452\t29854\tSARS-CoV-2_99\n""" in run_cmd.stdout


def test_plot_single_ref_chrom_ref():
    lib.plot_primers(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed",
    )
    run("rm -rf primer.html", cwd="./")


def test_plot_many_ref_chroms_ref():
    lib.plot_primers(data_dir / "many-ref-chroms/primer.bed")
    run("rm -rf primer.html", cwd="./")


def test_6to7_many_ref_chroms():
    scheme_bed_path = data_dir / "many-ref-chroms/scheme.bed"
    primer_bed_path = data_dir / "many-ref-chroms/primer.bed"
    reference_path = data_dir / "many-ref-chroms/reference.fasta"
    bed_str = lib.convert_scheme_bed_to_primer_bed(
        bed_path=scheme_bed_path, fasta_path=reference_path
    )
    with open(primer_bed_path) as fh:
        expected_bed_str = fh.read()
    assert bed_str == expected_bed_str


def test_invalid_duplicate_primers():
    with pytest.raises(ValueError):  # Also catches pydantic.ValidationError
        lib.validate(
            data_dir / "broken/duplicated-primers",
        )


def test_invalid_primer_bounds():
    with pytest.raises(ValueError):  # Also catches pydantic.ValidationError
        lib.validate(
            data_dir / "broken/primer-bounds",
        )


def test_invalid_amplicon_tiling():
    with pytest.raises(ValueError):  # Also catches pydantic.ValidationError
        lib.validate(
            data_dir / "broken/non-tiling",
        )


def test_format_primer_bed():
    """Sort BED into maximally compatible output order"""
    assert lib.format_primer_bed(data_dir / "unordered/primer.bed").strip() == (
        """MN908947.3	25	50	SARS-CoV-2_1_LEFT_1	1	+	AACAAACCAACCAACTTTCGATCTC
MN908947.3	408	431	SARS-CoV-2_1_RIGHT_1	1	-	CTTCTACTAAGCCACAAGTGCCA
MN908947.3	324	344	SARS-CoV-2_2_LEFT_1	2	+	TTTACAGGTTCGCGACGTGC
MN908947.3	705	727	SARS-CoV-2_2_RIGHT_1	2	-	ATAAGGATCAGTGCCAAGCTCG"""
    )


def test_invalid_missing_field():
    with pytest.raises(ValueError):  # Also catches pydantic.ValidationError
        lib.validate(data_dir / "broken/info-yml/missing-field")
        lib.validate(data_dir / "broken/info-yml/missing-field", full=True)


def test_invalid_extra_field():
    with pytest.raises(ValueError):  # Also catches pydantic.ValidationError
        lib.validate(data_dir / "broken/info-yml/extra-field")
        lib.validate(data_dir / "broken/info-yml/extra-field", full=True)


def test_subset():
    lib.subset(scheme_dir=data_dir / "many-ref-chroms", chrom="NC_038235.1")
    df = lib.parse_primer_bed("built/primer.bed")
    assert len(df) == 4
    run("rm -rf built", cwd="./")


def test_commented_bed():
    lib.validate(data_dir / "bed-comment")


def test_dev_scheme():
    # lib.validate(data_dir / "dev-scheme")
    lib.validate(data_dir / "dev-scheme", full=True)
