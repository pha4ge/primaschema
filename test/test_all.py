import shutil
import subprocess
from pathlib import Path

import pytest

import primaschema.cli as create_module
import primaschema.lib as lib
import primaschema.validate as validate_module
from primaschema.schema.info import PrimerScheme
from primaschema.util import sha256_checksum

data_dir = Path("test/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def test_hash_ref():
    assert (
        lib.hash_ref(
            "test/data/primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0/reference.fasta"
        )
        == "primaschema:b1acd7163146bf17"
    )


def test_checksum_case_normalisation():
    assert lib.hash_bed(
        data_dir / "primer-schemes/schemes/sars-cov-2/eden/2500/v1.0.0/primer.bed"
    ) == lib.hash_bed(data_dir / "different-case/eden.modified.primer.bed")


def test_hash_bed():
    lib.hash_bed(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    lib.hash_bed(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/scheme.bed"
    )


def test_build_index(tmp_path: Path):
    src = data_dir / "primer-schemes"
    dest = tmp_path / "primer-schemes"
    shutil.copytree(src, dest)
    lib.build_index(root_dir=dest, out_dir=tmp_path)
    index_path = tmp_path / "index.json"
    assert index_path.exists()


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


def test_calculate_intervals():
    all_intervals = lib.amplicon_intervals(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed"
    )
    assert "MN908947.3" in all_intervals
    intervals = all_intervals["MN908947.3"]
    assert "SARS-CoV-2_99" in intervals
    assert intervals["SARS-CoV-2_99"] == (29452, 29854)


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


def _copy_scheme(tmp_path: Path, rel_path: str) -> Path:
    src = data_dir / rel_path
    dest = tmp_path / rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)
    return dest


def _copy_scheme_to(src_rel_path: str, dest: Path) -> Path:
    src = data_dir / src_rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)
    return dest


def test_validate_autonormalize_primer_bed(tmp_path: Path):
    scheme_dir = _copy_scheme(
        tmp_path,
        "auto-normalisation/test/400/v2.0.0",
    )
    info_path = scheme_dir / "info.json"
    primer_path = scheme_dir / "primer.bed"
    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())

    assert sha256_checksum(primer_path) != primer_scheme.primer_file_sha256

    validate_module.validate(info_path, strict=True, fix=True)

    assert sha256_checksum(primer_path) == primer_scheme.primer_file_sha256


def test_validate_autonormalize_reference_fasta(tmp_path: Path):
    scheme_dir = _copy_scheme(
        tmp_path,
        "auto-normalisation/test/400/v2.0.0",
    )
    info_path = scheme_dir / "info.json"
    reference_path = scheme_dir / "reference.fasta"
    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())

    assert sha256_checksum(reference_path) != primer_scheme.reference_file_sha256

    validate_module.validate(info_path, strict=True, fix=True)

    assert sha256_checksum(reference_path) == primer_scheme.reference_file_sha256


def test_validate_all_aggregates_errors_create_cli(tmp_path: Path):
    bad1 = tmp_path / "bad1" / "999" / "v0.0.1"
    bad2 = tmp_path / "bad2" / "999" / "v0.0.2"
    _copy_scheme_to("auto-normalisation/test/400/v2.0.0", bad1)
    _copy_scheme_to("auto-normalisation/test/400/v2.0.0", bad2)

    with pytest.raises(ValueError) as exc:
        create_module.validate(
            path=tmp_path,
            all=True,
            additional_linkml=False,
            strict=True,
        )

    msg = str(exc.value)
    assert "bad1" in msg
    assert "bad2" in msg


def test_validate_all_aggregates_errors_module(tmp_path: Path):
    bad1 = tmp_path / "bad1" / "999" / "v0.0.1"
    bad2 = tmp_path / "bad2" / "999" / "v0.0.2"
    _copy_scheme_to("auto-normalisation/test/400/v2.0.0", bad1)
    _copy_scheme_to("auto-normalisation/test/400/v2.0.0", bad2)

    with pytest.raises(ValueError) as exc:
        validate_module.validate_all(tmp_path, additional_linkml=False, strict=True)

    msg = str(exc.value)
    assert "bad1" in msg
    assert "bad2" in msg


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


# def test_commented_bed():
#     lib.validate(data_dir / "bed-comment")


def test_dev_scheme():
    lib.validate(data_dir / "dev-scheme")
    lib.validate(data_dir / "dev-scheme", full=True)


def test_cli_create():
    run("mkdir -p built && rm -rf built/artic", cwd="./")
    run(
        "uv run primaschema create"
        " --name artic"
        " --amplicon-size 400"
        " --version v4.1.0"
        " --contributors 'ARTIC network'"
        " --target-organisms sars-cov-2"
        " --status DEPRECATED"
        " --bed-path test/data/dev-scheme/primer.bed"
        " --reference-path test/data/dev-scheme/reference.fasta"
        " --primer-schemes-path built",
        cwd="./",
    )
    assert Path("built/artic/400/v4.1.0/primer.bed").exists()
    assert Path("built/artic/400/v4.1.0/reference.fasta").exists()
    assert Path("built/artic/400/v4.1.0/info.json").exists()
    run("rm -rf built/artic", cwd="./")


def test_rebuild_syncs_metadata_from_path(tmp_path: Path):
    src = data_dir / "auto-normalisation/test/400/v2.0.0"
    dest = tmp_path / "artic-sars-cov-2" / "1200" / "v9.9.9"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)
    info_path = dest / "info.json"
    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())
    assert primer_scheme.name != "artic-sars-cov-2"
    assert primer_scheme.amplicon_size != 1200
    assert primer_scheme.version != "v9.9.9"

    from primaschema.cli import _rebuild_one

    _rebuild_one(info_path, sync_metadata=True)

    updated_scheme = PrimerScheme.model_validate_json(info_path.read_text())
    assert updated_scheme.name == "artic-sars-cov-2"
    assert updated_scheme.amplicon_size == 1200
    assert updated_scheme.version == "v9.9.9"
