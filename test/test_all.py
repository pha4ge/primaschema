import shutil
import subprocess
from pathlib import Path

import pytest

import primaschema.cli as create_module
import primaschema.get_scheme as get_scheme
import primaschema.lib as lib
import primaschema.validate as validate_module
from primaschema.schema.index import PrimerSchemeIndex
from primaschema.schema.info import PrimerScheme
from primaschema.util import sha256_checksum

data_dir = Path("test/data")


def run(cmd, cwd=data_dir):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def test_plot_single_ref_chrom_ref():
    lib.plot_primers(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed",
    )
    run("rm -rf primer.html", cwd="./")


def test_plot_many_ref_chroms_ref():
    lib.plot_primers(data_dir / "many-ref-chroms/primer.bed")
    run("rm -rf primer.html", cwd="./")


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

    assert sha256_checksum(primer_path) != primer_scheme.checksums.primer_sha256

    validate_module.validate(info_path, strict=True, fix=True)

    assert sha256_checksum(primer_path) == primer_scheme.checksums.primer_sha256


def test_validate_autonormalize_reference_fasta(tmp_path: Path):
    scheme_dir = _copy_scheme(
        tmp_path,
        "auto-normalisation/test/400/v2.0.0",
    )
    info_path = scheme_dir / "info.json"
    reference_path = scheme_dir / "reference.fasta"
    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())

    assert sha256_checksum(reference_path) != primer_scheme.checksums.reference_sha256

    validate_module.validate(info_path, strict=True, fix=True)

    assert sha256_checksum(reference_path) == primer_scheme.checksums.reference_sha256


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


@pytest.mark.network
def test_get_scheme(tmp_path: Path):
    psi = get_scheme.load_index(get_scheme.DEFAULT_INDEX_URL)
    schemes = get_scheme.resolve_schemes(
        index=psi,
        scheme_id="artic/400/v4.1.0",
        allow_multiple=False,
        all_schemes=False,
    )
    outputs = get_scheme.download_schemes(schemes=schemes, output=tmp_path)
    output_dir = outputs[0]
    assert (output_dir / "info.json").exists()
    assert (output_dir / "primer.bed").exists()
    assert (output_dir / "reference.fasta").exists()


def test_get_scheme_invalid_id():
    with pytest.raises(ValueError, match="expected format"):
        get_scheme.resolve_schemes(
            index=PrimerSchemeIndex(),
            scheme_id="artic/400/v1.0.0/extra",
            allow_multiple=False,
            all_schemes=False,
        )


@pytest.mark.network
def test_get_scheme_nonexistent(tmp_path: Path):
    psi = get_scheme.load_index(get_scheme.DEFAULT_INDEX_URL)
    with pytest.raises(ValueError, match="not found"):
        get_scheme.resolve_schemes(
            index=psi,
            scheme_id="nonexistent/999/v0.0.0",
            allow_multiple=False,
            all_schemes=False,
        )


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
