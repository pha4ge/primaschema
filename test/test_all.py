import shutil
import subprocess
from datetime import date
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
    """plot_primers completes without error for a single-chromosome scheme."""
    lib.plot_primers(
        data_dir / "primer-schemes/schemes/sars-cov-2/artic/400/v4.1.0/primer.bed",
    )
    run("rm -rf primer.html", cwd="./")


def test_plot_many_ref_chroms_ref():
    """plot_primers handles a scheme with multiple reference chromosomes."""
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
    """validate --fix rewrites an out-of-order primer.bed to match the stored checksum."""
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
    """validate --fix rewrites a malformatted reference.fasta to match the stored checksum."""
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
    """validate --all via CLI collects errors from multiple bad schemes into one message."""
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
    """validate_all() module function collects errors from multiple bad schemes into one message."""
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
    """primaschema create produces info.json, primer.bed, and reference.fasta on disk."""
    run("mkdir -p built && rm -rf built/artic", cwd="./")
    run(
        "uv run primaschema create"
        " --name artic"
        " --amplicon-size 400"
        " --version v4.1.0"
        " --contributors 'ARTIC network'"
        " --target-organisms sars-cov-2"
        " --status DEPRECATED"
        " --date-created 2020-09-04"
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
    """download_schemes fetches a real scheme from the default index."""
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
    """resolve_schemes raises ValueError for a scheme_id with too many path segments."""
    with pytest.raises(ValueError, match="expected format"):
        get_scheme.resolve_schemes(
            index=PrimerSchemeIndex(),
            scheme_id="artic/400/v1.0.0/extra",
            allow_multiple=False,
            all_schemes=False,
        )


@pytest.mark.network
def test_get_scheme_nonexistent(tmp_path: Path):
    """resolve_schemes raises ValueError when the scheme_id is not in the index."""
    psi = get_scheme.load_index(get_scheme.DEFAULT_INDEX_URL)
    with pytest.raises(ValueError, match="not found"):
        get_scheme.resolve_schemes(
            index=psi,
            scheme_id="nonexistent/999/v0.0.0",
            allow_multiple=False,
            all_schemes=False,
        )


def test_rebuild_syncs_metadata_from_path(tmp_path: Path):
    """rebuild --sync-metadata updates name, amplicon_size, and version from the directory path."""
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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _minimal_scheme(**kwargs) -> PrimerScheme:
    from primaschema.schema.info import Contributor, SchemeStatus, TargetOrganism

    defaults = dict(
        schema_version="1.0.0",
        name="test-scheme",
        amplicon_size=400,
        version="v1.0.0",
        contributors=[Contributor(name="Alice")],
        target_organisms=[TargetOrganism(common_name="SARS-CoV-2")],
        status=SchemeStatus.DRAFT,
    )
    defaults.update(kwargs)
    return PrimerScheme(**defaults)


# ---------------------------------------------------------------------------
# Unit tests — SchemeLicense
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "spdx",
    [
        "CC0-1.0",
        "CC-BY-4.0",
        "CC-BY-SA-4.0",
        "CC-BY-NC-4.0",
        "CC-BY-NC-SA-4.0",
        "CC-BY-ND-4.0",
        "CC-BY-NC-ND-4.0",
    ],
)
def test_scheme_license_all_values_accessible(spdx):
    """Every SPDX license string round-trips through SchemeLicense(value).value."""
    from primaschema.schema.info import SchemeLicense

    assert SchemeLicense(spdx).value == spdx


def test_license_footers_covers_all_licenses():
    """LICENSE_FOOTERS has an entry for every SchemeLicense member."""
    from primaschema.license_footers import LICENSE_FOOTERS
    from primaschema.schema.info import SchemeLicense

    for member in SchemeLicense:
        assert member in LICENSE_FOOTERS, f"Missing footer for {member}"


@pytest.mark.parametrize(
    "license,url_fragment",
    [
        ("CC0-1.0", "publicdomain/zero/1.0"),
        ("CC-BY-4.0", "licenses/by/4.0"),
        ("CC-BY-SA-4.0", "licenses/by-sa/4.0"),
        ("CC-BY-NC-4.0", "licenses/by-nc/4.0"),
        ("CC-BY-NC-SA-4.0", "licenses/by-nc-sa/4.0"),
        ("CC-BY-ND-4.0", "licenses/by-nd/4.0"),
        ("CC-BY-NC-ND-4.0", "licenses/by-nc-nd/4.0"),
    ],
)
def test_license_footer_contains_url(license, url_fragment):
    """Each license footer contains the canonical CC URL for that license."""
    from primaschema.license_footers import LICENSE_FOOTERS
    from primaschema.schema.info import SchemeLicense

    assert url_fragment in LICENSE_FOOTERS[SchemeLicense(license)]


# ---------------------------------------------------------------------------
# Unit tests — date fields
# ---------------------------------------------------------------------------


def test_primer_scheme_dates_optional():
    """PrimerScheme can be constructed without dates; both default to None."""
    ps = _minimal_scheme()
    assert ps.date_created is None
    assert ps.date_added is None


def test_primer_scheme_dates_accept_valid():
    """PrimerScheme accepts date objects for date_created and date_added."""
    ps = _minimal_scheme(date_created=date(2024, 1, 15), date_added=date(2024, 6, 1))
    assert ps.date_created == date(2024, 1, 15)
    assert ps.date_added == date(2024, 6, 1)


def test_cli_scheme_date_created_required():
    """CLIPrimerScheme raises ValidationError when date_created is omitted."""
    from pydantic import ValidationError
    from primaschema.cli import CLIPrimerScheme
    from primaschema.schema.info import Contributor, SchemeStatus, TargetOrganism

    with pytest.raises(ValidationError):
        CLIPrimerScheme(
            schema_version="1.0.0",
            name="test",
            amplicon_size=400,
            version="v1.0.0",
            status=SchemeStatus.DRAFT,
            contributors=[Contributor(name="Alice")],
            target_organisms=[TargetOrganism(common_name="SARS-CoV-2")],
        )


def test_cli_scheme_date_added_defaults_to_today():
    """CLIPrimerScheme sets date_added to today when not explicitly provided."""
    from unittest.mock import patch
    from primaschema.cli import CLIPrimerScheme
    from primaschema.schema.info import Contributor, SchemeStatus, TargetOrganism

    fixed = date(2026, 3, 30)
    with patch("primaschema.cli.date") as mock_date:
        mock_date.today.return_value = fixed
        mock_date.side_effect = lambda *a, **kw: date(*a, **kw)
        ps = CLIPrimerScheme(
            schema_version="1.0.0",
            name="test",
            amplicon_size=400,
            version="v1.0.0",
            status=SchemeStatus.DRAFT,
            contributors=[Contributor(name="Alice")],
            target_organisms=[TargetOrganism(common_name="SARS-CoV-2")],
            date_created=date(2024, 1, 1),
        )
    assert ps.date_added == fixed


# ---------------------------------------------------------------------------
# Integration tests — generate_readme
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "license,url_fragment",
    [
        ("CC0-1.0", "publicdomain/zero/1.0"),
        ("CC-BY-4.0", "licenses/by/4.0"),
        ("CC-BY-SA-4.0", "licenses/by-sa/4.0"),
        ("CC-BY-NC-4.0", "licenses/by-nc/4.0"),
        ("CC-BY-NC-SA-4.0", "licenses/by-nc-sa/4.0"),
        ("CC-BY-ND-4.0", "licenses/by-nd/4.0"),
        ("CC-BY-NC-ND-4.0", "licenses/by-nc-nd/4.0"),
    ],
)
def test_readme_license_footer_written(tmp_path, license, url_fragment):
    """generate_readme writes the correct license footer for each CC license."""
    from primaschema.cli import generate_readme
    from primaschema.schema.info import SchemeLicense

    ps = _minimal_scheme(license=SchemeLicense(license))
    generate_readme(tmp_path, ps)
    readme = (tmp_path / "README.md").read_text()
    assert (
        "------------------------------------------------------------------------"
        in readme
    )
    assert url_fragment in readme


def test_readme_no_footer_when_no_license(tmp_path):
    """generate_readme omits the license footer section when license is None."""
    from primaschema.cli import generate_readme

    ps = _minimal_scheme(license=None)
    generate_readme(tmp_path, ps)
    readme = (tmp_path / "README.md").read_text()
    assert (
        "------------------------------------------------------------------------"
        not in readme
    )


def test_readme_contains_dates_in_json(tmp_path):
    """generate_readme includes date_created and date_added in the embedded JSON block."""
    from primaschema.cli import generate_readme

    ps = _minimal_scheme(date_created=date(2024, 1, 15), date_added=date(2024, 6, 1))
    generate_readme(tmp_path, ps)
    readme = (tmp_path / "README.md").read_text()
    assert "2024-01-15" in readme
    assert "2024-06-01" in readme


def test_dates_round_trip():
    """Dates survive serialize → deserialize via serialize_primer_scheme_json."""
    from primaschema.util import serialize_primer_scheme_json

    ps = _minimal_scheme(date_created=date(2024, 1, 15), date_added=date(2024, 6, 1))
    json_bytes = serialize_primer_scheme_json(ps)
    restored = PrimerScheme.model_validate_json(json_bytes)
    assert restored.date_created == date(2024, 1, 15)
    assert restored.date_added == date(2024, 6, 1)


def test_dates_absent_when_none():
    """None dates are excluded from the serialized JSON (exclude_none=True)."""
    from primaschema.util import serialize_primer_scheme_json

    ps = _minimal_scheme()
    json_bytes = serialize_primer_scheme_json(ps)
    assert b"date_created" not in json_bytes
    assert b"date_added" not in json_bytes
