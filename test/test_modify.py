import shutil
from datetime import date
from pathlib import Path

from primaschema.cli import (
    add_contributor,
    add_target_organism,
    add_vendor,
    remove_contributor,
    remove_target_organism,
    update_date_added,
    update_date_created,
)
from primaschema.schema.info import (
    PrimerSchemeContributor,
    PrimerSchemeTargetOrganism,
    PrimerSchemeVendor,
)
from primaschema.schema.primer_scheme import PrimerScheme

data_dir = Path("test/data")
FIXTURE = "auto-normalisation/test/400/v2.0.0"


def _copy_scheme(tmp_path: Path, rel_path: str) -> Path:
    src = data_dir / rel_path
    dest = tmp_path / rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)
    return dest


def test_add_contributor_persists(tmp_path):
    """add_contributor appends a new contributor to info.json and persists it."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    add_contributor(
        info_path, PrimerSchemeContributor(primer_scheme_contributor_name="Alice")
    )
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert (
        len(after.primer_scheme_contributor)
        == len(before.primer_scheme_contributor) + 1
    )
    assert after.primer_scheme_contributor[-1].primer_scheme_contributor_name == "Alice"


def test_remove_contributor_persists(tmp_path):
    """remove_contributor removes the contributor at the given index from info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    remove_contributor(info_path, 0)
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert (
        len(after.primer_scheme_contributor)
        == len(before.primer_scheme_contributor) - 1
    )


def test_add_vendor_persists(tmp_path):
    """add_vendor appends a new vendor to info.json and persists it."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    add_vendor(info_path, PrimerSchemeVendor(primer_scheme_vendor_name="NewCo"))
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert len(after.primer_scheme_vendor) == len(before.primer_scheme_vendor) + 1
    assert after.primer_scheme_vendor[-1].primer_scheme_vendor_name == "NewCo"


def test_add_target_organism_persists(tmp_path):
    """add_target_organism appends a new organism to info.json and persists it."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    add_target_organism(
        info_path,
        PrimerSchemeTargetOrganism(primer_scheme_target_organism_name="Test virus"),
    )
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert (
        len(after.primer_scheme_target_organism)
        == len(before.primer_scheme_target_organism) + 1
    )
    assert (
        after.primer_scheme_target_organism[-1].primer_scheme_target_organism_name
        == "Test virus"
    )


def test_remove_target_organism_persists(tmp_path):
    """remove_target_organism removes the organism at the given index from info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    remove_target_organism(info_path, 0)
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert (
        len(after.primer_scheme_target_organism)
        == len(before.primer_scheme_target_organism) - 1
    )


def test_update_date_created_persists(tmp_path):
    """update_date_created writes the new date_created value to info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    update_date_created(info_path, date(2023, 6, 15))
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    assert ps.primer_scheme_creation_date == date(2023, 6, 15)


def test_update_date_added_persists(tmp_path):
    """update_date_added writes the new date_added value to info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    update_date_added(info_path, date(2024, 3, 1))
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    assert ps.primer_scheme_submission_date == date(2024, 3, 1)
