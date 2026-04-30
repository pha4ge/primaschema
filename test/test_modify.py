import shutil
from datetime import date
from pathlib import Path

from primaschema.cli import (
    add_contributor,
    add_tag,
    add_target_organism,
    add_vendor,
    remove_contributor,
    remove_tag,
    remove_target_organism,
    update_date_added,
    update_date_created,
)
from primaschema.schema.info import (
    Contributor,
    PrimerScheme,
    SchemeTag,
    TargetOrganism,
    Vendor,
)

data_dir = Path("test/data")
FIXTURE = "auto-normalisation/test/400/v2.0.0"


def _copy_scheme(tmp_path: Path, rel_path: str) -> Path:
    src = data_dir / rel_path
    dest = tmp_path / rel_path
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)
    return dest


def test_add_tag_persists(tmp_path):
    """add_tag must write the new tag to info.json (was silently dropped)."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    add_tag(info_path, SchemeTag.CLINICAL)
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    assert ps.tags is not None
    assert SchemeTag.CLINICAL in ps.tags


def test_remove_tag_persists(tmp_path):
    """remove_tag writes the deletion to info.json so the tag is absent on reload."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    add_tag(info_path, SchemeTag.CLINICAL)
    remove_tag(info_path, SchemeTag.CLINICAL)
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    assert ps.tags is None or SchemeTag.CLINICAL not in ps.tags


def test_add_contributor_persists(tmp_path):
    """add_contributor appends a new contributor to info.json and persists it."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    add_contributor(info_path, Contributor(name="Alice"))
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert len(after.contributors) == len(before.contributors) + 1
    assert after.contributors[-1].name == "Alice"


def test_remove_contributor_persists(tmp_path):
    """remove_contributor removes the contributor at the given index from info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    remove_contributor(info_path, 0)
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert len(after.contributors) == len(before.contributors) - 1


def test_add_vendor_persists(tmp_path):
    """add_vendor appends a new vendor to info.json and persists it."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    add_vendor(info_path, Vendor(organisation_name="NewCo"))
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert len(after.vendors) == len(before.vendors) + 1
    assert after.vendors[-1].organisation_name == "NewCo"


def test_add_target_organism_persists(tmp_path):
    """add_target_organism appends a new organism to info.json and persists it."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    add_target_organism(info_path, TargetOrganism(common_name="Test virus"))
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert len(after.target_organisms) == len(before.target_organisms) + 1
    assert after.target_organisms[-1].common_name == "Test virus"


def test_remove_target_organism_persists(tmp_path):
    """remove_target_organism removes the organism at the given index from info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    before = PrimerScheme.model_validate_json(info_path.read_text())
    remove_target_organism(info_path, 0)
    after = PrimerScheme.model_validate_json(info_path.read_text())
    assert len(after.target_organisms) == len(before.target_organisms) - 1


def test_update_date_created_persists(tmp_path):
    """update_date_created writes the new date_created value to info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    update_date_created(info_path, date(2023, 6, 15))
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    assert ps.date_created == date(2023, 6, 15)


def test_update_date_added_persists(tmp_path):
    """update_date_added writes the new date_added value to info.json."""
    info_path = _copy_scheme(tmp_path, FIXTURE) / "info.json"
    update_date_added(info_path, date(2024, 3, 1))
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    assert ps.date_added == date(2024, 3, 1)
