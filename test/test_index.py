from datetime import date
from pathlib import Path

import pytest

from primaschema.schema.index import (
    IndexPrimerScheme,
    PrimerSchemeIndex,
    create_index,
    update_index,
)
from primaschema.schema.primer_scheme import PrimerScheme

DATA_DIR = Path("test/data")


def _load_scheme() -> PrimerScheme:
    info_path = DATA_DIR / "auto-normalisation/test/400/v2.0.0/info.json"
    return PrimerScheme.model_validate_json(info_path.read_text())


def _clone_scheme(
    scheme: PrimerScheme,
    *,
    name: str | None = None,
    version: str | None = None,
    date_created: date | None = None,
    date_added: date | None = None,
    primer_sha256: str | None = None,
) -> PrimerScheme:
    data = scheme.model_dump()
    if name is not None:
        data["primer_scheme_name"] = name
    if version is not None:
        data["primer_scheme_version"] = version
    if date_created is not None:
        data["primer_scheme_creation_date"] = date_created
    if date_added is not None:
        data["primer_scheme_submission_date"] = date_added
    if primer_sha256 is not None:
        data["primer_scheme_checksums"]["primer_scheme_sha256"] = primer_sha256
    return PrimerScheme.model_validate(data)


def test_index_builds_urls_from_base_url():
    """update_index synthesises correct primer, reference and info URLs when given a base_url."""
    ps = _load_scheme()
    psi = PrimerSchemeIndex()
    base_url = "https://github.com/pha4ge/primer-schemes/main/v1b/schemes"

    update_index([ps], psi, base_url=base_url)

    data = psi.model_dump()
    entry = data["primerschemes"]["test"][400]["v2.0.0"]
    assert entry["primer_file_url"] == f"{base_url}/test/400/v2.0.0/primer.bed"
    assert entry["reference_file_url"] == f"{base_url}/test/400/v2.0.0/reference.fasta"
    assert entry["info_file_url"] == f"{base_url}/test/400/v2.0.0/info.json"


def test_index_preserves_urls_on_reload():
    """Serialised URL fields survive a model_dump_json → model_validate_json round-trip unchanged."""
    ps = _load_scheme()
    psi = PrimerSchemeIndex()
    base_url = "https://github.com/pha4ge/primer-schemes/main/v1b/schemes"

    update_index([ps], psi, base_url=base_url)
    dumped = psi.model_dump_json()

    reloaded = PrimerSchemeIndex.model_validate_json(dumped)
    entry = reloaded.primerschemes["test"][400]["v2.0.0"]

    assert entry.primer_file_url == f"{base_url}/test/400/v2.0.0/primer.bed"
    assert entry.reference_file_url == f"{base_url}/test/400/v2.0.0/reference.fasta"
    assert entry.info_file_url == f"{base_url}/test/400/v2.0.0/info.json"


def test_index_round_trip_stability():
    """A full serialize → deserialize cycle produces an identical model_dump (no data loss or mutation)."""
    ps = _load_scheme()
    psi = PrimerSchemeIndex()
    base_url = "https://github.com/pha4ge/primer-schemes/main/v1b/schemes"

    update_index([ps], psi, base_url=base_url)
    original = psi.model_dump()
    dumped = psi.model_dump_json()

    reloaded = PrimerSchemeIndex.model_validate_json(dumped)
    round_tripped = reloaded.model_dump()

    assert round_tripped == original


def test_index_flatten_and_get():
    """flatten() and get_schemes_from_index() correctly enumerate and filter indexed schemes."""
    ps = _load_scheme()
    ps2 = _clone_scheme(ps, version="v2.0.1")
    psi = PrimerSchemeIndex()
    base_url = "https://github.com/pha4ge/primer-schemes/main/v1b/schemes"

    update_index([ps, ps2], psi, base_url=base_url)

    flattened = psi.flatten()
    assert len(flattened) == 2
    versions = {entry.primer_scheme_version for entry in flattened}
    assert versions == {"v2.0.0", "v2.0.1"}

    all_for_name = psi.get_schemes_from_index("test")
    assert len(all_for_name) == 2

    filtered = psi.get_schemes_from_index("test", 400, "v2.0.1")
    assert len(filtered) == 1
    assert filtered[0].primer_scheme_version == "v2.0.1"


def test_index_includes_dates():
    """update_index copies date_created and date_added from PrimerScheme into the index entry."""
    ps = _clone_scheme(
        _load_scheme(),
        date_created=date(2024, 1, 15),
        date_added=date(2024, 6, 1),
    )
    psi = PrimerSchemeIndex()
    update_index([ps], psi)
    entry = psi.primerschemes["test"][400]["v2.0.0"]
    assert entry.primer_scheme_creation_date == date(2024, 1, 15)
    assert entry.primer_scheme_submission_date == date(2024, 6, 1)


def test_index_dates_round_trip():
    """Date fields survive a model_dump_json → model_validate_json cycle in the index."""
    ps = _clone_scheme(
        _load_scheme(),
        date_created=date(2024, 1, 15),
        date_added=date(2024, 6, 1),
    )
    psi = PrimerSchemeIndex()
    update_index([ps], psi)
    reloaded = PrimerSchemeIndex.model_validate_json(psi.model_dump_json())
    entry = reloaded.primerschemes["test"][400]["v2.0.0"]
    assert entry.primer_scheme_creation_date == date(2024, 1, 15)
    assert entry.primer_scheme_submission_date == date(2024, 6, 1)


def test_index_dates_absent_when_none():
    """IndexPrimerScheme.from_primer_scheme leaves date fields as None when the source scheme has no dates."""
    ps = _load_scheme()
    entry = IndexPrimerScheme.from_primer_scheme(ps)
    assert entry.primer_scheme_creation_date is None
    assert entry.primer_scheme_submission_date is None


def test_create_index_builds_full_structure():
    """create_index() builds an index from scratch, keyed by name/amplicon_size/version."""
    ps = _load_scheme()
    ps2 = _clone_scheme(ps, name="test2")

    psi = create_index([ps, ps2])

    assert set(psi.primerschemes) == {"test", "test2"}
    assert psi.primerschemes["test"][400]["v2.0.0"].primer_scheme_name == "test"
    assert psi.primerschemes["test2"][400]["v2.0.0"].primer_scheme_name == "test2"


def test_add_index_primer_scheme_strict_conflict_raises():
    """add_index_primer_scheme raises if an existing version's checksums differ and strict=True."""
    ps = _load_scheme()
    changed = _clone_scheme(ps, primer_sha256="0" * 64)
    psi = create_index([ps])

    with pytest.raises(ValueError, match="checksums have changed"):
        psi.add_index_primer_scheme(
            IndexPrimerScheme.from_primer_scheme(changed), strict=True
        )


def test_add_index_primer_scheme_non_strict_overwrites():
    """add_index_primer_scheme allows a checksum-changed overwrite when strict=False."""
    ps = _load_scheme()
    changed = _clone_scheme(ps, primer_sha256="0" * 64)
    psi = create_index([ps])

    psi.add_index_primer_scheme(
        IndexPrimerScheme.from_primer_scheme(changed), strict=False
    )

    entry = psi.primerschemes["test"][400]["v2.0.0"]
    assert entry.primer_scheme_checksums.primer_scheme_sha256 == "0" * 64


def test_remove_index_primer_scheme_removes_entry():
    """remove_index_primer_scheme removes a present entry and reports success."""
    ps = _load_scheme()
    psi = create_index([ps])
    entry = IndexPrimerScheme.from_primer_scheme(ps)

    assert psi.remove_index_primer_scheme(entry) is True
    assert psi.get_schemes_from_index("test") == []


def test_remove_index_primer_scheme_prunes_empty_branches():
    """Removing the only version under a name/amplicon_size prunes the now-empty branches."""
    ps = _load_scheme()
    psi = create_index([ps])
    entry = IndexPrimerScheme.from_primer_scheme(ps)

    psi.remove_index_primer_scheme(entry)

    assert "test" not in psi.primerschemes


def test_remove_index_primer_scheme_absent_returns_false():
    """remove_index_primer_scheme returns False for a name/amplicon_size/version not in the index."""
    ps = _load_scheme()
    psi = create_index([ps])
    other = IndexPrimerScheme.from_primer_scheme(_clone_scheme(ps, name="missing"))

    assert psi.remove_index_primer_scheme(other) is False


def test_source_commit_round_trip():
    """source_commit survives a serialize -> deserialize cycle and is excluded when unset."""
    psi = PrimerSchemeIndex()
    assert "source_commit" not in psi.model_dump_json(exclude_unset=True)

    psi.source_commit = "deadbeef"
    reloaded = PrimerSchemeIndex.model_validate_json(psi.model_dump_json())
    assert reloaded.source_commit == "deadbeef"


def _corrupt_index(**overrides) -> PrimerSchemeIndex:
    """Build an index whose single entry is filed under deliberately wrong keys."""
    ps = _load_scheme()
    psi = create_index([ps])
    entry = psi.flatten()[0]
    name = overrides.get("name", entry.primer_scheme_name)
    amplicon_size = overrides.get("amplicon_size", entry.amplicon_size)
    version = overrides.get("version", entry.primer_scheme_version)
    psi.primerschemes = {name: {amplicon_size: {version: entry}}}
    return psi


def test_validate_keys_accepts_a_well_formed_index():
    assert create_index([_load_scheme()]).validate_keys() == []


@pytest.mark.parametrize(
    "override,expected",
    [
        ({"name": "wrong-name"}, "primer_scheme_name"),
        ({"amplicon_size": 999}, "amplicon_size"),
        ({"version": "v9.9.9"}, "primer_scheme_version"),
    ],
)
def test_validate_keys_rejects_entry_filed_under_wrong_key(override, expected):
    """Each of the three key levels is checked, not just the ones a lookup filters on."""
    psi = _corrupt_index(**override)
    with pytest.raises(ValueError, match="Corrupted index") as excinfo:
        psi.validate_keys()
    assert expected in str(excinfo.value)


def test_validate_keys_warns_instead_of_raising_under_force(caplog):
    """--force leaves an escape hatch if a bad index is ever published."""
    psi = _corrupt_index(name="wrong-name")
    with caplog.at_level("WARNING"):
        mismatches = psi.validate_keys(force=True)
    assert mismatches
    assert "continuing due to --force" in caplog.text


def test_validate_keys_catches_corruption_a_name_only_lookup_would_miss():
    """A name-only lookup never filters on size/version, so those keys must still be checked."""
    psi = _corrupt_index(amplicon_size=999, version="v9.9.9")
    name = next(iter(psi.primerschemes))

    # The lookup itself returns the entry regardless; validate_keys is what catches it.
    assert len(psi.get_schemes_from_index(name)) == 1
    with pytest.raises(ValueError, match="Corrupted index"):
        psi.validate_keys()


def test_load_index_validates_keys(tmp_path):
    """load_index is the choke point, so --all and every other path get checked."""
    from primaschema import get_scheme

    index_path = tmp_path / "index.json"
    index_path.write_text(_corrupt_index(name="wrong-name").model_dump_json())

    with pytest.raises(ValueError, match="Corrupted index"):
        get_scheme.load_index(str(index_path))

    # force downgrades it to a warning and still returns a usable index
    assert get_scheme.load_index(str(index_path), force=True) is not None
