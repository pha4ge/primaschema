from datetime import date
from pathlib import Path

from primaschema.schema.index import IndexPrimerScheme, PrimerSchemeIndex, update_index
from primaschema.schema.info import PrimerScheme

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
) -> PrimerScheme:
    data = scheme.model_dump()
    if name is not None:
        data["name"] = name
    if version is not None:
        data["version"] = version
    if date_created is not None:
        data["date_created"] = date_created
    if date_added is not None:
        data["date_added"] = date_added
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
    versions = {entry.version for entry in flattened}
    assert versions == {"v2.0.0", "v2.0.1"}

    all_for_name = psi.get_schemes_from_index("test")
    assert len(all_for_name) == 2

    filtered = psi.get_schemes_from_index("test", 400, "v2.0.1")
    assert len(filtered) == 1
    assert filtered[0].version == "v2.0.1"


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
    assert entry.date_created == date(2024, 1, 15)
    assert entry.date_added == date(2024, 6, 1)


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
    assert entry.date_created == date(2024, 1, 15)
    assert entry.date_added == date(2024, 6, 1)


def test_index_dates_absent_when_none():
    """IndexPrimerScheme.from_primer_scheme leaves date fields as None when the source scheme has no dates."""
    ps = _load_scheme()
    entry = IndexPrimerScheme.from_primer_scheme(ps)
    assert entry.date_created is None
    assert entry.date_added is None
