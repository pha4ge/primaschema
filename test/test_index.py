from pathlib import Path

from primaschema.schema.index import PrimerSchemeIndex, update_index
from primaschema.schema.info import PrimerScheme

DATA_DIR = Path("test/data")


def _load_scheme() -> PrimerScheme:
    info_path = DATA_DIR / "auto-normalisation/test/400/v2.0.0/info.json"
    return PrimerScheme.model_validate_json(info_path.read_text())


def _clone_scheme(
    scheme: PrimerScheme, *, name: str | None = None, version: str | None = None
) -> PrimerScheme:
    data = scheme.model_dump()
    if name is not None:
        data["name"] = name
    if version is not None:
        data["version"] = version
    return PrimerScheme.model_validate(data)


def test_index_builds_urls_from_base_url():
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
