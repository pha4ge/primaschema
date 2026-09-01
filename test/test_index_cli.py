import gzip
import json
import shutil
import subprocess
from pathlib import Path

DATA_DIR = Path("test/data")
FIXTURES = DATA_DIR / "primer-schemes"


def run(cmd):
    return subprocess.run(
        cmd, cwd="./", shell=True, check=True, text=True, capture_output=True
    )


def test_cli_index_builds_from_fixtures(tmp_path: Path):
    """`primaschema index` builds index.json + .gz from a directory of real schemes."""
    run(
        f"uv run primaschema index --primer-schemes-path {FIXTURES}"
        " --base-url https://example.invalid/schemes"
        " --source-commit deadbeef123"
        f" --output-path {tmp_path}"
    )

    index_path = tmp_path / "index.json"
    gz_path = tmp_path / "index.json.gz"
    assert index_path.exists()
    assert gz_path.exists()

    data = json.loads(index_path.read_text())
    assert data["source_commit"] == "deadbeef123"

    schemes = data["primerschemes"]
    assert schemes["test-artic"]["400"]["v4.1.0"]
    assert schemes["test-eden"]["2500"]["v1.0.0"]
    assert schemes["test-midnight"]["1200"]["v1.0.0"]
    assert schemes["test-midnight"]["1200"]["v2.0.0"]

    entry = schemes["test-artic"]["400"]["v4.1.0"]
    assert (
        entry["primer_file_url"]
        == "https://example.invalid/schemes/test-artic/400/v4.1.0/primer.bed"
    )

    assert (
        gzip.decompress(gz_path.read_bytes()).decode("utf-8") == index_path.read_text()
    )


def test_cli_index_merges_with_existing_index_path(tmp_path: Path):
    """`primaschema index --index-path` merges newly found schemes into a pre-existing index."""
    subset_dir = tmp_path / "subset"
    shutil.copytree(FIXTURES / "test-artic", subset_dir / "test-artic")

    first_index = tmp_path / "first"
    first_index.mkdir()
    run(
        f"uv run primaschema index --primer-schemes-path {subset_dir}"
        f" --output-path {first_index}"
    )
    first_data = json.loads((first_index / "index.json").read_text())
    assert list(first_data["primerschemes"]) == ["test-artic"]

    merged_index = tmp_path / "merged"
    merged_index.mkdir()
    run(
        f"uv run primaschema index --primer-schemes-path {FIXTURES}"
        f" --index-path {first_index / 'index.json'}"
        f" --output-path {merged_index}"
    )
    merged_data = json.loads((merged_index / "index.json").read_text())
    assert set(merged_data["primerschemes"]) == {
        "test-artic",
        "test-eden",
        "test-midnight",
    }
