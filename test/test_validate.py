import json
import shutil
import subprocess
from pathlib import Path

import pytest

import primaschema.validate as validate_module
from primaschema.schema.info import PrimerScheme

DATA_DIR = Path("test/data")


def _copy_scheme(tmp_path: Path) -> Path:
    src = DATA_DIR / "auto-normalisation/test/400/v2.0.0"
    dest = tmp_path / "auto-normalisation/test/400/v2.0.0"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)
    return dest


def _run(cmd: list[str], cwd: Path = Path(".")) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=cwd, check=True, text=True, capture_output=True)


def test_validate_linkml_accepts_valid_info(tmp_path: Path):
    """LinkML validation passes for a valid scheme info.json."""
    scheme_dir = _copy_scheme(tmp_path)
    info_path = scheme_dir / "info.json"

    validate_module.validate_scheme_json_with_linkml(info_path)


def test_validate_linkml_rejects_missing_required_field(tmp_path: Path):
    """LinkML validation fails when a required info.json field is missing."""
    scheme_dir = _copy_scheme(tmp_path)
    info_path = scheme_dir / "info.json"

    data = json.loads(info_path.read_text())
    data.pop("name", None)
    info_path.write_text(json.dumps(data))

    with pytest.raises(ValueError):
        validate_module.validate_scheme_json_with_linkml(info_path)


def test_validate_with_linkml_succeeds(tmp_path: Path):
    """Top-level validate succeeds when LinkML checks are enabled on valid input."""
    scheme_dir = _copy_scheme(tmp_path)
    info_path = scheme_dir / "info.json"

    validate_module.validate(
        info_path,
        additional_linkml=True,
        strict=False,
        fix=True,
    )


def test_validate_name_mismatch_raises(tmp_path: Path):
    """validate_name raises when directory path metadata disagrees with info.json."""
    src = DATA_DIR / "auto-normalisation/test/400/v2.0.0"
    dest = tmp_path / "wrong-name/400/v2.0.0"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dest)

    info_path = dest / "info.json"

    with pytest.raises(ValueError):
        validate_module.validate_name(info_path)


def test_validate_readme_missing_raises(tmp_path: Path):
    """validate_readme raises when README.md is absent."""
    scheme_dir = _copy_scheme(tmp_path)
    readme_path = scheme_dir / "README.md"
    readme_path.unlink()

    info_path = scheme_dir / "info.json"

    with pytest.raises(FileNotFoundError):
        validate_module.validate_readme(info_path)


def test_validate_hashes_non_normalisable_mismatch_raises(tmp_path: Path):
    """validate_hashes raises for a non-normalisable primer.bed checksum mismatch."""
    scheme_dir = _copy_scheme(tmp_path)
    info_path = scheme_dir / "info.json"
    primer_path = scheme_dir / "primer.bed"

    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())

    primer_lines = primer_path.read_text().splitlines()
    header_lines = [line for line in primer_lines if line.startswith("#")]
    body_lines = [line for line in primer_lines if line and not line.startswith("#")]
    assert body_lines, "Expected at least one primer.bed record"

    updated_lines = header_lines + body_lines + [body_lines[0]]
    primer_path.write_text("\n".join(updated_lines) + "\n")

    with pytest.raises(ValueError):
        validate_module.validate_hashes(info_path, primer_scheme, fix=False)


def test_create_range_infos_all_pass_linkml(tmp_path: Path):
    """Create multiple schemes via CLI and ensure every generated info.json passes LinkML."""
    output_root = tmp_path / "generated-schemes"
    output_root.mkdir()

    bed_path = Path("test/data/dev-scheme/primer.bed")
    reference_path = Path("test/data/dev-scheme/reference.fasta")

    cases = [
        {
            "name": "batch-validate-a",
            "amplicon_size": "400",
            "version": "v1.0.0",
            "status": "DRAFT",
            "license": "CC0-1.0",
            "target_organism": "common_name=SARS-CoV-2,ncbi_tax_id=2697049",
        },
        {
            "name": "batch-validate-b",
            "amplicon_size": "500",
            "version": "v1.1.0",
            "status": "TESTED",
            "license": "CC-BY-4.0",
            "target_organism": "common_name=Influenza A,ncbi_tax_id=11320",
        },
        {
            "name": "batch-validate-c",
            "amplicon_size": "1200",
            "version": "v2.0.0",
            "status": "VALIDATED",
            "license": "CC-BY-SA-4.0",
            "target_organism": "common_name=Mpox virus,ncbi_tax_id=10244",
        },
        {
            "name": "batch-validate-d",
            "amplicon_size": "800",
            "version": "v3.0.0",
            "status": "DEPRECATED",
            "license": "CC-BY-NC-ND-4.0",
            "target_organism": "common_name=Dengue virus,ncbi_tax_id=12637",
        },
    ]

    for case in cases:
        _run(
            [
                "uv",
                "run",
                "primaschema",
                "create",
                "--name",
                case["name"],
                "--amplicon-size",
                case["amplicon_size"],
                "--version",
                case["version"],
                "--contributors",
                "name=Validation Bot,email=validation@example.org",
                "--target-organisms",
                case["target_organism"],
                "--status",
                case["status"],
                "--license",
                case["license"],
                "--date-created",
                "2024-01-01",
                "--bed-path",
                str(bed_path),
                "--reference-path",
                str(reference_path),
                "--primer-schemes-path",
                str(output_root),
            ]
        )

    info_paths = sorted(output_root.rglob("info.json"))
    assert len(info_paths) == len(cases)

    for info_path in info_paths:
        PrimerScheme.model_validate_json(info_path.read_text())
        validate_module.validate_scheme_json_with_linkml(info_path)
