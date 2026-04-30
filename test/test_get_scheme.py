import hashlib
import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch
from urllib.parse import urlparse

import dnaio
from primalbedtools.scheme import Scheme

from primaschema import METADATA_FILE_NAME, PRIMER_FILE_NAME, REFERENCE_FILE_NAME
from primaschema.get_scheme import (
    DownloadError,
    SanitisationMode,
    download_schemes,
    resolve_schemes,
)
from primaschema.schema.index import PrimerSchemeIndex, update_index
from primaschema.schema.info import PrimerScheme
from primaschema.util import serialize_fasta_records, serialize_primer_scheme_json

DATA_DIR = Path("test/data")
REPO_DIR = DATA_DIR / "dev-repo"
SCHEMES_DIR = REPO_DIR / "schemes"


def _load_scheme() -> PrimerScheme:
    info_path = SCHEMES_DIR / "dev-scheme/400/v4.1.0/info.json"
    return PrimerScheme.model_validate_json(info_path.read_text())


def _mutate_primer_bytes(primer_bytes: bytes) -> bytes:
    text = primer_bytes.decode("utf-8")
    lines = text.splitlines()
    first = lines[0].split("\t")
    first[-1] = first[-1][:-1] + ("A" if first[-1][-1] != "A" else "C")
    lines[0] = "\t".join(first)
    return ("\n".join(lines) + "\n").encode("utf-8")


def _build_index(*schemes: PrimerScheme) -> PrimerSchemeIndex:
    psi = PrimerSchemeIndex()
    update_index(list(schemes), psi, base_url="https://example.com/schemes")
    return psi


def _scheme_dir_from_url(url: str) -> Path:
    parsed = urlparse(url)
    parts = [part for part in Path(parsed.path).parts if part]
    if "schemes" in parts:
        parts = parts[parts.index("schemes") + 1 :]
    if len(parts) < 4:
        raise AssertionError(f"Unexpected URL path: {url}")
    name, size, version = parts[-4], parts[-3], parts[-2]

    fallback: list[Path] = []
    for info_path in SCHEMES_DIR.rglob("info.json"):
        data = json.loads(info_path.read_text())
        if data.get("name") != name or str(data.get("amplicon_size")) != size:
            continue
        if data.get("version") == version:
            return info_path.parent
        fallback.append(info_path)

    if fallback:
        return sorted(fallback)[0].parent

    raise AssertionError(f"No scheme found for {name}/{size}/{version}")


def _repo_bytes(*, primer_mutator=None):
    def _read(url: str, **_kwargs) -> bytes:
        scheme_dir = _scheme_dir_from_url(url)
        filename = Path(urlparse(url).path).name
        data = (scheme_dir / filename).read_bytes()
        if filename == PRIMER_FILE_NAME and primer_mutator:
            data = primer_mutator(data)
        return data

    return _read


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class TestDownloadSchemes(unittest.TestCase):
    def test_download_scheme_all_flag(self):
        # Downloads all schemes when all_schemes is True.
        ps = _load_scheme()
        ps2 = ps.model_copy(update={"version": "v4.1.1"})
        psi = _build_index(ps, ps2)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id=None,
                    allow_multiple=False,
                    all_schemes=True,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                )

            self.assertEqual(len(outputs), 2)

    def test_download_scheme_sanitise_raw(self):
        # RAW sanitisation writes bytes exactly as downloaded.
        ps = _load_scheme()
        psi = _build_index(ps)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id="artic/400/v4.1.0",
                    allow_multiple=False,
                    all_schemes=False,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                )

            scheme_dir = _scheme_dir_from_url(
                "https://example.com/schemes/artic/400/v4.1.0/info.json"
            )
            output_dir = outputs[0]
            self.assertEqual(
                (output_dir / METADATA_FILE_NAME).read_bytes(),
                (scheme_dir / METADATA_FILE_NAME).read_bytes(),
            )
            self.assertEqual(
                (output_dir / PRIMER_FILE_NAME).read_bytes(),
                (scheme_dir / PRIMER_FILE_NAME).read_bytes(),
            )
            self.assertEqual(
                (output_dir / REFERENCE_FILE_NAME).read_bytes(),
                (scheme_dir / REFERENCE_FILE_NAME).read_bytes(),
            )

    def test_download_scheme_sanitise_canonical(self):
        # CANONICAL sanitisation writes serialised model outputs.
        ps = _load_scheme()
        psi = _build_index(ps)
        entry = psi.primerschemes["artic"][400]["v4.1.0"]

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            scheme_dir = _scheme_dir_from_url(
                "https://example.com/schemes/artic/400/v4.1.0/info.json"
            )
            info_obj = PrimerScheme.model_validate_json(
                (scheme_dir / METADATA_FILE_NAME).read_text()
            )
            expected_info = serialize_primer_scheme_json(info_obj)

            primer_text = (scheme_dir / PRIMER_FILE_NAME).read_text()
            scheme_bed = Scheme.from_str(primer_text)
            scheme_bed.sort_bedlines()
            primer_out = scheme_bed.to_str()
            if not primer_out.endswith("\n"):
                primer_out = f"{primer_out}\n"
            expected_primer = primer_out.encode("utf-8")

            fasta_bytes = (scheme_dir / REFERENCE_FILE_NAME).read_bytes()
            with dnaio.open(BytesIO(fasta_bytes), fileformat="fasta") as reader:
                records = list(reader)
            expected_ref = serialize_fasta_records(records)

            entry.checksums.primer_sha256 = _sha256_bytes(expected_primer)
            entry.checksums.reference_sha256 = _sha256_bytes(expected_ref)

            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id="artic/400/v4.1.0",
                    allow_multiple=False,
                    all_schemes=False,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                    sanitisation=SanitisationMode.CANONICAL,
                )

            output_dir = outputs[0]

            self.assertEqual(
                (output_dir / METADATA_FILE_NAME).read_bytes(), expected_info
            )
            self.assertEqual(
                (output_dir / PRIMER_FILE_NAME).read_bytes(), expected_primer
            )
            self.assertEqual(
                (output_dir / REFERENCE_FILE_NAME).read_bytes(), expected_ref
            )

    def test_download_scheme_happy_path(self):
        # Downloads a single scheme successfully with valid data.
        ps = _load_scheme()
        psi = _build_index(ps)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id="artic/400/v4.1.0",
                    allow_multiple=False,
                    all_schemes=False,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                )

            self.assertEqual(len(outputs), 1)
            output_dir = outputs[0]
            self.assertTrue((output_dir / PRIMER_FILE_NAME).exists())
            self.assertTrue((output_dir / REFERENCE_FILE_NAME).exists())

    def test_download_scheme_checksum_mismatch_fails(self):
        # Fails when primer.bed checksum does not match and force is False.
        ps = _load_scheme()
        psi = _build_index(ps)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(primer_mutator=_mutate_primer_bytes),
            ):
                with self.assertRaises(DownloadError):
                    schemes = resolve_schemes(
                        index=psi,
                        scheme_id="artic/400/v4.1.0",
                        allow_multiple=False,
                        all_schemes=False,
                    )
                    download_schemes(
                        schemes=schemes,
                        output=tmp_path,
                    )

    def test_download_scheme_checksum_mismatch_force(self):
        # Allows checksum mismatch when force is True.
        ps = _load_scheme()
        psi = _build_index(ps)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(primer_mutator=_mutate_primer_bytes),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id="artic/400/v4.1.0",
                    allow_multiple=False,
                    all_schemes=False,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                    force=True,
                )

            self.assertEqual(len(outputs), 1)

    def test_download_scheme_missing_checksums_strict(self):
        # Fails when checksums are missing from the index and force is False.
        ps = _load_scheme()
        psi = _build_index(ps)
        entry = psi.primerschemes["artic"][400]["v4.1.0"]
        entry.checksums = None

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                with self.assertRaises(DownloadError):
                    schemes = resolve_schemes(
                        index=psi,
                        scheme_id="artic/400/v4.1.0",
                        allow_multiple=False,
                        all_schemes=False,
                    )
                    download_schemes(
                        schemes=schemes,
                        output=tmp_path,
                    )

    def test_download_scheme_missing_checksums_force(self):
        # Allows missing checksums when force is True.
        ps = _load_scheme()
        psi = _build_index(ps)
        entry = psi.primerschemes["artic"][400]["v4.1.0"]
        entry.checksums = None

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id="artic/400/v4.1.0",
                    allow_multiple=False,
                    all_schemes=False,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                    force=True,
                )

            self.assertEqual(len(outputs), 1)

    def test_download_scheme_allow_multiple(self):
        # Downloads all matches when allow_multiple is True.
        ps = _load_scheme()
        ps2 = ps.model_copy(update={"version": "v4.1.1"})
        psi = _build_index(ps, ps2)

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            with patch(
                "primaschema.get_scheme._download_bytes",
                side_effect=_repo_bytes(),
            ):
                schemes = resolve_schemes(
                    index=psi,
                    scheme_id="artic/400",
                    allow_multiple=True,
                    all_schemes=False,
                )
                outputs = download_schemes(
                    schemes=schemes,
                    output=tmp_path,
                )

            self.assertEqual(len(outputs), 2)
            versions = {p.name for p in outputs}
            self.assertEqual(versions, {"v4.1.0", "v4.1.1"})

    def test_download_scheme_multiple_without_flag(self):
        # Raises when multiple matches are found but allow_multiple is False.
        ps = _load_scheme()
        ps2 = ps.model_copy(update={"version": "v4.1.1"})
        psi = _build_index(ps, ps2)

        with self.assertRaises(ValueError):
            resolve_schemes(
                index=psi,
                scheme_id="artic/400",
                allow_multiple=False,
                all_schemes=False,
            )
