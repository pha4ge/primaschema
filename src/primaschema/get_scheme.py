import gzip
import hashlib
import logging
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from enum import Enum
from io import BytesIO
from pathlib import Path
from urllib.parse import urlparse

import dnaio
import httpx
from primalbedtools.scheme import Scheme
from primalbedtools.validate import validate_ref_and_bed

from primaschema import (
    DEFAULT_INDEX_URL,
    METADATA_FILE_NAME,
    PRIMER_FILE_NAME,
    REFERENCE_FILE_NAME,
)
from primaschema.schema.index import IndexPrimerScheme, PrimerSchemeIndex
from primaschema.schema.info import PrimerScheme
from primaschema.util import serialize_fasta_records, serialize_primer_scheme_json

logger = logging.getLogger(__name__)

MAX_DOWNLOAD_MB = int(os.getenv("PRIMASCHEMA_MAX_DOWNLOAD_MB", "100"))
MAX_DOWNLOAD_BYTES = MAX_DOWNLOAD_MB * 1024 * 1024
MAX_DOWNLOAD_WORKERS = int(os.getenv("PRIMASCHEMA_DOWNLOAD_WORKERS", "4"))


class DownloadError(RuntimeError):
    pass


class SanitisationMode(str, Enum):
    RAW = "raw"
    CANONICAL = "canonical"


def _ensure_https(url: str) -> None:
    parsed = urlparse(url)
    if parsed.scheme != "https":
        raise ValueError(f"Only https URLs are allowed: {url}")


def _read_index_bytes(source: str) -> bytes:
    if source.startswith("http"):
        logger.info(f"Fetching index from URL {source}")
        return _download_bytes(source)

    index_path = Path(source)
    logger.info(f"Loading index from path {index_path}")
    return index_path.read_bytes()


def _download_bytes(url: str) -> bytes:
    _ensure_https(url)
    with httpx.stream("GET", url, follow_redirects=True) as resp:
        resp.raise_for_status()
        content_length = resp.headers.get("Content-Length")
        if content_length is not None:
            try:
                if int(content_length) > MAX_DOWNLOAD_BYTES:
                    raise DownloadError(
                        f"{url} exceeds max size of {MAX_DOWNLOAD_MB} MB"
                    )
            except ValueError:
                pass
        total = 0
        chunks: list[bytes] = []
        for chunk in resp.iter_bytes():
            total += len(chunk)
            if total > MAX_DOWNLOAD_BYTES:
                raise DownloadError(f"{url} exceeds max size of {MAX_DOWNLOAD_MB} MB")
            chunks.append(chunk)
    return b"".join(chunks)


def _get_and_validate_info_json(url: str) -> tuple[PrimerScheme, bytes]:
    if not url:
        raise DownloadError("Missing URL for info.json in index")
    info_bytes = _download_bytes(url)
    try:
        scheme = PrimerScheme.model_validate_json(info_bytes)
    except Exception as exc:
        raise DownloadError(f"Invalid info.json: {exc}") from exc
    return scheme, info_bytes


def _verify_checksum(path: Path, expected: str, force: bool) -> None:
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        if force:
            logger.warning(
                f"Checksum mismatch for {path.name}: expected {expected}, got {actual}; continuing due to --force"
            )
            return
        raise DownloadError(
            f"Checksum mismatch for {path.name}: expected {expected}, got {actual}"
        )


def _validate_primer_bed_bytes(bed_bytes: bytes) -> Scheme:
    try:
        scheme = Scheme.from_str(bed_bytes.decode("utf-8"))
        scheme.sort_bedlines()
        return scheme
    except Exception as exc:
        raise DownloadError(f"Invalid primer.bed: {exc}") from exc


def _get_and_validate_primer_bed(url: str) -> tuple[Scheme, bytes]:
    if not url:
        raise DownloadError("Missing URL for primer.bed in index")
    bed_bytes = _download_bytes(url)
    scheme = _validate_primer_bed_bytes(bed_bytes)
    return scheme, bed_bytes


def _get_and_validate_reference_fasta(
    url: str,
) -> tuple[list[dnaio.SequenceRecord], bytes]:
    if not url:
        raise DownloadError("Missing URL for reference.fasta in index")
    fasta_bytes = _download_bytes(url)
    records: list[dnaio.SequenceRecord] = []
    try:
        count = 0
        with dnaio.open(BytesIO(fasta_bytes), fileformat="fasta") as reader:
            for record in reader:
                count += 1
                if not record.sequence:
                    raise DownloadError(
                        f"reference.fasta contains empty sequence for {record.id}"
                    )
                records.append(record)
        if count == 0:
            raise DownloadError("reference.fasta contains no records")
    except DownloadError:
        raise
    except Exception as exc:
        raise DownloadError(f"Invalid reference.fasta: {exc}") from exc
    return records, fasta_bytes


def load_index(source: str | None = None) -> PrimerSchemeIndex:
    if source is None:
        source = DEFAULT_INDEX_URL
    raw = _read_index_bytes(source)
    data = gzip.decompress(raw) if source.endswith(".gz") else raw
    return PrimerSchemeIndex.model_validate_json(data)


def _require_checksums(index: IndexPrimerScheme, force: bool) -> None:
    if index.checksums is None:
        if force:
            logger.warning("Missing checksums in index; continuing due to --force")
            return
        raise DownloadError("Checksums are required in the index")
    if not index.checksums.primer_sha256 or not index.checksums.reference_sha256:
        if force:
            logger.warning("Incomplete checksums in index; continuing due to --force")
            return
        raise DownloadError("Checksums are required in the index")


def _resolve_schemes(
    index: PrimerSchemeIndex, scheme_id: str, allow_multiple: bool
) -> list[IndexPrimerScheme]:
    parts = scheme_id.strip("/").split("/")
    if len(parts) not in (1, 2, 3):
        raise ValueError(
            f"Invalid scheme_id '{scheme_id}': expected format name[/amplicon_size[/version]]"
        )
    name = parts[0]
    amplicon_size = parts[1] if len(parts) > 1 else None
    version = parts[2] if len(parts) > 2 else None
    schemes = index.get_schemes_from_index(name, amplicon_size, version)
    if len(parts) < 3 and not allow_multiple and len(schemes) != 1:
        raise ValueError(
            f"Scheme id '{scheme_id}' matched {len(schemes)} schemes; rerun with --allow-multiple to download all"
        )
    return schemes


def resolve_schemes(
    index: PrimerSchemeIndex,
    scheme_id: str | None,
    allow_multiple: bool,
    all_schemes: bool,
) -> list[IndexPrimerScheme]:
    if all_schemes:
        schemes = index.flatten()
        if not schemes:
            raise ValueError("No schemes found in the index")
        return schemes

    if not scheme_id:
        raise ValueError("scheme_id is required unless --all is set")
    schemes = _resolve_schemes(index, scheme_id, allow_multiple=allow_multiple)
    if not schemes:
        raise ValueError(f"Scheme ({scheme_id}) not found in the index")

    if not allow_multiple and len(schemes) > 1:
        raise ValueError(
            f"Scheme id '{scheme_id}' matched {len(schemes)} schemes; rerun with --allow-multiple to download all"
        )
    return schemes


def _download_scheme_entry(
    scheme: IndexPrimerScheme,
    output: Path,
    strict: bool,
    force: bool,
    sanitisation: SanitisationMode,
    requested_id: str | None = None,
) -> Path:
    if requested_id and scheme.relative_path != requested_id:
        msg = f"Index entry path {scheme.relative_path} does not match requested {requested_id}"
        if strict:
            raise DownloadError(msg)
        logger.warning(msg)

    _require_checksums(scheme, force=force)

    output_dir = output / scheme.name / str(scheme.amplicon_size) / scheme.version
    if output_dir.exists():
        msg = f"Output directory already exists: {output_dir}"
        if strict:
            raise DownloadError(msg)
        logger.warning(msg)

    tmp_dir_path = None
    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir_path = Path(tmp_dir)
            tmp_scheme_dir = (
                tmp_dir_path / scheme.name / str(scheme.amplicon_size) / scheme.version
            )
            tmp_scheme_dir.mkdir(parents=True, exist_ok=True)

            info_scheme, info_bytes = _get_and_validate_info_json(scheme.info_file_url)
            if sanitisation == SanitisationMode.RAW:
                (tmp_scheme_dir / METADATA_FILE_NAME).write_bytes(info_bytes)
            else:
                (tmp_scheme_dir / METADATA_FILE_NAME).write_bytes(
                    serialize_primer_scheme_json(info_scheme)
                )

            scheme_bed, bed_bytes = _get_and_validate_primer_bed(scheme.primer_file_url)
            if sanitisation == SanitisationMode.RAW:
                (tmp_scheme_dir / PRIMER_FILE_NAME).write_bytes(bed_bytes)
            else:
                (tmp_scheme_dir / PRIMER_FILE_NAME).write_text(scheme_bed.to_str())

            ref_records, ref_bytes = _get_and_validate_reference_fasta(
                scheme.reference_file_url
            )
            if sanitisation == SanitisationMode.RAW:
                (tmp_scheme_dir / REFERENCE_FILE_NAME).write_bytes(ref_bytes)
            else:
                (tmp_scheme_dir / REFERENCE_FILE_NAME).write_bytes(
                    serialize_fasta_records(ref_records)
                )

            validate_ref_and_bed(
                scheme_bed.bedlines, str(tmp_scheme_dir / REFERENCE_FILE_NAME)
            )

            if scheme.checksums:
                if scheme.checksums.primer_sha256:
                    _verify_checksum(
                        tmp_scheme_dir / PRIMER_FILE_NAME,
                        scheme.checksums.primer_sha256,
                        force=force,
                    )
                if scheme.checksums.reference_sha256:
                    _verify_checksum(
                        tmp_scheme_dir / REFERENCE_FILE_NAME,
                        scheme.checksums.reference_sha256,
                        force=force,
                    )

            output_dir.mkdir(parents=True, exist_ok=True)
            for filename in (METADATA_FILE_NAME, PRIMER_FILE_NAME, REFERENCE_FILE_NAME):
                shutil.copy2(tmp_scheme_dir / filename, output_dir / filename)
    except Exception:
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Cleaning up temporary download dir {tmp_dir_path}")
        raise

    logger.info(f"Downloaded scheme {scheme.relative_path} to {output_dir}")
    return output_dir


def download_schemes(
    schemes: list[IndexPrimerScheme],
    output: Path = Path("."),
    strict: bool = False,
    force: bool = False,
    sanitisation: SanitisationMode = SanitisationMode.RAW,
    workers: int | None = None,
) -> list[Path]:
    if len(schemes) == 1:
        return [
            _download_scheme_entry(
                schemes[0],
                output=output,
                strict=strict,
                force=force,
                sanitisation=sanitisation,
                requested_id=schemes[0].relative_path,
            )
        ]

    max_workers = workers or MAX_DOWNLOAD_WORKERS
    results: list[Path] = []
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(
                _download_scheme_entry,
                scheme,
                output,
                strict,
                force,
                sanitisation,
                scheme.relative_path,
            ): scheme.relative_path
            for scheme in schemes
        }
        for future in as_completed(future_map):
            scheme_path = future_map[future]
            try:
                results.append(future.result())
            except Exception as exc:
                errors.append(f"{scheme_path}: {exc}")

    if errors:
        raise DownloadError(
            f"Failed to download {len(errors)} scheme(s): {', '.join(errors)}"
        )
    return results
