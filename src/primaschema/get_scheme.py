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
DEFAULT_MAX_WORKERS = 4
MAX_DOWNLOAD_WORKERS = os.getenv(
    "PRIMASCHEMA_DOWNLOAD_WORKERS", str(DEFAULT_MAX_WORKERS)
)
DEFAULT_HTTP_TIMEOUT_SECONDS = float(os.getenv("PRIMASCHEMA_HTTP_TIMEOUT", "30"))
MIN_DOWNLOAD_WORKERS = 1
MAX_DOWNLOAD_WORKERS_LIMIT = 32


class DownloadError(RuntimeError):
    pass


class SanitisationMode(str, Enum):
    RAW = "raw"
    CANONICAL = "canonical"


def _ensure_https(url: str) -> None:
    """Ensure a URL uses HTTPS.

    Args:
        url: URL to validate.

    Raises:
        ValueError: If the URL scheme is not HTTPS.
    """
    parsed = urlparse(url)
    if parsed.scheme != "https":
        raise ValueError(f"Only https URLs are allowed: {url}")


def _read_index_bytes(source: str) -> bytes:
    """Read index bytes from a URL or local path.

    Args:
        source: URL or filesystem path for the index.

    Returns:
        Raw bytes of the index file.
    """
    if source.startswith("http"):
        logger.info(f"Fetching index from URL {source}")
        return _download_bytes(source)

    index_path = Path(source)
    logger.info(f"Loading index from path {index_path}")
    return index_path.read_bytes()


def _resolve_timeout(timeout: float | None) -> httpx.Timeout:
    """Resolve an HTTP timeout value to an httpx.Timeout.

    Args:
        timeout: Timeout in seconds, or None to use default.

    Returns:
        An httpx.Timeout instance.
    """
    if timeout is None:
        timeout = DEFAULT_HTTP_TIMEOUT_SECONDS
    return httpx.Timeout(timeout)


def _resolve_workers(workers: int | None) -> int:
    """Clamp worker count to a safe range.

    Args:
        workers: Requested worker count, or None to use default.

    Returns:
        A validated worker count within the allowed range.
    """
    if workers is None:
        try:
            workers = int(MAX_DOWNLOAD_WORKERS)
        except ValueError:
            logger.warning(
                f"Invalid PRIMASCHEMA_DOWNLOAD_WORKERS value '{MAX_DOWNLOAD_WORKERS}', using default {DEFAULT_MAX_WORKERS}"
            )
            workers = DEFAULT_MAX_WORKERS

    if workers < MIN_DOWNLOAD_WORKERS:
        logger.warning(
            f"Invalid worker count {workers}, using minimum {MIN_DOWNLOAD_WORKERS}"
        )
        return MIN_DOWNLOAD_WORKERS
    if workers > MAX_DOWNLOAD_WORKERS_LIMIT:
        logger.warning(
            f"Worker count {workers} exceeds limit {MAX_DOWNLOAD_WORKERS_LIMIT}, clamping"
        )
        return MAX_DOWNLOAD_WORKERS_LIMIT
    return workers


def _download_bytes(url: str, *, timeout: float | None = None) -> bytes:
    """Download bytes from a URL with size and timeout limits.

    Args:
        url: HTTPS URL to download.
        timeout: Timeout in seconds, or None to use default.

    Returns:
        Downloaded bytes.

    Raises:
        ValueError: If the URL scheme is not HTTPS.
        DownloadError: If the download exceeds size limits.
        httpx.HTTPError: If the request fails.
    """
    _ensure_https(url)
    with httpx.stream(
        "GET",
        url,
        follow_redirects=True,
        timeout=_resolve_timeout(timeout),
    ) as resp:
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


def _get_and_validate_info_json(
    url: str, *, timeout: float | None = None
) -> tuple[PrimerScheme, bytes]:
    """Fetch and validate info.json.

    Args:
        url: HTTPS URL to info.json.
        timeout: Timeout in seconds, or None to use default.

    Returns:
        A tuple of (PrimerScheme, raw bytes).

    Raises:
        DownloadError: If the URL is missing or JSON is invalid.
    """
    if not url:
        raise DownloadError("Missing URL for info.json in index")
    info_bytes = _download_bytes(url, timeout=timeout)
    try:
        scheme = PrimerScheme.model_validate_json(info_bytes)
    except Exception as exc:
        raise DownloadError(f"Invalid info.json: {exc}") from exc
    return scheme, info_bytes


def _verify_checksum(path: Path, expected: str, force: bool) -> None:
    """Verify a file's SHA256 checksum.

    Args:
        path: File path to hash.
        expected: Expected hex digest.
        force: If True, log warning instead of raising on mismatch.

    Raises:
        DownloadError: If checksum mismatches and force is False.
    """
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
    """Parse and validate primer.bed bytes into a Scheme.

    Args:
        bed_bytes: Raw primer.bed bytes.

    Returns:
        Parsed Scheme with sorted bedlines.

    Raises:
        DownloadError: If parsing fails.
    """
    try:
        scheme = Scheme.from_str(bed_bytes.decode("utf-8"))
        scheme.sort_bedlines()
        return scheme
    except Exception as exc:
        raise DownloadError(f"Invalid primer.bed: {exc}") from exc


def _get_and_validate_primer_bed(
    url: str, *, timeout: float | None = None
) -> tuple[Scheme, bytes]:
    """Fetch and validate primer.bed.

    Args:
        url: HTTPS URL to primer.bed.
        timeout: Timeout in seconds, or None to use default.

    Returns:
        A tuple of (Scheme, raw bytes).

    Raises:
        DownloadError: If URL is missing or file is invalid.
    """
    if not url:
        raise DownloadError("Missing URL for primer.bed in index")
    bed_bytes = _download_bytes(url, timeout=timeout)
    scheme = _validate_primer_bed_bytes(bed_bytes)
    return scheme, bed_bytes


def _get_and_validate_reference_fasta(
    url: str, *, timeout: float | None = None
) -> tuple[list[dnaio.SequenceRecord], bytes]:
    """Fetch and validate reference.fasta.

    Args:
        url: HTTPS URL to reference.fasta.
        timeout: Timeout in seconds, or None to use default.

    Returns:
        A tuple of (records, raw bytes).

    Raises:
        DownloadError: If URL is missing or FASTA is invalid.
    """
    if not url:
        raise DownloadError("Missing URL for reference.fasta in index")
    fasta_bytes = _download_bytes(url, timeout=timeout)
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


def load_index(source: str, *, timeout: float | None = None) -> PrimerSchemeIndex:
    """Load a PrimerSchemeIndex from a URL or local file.

    Args:
        source: URL or path to index JSON (optionally gzipped).
        timeout: Timeout in seconds for URL fetches.

    Returns:
        Parsed PrimerSchemeIndex.
    """
    if source.startswith("http"):
        raw = _download_bytes(source, timeout=timeout)
    else:
        raw = _read_index_bytes(source)
    data = gzip.decompress(raw) if source.endswith(".gz") else raw
    return PrimerSchemeIndex.model_validate_json(data)


def _require_checksums(index: IndexPrimerScheme, force: bool) -> None:
    """Ensure checksums exist in the index entry.

    Args:
        index: Index entry to inspect.
        force: If True, allow missing checksums with warning.

    Raises:
        DownloadError: If checksums are missing and force is False.
    """
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
    """Resolve a scheme identifier to index entries.

    Args:
        index: PrimerSchemeIndex to query.
        scheme_id: Scheme id in name[/amplicon_size[/version]] format.
        allow_multiple: Whether multiple matches are allowed.

    Returns:
        Matching IndexPrimerScheme entries.

    Raises:
        ValueError: If scheme_id format is invalid or matches are ambiguous.
    """
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
    """Resolve schemes for a download request.

    Args:
        index: PrimerSchemeIndex to query.
        scheme_id: Scheme id or None when all_schemes is True.
        allow_multiple: Whether multiple matches are allowed.
        all_schemes: If True, return all schemes from the index.

    Returns:
        List of IndexPrimerScheme entries.

    Raises:
        ValueError: If input is invalid or no schemes are found.
    """
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
    timeout: float | None,
    check_output_exists: bool,
    suppress_log: bool,
    requested_id: str | None = None,
) -> Path:
    """Download and validate a single scheme into a target root directory.

    Args:
        scheme: Index entry to download.
        output: Root output directory.
        strict: Fail on mismatches or existing outputs.
        force: Allow missing or mismatched checksums.
        sanitisation: RAW or CANONICAL output mode.
        timeout: Timeout in seconds for downloads.
        check_output_exists: Whether to check existing output dir.
        requested_id: Optional requested scheme id for validation.

    Returns:
        Path to the downloaded scheme directory.

    Raises:
        DownloadError: If validation or download fails.
    """
    if requested_id and scheme.relative_path != requested_id:
        msg = f"Index entry path {scheme.relative_path} does not match requested {requested_id}"
        if strict:
            raise DownloadError(msg)
        logger.warning(msg)

    _require_checksums(scheme, force=force)

    output_dir = output / scheme.name / str(scheme.amplicon_size) / scheme.version
    if check_output_exists and output_dir.exists():
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

            info_scheme, info_bytes = _get_and_validate_info_json(
                scheme.info_file_url,
                timeout=timeout,
            )
            if sanitisation == SanitisationMode.RAW:
                (tmp_scheme_dir / METADATA_FILE_NAME).write_bytes(info_bytes)
            else:
                (tmp_scheme_dir / METADATA_FILE_NAME).write_bytes(
                    serialize_primer_scheme_json(info_scheme)
                )

            scheme_bed, bed_bytes = _get_and_validate_primer_bed(
                scheme.primer_file_url,
                timeout=timeout,
            )
            if sanitisation == SanitisationMode.RAW:
                (tmp_scheme_dir / PRIMER_FILE_NAME).write_bytes(bed_bytes)
            else:
                (tmp_scheme_dir / PRIMER_FILE_NAME).write_text(scheme_bed.to_str())

            ref_records, ref_bytes = _get_and_validate_reference_fasta(
                scheme.reference_file_url,
                timeout=timeout,
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

    logger.info(f"Downloaded scheme {scheme.relative_path}")
    return output_dir


def download_schemes(
    schemes: list[IndexPrimerScheme],
    output: Path = Path("."),
    strict: bool = False,
    force: bool = False,
    sanitisation: SanitisationMode = SanitisationMode.RAW,
    timeout: float | None = None,
    workers: int | None = None,
) -> list[Path]:
    """Download multiple schemes with all-or-nothing semantics.

    Schemes are downloaded into a staging directory first. If all succeed,
    outputs are copied to the final destination. If any fail, nothing is written.

    Args:
        schemes: List of schemes to download.
        output: Root output directory.
        strict: Fail on mismatches or existing outputs.
        force: Allow missing or mismatched checksums.
        sanitisation: RAW or CANONICAL output mode.
        timeout: Timeout in seconds for downloads.
        workers: Max worker count for parallel downloads.

    Returns:
        List of final output directories.

    Raises:
        DownloadError: If any scheme fails to download or validate.
        ValueError: If no schemes are provided.
    """
    if not schemes:
        raise ValueError("No schemes provided")

    for scheme in schemes:
        output_dir = output / scheme.name / str(scheme.amplicon_size) / scheme.version
        if output_dir.exists():
            msg = f"Output directory already exists: {output_dir}"
            if strict:
                raise DownloadError(msg)
            logger.warning(msg)

    max_workers = _resolve_workers(workers)
    results: list[Path] = []
    errors: list[str] = []

    with tempfile.TemporaryDirectory() as staging_dir:
        staging_root = Path(staging_dir)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(
                    _download_scheme_entry,
                    scheme,
                    staging_root,
                    strict,
                    force,
                    sanitisation,
                    timeout,
                    False,
                    True,
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

        final_outputs: list[Path] = []
        for scheme in schemes:
            source_dir = (
                staging_root / scheme.name / str(scheme.amplicon_size) / scheme.version
            )
            dest_dir = output / scheme.name / str(scheme.amplicon_size) / scheme.version
            dest_dir.mkdir(parents=True, exist_ok=True)
            for filename in (METADATA_FILE_NAME, PRIMER_FILE_NAME, REFERENCE_FILE_NAME):
                shutil.copy2(source_dir / filename, dest_dir / filename)
            final_outputs.append(dest_dir)

        logger.info(f"Downloaded {len(final_outputs)} scheme(s) to {output}")
        return final_outputs
