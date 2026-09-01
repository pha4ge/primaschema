import enum
import logging
import tempfile
from pathlib import Path

import linkml.validator
from primalbedtools.bedfiles import BedLineParser, sort_bedlines
from primalbedtools.scheme import Scheme
from primalbedtools.validate import validate_ref_and_bed
from pydantic_core import from_json

from primaschema import (
    METADATA_FILE_NAME,
    PRIMER_FILE_NAME,
    README_FILE_NAME,
    REFERENCE_FILE_NAME,
    SCHEMA_DIR,
)
from primaschema.schema.primer_scheme import (
    PrimerScheme,
    check_primer_scheme_identifier,
)
from primaschema.util import read_fasta_records, sha256_checksum, write_fasta_records

logger = logging.getLogger(__name__)


class ValidationEngine(enum.Enum):
    PYDANTIC = "pydantic"
    LINKML = "linkml"


def validate_scheme_json_with_pydantic(info_path: Path) -> PrimerScheme:
    """Validate an info.json file using the Pydantic model.

    Args:
        info_path: Path to info.json.

    Returns:
        Parsed PrimerScheme instance.
    """
    # Use the derived pydantic model
    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())
    logger.info(f"Validated {info_path} with pydantic")

    # TODO Handle the rules not encoded in base model.

    return primer_scheme


def validate_scheme_json_with_linkml(info_path: Path) -> None:
    """Validate an info.json file using the LinkML schema.

    Args:
        info_path: Path to info.json.

    Raises:
        ValueError: If validation fails.
    """
    # Use the linkml.validator to validate the output
    report = linkml.validator.validate(
        from_json(info_path.read_text()), str(SCHEMA_DIR / "info.yml"), "PrimerScheme"
    )
    if report.results:
        msg = ""
        for result in report.results:
            msg += f"{result.message}\n"
        raise ValueError(msg)


def validate_primer_bed(infopath: Path, strict: bool = False) -> Scheme:
    """Validate the adjacent primer.bed file using primalbedtools.

    Args:
        infopath: Path to info.json metadata file.
        strict: If True, fail on non-canonical primer ordering.

    Returns:
        Parsed Scheme object.

    Raises:
        ValueError: If primer.bed is invalid or order differs in strict mode.
    """

    primer_path = infopath.parent / PRIMER_FILE_NAME
    primer_txt = primer_path.read_text()

    # Round trip the test
    scheme = Scheme.from_str(primer_txt)
    scheme.sort_bedlines()

    if strict and primer_txt != scheme.to_str():
        raise ValueError(
            f"Change detected for {infopath}: {PRIMER_FILE_NAME} has changed order."
        )

    return scheme


def validate_name(infopath: Path, primer_scheme: PrimerScheme | None = None):
    """Validate scheme name, amplicon size, and version against the path.

    Args:
        infopath: Path to info.json.
        primer_scheme: Optional PrimerScheme instance to avoid re-parsing.

    Raises:
        ValueError: If any path component mismatches info.json.
    """
    """
    Validate the schemename, ampliconsize, and schemeversion in the path, README.md, and info.json
    :raises ValueError: If a mismatch is found
    :raises FileNotFoundError: If the README.md does not exist
    """
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    scheme_subpath = (
        Path(primer_scheme.primer_scheme_name)
        / str(primer_scheme.amplicon_size)
        / primer_scheme.primer_scheme_version
    )

    # Check the info version matches path version
    version_path = infopath.parent.name
    if primer_scheme.primer_scheme_version != version_path:
        raise ValueError(
            f"Version mismatch for {scheme_subpath}: info ({primer_scheme.primer_scheme_version}) != path ({version_path})"
        )

    # Check the amplicon size matches the schemepath
    ampliconsize_path = infopath.parent.parent.name
    if str(primer_scheme.amplicon_size) != ampliconsize_path:
        raise ValueError(
            f"Amplicon size mismatch for {scheme_subpath}: info ({primer_scheme.amplicon_size}) != path ({ampliconsize_path})"
        )

    # Check the schemepath matches the path
    schemeid_path = infopath.parent.parent.parent.name
    if primer_scheme.primer_scheme_name != schemeid_path:
        raise ValueError(
            f"Name mismatch for {scheme_subpath}: info ({primer_scheme.primer_scheme_name}) != path ({schemeid_path})"
        )


def validate_identifier(infopath: Path, primer_scheme: PrimerScheme | None = None):
    """Validate a stored primer_scheme_identifier against the other fields.

    primer_scheme_identifier is computed and self-healed by the PrimerScheme
    model itself, so a wrong value in the parsed object would already be
    silently corrected. To catch a genuinely wrong value on disk, this reads
    the file's raw JSON directly instead of relying on the parsed object.

    Args:
        infopath: Path to info.json.
        primer_scheme: Optional PrimerScheme instance to avoid re-parsing.

    Raises:
        ValueError: If the file provides a primer_scheme_identifier that
            disagrees with the one computed from name/amplicon_size/version.
    """
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    raw = from_json(infopath.read_text())
    check_primer_scheme_identifier(
        raw.get("primer_scheme_identifier"),
        primer_scheme.primer_scheme_name,
        primer_scheme.amplicon_size,
        primer_scheme.primer_scheme_version,
        source=str(infopath),
    )


def validate_readme(infopath: Path, primer_scheme: PrimerScheme | None = None):
    """Validate README.md exists and contains scheme metadata.

    Args:
        infopath: Path to info.json.
        primer_scheme: Optional PrimerScheme instance to avoid re-parsing.

    Raises:
        FileNotFoundError: If README.md is missing.
        ValueError: If README.md content does not match info.json.
    """
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    # Check the README.md
    readme = infopath.parent / README_FILE_NAME
    if not readme.exists():
        raise FileNotFoundError(f"{readme} does not exist")

    # Check the readme has been updated
    readme_text = readme.read_text()
    if readme_text.find(primer_scheme.primer_scheme_name) == -1:
        raise ValueError(
            f"Scheme name ({primer_scheme.primer_scheme_name}) not found in {README_FILE_NAME}"
        )
    if readme_text.find(str(primer_scheme.amplicon_size)) == -1:
        raise ValueError(
            f"Amplicon size ({primer_scheme.amplicon_size}) not found in {README_FILE_NAME}"
        )
    if readme_text.find(primer_scheme.primer_scheme_version) == -1:
        raise ValueError(
            f"Scheme version ({primer_scheme.primer_scheme_version}) not found in {README_FILE_NAME}"
        )


def validate_hashes(
    infopath: Path,
    primer_scheme: PrimerScheme | None = None,
    fix: bool = False,
):
    """Validate primer.bed and reference.fasta checksums against info.json.

    Args:
        infopath: Path to info.json.
        primer_scheme: Optional PrimerScheme instance to avoid re-parsing.
        fix: If True, rewrite files when normalised hashes match expected.

    Raises:
        ValueError: If hashes mismatch and cannot be normalised (or fix is False).
    """
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    scheme_subpath = (
        Path(primer_scheme.primer_scheme_name)
        / str(primer_scheme.amplicon_size)
        / primer_scheme.primer_scheme_version
    )

    if not primer_scheme.primer_scheme_checksums:
        logger.warning(f"No checksums found for {scheme_subpath}")
        return

    # Check sha256 hash bedfile
    primer_path = infopath.parent / PRIMER_FILE_NAME
    primer_sha = sha256_checksum(primer_path)
    if primer_sha != primer_scheme.primer_scheme_checksums.primer_scheme_sha256:
        logger.warning(
            f"primer.bed sha256 mismatch for {scheme_subpath}. Attempting to normalise and recheck."
        )
        reformatted_sha = None
        header = None
        bedlines = None
        try:
            header, bedlines = BedLineParser.from_file(primer_path)
            bedlines = sort_bedlines(bedlines)
            with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
                tmp_path = Path(tmp.name)
            try:
                BedLineParser.to_file(tmp_path, header, bedlines)
                reformatted_sha = sha256_checksum(tmp_path)
            finally:
                tmp_path.unlink(missing_ok=True)
        except Exception as exc:
            logger.warning(
                f"Failed to normalise primer.bed for {scheme_subpath}: {exc}"
            )

        if (
            reformatted_sha
            == primer_scheme.primer_scheme_checksums.primer_scheme_sha256
        ):
            if fix:
                BedLineParser.to_file(primer_path, header, bedlines)
                primer_sha = reformatted_sha
                logger.warning(
                    f"primer.bed reformatted to match expected sha256 for {scheme_subpath}."
                )
            else:
                raise ValueError(
                    f"{PRIMER_FILE_NAME} sha256 mismatch for {scheme_subpath} (normalisable; rerun with --fix)"
                )
        else:
            raise ValueError(
                f"{PRIMER_FILE_NAME} sha256 ({primer_sha} != info sha256 ({primer_scheme.primer_scheme_checksums.primer_scheme_sha256}): {scheme_subpath}"
            )

    # Check sha256 hash ref
    reference_path = infopath.parent / REFERENCE_FILE_NAME
    reference_sha = sha256_checksum(reference_path)
    if reference_sha != primer_scheme.primer_scheme_checksums.reference_sequence_sha256:
        logger.warning(
            f"reference.fasta sha256 mismatch for {scheme_subpath}. Attempting to normalise and recheck."
        )
        reformatted_sha = None
        try:
            reference_records = read_fasta_records(reference_path)
            with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
                tmp_path = Path(tmp.name)
            try:
                write_fasta_records(tmp_path, reference_records)
                reformatted_sha = sha256_checksum(tmp_path)
            finally:
                tmp_path.unlink(missing_ok=True)
        except Exception as exc:
            logger.warning(
                f"Failed to normalise reference.fasta for {scheme_subpath}: {exc}"
            )

        if (
            reformatted_sha
            == primer_scheme.primer_scheme_checksums.reference_sequence_sha256
        ):
            if fix:
                write_fasta_records(reference_path, reference_records)
                reference_sha = reformatted_sha
                logger.warning(
                    f"reference.fasta reformatted to match expected sha256 for {scheme_subpath}."
                )
            else:
                raise ValueError(
                    f"{REFERENCE_FILE_NAME} sha256 mismatch for {scheme_subpath} (normalisable; rerun with --fix)"
                )
        else:
            raise ValueError(
                f"{REFERENCE_FILE_NAME} sha256 ({reference_sha}) != info sha256 ({primer_scheme.primer_scheme_checksums.reference_sequence_sha256}): {scheme_subpath}"
            )


def validate(
    infopath: Path,
    primer_scheme: PrimerScheme | None = None,
    additional_linkml: bool = False,
    strict: bool = False,
    fix: bool = False,
):
    """Validate a scheme directory (info.json, primer.bed, reference.fasta).

    Args:
        infopath: Path to info.json.
        primer_scheme: Optional PrimerScheme instance to avoid re-parsing.
        additional_linkml: If True, run LinkML validation.
        strict: If True, enforce strict primer.bed ordering.
        fix: If True, rewrite files when normalised hashes match expected.

    Raises:
        ValueError: If any validation step fails.
    """
    logger.debug(f"Validating {'strict' if strict else ''} {infopath}")
    if additional_linkml:
        logger.debug(f"Validated with LinkML: {infopath}")
        validate_scheme_json_with_linkml(infopath)

    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())
    validate_name(infopath, primer_scheme)
    validate_identifier(infopath, primer_scheme)
    logger.debug(f"Validated with Pydantic: {infopath}")

    # Validate primer + ref
    try:
        scheme = validate_primer_bed(infopath, strict)
    except ValueError as exc:
        if strict and "primer.bed has changed order" in str(exc):
            logger.warning(
                f"primer.bed order differs for {infopath}; attempting normalisation if hashes mismatch."
            )
            scheme = validate_primer_bed(infopath, strict=False)
        else:
            raise
    validate_ref_and_bed(scheme.bedlines, str(infopath.parent / REFERENCE_FILE_NAME))
    logger.debug(f"Validated primer.bed files:  {infopath}")

    # Validate hashes
    validate_hashes(infopath, primer_scheme, fix=fix)
    validate_readme(infopath, primer_scheme)
    logger.debug(f"Validated hashes and README:  {infopath}")


def validate_all(
    primer_schemes_path: Path, additional_linkml: bool = False, strict: bool = True
):
    """Validate all schemes under a root directory.

    Args:
        primer_schemes_path: Root path to search for info.json files.
        additional_linkml: If True, run LinkML validation.
        strict: If True, enforce strict primer.bed ordering.

    Raises:
        ValueError: If any scheme fails validation.
    """

    errors: list[str] = []
    for schemeinfo in primer_schemes_path.rglob(f"*/{METADATA_FILE_NAME}"):
        try:
            validate(schemeinfo, None, additional_linkml, strict)
        except Exception as exc:
            logger.error(f"Validation failed for {schemeinfo}: {exc}")
            errors.append(f"{schemeinfo}: {exc}")
    if errors:
        raise ValueError(
            f"Validation failed for {len(errors)} scheme(s):\n" + "\n".join(errors)
        )
