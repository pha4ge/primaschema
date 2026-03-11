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
    REFERENCE_FILE_NAME,
    SCHEMA_DIR,
)
from primaschema.schema.info import PrimerScheme
from primaschema.util import read_fasta_records, sha256_checksum, write_fasta_records

logger = logging.getLogger(__name__)


class ValidationEngine(enum.Enum):
    PYDANTIC = "pydantic"
    LINKML = "linkml"


def validate_scheme_json_with_pydantic(info_path: Path) -> PrimerScheme:
    # Use the derived pydantic model
    primer_scheme = PrimerScheme.model_validate_json(info_path.read_text())
    logger.info(f"Validated {info_path} with pydantic")

    # TODO Handle the rules not encoded in base model.

    return primer_scheme


def validate_scheme_json_with_linkml(info_path: Path) -> None:
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
    """
    Parses the adjacent primer.bed file with PBT

    :param infopath: The path to the info.json metadata file
    :type infopath: Path
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
    """
    Validate the schemename, ampliconsize, and schemeversion in the path, README.md, and info.json
    :raises ValueError: If a mismatch is found
    :raises FileNotFoundError: If the README.md does not exist
    """
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    scheme_subpath = (
        Path(primer_scheme.name)
        / str(primer_scheme.amplicon_size)
        / primer_scheme.version
    )

    # Check the info version matches path version
    version_path = infopath.parent.name
    if primer_scheme.version != version_path:
        raise ValueError(
            f"Version mismatch for {scheme_subpath}: info ({primer_scheme.version}) != path ({version_path})"
        )

    # Check the amplicon size matches the schemepath
    ampliconsize_path = infopath.parent.parent.name
    if str(primer_scheme.amplicon_size) != ampliconsize_path:
        raise ValueError(
            f"Amplicon size mismatch for {scheme_subpath}: info ({primer_scheme.amplicon_size}) != path ({ampliconsize_path})"
        )

    # Check the schemepath matches the path
    schemeid_path = infopath.parent.parent.parent.name
    if primer_scheme.name != schemeid_path:
        raise ValueError(
            f"Name mismatch for {scheme_subpath}: info ({primer_scheme.name}) != path ({schemeid_path})"
        )


def validate_readme(infopath: Path, primer_scheme: PrimerScheme | None = None):
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    scheme_subpath = (
        Path(primer_scheme.name)
        / str(primer_scheme.amplicon_size)
        / primer_scheme.version
    )

    # Check the ReadME.md
    readme = infopath.parent / "README.md"
    if not readme.exists():
        raise FileNotFoundError(f"{readme} does not exist")

    # Check the readme has been updated
    readme = readme.read_text()
    if readme.find(primer_scheme.name) == -1:
        raise ValueError(
            f"Scheme name ({primer_scheme.name}) not found in {readme}: {scheme_subpath}"
        )
    if readme.find(str(primer_scheme.amplicon_size)) == -1:
        raise ValueError(
            f"Amplicon size ({primer_scheme.amplicon_size}) not found in {readme}: {scheme_subpath}"
        )
    if readme.find(primer_scheme.version) == -1:
        raise ValueError(
            f"Scheme version ({primer_scheme.version}) not found in {readme}: {scheme_subpath}"
        )


def validate_hashes(
    infopath: Path,
    primer_scheme: PrimerScheme | None = None,
    edit_inplace: bool = False,
):
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())

    scheme_subpath = (
        Path(primer_scheme.name)
        / str(primer_scheme.amplicon_size)
        / primer_scheme.version
    )

    # Check sha256 hash bedfile
    primer_path = infopath.parent / PRIMER_FILE_NAME
    primer_sha = sha256_checksum(primer_path)
    if primer_sha != primer_scheme.primer_file_sha256:
        logger.warning(
            f"primer.bed sha256 mismatch for {scheme_subpath}. Attempting to normalize and recheck."
        )
        reformatted_sha = None
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
                f"Failed to normalize primer.bed for {scheme_subpath}: {exc}"
            )

        if reformatted_sha == primer_scheme.primer_file_sha256:
            if edit_inplace:
                BedLineParser.to_file(primer_path, header, bedlines)
            primer_sha = reformatted_sha
            logger.warning(
                f"primer.bed reformatted to match expected sha256 for {scheme_subpath}."
            )
        else:
            raise ValueError(
                f"{PRIMER_FILE_NAME} sha256 ({primer_sha} != info sha256 ({primer_scheme.primer_file_sha256}): {scheme_subpath}"
            )

    # Check sha256 hash ref
    reference_path = infopath.parent / REFERENCE_FILE_NAME
    reference_sha = sha256_checksum(reference_path)
    if reference_sha != primer_scheme.reference_file_sha256:
        logger.warning(
            f"reference.fasta sha256 mismatch for {scheme_subpath}. Attempting to normalize and recheck."
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
                f"Failed to normalize reference.fasta for {scheme_subpath}: {exc}"
            )

        if reformatted_sha == primer_scheme.reference_file_sha256:
            if edit_inplace:
                write_fasta_records(reference_path, reference_records)
            reference_sha = reformatted_sha
            logger.warning(
                f"reference.fasta reformatted to match expected sha256 for {scheme_subpath}."
            )
        else:
            raise ValueError(
                f"{REFERENCE_FILE_NAME} sha256 ({reference_sha} != info sha256 ({primer_scheme.reference_file_sha256}): {scheme_subpath}"
            )

    # TODO Check primaschema hashes.


def validate(
    infopath: Path,
    primer_scheme: PrimerScheme | None = None,
    additional_linkml: bool = False,
    strict: bool = False,
    edit_inplace: bool = False,
):
    logger.debug(f"Validating {'strict' if strict else ''} {infopath}")
    if additional_linkml:
        logger.debug(f"Validated with LinkML: {infopath}")
        validate_scheme_json_with_linkml(infopath)

    if primer_scheme is None:
        primer_scheme = PrimerScheme.model_validate_json(infopath.read_text())
    validate_name(infopath, primer_scheme)
    logger.debug(f"Validated with Pydantic: {infopath}")

    # Validate primer + ref
    try:
        scheme = validate_primer_bed(infopath, strict)
    except ValueError as exc:
        if strict and "primer.bed has changed order" in str(exc):
            logger.warning(
                f"primer.bed order differs for {infopath}; attempting normalization if hashes mismatch."
            )
            scheme = validate_primer_bed(infopath, strict=False)
        else:
            raise
    validate_ref_and_bed(scheme.bedlines, str(infopath.parent / REFERENCE_FILE_NAME))
    logger.debug(f"Validated primer.bed files:  {infopath}")

    # Validate hashes
    validate_hashes(infopath, primer_scheme, edit_inplace=edit_inplace)
    validate_readme(infopath, primer_scheme)
    logger.debug(f"Validated hashes and README:  {infopath}")


def validate_all(
    primer_schemes_path: Path, additional_linkml: bool = False, strict: bool = True
):
    """
    Recursively searches through the primer_schemes_path for {METADATA_FILE_NAME}
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
