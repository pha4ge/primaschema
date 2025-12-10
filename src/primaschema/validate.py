from primaschema.schema.info import PrimerScheme
from pathlib import Path
import linkml.validator

from primaschema import SCHEMA_DIR
from pydantic_core import from_json


from primaschema import (
    METADATA_FILE_NAME,
    REFERENCE_FILE_NAME,
    PRIMER_FILE_NAME,
    logger,
)
from primaschema.util import sha256_checksum
import enum


class ValidationEngine(enum.Enum):
    PYDANTIC = "pydantic"
    LINKML = "linkml"


def validate_scheme_json_with_pydantic(info_path: Path) -> PrimerScheme:
    # Use the derived pydantic model
    primerscheme = PrimerScheme.model_validate_json(info_path.read_text())

    # TODO Handle the rules not encoded in base model.

    return primerscheme


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


def validate_name(infopath: Path, primerscheme: PrimerScheme | None = None):
    """
    Validate the schemename, ampliconsize, and schemeversion in the path, README.md, and info.json
    :raises ValueError: If a mismatch is found
    :raises FileNotFoundError: If the README.md does not exist
    """
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primerscheme is None:
        primerscheme = PrimerScheme.model_validate_json(infopath.read_text())

    primerscheme_derived_scheme_path = (
        Path(primerscheme.name) / str(primerscheme.amplicon_size) / primerscheme.version
    )

    # Check the info version matches path version
    version_path = infopath.parent.name
    if primerscheme.version != version_path:
        raise ValueError(
            f"Version mismatch for {primerscheme_derived_scheme_path}: info ({primerscheme.version}) != path ({version_path})"
        )

    # Check the amplicon size matches the schemepath
    ampliconsize_path = infopath.parent.parent.name
    if str(primerscheme.amplicon_size) != ampliconsize_path:
        raise ValueError(
            f"Amplicon size mismatch for {primerscheme_derived_scheme_path}: info ({primerscheme.amplicon_size}) != path ({ampliconsize_path})"
        )

    # Check the schemepath matches the path
    schemeid_path = infopath.parent.parent.parent.name
    if primerscheme.name != schemeid_path:
        raise ValueError(
            f"Name mismatch for {primerscheme_derived_scheme_path}: info ({primerscheme.name}) != path ({schemeid_path})"
        )


def validate_readme(infopath: Path, primerscheme: PrimerScheme | None = None):
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primerscheme is None:
        primerscheme = PrimerScheme.model_validate_json(infopath.read_text())

    primerscheme_derived_scheme_path = (
        Path(primerscheme.name) / str(primerscheme.amplicon_size) / primerscheme.version
    )

    # Check the ReadME.md
    readme = infopath.parent / "README.md"
    if not readme.exists():
        raise FileNotFoundError(f"{readme} does not exist")

    # Check the readme has been updated
    readme = readme.read_text()
    if readme.find(primerscheme.name) == -1:
        raise ValueError(
            f"Scheme name ({primerscheme.name}) not found in {readme}: {primerscheme_derived_scheme_path}"
        )
    if readme.find(str(primerscheme.amplicon_size)) == -1:
        raise ValueError(
            f"Amplicon size ({primerscheme.amplicon_size}) not found in {readme}: {primerscheme_derived_scheme_path}"
        )
    if readme.find(primerscheme.version) == -1:
        raise ValueError(
            f"Scheme version ({primerscheme.version}) not found in {readme}: {primerscheme_derived_scheme_path}"
        )


def validate_hashes(infopath: Path, primerscheme: PrimerScheme | None = None):
    # Use provided PrimerScheme (prevent duplication) or read in for validation
    if primerscheme is None:
        primerscheme = PrimerScheme.model_validate_json(infopath.read_text())

    primerscheme_derived_scheme_path = (
        Path(primerscheme.name) / str(primerscheme.amplicon_size) / primerscheme.version
    )

    # Check sha256 hash bedfile
    primer_path = infopath.parent / PRIMER_FILE_NAME
    primer_sha = sha256_checksum(primer_path)
    if primer_sha != primerscheme.primer_file_sha256:
        raise ValueError(
            f"File sha256 ({primer_sha} != info sha256 ({primerscheme.primer_file_sha256}): {primerscheme_derived_scheme_path}"
        )

    # Check sha256 hash ref
    reference_path = infopath.parent / REFERENCE_FILE_NAME
    reference_sha = sha256_checksum(reference_path)
    if reference_sha != primerscheme.reference_file_sha256:
        raise ValueError(
            f"File sha256 ({reference_sha} != info sha256 ({primerscheme.reference_file_sha256}): {primerscheme_derived_scheme_path}"
        )

    # TODO Check primaschema hashes.


def validate(
    infopath: Path,
    primerscheme: PrimerScheme | None = None,
    additional_linkml: bool = False,
):
    logger.debug(f"Validating {infopath}")
    if additional_linkml:
        logger.debug(f"Validating with LinkML: {infopath}")
        validate_scheme_json_with_linkml(infopath)
    


    if primerscheme is None:
        primerscheme = PrimerScheme.model_validate_json(infopath.read_text())
    validate_name(infopath, primerscheme)
    validate_hashes(infopath, primerscheme)
    validate_readme(infopath, primerscheme)


def validate_all(primerschemes_repo: Path, additional_linkml: bool = False):
    """
    Recursively searches through the primerschemes_repo for {METADATA_FILE_NAME}
    """

    for schemeinfo in primerschemes_repo.rglob(f"*/{METADATA_FILE_NAME}"):
        validate(schemeinfo, None, additional_linkml)
