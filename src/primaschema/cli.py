import json
import logging
import pathlib
import shutil
import sys
import tempfile
from typing import Annotated, Any, List, Optional

from cyclopts import App, Parameter, validators
from primalbedtools.bedfiles import BedLineParser, sort_bedlines
from primalbedtools.validate import validate_ref_and_bed
from pydantic import BeforeValidator, field_validator, model_validator
from rich.console import Console
from rich.traceback import install as install_rich_traceback

from primaschema import (
    METADATA_FILE_NAME,
    PRIMER_FILE_NAME,
    REFERENCE_FILE_NAME,
)
from primaschema.lib import plot_primers
from primaschema.schema.info import (
    Algorithm,
    Contributor,
    PrimerScheme,
    SchemeLicense,
    SchemeStatus,
    SchemeTag,
    TargetOrganism,
    Vendor,
)
from primaschema.schema.info import (
    version as SCHEMA_VERSION,
)
from primaschema.schema.index import (
    PrimerSchemeIndex,
    update_index,
)
from primaschema.setup_logging import LogLevel, configure_logging
from primaschema.util import (
    find_all_info_json,
    primaschema_bed_hash,
    primaschema_ref_hash,
    read_fasta_records,
    sha256_checksum,
    write_fasta_records,
)
from primaschema.validate import validate as validate_scheme

logger = logging.getLogger(__name__)


LICENSE_TXT_CC_BY_SA_4_0 = """\n\n------------------------------------------------------------------------

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/)

![](https://i.creativecommons.org/l/by-sa/4.0/88x31.png)"""

# Patch PrimerScheme to fix cyclopts issue with string defaults for Enums
# See https://github.com/pha4ge/primaschema/issues/new
if PrimerScheme.model_fields["license"].default == "CC-BY-SA-4.0":
    PrimerScheme.model_fields["license"].default = SchemeLicense.CC_BY_SA_4FULL_STOP0

# Add rich formatted errors
error_console = Console(stderr=True)
install_rich_traceback(console=error_console)

# Create the apps
app = App(
    name="primaschema",
    version_flags="--show-version",
    error_console=error_console,
    default_parameter=Parameter(
        show_default=True,
    ),
)


@app.meta.default
def cli_launcher(
    *tokens: Annotated[str, Parameter(show=False, allow_leading_hyphen=True)],
    log_level: Annotated[
        LogLevel | None,
        Parameter(name=["--log-level", "-l"], show_default=True),
    ] = LogLevel.INFO,
):
    configure_logging(log_level=log_level)
    app(tokens)


modify_app = App(name="modify", help="Modify fields of an existing primer scheme")
app.command(modify_app)


def parse_contributor_single(v: Any) -> Contributor:
    """Parses a single contributor from various input formats.

    Args:
        v (Any): The input value to parse. Can be a Contributor object, a dictionary,
            or a string. If a string, it can be a JSON object, a comma-separated
            key-value string (e.g., "name=John,email=john@example.com"), or
            simply the name of the contributor.

    Returns:
        Contributor: A Contributor object parsed from the input.

    Raises:
        ValueError: If the input cannot be parsed into a Contributor.
    """
    if isinstance(v, Contributor):
        return v
    if isinstance(v, dict):
        return Contributor(**v)
    if isinstance(v, str):
        # Try JSON first
        try:
            data = json.loads(v)
            if isinstance(data, dict):
                return Contributor(**data)
        except json.JSONDecodeError:
            pass

        # Key-value parsing
        if "=" in v:
            parts = {}
            for part in v.split(","):
                if "=" in part:
                    key, val = part.split("=", 1)
                    parts[key.strip()] = val.strip()
            return Contributor(**parts)

        # Fallback to name only
        return Contributor(name=v)
    raise ValueError(f"Cannot parse contributor: {v}")


def parse_contributors_pydantic(v: Any) -> List[Contributor]:
    if isinstance(v, list):
        return [parse_contributor_single(x) for x in v]
    return v


def parse_vendor_single(v: Any) -> Vendor:
    """Parses a single vendor from various input formats."""
    if isinstance(v, Vendor):
        return v
    if isinstance(v, dict):
        return Vendor(**v)
    if isinstance(v, str):
        # Try JSON first
        try:
            data = json.loads(v)
            if isinstance(data, dict):
                return Vendor(**data)
        except json.JSONDecodeError:
            pass
        # Key-value parsing
        if "=" in v:
            parts = {}
            for part in v.split(","):
                if "=" in part:
                    key, val = part.split("=", 1)
                    parts[key.strip()] = val.strip()
            return Vendor(**parts)

        # Fallback to organisation_name only
        return Vendor(organisation_name=v)
    raise ValueError(f"Cannot parse vendor: {v}")


def _save_and_rebuild_readme(
    info_path: pathlib.Path, primer_scheme: PrimerScheme, rebuild_plot: bool = False
):
    """Saves the PrimerScheme to info.json and rebuilds the README."""
    # Save info.json
    logger.debug(f"Writing info.json to {info_path}")
    with open(info_path, "w") as f:
        f.write(
            primer_scheme.model_dump_json(
                exclude_unset=True, exclude_none=True, indent=4
            )
        )

    # Regenerate README
    scheme_dir = info_path.parent
    logger.debug(f"Regenerating README.md in {scheme_dir}")
    generate_readme(scheme_dir, primer_scheme)

    if rebuild_plot:
        logger.debug(f"Ensuring plot output directory in {scheme_dir / 'work'}")
        (scheme_dir / "work").mkdir(exist_ok=True)
        logger.debug(f"Rendering primer plot to {scheme_dir / 'work' / 'primer.svg'}")
        plot_primers(scheme_dir / PRIMER_FILE_NAME, scheme_dir / "work" / "primer.svg")


def create_status_badge(primer_scheme: PrimerScheme) -> str:
    """
    Create a badge for the README.md file
    """
    match primer_scheme.status:
        case SchemeStatus.VALIDATED:
            color = "green"
        case SchemeStatus.WITHDRAWN | SchemeStatus.DEPRECATED:
            color = "red"
        case _:
            color = "blue"

    return f"![Generic badge](https://img.shields.io/badge/STATUS-{primer_scheme.status}-{color}.svg)"


def generate_readme(path: pathlib.Path, primer_scheme: PrimerScheme):
    """
    Generate the README.md file for a primer scheme

    :param path: The path to the scheme directory
    :type path: pathlib.Path
    :param info: The scheme information
    :type info: Info
    :param pngs: The list of PNG files
    :type pngs: list[pathlib.Path]
    """

    with open(path / "README.md", "w") as readme:
        readme.write(
            f"# {primer_scheme.name} {primer_scheme.amplicon_size}bp {primer_scheme.version}\n\n"
        )
        # Add the status badge
        readme.write(f"{create_status_badge(primer_scheme)}\n\n")

        # Add citation if present
        if primer_scheme.citations and primer_scheme.citations is not None:
            for cit in primer_scheme.citations:
                readme.write(f"> If you use this scheme please cite: {cit}\n\n")

        readme.write(
            f"[primalscheme labs](https://labs.primalscheme.com/detail/{primer_scheme.name}/{primer_scheme.amplicon_size}/{primer_scheme.version})\n\n"
        )

        if primer_scheme.notes and primer_scheme.notes is not None:
            readme.write("## Notes\n\n")
            for note in primer_scheme.notes:
                readme.write(note + "\n\n")

        readme.write("## Metadata\n\n")
        if primer_scheme.target_organisms:
            readme.write("**Target Organisms:**\n")
            for to in primer_scheme.target_organisms:
                to_str = f"- {to.common_name or ''}"
                if to.ncbi_tax_id:
                    to_str += f" (Tax ID: {to.ncbi_tax_id})"
                readme.write(f"{to_str}\n")
            readme.write("\n")

        if primer_scheme.derived_from:
            readme.write(f"**Derived from:** {primer_scheme.derived_from}\n\n")

        if primer_scheme.tags:
            readme.write(f"**Tags:** {', '.join(primer_scheme.tags)}\n\n")

        if primer_scheme.contributors:
            readme.write("## Contributors\n\n")
            for contributor in primer_scheme.contributors:
                contrib_str = f"- {contributor.name}"
                if contributor.email:
                    contrib_str += f" <{contributor.email}>"
                if contributor.orcid_id:
                    contrib_str += f" (ORCID: {contributor.orcid_id})"
                readme.write(f"{contrib_str}\n")
            readme.write("\n")

        if primer_scheme.vendors:
            readme.write("## Vendors\n\n")
            for vendor in primer_scheme.vendors:
                vendor_str = f"- {vendor.organisation_name}"
                if vendor.kit_name:
                    vendor_str += f": {vendor.kit_name}"
                if vendor.home_page:
                    vendor_str += f" ([Website]({vendor.home_page}))"
                readme.write(f"{vendor_str}\n")
            readme.write("\n")

        readme.write("## Overviews\n\n")
        readme.write(
            '<div style="width: 100%;"><img src="work/primer.svg" style="width: 100%;" alt="Click to see the source"></div>\n\n'
        )

        readme.write("## Details\n\n")

        # Write the details into the readme
        readme.write(
            f"""```json\n{primer_scheme.model_dump_json(indent=4, exclude_unset=True, exclude_none=True)}\n```\n\n"""
        )

        if primer_scheme.license == SchemeLicense.CC_BY_SA_4FULL_STOP0:
            readme.write(LICENSE_TXT_CC_BY_SA_4_0)


def parse_algorithm(v: Any) -> Optional[Algorithm]:
    if v is None:
        return None
    if isinstance(v, Algorithm):
        return v
    if isinstance(v, dict):
        return Algorithm(**v)
    if isinstance(v, str):
        if ":" in v:
            name, version = v.split(":", 1)
            return Algorithm(name=name, version=version)
        return Algorithm(name=v)
    raise ValueError(f"Cannot parse algorithm: {v}")


def parse_target_organism_single(v: Any) -> TargetOrganism:
    if isinstance(v, TargetOrganism):
        return v
    if isinstance(v, dict):
        return TargetOrganism(**v)
    if isinstance(v, str):
        # Try JSON first
        try:
            data = json.loads(v)
            if isinstance(data, dict):
                return TargetOrganism(**data)
        except json.JSONDecodeError:
            pass

        # Key-value parsing
        if "=" in v:
            parts = {}
            for part in v.split(","):
                if "=" in part:
                    key, val = part.split("=", 1)
                    parts[key.strip()] = val.strip()
            return TargetOrganism(**parts)

        # If it looks like an int, assume it's a tax id
        if v.isdigit():
            return TargetOrganism(ncbi_tax_id=v)

        # Otherwise assume common name
        return TargetOrganism(common_name=v)
    raise ValueError(f"Cannot parse target organism: {v}")


def parse_target_organisms_pydantic(v: Any) -> List[TargetOrganism]:
    if isinstance(v, list):
        return [parse_target_organism_single(x) for x in v]
    if isinstance(v, (str, dict, TargetOrganism)):
        return [parse_target_organism_single(v)]
    return v


def parse_vendors_pydantic(v: Any) -> Optional[List[Vendor]]:
    if v is None:
        return None
    if isinstance(v, list):
        return [parse_vendor_single(x) for x in v]
    if isinstance(v, (str, dict, Vendor)):
        return [parse_vendor_single(v)]
    return v


class CLIPrimerScheme(PrimerScheme):
    schema_version: Annotated[str, Parameter(parse=False)] = SCHEMA_VERSION
    contributors: Annotated[  # type: ignore
        List[Contributor], BeforeValidator(parse_contributors_pydantic)
    ]
    target_organisms: Annotated[  # type: ignore
        List[TargetOrganism], BeforeValidator(parse_target_organisms_pydantic)
    ]
    vendors: Annotated[
        Optional[List[Vendor]], BeforeValidator(parse_vendors_pydantic)
    ] = None
    algorithm: Annotated[Optional[Algorithm], Parameter(parse=False)] = None
    # Don't expose the checksums to cli
    primer_checksum: Annotated[str | None, Parameter(parse=False)] = None
    primer_file_sha256: Annotated[str | None, Parameter(parse=False)] = None
    reference_checksum: Annotated[str | None, Parameter(parse=False)] = None
    reference_file_sha256: Annotated[str | None, Parameter(parse=False)] = None

    @field_validator("target_organisms")
    def validate_target_organisms(cls, v):
        for to in v:
            if not to.common_name and not to.ncbi_tax_id:
                raise ValueError(
                    "TargetOrganism must have at least one of 'common_name' or 'ncbi_tax_id'"
                )
        return v

    @model_validator(mode="before")
    @classmethod
    def uppercase_enums(cls, data: Any) -> Any:
        if isinstance(data, dict):
            # Uppercase status if it's a string
            if "status" in data and isinstance(data["status"], str):
                data["status"] = data["status"].upper()

            # Uppercase license if it's a string
            if "license" in data and isinstance(data["license"], str):
                data["license"] = data["license"].upper()

            # Uppercase tags if it's a list of strings
            if "tags" in data and isinstance(data["tags"], list):
                data["tags"] = [
                    t.upper() if isinstance(t, str) else t for t in data["tags"]
                ]
        return data


@app.command
def create(
    cli_ps: Annotated[
        CLIPrimerScheme,
        Parameter(name="*"),
    ],
    bed_path: Annotated[
        pathlib.Path,
        Parameter(
            validator=validators.Path(exists=True, file_okay=True),
            help="The path to the corresponding primer.bed file",
        ),
    ],
    reference_path: Annotated[
        pathlib.Path,
        Parameter(
            validator=validators.Path(exists=True, file_okay=True),
            help="The path to the corresponding reference.fasta file",
        ),
    ],
    primer_schemes_path: Annotated[
        pathlib.Path,
        Parameter(
            env_var="PRIMER_SCHEMES_PATH",
            validator=validators.Path(exists=True, dir_okay=True, file_okay=False),
            help="The path to the primer schemes directory. Will use the ENV VAR PRIMER_SCHEMES_PATH",
        ),
    ],
    algorithm: Annotated[
        Optional[str],
        Parameter(
            help="The algorithm used to generate the scheme (e.g. primalscheme:3.0.3)"
        ),
    ] = None,
):
    """Create a new primer scheme definition"""
    # Parse algorithm if provided
    if algorithm:
        cli_ps.algorithm = parse_algorithm(algorithm)
        logger.debug(f"Parsed algorithm '{algorithm}' -> Algorithm({cli_ps.algorithm})")

    # Convert to base PrimerScheme to ensure strict adherence to the schema
    ps = PrimerScheme.model_validate(cli_ps.model_dump())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    _headers, bedlines = BedLineParser.from_file(str(bed_path))
    bedlines = sort_bedlines(bedlines)
    logger.debug(f"Loaded and sorted bedlines from {bed_path}")

    # Create a directory to store the new scheme in.
    output_dir = primer_schemes_path / ps.name / str(ps.amplicon_size) / ps.version
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")

    logger.debug(f"Creating scheme at {output_dir}")

    # Use a tmp dir to ensure atomic
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = pathlib.Path(tmp_dir)
        tmp_version_level = tmp_path / ps.version
        tmp_version_level.mkdir()
        logger.debug("Created tmp dir")

        # Move / Write the bedfile
        BedLineParser.to_file(tmp_version_level / PRIMER_FILE_NAME, _headers, bedlines)
        # Parse ref
        reference_records = read_fasta_records(reference_path)
        write_fasta_records(tmp_version_level / REFERENCE_FILE_NAME, reference_records)

        # Validate the bed and ref files files
        validate_ref_and_bed(
            bedlines, str((tmp_version_level / REFERENCE_FILE_NAME).absolute())
        )
        logger.debug(
            f"Generated validated {PRIMER_FILE_NAME} and {REFERENCE_FILE_NAME}"
        )

        # Generate hashes of the files
        ps.primer_file_sha256 = sha256_checksum(tmp_version_level / PRIMER_FILE_NAME)
        ps.reference_file_sha256 = sha256_checksum(
            tmp_version_level / REFERENCE_FILE_NAME
        )
        logger.debug(
            f"Generated sha256 hashes for {PRIMER_FILE_NAME} ({ps.primer_file_sha256})"
            f" and {REFERENCE_FILE_NAME} ({ps.reference_file_sha256})"
        )

        # add primaschema hashes
        ps.primer_checksum = primaschema_bed_hash(None, bedlines)
        ps.reference_checksum = primaschema_ref_hash(None, reference_records)
        logger.debug(
            f"Generated primaschema hashes for {PRIMER_FILE_NAME} ({ps.primer_checksum})"
            f" and {REFERENCE_FILE_NAME} ({ps.reference_checksum})"
        )

        # Write info.json to tmp
        _save_and_rebuild_readme(tmp_version_level / METADATA_FILE_NAME, ps, True)
        # if all valid copy the tmp_version_level to output_dir
        shutil.copytree(tmp_version_level, output_dir)
        logger.debug(f"Copied tmp dir -> {output_dir}")
    # log
    logger.info(f"Created scheme {scheme_label} at {output_dir}")


@modify_app.command
def add_contributor(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    contributor: Annotated[
        Contributor,
        Parameter(name="*", converter=parse_contributor_single),
    ],
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Add a contributor to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx is not None:
        logger.debug(f"Inserting contributor at idx={idx}: {contributor}")
        ps.contributors.insert(idx, contributor)
        actual_idx = idx
    else:
        logger.debug(f"Appending contributor: {contributor}")
        ps.contributors.append(contributor)
        actual_idx = len(ps.contributors) - 1
    _save_and_rebuild_readme(info_path, ps)
    logger.info(
        f"Updated contributors for {scheme_label}: added {contributor} at idx {actual_idx}"
    )


@modify_app.command
def remove_contributor(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
):
    """Remove a contributor by index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx >= len(ps.contributors):
        logger.error(
            f"Index {idx} out of range for contributors in {info_path}. "
            f"Valid range is 0..{len(ps.contributors) - 1}."
        )
        sys.exit(1)
    removed = ps.contributors[idx]
    logger.debug(f"Removing contributor at idx={idx}: {removed}")
    ps.contributors.pop(idx)
    _save_and_rebuild_readme(info_path, ps)
    logger.info(
        f"Updated contributors for {scheme_label}: removed {removed} at idx {idx}"
    )


@modify_app.command
def update_contributor(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
    contributor: Annotated[
        Contributor,
        Parameter(name="*", converter=parse_contributor_single),
    ],
):
    """Update a contributor at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx >= len(ps.contributors):
        logger.error(
            f"Index {idx} out of range for contributors in {info_path}. "
            f"Valid range is 0..{len(ps.contributors) - 1}."
        )
        sys.exit(1)
    previous = ps.contributors[idx]
    logger.debug(f"Updating contributor at idx={idx}: {previous} -> {contributor}")
    ps.contributors[idx] = contributor
    _save_and_rebuild_readme(info_path, ps)
    logger.info(
        f"Updated contributors for {scheme_label}: idx {idx} {previous} -> {contributor}"
    )


@modify_app.command
def add_vendor(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    vendor: Annotated[
        Vendor,
        Parameter(name="*", converter=parse_vendor_single),
    ],
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Add a vendor to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if ps.vendors is None:
        logger.debug("Initializing vendors list")
        ps.vendors = []
    if idx is not None:
        logger.debug(f"Inserting vendor at idx={idx}: {vendor}")
        ps.vendors.insert(idx, vendor)
        actual_idx = idx
    else:
        logger.debug(f"Appending vendor: {vendor}")
        ps.vendors.append(vendor)
        actual_idx = len(ps.vendors) - 1
    _save_and_rebuild_readme(info_path, ps)
    logger.info(
        f"Updated vendors for {scheme_label}: added {vendor} at idx {actual_idx}"
    )


@modify_app.command
def remove_vendor(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
):
    """Remove a vendor by index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if not ps.vendors or idx >= len(ps.vendors):
        max_idx = len(ps.vendors) - 1 if ps.vendors else -1
        logger.error(
            f"Index {idx} out of range for vendors in {info_path}. "
            f"Valid range is 0..{max_idx}."
        )
        sys.exit(1)
    removed = ps.vendors[idx]
    logger.debug(f"Removing vendor at idx={idx}: {removed}")
    ps.vendors.pop(idx)
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated vendors for {scheme_label}: removed {removed} at idx {idx}")


@modify_app.command
def update_vendor(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
    vendor: Annotated[
        Vendor,
        Parameter(name="*", converter=parse_vendor_single),
    ],
):
    """Update a vendor at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if not ps.vendors or idx >= len(ps.vendors):
        max_idx = len(ps.vendors) - 1 if ps.vendors else -1
        logger.error(
            f"Index {idx} out of range for vendors in {info_path}. "
            f"Valid range is 0..{max_idx}."
        )
        sys.exit(1)
    previous = ps.vendors[idx]
    logger.debug(f"Updating vendor at idx={idx}: {previous} -> {vendor}")
    ps.vendors[idx] = vendor
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated vendors for {scheme_label}: idx {idx} {previous} -> {vendor}")


@modify_app.command
def add_tag(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    tag: SchemeTag,
):
    """Add a tag to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if ps.tags is None:
        logger.debug("Initializing tags list")
        ps.tags = []
    if tag not in ps.tags:
        logger.debug(f"Adding tag: {tag}")
        ps.tags.append(tag)
        _save_and_rebuild_readme(info_path, ps)
        logger.info(f"Updated tags for {scheme_label}: added {tag}")
        return
    logger.info(f"No change for tags on {scheme_label}: {tag} already present")


@modify_app.command
def remove_tag(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    tag: SchemeTag,
):
    """Remove a tag from the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if ps.tags and tag in ps.tags:
        logger.debug(f"Removing tag: {tag}")
        ps.tags.remove(tag)
        _save_and_rebuild_readme(info_path, ps)
        logger.info(f"Updated tags for {scheme_label}: removed {tag}")
        return
    logger.info(f"No change for tags on {scheme_label}: {tag} not present")


@modify_app.command
def update_license(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    license: SchemeLicense,
):
    """Update the scheme license."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    previous = ps.license
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    logger.debug(f"Updating license: {previous} -> {license}")
    ps.license = license
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated license for {scheme_label}: {previous} -> {license}")


@modify_app.command
def update_status(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    status: SchemeStatus,
):
    """Update the scheme status."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    previous = ps.status
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    logger.debug(f"Updating status: {previous} -> {status}")
    ps.status = status
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated status for {scheme_label}: {previous} -> {status}")


@modify_app.command
def remove_target_organism(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
):
    """Remove a target organism by index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx >= len(ps.target_organisms):
        logger.error(
            f"Index {idx} out of range for target_organisms in {info_path}. "
            f"Valid range is 0..{len(ps.target_organisms) - 1}."
        )
        sys.exit(1)
    removed = ps.target_organisms[idx]
    logger.debug(f"Removing target_organism at idx={idx}: {removed}")
    ps.target_organisms.pop(idx)
    _save_and_rebuild_readme(info_path, ps)
    logger.info(
        f"Updated target_organisms for {scheme_label}: removed {removed} at idx {idx}"
    )


@modify_app.command
def add_target_organism(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    target_organism: Annotated[Optional[TargetOrganism], Parameter(name="*")] = None,
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Adds a target organism at a specific index."""
    if target_organism is None:
        target_organism = TargetOrganism()

    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")

    # append
    if idx is None:
        idx = len(ps.target_organisms)

    logger.debug(f"Adding target_organism at idx={idx}: {target_organism}")
    ps.target_organisms.insert(idx, target_organism)
    _save_and_rebuild_readme(info_path, ps)
    logger.info(
        f"Updated target_organisms for {scheme_label}: added {target_organism} at idx {idx}"
    )


@modify_app.command
def update_algorithm(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    algorithm: Algorithm,
):
    """Update the algorithm."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    previous = ps.algorithm
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    logger.debug(f"Updating algorithm: {previous} -> {algorithm}")
    ps.algorithm = algorithm
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated algorithm for {scheme_label}: {previous} -> {algorithm}")


# Index commands
@app.command
def index(
    primer_schemes_path: Annotated[
        pathlib.Path,
        Parameter(
            env_var="PRIMER_SCHEMES_PATH",
            validator=validators.Path(exists=True, dir_okay=True, file_okay=False),
            help="The path to the primer schemes directory. Will use the ENV VAR PRIMER_SCHEMES_PATH",
        ),
    ],
    index_path: Optional[pathlib.Path] = None,
    base_url: str = "",
):
    """Build a JSON index of all primer schemes in a directory"""
    # Read in current index
    if index_path is not None:
        psi = PrimerSchemeIndex.model_validate_json(index_path.read_text())
    else:
        psi = PrimerSchemeIndex()

    # find all primer schemes
    ps = []
    for ps_info in find_all_info_json(primer_schemes_path):
        logger.debug(f"found {ps_info}")
        ps.append(PrimerScheme.model_validate_json(ps_info.read_text()))

    update_index(ps, psi, base_url=base_url)

    # Ensure schemes is marked as set for exclude_unset=True
    psi.primerschemes = psi.primerschemes

    print(psi.model_dump_json(exclude_unset=True, exclude_none=True))


# Validate commands
@app.command
def validate(
    path: Annotated[
        pathlib.Path,
        Parameter(
            env_var="PRIMER_SCHEMES_PATH",
            validator=validators.Path(exists=True),
            help="Path to an info.json file, or a directory of schemes when using --all",
        ),
    ],
    all: bool = False,
    additional_linkml: bool = False,
    strict: bool = True,
    fix: Annotated[
        bool,
        Parameter(
            name="--fix",
            help="Normalise primer.bed and reference.fasta in place if they differ only by formatting",
        ),
    ] = False,
):
    """Validate primer scheme definitions"""
    if all:
        logger.debug(f"Validating all schemes under {path}")
        errors: list[str] = []
        for info_path in find_all_info_json(path):
            ps = PrimerScheme.model_validate_json(info_path.read_text())
            scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
            logger.debug(f"Validating scheme {scheme_label} from {info_path}")
            try:
                validate_scheme(
                    info_path,
                    ps,
                    additional_linkml,
                    strict,
                    fix=fix,
                )
                logger.info(f"Validated scheme {scheme_label}")
            except Exception as exc:
                logger.error(f"Validation failed for {info_path}: {exc}")
                errors.append(f"{info_path}: {exc}")
        if errors:
            raise ValueError(
                f"Validation failed for {len(errors)} scheme(s):\n" + "\n".join(errors)
            )
    else:
        ps = PrimerScheme.model_validate_json(path.read_text())
        scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
        logger.debug(f"Validating scheme {scheme_label} from {path}")
        validate_scheme(
            path,
            ps,
            additional_linkml,
            strict,
            fix=fix,
        )
        logger.info(f"Validated scheme {scheme_label}")


def _sync_metadata_from_path(
    primer_scheme: PrimerScheme, info_path: pathlib.Path
) -> bool:
    scheme_dir = info_path.parent
    version = scheme_dir.name
    amplicon_size_raw = scheme_dir.parent.name
    name = scheme_dir.parent.parent.name

    try:
        amplicon_size = int(amplicon_size_raw)
    except ValueError as exc:
        raise ValueError(
            f"Invalid amplicon size in path {scheme_dir}: {amplicon_size_raw}"
        ) from exc

    changed = False
    if primer_scheme.name != name:
        logger.debug(
            f"Syncing scheme name from {primer_scheme.name} to {name} for {info_path}"
        )
        primer_scheme.name = name
        changed = True
    if primer_scheme.amplicon_size != amplicon_size:
        logger.debug(
            f"Syncing amplicon_size from {primer_scheme.amplicon_size} to {amplicon_size} for {info_path}"
        )
        primer_scheme.amplicon_size = amplicon_size
        changed = True
    if primer_scheme.version != version:
        logger.debug(
            f"Syncing version from {primer_scheme.version} to {version} for {info_path}"
        )
        primer_scheme.version = version
        changed = True
    return changed


def _rebuild_one(
    info_path: pathlib.Path,
    reformat_primer_bed: bool = False,
    sync_metadata: bool = True,
) -> str:
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if sync_metadata:
        if _sync_metadata_from_path(ps, info_path):
            logger.debug(f"Synced scheme metadata from path for {info_path}")
    scheme_label = f"{ps.name}/{ps.amplicon_size}/{ps.version}"
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    _h, bls = BedLineParser.from_file(info_path.parent / PRIMER_FILE_NAME)
    logger.debug(f"Loaded bedlines from {info_path.parent / PRIMER_FILE_NAME}")
    if reformat_primer_bed:
        logger.debug("Sorting bedlines for reformat_primer_bed")
        bls = sort_bedlines(bls)
        BedLineParser.to_file(info_path.parent / PRIMER_FILE_NAME, _h, bls)
        logger.debug(f"Wrote sorted bedlines to {info_path.parent / PRIMER_FILE_NAME}")
    logger.debug("Validating primer.bed against reference.fasta")
    validate_ref_and_bed(bls, str((info_path.parent / REFERENCE_FILE_NAME).absolute()))
    logger.debug("Computing sha256 and primaschema hashes")
    ps.primer_file_sha256 = sha256_checksum(info_path.parent / PRIMER_FILE_NAME)
    ps.reference_file_sha256 = sha256_checksum(info_path.parent / REFERENCE_FILE_NAME)
    ps.primer_checksum = primaschema_bed_hash(None, bls)
    ps.reference_checksum = primaschema_ref_hash(
        info_path.parent / REFERENCE_FILE_NAME, None
    )
    _save_and_rebuild_readme(info_path, ps)
    return scheme_label


@app.command
def rebuild(
    path: Annotated[
        pathlib.Path,
        Parameter(
            validator=validators.Path(exists=True),
            help="Path to an info.json file, or a directory of schemes when using --all",
        ),
    ],
    all: bool = False,
    reformat_primer_bed: bool = False,
    sync_metadata: Annotated[
        bool,
        Parameter(
            name="--sync-metadata",
            help="Sync name/amplicon_size/version from the scheme path",
        ),
    ] = True,
):
    """Rebuild and normalise primer scheme metadata"""
    if all:
        for info_path in find_all_info_json(path):
            scheme_label = _rebuild_one(
                info_path,
                reformat_primer_bed=reformat_primer_bed,
                sync_metadata=sync_metadata,
            )
            logger.info(f"Rebuilt scheme {scheme_label}")
    else:
        scheme_label = _rebuild_one(
            path,
            reformat_primer_bed=reformat_primer_bed,
            sync_metadata=sync_metadata,
        )
        logger.info(f"Rebuilt scheme {scheme_label}")


def main():
    app.meta()


if __name__ == "__main__":
    main()
