from primaschema.schema.info import (
    PrimerScheme,
    Vendor,
    Contributor,
    SchemeLicense,
    SchemeStatus,
    SchemeTag,
)
from typing import Annotated, Any, List, Optional
import json
import pathlib
import shutil
import sys
from primalbedtools.validate import validate
from primalbedtools.bedfiles import BedLineParser, sort_bedlines

from cyclopts import App, Parameter, validators
from pydantic import BeforeValidator, model_validator

from rich.console import Console
from rich.traceback import install as install_rich_traceback
from primaschema import METADATA_FILE_NAME, logger
from primaschema.schema.manifest import (
    PrimerSchemeIndex,
    ManifestPrimerScheme,
    update_index,
)
from primaschema.cli import configure_logging
from primaschema.util import sha256_checksum, find_all_info_json
from primaschema.validate import validate_all

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
    name="primaschema", version_flags="--show-version", error_console=error_console
)
modify_app = App(name="modify")
app.command(modify_app)

index_app = App(name="index")
app.command(index_app)

validate_app = App(name="validate")
app.command(validate_app)


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


def _save_and_regenerate(info_path: pathlib.Path, ps: PrimerScheme):
    """Saves the PrimerScheme to info.json and regenerates the README."""
    # Save info.json
    with open(info_path, "w") as f:
        f.write(ps.model_dump_json(exclude_unset=True, exclude_none=True, indent=1))

    # Regenerate README
    scheme_dir = info_path.parent
    generate_readme(scheme_dir, ps, [])
    print(f"Updated {info_path} and regenerated README.md")


def create_status_badge(primerscheme: PrimerScheme) -> str:
    """
    Create a badge for the README.md file
    """
    match primerscheme.status:
        case SchemeStatus.VALIDATED:
            color = "green"
        case SchemeStatus.WITHDRAWN | SchemeStatus.DEPRECATED:
            color = "red"
        case _:
            color = "blue"

    return f"[![Generic badge](https://img.shields.io/badge/STATUS-{primerscheme.status}-{color}.svg)]"


def generate_readme(
    path: pathlib.Path, primerscheme: PrimerScheme, pngs: list[pathlib.Path]
):
    """
    Generate the README.md file for a scheme

    :param path: The path to the scheme directory
    :type path: pathlib.Path
    :param info: The scheme information
    :type info: Info
    :param pngs: The list of PNG files
    :type pngs: list[pathlib.Path]
    """

    with open(path / "README.md", "w") as readme:
        readme.write(
            f"# {primerscheme.name} {primerscheme.amplicon_size}bp {primerscheme.version}\n\n"
        )
        # Add the status badge
        readme.write(f"{create_status_badge(primerscheme)}\n\n")

        # Add citation if present
        if primerscheme.citations and primerscheme.citations is not None:
            for cit in primerscheme.citations:
                readme.write(f"> If you use this scheme please cite: {cit}\n\n")

        readme.write(
            f"[primalscheme labs](https://labs.primalscheme.com/detail/{primerscheme.name}/{primerscheme.amplicon_size}/{primerscheme.version})\n\n"
        )

        if primerscheme.notes and primerscheme.notes is not None:
            readme.write("## Notes\n\n")
            for note in primerscheme.notes:
                readme.write(note + "\n\n")

        readme.write("## Metadata\n\n")
        readme.write(f"**Organism:** {primerscheme.organism}\n\n")

        if primerscheme.derived_from:
            readme.write(f"**Derived from:** {primerscheme.derived_from}\n\n")

        if primerscheme.tags:
            readme.write(f"**Tags:** {', '.join(primerscheme.tags)}\n\n")

        if primerscheme.contributors:
            readme.write("## Contributors\n\n")
            for contributor in primerscheme.contributors:
                contrib_str = f"- {contributor.name}"
                if contributor.email:
                    contrib_str += f" <{contributor.email}>"
                if contributor.orcid:
                    contrib_str += f" (ORCID: {contributor.orcid})"
                readme.write(f"{contrib_str}\n")
            readme.write("\n")

        if primerscheme.vendors:
            readme.write("## Vendors\n\n")
            for vendor in primerscheme.vendors:
                vendor_str = f"- {vendor.organisation_name}"
                if vendor.kit_name:
                    vendor_str += f": {vendor.kit_name}"
                if vendor.home_page:
                    vendor_str += f" ([Website]({vendor.home_page}))"
                readme.write(f"{vendor_str}\n")
            readme.write("\n")

        if pngs:
            readme.write("## Overviews\n\n")
            for png in pngs:
                readme.write(f"![{png.name}](work/{png.name})\n\n")

        readme.write("## Details\n\n")

        # Write the details into the readme
        readme.write(
            f"""```json\n{primerscheme.model_dump_json(indent=4, exclude_unset=True, exclude_none=True)}\n```\n\n"""
        )

        if primerscheme.license == SchemeLicense.CC_BY_SA_4FULL_STOP0:
            readme.write(LICENSE_TXT_CC_BY_SA_4_0)


class CLIPrimerScheme(PrimerScheme):
    contributors: Annotated[
        List[Contributor], BeforeValidator(parse_contributors_pydantic)
    ]
    # Don't expose the checksums to cli
    primer_checksum: Annotated[str | None, Parameter(parse=False)] = None
    primer_file_sha256: Annotated[str | None, Parameter(parse=False)] = None
    reference_checksum: Annotated[str | None, Parameter(parse=False)] = None
    reference_file_sha256: Annotated[str | None, Parameter(parse=False)] = None

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
            help="The path to the primerschemes directory. Will use the ENV VAR PRIMER_SCHEMES_PATH",
        ),
    ],
):
    # Convert to base PrimerScheme to ensure strict adherence to the schema
    ps = PrimerScheme.model_validate(cli_ps.model_dump())

    # Validate the bed and ref files files
    validate(str(bed_path), str(reference_path.absolute()))

    # Create a directory to store the new scheme in.
    output_dir = primer_schemes_path / ps.name / str(ps.amplicon_size) / ps.version
    if output_dir.exists():
        print(f"Output directory already exists: {output_dir}", file=sys.stderr)
        raise ValueError(f"Output directory already exists: {output_dir}")

    output_dir.mkdir(parents=True)
    output_bed = output_dir / "primer.bed"
    output_ref = output_dir / "reference.fasta"

    # Move / Write the bedfile
    _headers, bedlines = BedLineParser.from_file(str(bed_path))
    bedlines = sort_bedlines(bedlines)
    BedLineParser.to_file(output_bed, _headers, bedlines)
    # Copy ref
    shutil.copy(reference_path, output_ref)

    # Generate hashes of the files
    ps.primer_file_sha256 = sha256_checksum(output_bed)
    ps.reference_file_sha256 = sha256_checksum(output_ref)

    # TODO add primaschema hashes

    # Write info.json
    with open(output_dir / METADATA_FILE_NAME, "w") as f:
        data = ps.model_dump_json(exclude_unset=True, exclude_none=True, indent=1)
        f.write(data)

    # Write the readme
    generate_readme(output_dir, ps, [])  # TODO auto generate pngs for the scheme


@modify_app.command
def add_contributor(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    contributor: Annotated[
        Contributor,
        Parameter(name="*", converter=parse_contributor_single),
    ],
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Add a contributor to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if idx is not None:
        ps.contributors.insert(idx, contributor)
    else:
        ps.contributors.append(contributor)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def remove_contributor(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
):
    """Remove a contributor by index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if idx >= len(ps.contributors):
        print(
            f"Index {idx} out of range. Max index is {len(ps.contributors) - 1}",
            file=sys.stderr,
        )
        sys.exit(1)
    ps.contributors.pop(idx)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def update_contributor(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
    contributor: Annotated[
        Contributor,
        Parameter(name="*", converter=parse_contributor_single),
    ],
):
    """Update a contributor at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if idx >= len(ps.contributors):
        print(
            f"Index {idx} out of range. Max index is {len(ps.contributors) - 1}",
            file=sys.stderr,
        )
        sys.exit(1)
    ps.contributors[idx] = contributor
    _save_and_regenerate(info_path, ps)


@modify_app.command
def add_vendor(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    vendor: Annotated[
        Vendor,
        Parameter(name="*", converter=parse_vendor_single),
    ],
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Add a vendor to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if ps.vendors is None:
        ps.vendors = []
    if idx is not None:
        ps.vendors.insert(idx, vendor)
    else:
        ps.vendors.append(vendor)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def remove_vendor(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
):
    """Remove a vendor by index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if not ps.vendors or idx >= len(ps.vendors):
        print(
            f"Index {idx} out of range.",
            file=sys.stderr,
        )
        sys.exit(1)
    ps.vendors.pop(idx)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def update_vendor(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
    vendor: Annotated[
        Vendor,
        Parameter(name="*", converter=parse_vendor_single),
    ],
):
    """Update a vendor at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if not ps.vendors or idx >= len(ps.vendors):
        print(
            f"Index {idx} out of range.",
            file=sys.stderr,
        )
        sys.exit(1)
    ps.vendors[idx] = vendor
    _save_and_regenerate(info_path, ps)


@modify_app.command
def add_tag(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    tag: SchemeTag,
):
    """Add a tag to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if ps.tags is None:
        ps.tags = []
    if tag not in ps.tags:
        ps.tags.append(tag)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def remove_tag(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    tag: SchemeTag,
):
    """Remove a tag from the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if ps.tags and tag in ps.tags:
        ps.tags.remove(tag)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def update_license(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    license: SchemeLicense,
):
    """Update the scheme license."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    ps.license = license
    _save_and_regenerate(info_path, ps)


@modify_app.command
def update_status(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    status: SchemeStatus,
):
    """Update the scheme status."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    ps.status = status
    _save_and_regenerate(info_path, ps)


# Index commands
@app.command
def build_index(
    primer_schemes_path: Annotated[
        pathlib.Path,
        Parameter(
            env_var="PRIMER_SCHEMES_PATH",
            validator=validators.Path(exists=True, dir_okay=True, file_okay=False),
            help="The path to the primerschemes directory. Will use the ENV VAR PRIMER_SCHEMES_PATH",
        ),
    ],
    manifest_path: Optional[pathlib.Path] = None,
    base_url: str = "",
):
    # Set up logging
    configure_logging(debug=True)

    # Read in current manifest
    if manifest_path is not None:
        psi = PrimerSchemeIndex.model_validate_json(manifest_path.read_text())
    else:
        psi = PrimerSchemeIndex()

    # find all primerschemes
    ps = []
    for ps_info in find_all_info_json(primer_schemes_path):
        logger.debug(f"found {ps_info}")
        ps.append(PrimerScheme.model_validate_json(ps_info.read_text()))

    update_index(ps, psi, base_url=base_url)

    # Ensure schemes is marked as set for exclude_unset=True
    psi.primerschemes = psi.primerschemes

    print(psi.model_dump_json(exclude_unset=True, exclude_none=True))


# Validate commands
@validate_app.command
def all(
    primer_schemes_path: Annotated[
        pathlib.Path,
        Parameter(
            env_var="PRIMER_SCHEMES_PATH",
            validator=validators.Path(exists=True, dir_okay=True, file_okay=False),
            help="The path to the primerschemes directory. Will use the ENV VAR PRIMER_SCHEMES_PATH",
        ),
    ],
    additional_linkml: bool = False,
):
    configure_logging(debug=True)
    validate_all(primer_schemes_path, additional_linkml)


app()
