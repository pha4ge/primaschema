from primaschema.schema.info import (
    PrimerScheme,
    Vendor,
    Contributor,
    SchemeLicense,
    SchemeStatus,
    SchemeTag,
    Algorithm,
    TargetOrganism,
    version as SCHEMA_VERSION,
)
import tempfile
from typing import Annotated, Any, List, Optional
import json
import pathlib
import shutil
import sys
from primalbedtools.validate import validate_ref_and_bed
from primalbedtools.bedfiles import BedLineParser, sort_bedlines
from Bio import SeqIO
from cyclopts import App, Parameter, validators
from pydantic import BeforeValidator, model_validator, field_validator

from rich.console import Console
from rich.traceback import install as install_rich_traceback
from primaschema import (
    METADATA_FILE_NAME,
    logger,
    PRIMER_FILE_NAME,
    REFERENCE_FILE_NAME,
)
from primaschema.schema.manifest import (
    PrimerSchemeIndex,
    ManifestPrimerScheme,
    update_index,
)
from primaschema.cli import configure_logging
from primaschema.util import sha256_checksum, find_all_info_json, primaschema_bed_hash, primaschema_ref_hash
from primaschema.validate import validate_all
from primaschema.lib import plot_primers


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


def _save_and_regenerate(
    info_path: pathlib.Path, ps: PrimerScheme, regenerate_plot: bool = False
):
    """Saves the PrimerScheme to info.json and regenerates the README."""
    # Save info.json
    with open(info_path, "w") as f:
        f.write(ps.model_dump_json(exclude_unset=True, exclude_none=True, indent=4))

    # Regenerate README
    scheme_dir = info_path.parent
    generate_readme(scheme_dir, ps)
    logger.debug(f"Scheme saved and regenerated README.md ({info_path})")

    if regenerate_plot:
        (scheme_dir / "work").mkdir(exist_ok=True)
        plot_primers(scheme_dir / PRIMER_FILE_NAME, scheme_dir / "work" / "primer.svg")


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


def generate_readme(path: pathlib.Path, primerscheme: PrimerScheme):
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
        if primerscheme.target_organisms:
            readme.write("**Target Organisms:**\n")
            for to in primerscheme.target_organisms:
                to_str = f"- {to.common_name or ''}"
                if to.ncbi_tax_id:
                    to_str += f" (Tax ID: {to.ncbi_tax_id})"
                readme.write(f"{to_str}\n")
            readme.write("\n")

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

        readme.write("## Overviews\n\n")
        readme.write(
            '<div style="width: 100%;"><img src="work/primer.svg" style="width: 100%;" alt="Click to see the source"></div>\n\n'
        )

        readme.write("## Details\n\n")

        # Write the details into the readme
        readme.write(
            f"""```json\n{primerscheme.model_dump_json(indent=4, exclude_unset=True, exclude_none=True)}\n```\n\n"""
        )

        if primerscheme.license == SchemeLicense.CC_BY_SA_4FULL_STOP0:
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
            help="The path to the primerschemes directory. Will use the ENV VAR PRIMER_SCHEMES_PATH",
        ),
    ],
    algorithm: Annotated[
        Optional[str],
        Parameter(
            help="The algorithm used to generate the scheme (e.g. primalscheme:3.0.3)"
        ),
    ] = None,
):
    # Parse algorithm if provided
    if algorithm:
        cli_ps.algorithm = parse_algorithm(algorithm)

    # Convert to base PrimerScheme to ensure strict adherence to the schema
    ps = PrimerScheme.model_validate(cli_ps.model_dump())
    _headers, bedlines = BedLineParser.from_file(str(bed_path))
    bedlines = sort_bedlines(bedlines)

    # Create a directory to store the new scheme in.
    output_dir = primer_schemes_path / ps.name / str(ps.amplicon_size) / ps.version
    if output_dir.exists():
        print(f"Output directory already exists: {output_dir}", file=sys.stderr)
        raise ValueError(f"Output directory already exists: {output_dir}")

    # Use a tmp dir to ensure atomic
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = pathlib.Path(tmp_dir)
        tmp_version_level = tmp_path / ps.version
        tmp_version_level.mkdir()

        # Move / Write the bedfile
        BedLineParser.to_file(tmp_version_level / PRIMER_FILE_NAME, _headers, bedlines)
        # Parse ref
        reference_records = list(SeqIO.parse(reference_path, "fasta"))
        with open(tmp_version_level / REFERENCE_FILE_NAME, "w") as ref_file:
            SeqIO.write(reference_records, ref_file, "fasta")
        # Validate the bed and ref files files
        validate_ref_and_bed(
            bedlines, str((tmp_version_level / REFERENCE_FILE_NAME).absolute())
        )

        # Generate hashes of the files
        ps.primer_file_sha256 = sha256_checksum(tmp_version_level / PRIMER_FILE_NAME)
        ps.reference_file_sha256 = sha256_checksum(
            tmp_version_level / REFERENCE_FILE_NAME
        )

        # add primaschema hashes
        ps.primer_checksum = primaschema_bed_hash(None, bedlines)
        ps.reference_checksum = primaschema_ref_hash(None, reference_records)

        # Write info.json to tmp
        _save_and_regenerate(tmp_version_level / METADATA_FILE_NAME, ps, True)
        # if all valid copy the tmp_version_level to output_dir
        shutil.copytree(tmp_version_level, output_dir)


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



@modify_app.command
def remove_target_organism(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    idx: Annotated[int, Parameter(validator=validators.Number(gte=0))],
):
    """Remove a target organism by index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    if idx >= len(ps.target_organisms):
        logger.warning(
            f"Index {idx} out of range. Max index is {len(ps.target_organisms) - 1}",
        )
        sys.exit(1)
    ps.target_organisms.pop(idx)
    _save_and_regenerate(info_path, ps)


@modify_app.command
def add_target_organism(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    target_organism: TargetOrganism,
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,

):
    """Adds a target organism at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())

    # append
    if idx is None:
        idx = len(ps.target_organisms)

    ps.target_organisms.insert(idx, target_organism) 
    _save_and_regenerate(info_path, ps)


@modify_app.command
def update_algorithm(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    algorithm: 
        Algorithm,
):
    """Update the algorithm."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    ps.algorithm = algorithm
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
    configure_logging(debug=False)

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
    strict: bool = True,
):
    configure_logging(debug=True)
    validate_all(primer_schemes_path, additional_linkml, strict)


@app.command
def regenerate(
    info_path: Annotated[
        pathlib.Path, Parameter(validator=validators.Path(exists=True, file_okay=True))
    ],
    reformat_primer_bed: bool = False,
):
    """
    Regenerates the metadata.
    - (optionally) reformats the primer.bed file
    - updates hashes
    - updates readme
    """
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    _h, bls = BedLineParser.from_file(info_path.parent / PRIMER_FILE_NAME)

    # Read in the primer.bed
    if reformat_primer_bed:
        bls = sort_bedlines(bls)
        BedLineParser.to_file(info_path.parent / PRIMER_FILE_NAME, _h, bls)

    validate_ref_and_bed(bls, str((info_path.parent / REFERENCE_FILE_NAME).absolute()))

    # Regenerate the hashes
    ps.primer_file_sha256 = sha256_checksum(info_path.parent / PRIMER_FILE_NAME)
    ps.reference_file_sha256 = sha256_checksum(info_path.parent / REFERENCE_FILE_NAME)
    ps.primer_checksum = primaschema_bed_hash(None, bls)
    ps.reference_checksum = primaschema_ref_hash(info_path.parent / REFERENCE_FILE_NAME, None)

    _save_and_regenerate(info_path, ps)



if __name__ == "__main__":
    app()
