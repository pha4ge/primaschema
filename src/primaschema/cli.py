import csv
import gzip
import json
import logging
import pathlib
import shutil
import tempfile
from datetime import date
from typing import Annotated, Any, List, Literal, Optional

from cyclopts import App, Parameter, validators
from cyclopts.utils import default_name_transform
from primalbedtools.bedfiles import BedLineParser, sort_bedlines
from primalbedtools.validate import validate_ref_and_bed
from pydantic import (
    BeforeValidator,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)
from rich.console import Console
from rich.traceback import install as install_rich_traceback

from primaschema import (
    DEFAULT_INDEX_URL,
    INDEX_FILE_NAME,
    METADATA_FILE_NAME,
    PRIMER_FILE_NAME,
    README_FILE_NAME,
    REFERENCE_FILE_NAME,
)
from primaschema.get_scheme import (
    DEFAULT_HTTP_TIMEOUT_SECONDS,
    DownloadError,
    SanitisationMode,
    download_schemes,
    load_index,
    resolve_schemes,
)
from primaschema.lib import plot_primers
from primaschema.license_footers import LICENSE_FOOTERS
from primaschema.schema.flatten import CSV_FIELDNAMES, flatten_scheme, unflatten_scheme
from primaschema.schema.index import (
    PrimerSchemeIndex,
    update_index,
)
from primaschema.schema.info import (
    PrimerSchemeGenerator,
    PrimerSchemeChecksums,
    PrimerSchemeContributor,
    PrimerSchemeApplication,
    PrimerSchemeLicense,
    PrimerSchemeScope,
    PrimerSchemeDevelopmentStatus,
    PrimerSchemeTargetOrganism,
    PrimerSchemeVendor,
)
from primaschema.schema.info import (
    version as SCHEMA_VERSION,
)
from primaschema.schema.primer_scheme import PrimerScheme
from primaschema.setup_logging import LogLevel, configure_logging
from primaschema.util import (
    ensure_path_within,
    find_all_info_json,
    read_fasta_records,
    serialize_primer_scheme_json,
    sha256_checksum,
    write_fasta_records,
)
from primaschema.validate import validate as validate_scheme

logger = logging.getLogger(__name__)

# Literal type built from PrimerSchemeLicense values so Cyclopts displays proper SPDX strings
_LicenseLiteral = Literal[
    "CC0-1.0",
    "CC-BY-4.0",
    "CC-BY-SA-4.0",
    "CC-BY-NC-4.0",
    "CC-BY-NC-SA-4.0",
    "CC-BY-ND-4.0",
    "CC-BY-NC-ND-4.0",
]

# Patch PrimerScheme to fix cyclopts issue with string defaults for Enums
# See https://github.com/pha4ge/primaschema/issues/new
if isinstance(PrimerScheme.model_fields["primer_scheme_license"].default, str):
    PrimerScheme.model_fields["primer_scheme_license"].default = PrimerSchemeLicense(
        PrimerScheme.model_fields["primer_scheme_license"].default
    )

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


# Errors we raise ourselves for expected, user-facing problems (bad input,
# failed validation, a mismatched checksum/identifier, a missing file) —
# as opposed to AttributeError/TypeError/KeyError/etc., which indicate an
# actual bug and should always show a full traceback so it can be diagnosed.
_EXPECTED_EXCEPTIONS = (ValueError, DownloadError, FileNotFoundError, ValidationError)


@app.meta.default
def cli_launcher(
    *tokens: Annotated[str, Parameter(show=False, allow_leading_hyphen=True)],
    log_level: Annotated[
        LogLevel | None,
        Parameter(name=["--log-level", "-l"], show_default=True),
    ] = LogLevel.INFO,
):
    configure_logging(log_level=log_level)
    try:
        app(tokens)
    except _EXPECTED_EXCEPTIONS as exc:
        if log_level == LogLevel.DEBUG:
            raise
        error_console.print(f"[bold red]Error:[/bold red] {exc}")
        raise SystemExit(1) from None


modify_app = App(name="modify", help="Modify fields of an existing primer scheme")
app.command(modify_app)


def _expand_short_keys(
    parts: dict[str, str], model_cls: type, prefixes: tuple[str, ...]
) -> dict[str, str]:
    """Expand short `key=value` keys to a model's full field name.

    Lets CLI input use short keys (e.g. `name=`) instead of the model's full,
    schema-compatible field name (e.g. `primer_scheme_contributor_name=`).
    A key is only rewritten if `prefix + key` matches a real field on
    `model_cls`; anything else (including keys that are already full field
    names) passes through unchanged so unknown keys still fail loudly.
    """
    fields = model_cls.model_fields
    expanded = {}
    for key, value in parts.items():
        if key not in fields:
            for prefix in prefixes:
                candidate = f"{prefix}{key}"
                if candidate in fields:
                    key = candidate
                    break
        expanded[key] = value
    return expanded


_CONTRIBUTOR_KEY_PREFIXES = ("primer_scheme_contributor_", "primer_scheme_")
_VENDOR_KEY_PREFIXES = ("primer_scheme_vendor_", "primer_scheme_")
_TARGET_ORGANISM_KEY_PREFIXES = (
    "primer_scheme_target_organism_",
    "primer_scheme_",
)


def parse_contributor_single(v: Any) -> PrimerSchemeContributor:
    """Parses a single contributor from various input formats.

    Args:
        v (Any): The input value to parse. Can be a PrimerSchemeContributor object, a dictionary,
            or a string. If a string, it can be a JSON object, a comma-separated
            key-value string (e.g., "name=John,email=john@example.com" or the full
            "primer_scheme_contributor_name=John,primer_scheme_contributor_email=john@example.com"),
            or simply the name of the contributor.

    Returns:
        PrimerSchemeContributor: A PrimerSchemeContributor object parsed from the input.

    Raises:
        ValueError: If the input cannot be parsed into a PrimerSchemeContributor.
    """
    if isinstance(v, PrimerSchemeContributor):
        return v
    if isinstance(v, dict):
        return PrimerSchemeContributor(
            **_expand_short_keys(v, PrimerSchemeContributor, _CONTRIBUTOR_KEY_PREFIXES)
        )
    if isinstance(v, str):
        # Try JSON first
        try:
            data = json.loads(v)
            if isinstance(data, dict):
                return PrimerSchemeContributor(
                    **_expand_short_keys(
                        data, PrimerSchemeContributor, _CONTRIBUTOR_KEY_PREFIXES
                    )
                )
        except json.JSONDecodeError:
            pass

        # Key-value parsing
        if "=" in v:
            parts = {}
            for part in v.split(","):
                if "=" in part:
                    key, val = part.split("=", 1)
                    parts[key.strip()] = val.strip()
            return PrimerSchemeContributor(
                **_expand_short_keys(
                    parts, PrimerSchemeContributor, _CONTRIBUTOR_KEY_PREFIXES
                )
            )

        # Fallback to name only
        return PrimerSchemeContributor(primer_scheme_contributor_name=v)
    raise ValueError(f"Cannot parse contributor: {v}")


def parse_contributors_pydantic(v: Any) -> List[PrimerSchemeContributor]:
    if isinstance(v, list):
        return [parse_contributor_single(x) for x in v]
    return v


def parse_vendor_single(v: Any) -> PrimerSchemeVendor:
    """Parses a single vendor from various input formats."""
    if isinstance(v, PrimerSchemeVendor):
        return v
    if isinstance(v, dict):
        return PrimerSchemeVendor(
            **_expand_short_keys(v, PrimerSchemeVendor, _VENDOR_KEY_PREFIXES)
        )
    if isinstance(v, str):
        # Try JSON first
        try:
            data = json.loads(v)
            if isinstance(data, dict):
                return PrimerSchemeVendor(
                    **_expand_short_keys(data, PrimerSchemeVendor, _VENDOR_KEY_PREFIXES)
                )
        except json.JSONDecodeError:
            pass
        # Key-value parsing
        if "=" in v:
            parts = {}
            for part in v.split(","):
                if "=" in part:
                    key, val = part.split("=", 1)
                    parts[key.strip()] = val.strip()
            return PrimerSchemeVendor(
                **_expand_short_keys(parts, PrimerSchemeVendor, _VENDOR_KEY_PREFIXES)
            )

        # Fallback to organisation_name only
        return PrimerSchemeVendor(primer_scheme_vendor_name=v)
    raise ValueError(f"Cannot parse vendor: {v}")


def _save_and_rebuild_readme(
    info_path: pathlib.Path, primer_scheme: PrimerScheme, rebuild_plot: bool = False
):
    """Saves the PrimerScheme to info.json and rebuilds the README."""
    # Save info.json
    logger.debug(f"Writing info.json to {info_path}")
    info_bytes = serialize_primer_scheme_json(primer_scheme)
    info_path.write_bytes(info_bytes)

    # Regenerate README
    scheme_dir = info_path.parent
    logger.debug(f"Regenerating README.md in {scheme_dir}")
    generate_readme(scheme_dir, primer_scheme)

    if rebuild_plot:
        logger.debug(f"Ensuring plot output directory in {scheme_dir / 'assets'}")
        (scheme_dir / "assets").mkdir(exist_ok=True)
        logger.debug(f"Rendering primer plot to {scheme_dir / 'assets' / 'primer.svg'}")
        plot_primers(
            scheme_dir / PRIMER_FILE_NAME, scheme_dir / "assets" / "primer.svg"
        )


def create_status_badge(primer_scheme: PrimerScheme) -> str:
    """
    Create a badge for the README.md file
    """
    match primer_scheme.primer_scheme_development_status:
        case PrimerSchemeDevelopmentStatus.VALIDATED:
            color = "green"
        case (
            PrimerSchemeDevelopmentStatus.WITHDRAWN
            | PrimerSchemeDevelopmentStatus.DEPRECATED
        ):
            color = "red"
        case _:
            color = "blue"

    return f"![Generic badge](https://img.shields.io/badge/STATUS-{primer_scheme.primer_scheme_development_status}-{color}.svg)"


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

    with open(path / README_FILE_NAME, "w", encoding="utf-8") as readme:
        readme.write(
            f"# {primer_scheme.primer_scheme_name} {primer_scheme.amplicon_size}bp {primer_scheme.primer_scheme_version}\n\n"
        )
        # Add the status badge
        readme.write(f"{create_status_badge(primer_scheme)}\n\n")

        # Add citation if present
        if primer_scheme.citation and primer_scheme.citation is not None:
            for cit in primer_scheme.citation:
                readme.write(f"> If you use this scheme please cite: {cit}\n\n")

        if (
            primer_scheme.primer_scheme_details
            and primer_scheme.primer_scheme_details is not None
        ):
            readme.write("## Notes\n\n")
            for note in primer_scheme.primer_scheme_details:
                readme.write(note + "\n\n")

        readme.write("## Metadata\n\n")
        if primer_scheme.primer_scheme_target_organism:
            readme.write("**Target Organisms:**\n")
            for to in primer_scheme.primer_scheme_target_organism:
                to_str = f"- {to.primer_scheme_target_organism_name or ''}"
                if to.primer_scheme_target_organism_ncbi_taxon_id:
                    to_str += (
                        f" (Tax ID: {to.primer_scheme_target_organism_ncbi_taxon_id})"
                    )
                readme.write(f"{to_str}\n")
            readme.write("\n")

        if primer_scheme.primer_scheme_derived_from:
            readme.write(
                f"**Derived from:** {primer_scheme.primer_scheme_derived_from}\n\n"
            )

        if primer_scheme.primer_scheme_contributor:
            readme.write("## Contributors\n\n")
            for contributor in primer_scheme.primer_scheme_contributor:
                contrib_str = f"- {contributor.primer_scheme_contributor_name}"
                if contributor.primer_scheme_contributor_email:
                    contrib_str += f" <{contributor.primer_scheme_contributor_email}>"
                if contributor.primer_scheme_contributor_orcid:
                    contrib_str += (
                        f" (ORCID: {contributor.primer_scheme_contributor_orcid})"
                    )
                readme.write(f"{contrib_str}\n")
            readme.write("\n")

        if primer_scheme.primer_scheme_vendor:
            readme.write("## Vendors\n\n")
            for vendor in primer_scheme.primer_scheme_vendor:
                vendor_str = f"- {vendor.primer_scheme_vendor_name}"
                if vendor.primer_scheme_vendor_kit_name:
                    vendor_str += f": {vendor.primer_scheme_vendor_kit_name}"
                if vendor.primer_scheme_vendor_url:
                    vendor_str += f" ([Website]({vendor.primer_scheme_vendor_url}))"
                readme.write(f"{vendor_str}\n")
            readme.write("\n")

        readme.write("## Overviews\n\n")
        readme.write(
            '<div style="width: 100%;"><img src="assets/primer.svg" style="width: 100%;" alt="Click to see the source"></div>\n\n'
        )

        readme.write("## Details\n\n")

        # Write the details into the readme
        details_json = serialize_primer_scheme_json(primer_scheme).decode("utf-8")
        readme.write(f"""```json\n{details_json}\n```\n\n""")

        if primer_scheme.primer_scheme_license and (
            footer := LICENSE_FOOTERS.get(primer_scheme.primer_scheme_license)
        ):
            readme.write(footer)


def parse_algorithm(v: Any) -> Optional[PrimerSchemeGenerator]:
    if v is None:
        return None
    if isinstance(v, PrimerSchemeGenerator):
        return v
    if isinstance(v, dict):
        return PrimerSchemeGenerator(**v)
    if isinstance(v, str):
        if ":" in v:
            name, version = v.split(":", 1)
            return PrimerSchemeGenerator(
                primer_scheme_generator_name=name,
                primer_scheme_generator_version=version,
            )
        return PrimerSchemeGenerator(primer_scheme_generator_name=v)
    raise ValueError(f"Cannot parse algorithm: {v}")


def parse_target_organism_single(v: Any) -> PrimerSchemeTargetOrganism:
    if isinstance(v, PrimerSchemeTargetOrganism):
        return v
    if isinstance(v, dict):
        return PrimerSchemeTargetOrganism(
            **_expand_short_keys(
                v, PrimerSchemeTargetOrganism, _TARGET_ORGANISM_KEY_PREFIXES
            )
        )
    if isinstance(v, str):
        # Try JSON first
        try:
            data = json.loads(v)
            if isinstance(data, dict):
                return PrimerSchemeTargetOrganism(
                    **_expand_short_keys(
                        data, PrimerSchemeTargetOrganism, _TARGET_ORGANISM_KEY_PREFIXES
                    )
                )
        except json.JSONDecodeError:
            pass

        # Key-value parsing
        if "=" in v:
            parts = {}
            for part in v.split(","):
                if "=" in part:
                    key, val = part.split("=", 1)
                    parts[key.strip()] = val.strip()
            return PrimerSchemeTargetOrganism(
                **_expand_short_keys(
                    parts, PrimerSchemeTargetOrganism, _TARGET_ORGANISM_KEY_PREFIXES
                )
            )

        # If it looks like an int, assume it's a tax id
        if v.isdigit():
            return PrimerSchemeTargetOrganism(
                primer_scheme_target_organism_ncbi_taxon_id=v
            )

        # Otherwise assume common name
        return PrimerSchemeTargetOrganism(primer_scheme_target_organism_name=v)
    raise ValueError(f"Cannot parse target organism: {v}")


def parse_target_organisms_pydantic(v: Any) -> List[PrimerSchemeTargetOrganism]:
    if isinstance(v, list):
        return [parse_target_organism_single(x) for x in v]
    if isinstance(v, (str, dict, PrimerSchemeTargetOrganism)):
        return [parse_target_organism_single(v)]
    return v


def parse_vendors_pydantic(v: Any) -> Optional[List[PrimerSchemeVendor]]:
    if v is None:
        return None
    if isinstance(v, list):
        return [parse_vendor_single(x) for x in v]
    if isinstance(v, (str, dict, PrimerSchemeVendor)):
        return [parse_vendor_single(v)]
    return v


def _normalize_license(v: Any) -> Any:
    """Case-insensitive match against valid SPDX values; pass through PrimerSchemeLicense instances."""
    if isinstance(v, PrimerSchemeLicense):
        return v.value
    if isinstance(v, str):
        for member in PrimerSchemeLicense:
            if v.lower() == member.value.lower():
                return member.value
    return v


_FIELD_PREFIX = "primer_scheme_"


def _strip_field_prefix_name_transform(name: str) -> str:
    """Strip the redundant 'primer_scheme_' prefix from CLI flag names.

    The Pydantic model retains the full field names (e.g. `primer_scheme_name`)
    for schema/data-file compatibility, but typing `--primer-scheme-name` on
    every flag is redundant, so the CLI drops the prefix (e.g. `--name`).
    """
    if name.startswith(_FIELD_PREFIX):
        name = name[len(_FIELD_PREFIX) :]
    return default_name_transform(name)


class CLIPrimerScheme(PrimerScheme):
    schema_version: Annotated[str, Parameter(parse=False)] = SCHEMA_VERSION
    primer_scheme_identifier: Annotated[Optional[str], Parameter(parse=False)] = None
    primer_scheme_contributor: Annotated[  # type: ignore
        List[PrimerSchemeContributor],
        BeforeValidator(parse_contributors_pydantic),
        Parameter(
            help="Individuals, organisations, or institutions that have contributed to the development. e.g. `name=Alice Smith,email=alice@example.org,orcid=0000-0001-2345-6789`"
        ),
    ]
    primer_scheme_target_organism: Annotated[  # type: ignore
        List[PrimerSchemeTargetOrganism],
        BeforeValidator(parse_target_organisms_pydantic),
        Parameter(
            help="The organism(s) targeted by this primer scheme. e.g. `name=SARS-CoV-2,ncbi_taxon_id=2697049`"
        ),
    ]
    primer_scheme_vendor: Annotated[
        Optional[List[PrimerSchemeVendor]],
        BeforeValidator(parse_vendors_pydantic),
        Parameter(
            help="Vendors where one can purchase the primers or a kit containing them. e.g. `name=IDT,kit_name=10011442,url=https://example.com`"
        ),
    ] = None
    primer_scheme_generator: Annotated[
        Optional[PrimerSchemeGenerator], Parameter(parse=False)
    ] = None
    # Don't expose the checksums to cli
    primer_scheme_checksums: Annotated[
        PrimerSchemeChecksums | None, Parameter(parse=False)
    ] = None
    # Override with Literal so Cyclopts displays proper SPDX strings instead of mangled enum names
    primer_scheme_license: Annotated[  # type: ignore
        Optional[_LicenseLiteral],
        BeforeValidator(_normalize_license),
    ] = PrimerSchemeLicense.CC_BY_SA_4FULL_STOP0.value
    # show_default=False: cyclopts' Enum default-rendering assumes a non-None
    # member (`default.name`), which crashes --help for an Optional[Enum]
    # whose default is None. Neither field has a meaningful default to show.
    primer_scheme_application: Annotated[
        Optional[PrimerSchemeApplication], Parameter(show_default=False)
    ] = None
    primer_scheme_scope: Annotated[
        Optional[PrimerSchemeScope], Parameter(show_default=False)
    ] = None
    primer_scheme_creation_date: Annotated[
        date,
        Parameter(help="Date the primer scheme was originally created by its authors"),
    ]
    primer_scheme_submission_date: Annotated[
        date,
        Parameter(help="Date the scheme was added to this registry [default: today]"),
    ] = Field(default_factory=date.today)

    @field_validator("primer_scheme_target_organism")
    def validate_target_organisms(cls, v):
        for to in v:
            if (
                not to.primer_scheme_target_organism_name
                and not to.primer_scheme_target_organism_ncbi_taxon_id
            ):
                raise ValueError(
                    "PrimerSchemeTargetOrganism must have at least one of "
                    "'primer_scheme_target_organism_name' or 'primer_scheme_target_organism_ncbi_taxon_id'"
                )
        return v

    @model_validator(mode="before")
    @classmethod
    def uppercase_enums(cls, data: Any) -> Any:
        if isinstance(data, dict):
            # Uppercase status if it's a string
            if "primer_scheme_development_status" in data and isinstance(
                data["primer_scheme_development_status"], str
            ):
                data["primer_scheme_development_status"] = data[
                    "primer_scheme_development_status"
                ].upper()
        return data


# Keep the full `--primer-scheme-*` flag as a working alias alongside the
# shortened name produced by _strip_field_prefix_name_transform, so existing
# scripts/docs written against the full flag names keep working.
for _field_name, _field_info in CLIPrimerScheme.model_fields.items():
    if _field_name.startswith(_FIELD_PREFIX):
        _field_info.metadata.append(
            Parameter(alias="--" + _field_name.replace("_", "-"))
        )
del _field_name, _field_info


@app.command
def create(
    cli_ps: Annotated[
        CLIPrimerScheme,
        Parameter(name="*", name_transform=_strip_field_prefix_name_transform),
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
        cli_ps.primer_scheme_generator = parse_algorithm(algorithm)
        logger.debug(
            f"Parsed algorithm '{algorithm}' -> PrimerSchemeGenerator({cli_ps.primer_scheme_generator})"
        )

    # Convert to base PrimerScheme to ensure strict adherence to the schema
    ps = PrimerScheme.model_validate(cli_ps.model_dump())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    _headers, bedlines = BedLineParser.from_file(str(bed_path))
    bedlines = sort_bedlines(bedlines)
    logger.debug(f"Loaded and sorted bedlines from {bed_path}")

    # Create a directory to store the new scheme in.
    output_dir = (
        primer_schemes_path
        / ps.primer_scheme_name
        / str(ps.amplicon_size)
        / ps.primer_scheme_version
    )
    # Defense-in-depth: name/version are pattern-constrained, but this
    # confirms the actual resolved path stays under primer_schemes_path
    # regardless of what produced it.
    output_dir = ensure_path_within(primer_schemes_path, output_dir)
    if output_dir.exists():
        raise ValueError(f"Output directory already exists: {output_dir}")

    logger.debug(f"Creating scheme at {output_dir}")

    # Use a tmp dir to ensure atomic
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = pathlib.Path(tmp_dir)
        tmp_version_level = tmp_path / ps.primer_scheme_version
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

        # Generate checksums
        ps.primer_scheme_checksums = PrimerSchemeChecksums(
            primer_scheme_sha256=sha256_checksum(tmp_version_level / PRIMER_FILE_NAME),
            reference_sequence_sha256=sha256_checksum(
                tmp_version_level / REFERENCE_FILE_NAME
            ),
        )
        logger.debug(
            f"Generated checksums for {PRIMER_FILE_NAME} ({ps.primer_scheme_checksums.primer_scheme_sha256})"
            f" and {REFERENCE_FILE_NAME} ({ps.primer_scheme_checksums.reference_sequence_sha256})"
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
        PrimerSchemeContributor,
        Parameter(name="*", converter=parse_contributor_single),
    ],
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Add a contributor to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx is not None:
        logger.debug(f"Inserting contributor at idx={idx}: {contributor}")
        ps.primer_scheme_contributor = [
            *ps.primer_scheme_contributor[:idx],
            contributor,
            *ps.primer_scheme_contributor[idx:],
        ]
        actual_idx = idx
    else:
        logger.debug(f"Appending contributor: {contributor}")
        ps.primer_scheme_contributor = [*ps.primer_scheme_contributor, contributor]
        actual_idx = len(ps.primer_scheme_contributor) - 1
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
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx >= len(ps.primer_scheme_contributor):
        raise ValueError(
            f"Index {idx} out of range for contributors in {info_path}. "
            f"Valid range is 0..{len(ps.primer_scheme_contributor) - 1}."
        )
    if len(ps.primer_scheme_contributor) == 1:
        raise ValueError(
            f"Cannot remove the only contributor from {scheme_label}. "
            "At least one contributor is required."
        )
    removed = ps.primer_scheme_contributor[idx]
    logger.debug(f"Removing contributor at idx={idx}: {removed}")
    ps.primer_scheme_contributor = [
        c for i, c in enumerate(ps.primer_scheme_contributor) if i != idx
    ]
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
        PrimerSchemeContributor,
        Parameter(name="*", converter=parse_contributor_single),
    ],
):
    """Update a contributor at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if idx >= len(ps.primer_scheme_contributor):
        raise ValueError(
            f"Index {idx} out of range for contributors in {info_path}. "
            f"Valid range is 0..{len(ps.primer_scheme_contributor) - 1}."
        )
    previous = ps.primer_scheme_contributor[idx]
    logger.debug(f"Updating contributor at idx={idx}: {previous} -> {contributor}")
    ps.primer_scheme_contributor[idx] = contributor
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
        PrimerSchemeVendor,
        Parameter(name="*", converter=parse_vendor_single),
    ],
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Add a vendor to the scheme."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if ps.primer_scheme_vendor is None:
        ps.primer_scheme_vendor = []
    if idx is not None:
        logger.debug(f"Inserting vendor at idx={idx}: {vendor}")
        ps.primer_scheme_vendor = [
            *ps.primer_scheme_vendor[:idx],
            vendor,
            *ps.primer_scheme_vendor[idx:],
        ]
        actual_idx = idx
    else:
        logger.debug(f"Appending vendor: {vendor}")
        ps.primer_scheme_vendor = [*ps.primer_scheme_vendor, vendor]
        actual_idx = len(ps.primer_scheme_vendor) - 1
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
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if not ps.primer_scheme_vendor or idx >= len(ps.primer_scheme_vendor):
        max_idx = len(ps.primer_scheme_vendor) - 1 if ps.primer_scheme_vendor else -1
        raise ValueError(
            f"Index {idx} out of range for vendors in {info_path}. "
            f"Valid range is 0..{max_idx}."
        )
    removed = ps.primer_scheme_vendor[idx]
    logger.debug(f"Removing vendor at idx={idx}: {removed}")
    ps.primer_scheme_vendor.pop(idx)
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
        PrimerSchemeVendor,
        Parameter(name="*", converter=parse_vendor_single),
    ],
):
    """Update a vendor at a specific index."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if not ps.primer_scheme_vendor or idx >= len(ps.primer_scheme_vendor):
        max_idx = len(ps.primer_scheme_vendor) - 1 if ps.primer_scheme_vendor else -1
        raise ValueError(
            f"Index {idx} out of range for vendors in {info_path}. "
            f"Valid range is 0..{max_idx}."
        )
    previous = ps.primer_scheme_vendor[idx]
    logger.debug(f"Updating vendor at idx={idx}: {previous} -> {vendor}")
    ps.primer_scheme_vendor[idx] = vendor
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated vendors for {scheme_label}: idx {idx} {previous} -> {vendor}")


@modify_app.command
def update_license(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    license: Annotated[
        _LicenseLiteral,
        Parameter(converter=_normalize_license),
    ],
):
    """Update the scheme license."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    previous = ps.primer_scheme_license
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    logger.debug(f"Updating license: {previous} -> {license}")
    ps.primer_scheme_license = license
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated license for {scheme_label}: {previous} -> {license}")


@modify_app.command
def update_status(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    status: PrimerSchemeDevelopmentStatus,
):
    """Update the scheme status."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    previous = ps.primer_scheme_development_status
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    logger.debug(f"Updating status: {previous} -> {status}")
    ps.primer_scheme_development_status = status
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated status for {scheme_label}: {previous} -> {status}")


@modify_app.command
def update_date_created(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    date_created: date,
):
    """Update the date the primer scheme was originally created."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    previous = ps.primer_scheme_creation_date
    ps.primer_scheme_creation_date = date_created
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated date_created: {previous} -> {date_created}")


@modify_app.command
def update_date_added(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    date_added: date,
):
    """Update the date the scheme was added to the registry."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    previous = ps.primer_scheme_submission_date
    ps.primer_scheme_submission_date = date_added
    _save_and_rebuild_readme(info_path, ps)
    logger.info(f"Updated date_added: {previous} -> {date_added}")


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
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    if len(ps.primer_scheme_target_organism) == 1:
        raise ValueError(
            f"Cannot remove the only target organism from {scheme_label}. "
            "At least one target organism is required."
        )
    if idx >= len(ps.primer_scheme_target_organism):
        raise ValueError(
            f"Index {idx} out of range for target_organisms in {info_path}. "
            f"Valid range is 0..{len(ps.primer_scheme_target_organism) - 1}."
        )
    removed = ps.primer_scheme_target_organism[idx]
    logger.debug(f"Removing target_organism at idx={idx}: {removed}")
    ps.primer_scheme_target_organism = [
        to for i, to in enumerate(ps.primer_scheme_target_organism) if i != idx
    ]
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
    target_organism: Annotated[
        Optional[PrimerSchemeTargetOrganism], Parameter(name="*")
    ] = None,
    idx: Annotated[None | int, Parameter(validator=validators.Number(gte=0))] = None,
):
    """Adds a target organism at a specific index."""
    if target_organism is None:
        target_organism = PrimerSchemeTargetOrganism()

    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")

    # append
    if idx is None:
        idx = len(ps.primer_scheme_target_organism)

    logger.debug(f"Adding target_organism at idx={idx}: {target_organism}")
    ps.primer_scheme_target_organism = [
        *ps.primer_scheme_target_organism[:idx],
        target_organism,
        *ps.primer_scheme_target_organism[idx:],
    ]
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
    algorithm: PrimerSchemeGenerator,
):
    """Update the algorithm."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
    previous = ps.primer_scheme_generator
    logger.debug(f"Loaded scheme {scheme_label} from {info_path}")
    logger.debug(f"Updating algorithm: {previous} -> {algorithm}")
    ps.primer_scheme_generator = algorithm
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
    base_url: Annotated[
        str,
        Parameter(
            help="The URL source at which the primer schemes can be found. i.e `https://github.com/pha4ge/primer-schemes/main/v1b/schemes`",
        ),
    ] = "",
    output_path: Annotated[
        pathlib.Path,
        Parameter(
            validator=validators.Path(exists=True, dir_okay=True, file_okay=False),
            help=f"The directory to write the {INDEX_FILE_NAME} and {INDEX_FILE_NAME}.gz",
        ),
    ] = pathlib.Path("."),
    source_commit: Annotated[
        Optional[str],
        Parameter(
            help="Commit SHA of the primer schemes repository this index is built from, e.g. $GITHUB_SHA",
        ),
    ] = None,
):
    """Build a JSON index of all primer schemes in a directory"""
    # Read in current index
    if index_path is not None:
        psi = PrimerSchemeIndex.model_validate_json(index_path.read_text())
    else:
        psi = PrimerSchemeIndex()

    if source_commit:
        psi.source_commit = source_commit

    # Sanitise the base_url
    base_url = base_url.strip("/")

    # find all primer schemes
    ps = []
    for ps_info in find_all_info_json(primer_schemes_path):
        logger.debug(f"found {ps_info}")
        ps.append(PrimerScheme.model_validate_json(ps_info.read_text()))
    update_index(ps, psi, base_url=base_url)

    # Ensure schemes is marked as set for exclude_unset=True
    psi.primerschemes = psi.primerschemes

    index_str = psi.model_dump_json(
        exclude_unset=True, exclude_none=True, ensure_ascii=True
    )

    # Write out the text and compressed index
    (output_path / INDEX_FILE_NAME).write_text(index_str)
    (output_path / (INDEX_FILE_NAME + ".gz")).write_bytes(
        gzip.compress(index_str.encode("utf-8"))
    )
    logger.debug(f"wrote {INDEX_FILE_NAME} to `{output_path}`")


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
            scheme_label = (
                f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
            )
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
        scheme_label = (
            f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
        )
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
    if primer_scheme.primer_scheme_name != name:
        logger.debug(
            f"Syncing scheme name from {primer_scheme.primer_scheme_name} to {name} for {info_path}"
        )
        primer_scheme.primer_scheme_name = name
        changed = True
    if primer_scheme.amplicon_size != amplicon_size:
        logger.debug(
            f"Syncing amplicon_size from {primer_scheme.amplicon_size} to {amplicon_size} for {info_path}"
        )
        primer_scheme.amplicon_size = amplicon_size
        changed = True
    if primer_scheme.primer_scheme_version != version:
        logger.debug(
            f"Syncing version from {primer_scheme.primer_scheme_version} to {version} for {info_path}"
        )
        primer_scheme.primer_scheme_version = version
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
    scheme_label = (
        f"{ps.primer_scheme_name}/{ps.amplicon_size}/{ps.primer_scheme_version}"
    )
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
    logger.debug("Computing sha256 checksums")
    ps.primer_scheme_checksums = PrimerSchemeChecksums(
        primer_scheme_sha256=sha256_checksum(info_path.parent / PRIMER_FILE_NAME),
        reference_sequence_sha256=sha256_checksum(
            info_path.parent / REFERENCE_FILE_NAME
        ),
    )
    _save_and_rebuild_readme(info_path, ps, rebuild_plot=True)
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


@app.command
def get(
    scheme_id: Annotated[
        Optional[str],
        Parameter(
            help="Scheme identifier, e.g. artic/400/v5.4.2 (required unless --all)"
        ),
    ] = None,
    output: Annotated[
        pathlib.Path,
        Parameter(name=["--output", "-o"], help="Output directory"),
    ] = pathlib.Path("."),
    index: Annotated[
        str,
        Parameter(
            help=f"Path or URL to an {METADATA_FILE_NAME}",
        ),
    ] = DEFAULT_INDEX_URL,
    strict: Annotated[
        bool,
        Parameter(
            name="--strict",
            help="Fail on any index mismatch or pre-existing output directory",
        ),
    ] = False,
    force: Annotated[
        bool,
        Parameter(
            name="--force",
            help="Allow missing or mismatched checksums",
        ),
    ] = False,
    allow_multiple: Annotated[
        bool,
        Parameter(
            name="--allow-multiple",
            help="Allow partial scheme_id and download all matches in parallel",
        ),
    ] = False,
    sanitisation: Annotated[
        SanitisationMode,
        Parameter(
            name="--sanitise",
            help="Sanitisation mode for downloaded files",
        ),
    ] = SanitisationMode.RAW,
    timeout: Annotated[
        float,
        Parameter(
            name="--timeout",
            help="HTTP timeout in seconds",
        ),
    ] = DEFAULT_HTTP_TIMEOUT_SECONDS,
    all_schemes: Annotated[
        bool,
        Parameter(
            name="--all",
            help="Download all schemes in the index",
        ),
    ] = False,
):
    """Download a primer scheme by identifier"""
    psi = load_index(index, timeout=timeout, force=force)
    schemes = resolve_schemes(
        index=psi,
        scheme_id=scheme_id,
        allow_multiple=allow_multiple,
        all_schemes=all_schemes,
    )
    output_dirs = download_schemes(
        schemes=schemes,
        output=output,
        strict=strict,
        force=force,
        sanitisation=sanitisation,
        timeout=timeout,
    )
    if len(output_dirs) == 1:
        logger.info(f"Scheme files written to {output_dirs[0]}")
    else:
        logger.info(f"Scheme files written to {len(output_dirs)} directories")


@app.command
def flatten(
    info_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    output_path: Annotated[
        pathlib.Path,
        Parameter(help="Path to write the single-row CSV to"),
    ] = pathlib.Path("scheme.csv"),
):
    """Flatten a single info.json into a one-row, semicolon-delimited CSV."""
    ps = PrimerScheme.model_validate_json(info_path.read_text())
    row = flatten_scheme(ps)
    with output_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerow(row)
    logger.info(f"Flattened {info_path} to {output_path}")


@app.command
def unflatten(
    csv_path: Annotated[
        pathlib.Path,
        Parameter(validator=validators.Path(exists=True, file_okay=True)),
    ],
    output_path: Annotated[
        pathlib.Path,
        Parameter(help="Path to write the reconstructed info.json to"),
    ],
):
    """Reconstruct an info.json from a one-row CSV produced by `flatten`."""
    with csv_path.open(newline="") as f:
        row = next(csv.DictReader(f))
    ps = unflatten_scheme(row)
    output_path.write_bytes(serialize_primer_scheme_json(ps))
    logger.info(f"Unflattened {csv_path} to {output_path}")


def main():
    app.meta()


if __name__ == "__main__":
    main()
