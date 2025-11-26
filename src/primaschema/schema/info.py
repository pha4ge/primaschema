from __future__ import annotations

import re
from enum import Enum
from typing import Any, ClassVar, Optional

from pydantic import BaseModel, ConfigDict, Field, RootModel, field_validator


metamodel_version = "None"
version = "1.0.0-alpha"


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        validate_assignment=True,
        validate_default=True,
        extra="forbid",
        arbitrary_types_allowed=True,
        use_enum_values=True,
        strict=False,
    )
    pass


class LinkMLMeta(RootModel):
    root: dict[str, Any] = {}
    model_config = ConfigDict(frozen=True)

    def __getattr__(self, key: str):
        return getattr(self.root, key)

    def __getitem__(self, key: str):
        return self.root[key]

    def __setitem__(self, key: str, value):
        self.root[key] = value

    def __contains__(self, key: str) -> bool:
        return key in self.root


linkml_meta = LinkMLMeta(
    {
        "default_curi_maps": ["semweb_context"],
        "default_prefix": "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/",
        "default_range": "string",
        "description": "Data model for tiling primer scheme definitions",
        "id": "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme",
        "imports": ["linkml:types"],
        "name": "primer-scheme",
        "prefixes": {
            "GENEPIO": {
                "prefix_prefix": "GENEPIO",
                "prefix_reference": "http://purl.obolibrary.org/obo/GENEPIO_",
            },
            "IAO": {
                "prefix_prefix": "IAO",
                "prefix_reference": "http://purl.obolibrary.org/obo/IAO_",
            },
            "ORCID": {
                "prefix_prefix": "ORCID",
                "prefix_reference": "http://identifiers.org/orcid/",
            },
            "linkml": {
                "prefix_prefix": "linkml",
                "prefix_reference": "https://w3id.org/linkml/",
            },
            "schema": {
                "prefix_prefix": "schema",
                "prefix_reference": "http://schema.org/",
            },
        },
        "source_file": "/Users/bede/Research/git/primaschema/src/primaschema/schema/info.yml",
    }
)


class SchemeStatus(str, Enum):
    """
    Status of this amplicon primer scheme
    """

    DRAFT = "DRAFT"
    TESTED = "TESTED"
    VALIDATED = "VALIDATED"
    DEPRECATED = "DEPRECATED"
    WITHDRAWN = "WITHDRAWN"


class PrimerScheme(ConfiguredBaseModel):
    """
    A tiled amplicon PCR primer scheme definition
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {
            "from_schema": "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme",
            "tree_root": True,
        }
    )

    schema_version: str = Field(
        default=...,
        description="""The version of the schema used to create this scheme definition""",
        json_schema_extra={
            "linkml_meta": {"alias": "schema_version", "domain_of": ["PrimerScheme"]}
        },
    )
    name: str = Field(
        default=...,
        description="""The canonical name of the primer scheme (lowercase)""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "name",
                "domain_of": ["PrimerScheme", "Contributor", "Mask"],
                "slot_uri": "GENEPIO:0001456",
            }
        },
    )
    amplicon_size: int = Field(
        default=...,
        description="""The length (in base pairs) of an amplicon in the primer scheme""",
        ge=1,
        json_schema_extra={
            "linkml_meta": {
                "alias": "amplicon_size",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0001449",
            }
        },
    )
    contributors: list[Contributor] = Field(
        default=...,
        json_schema_extra={
            "linkml_meta": {"alias": "contributors", "domain_of": ["PrimerScheme"]}
        },
    )
    version: str = Field(
        default=...,
        json_schema_extra={
            "linkml_meta": {"alias": "version", "domain_of": ["PrimerScheme"]}
        },
    )
    organism: str = Field(
        default=...,
        description="""The organism against which this primer scheme is targeted. Lowercase, e.g. sars-cov-2""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "organism",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100682",
            }
        },
    )
    source_url: Optional[str] = Field(
        default=None,
        description="""Source URL of primer scheme BED file, if available, e.g. GitHub repository URL""",
        json_schema_extra={
            "linkml_meta": {"alias": "source_url", "domain_of": ["PrimerScheme"]}
        },
    )
    definition_url: Optional[str] = Field(
        default=None,
        description="""GitHub URL of PHA4GE compatible primer scheme scheme definition""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "definition_url",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100683",
            }
        },
    )
    aliases: Optional[list[str]] = Field(
        default=None,
        description="""Aliases for primer scheme name""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "aliases",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100670",
            }
        },
    )
    license: Optional[str] = Field(
        default=None,
        description="""License under which the primer scheme is distributed""",
        json_schema_extra={
            "linkml_meta": {"alias": "license", "domain_of": ["PrimerScheme"]}
        },
    )
    status: SchemeStatus = Field(
        default=...,
        description="""The status of this primer scheme (e.g. published, deprecated)""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "status",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100681",
            }
        },
    )
    derived_from: Optional[str] = Field(
        default=None,
        description="""Canonical name of the primer scheme from which this scheme was derived""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "derived_from",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100671",
            }
        },
    )
    citations: Optional[list[str]] = Field(
        default=None,
        description="""URLs of publications describing the scheme (DOIs preferred when available)""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "citations",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "IAO:0000301",
            }
        },
    )
    notes: Optional[list[str]] = Field(
        default=None,
        description="""Notes about the amplicon primer scheme""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "notes",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100672",
            }
        },
    )
    vendors: Optional[list[Vendor]] = Field(
        default=None,
        description="""Vendors where one can purchase the primers described in the amplicon scheme or a kit containing these primers""",
        json_schema_extra={
            "linkml_meta": {"alias": "vendors", "domain_of": ["PrimerScheme"]}
        },
    )
    masks: Optional[list[Mask]] = Field(
        default=None,
        description="""Regions of the reference genome that should be masked out with N""",
        json_schema_extra={
            "linkml_meta": {"alias": "masks", "domain_of": ["PrimerScheme"]}
        },
    )
    primer_checksum: Optional[str] = Field(
        default=None,
        description="""Checksum for the primer scheme BED file, in format checksum_type:checksum, where checksum_type is lowercase name of checksum generator e.g. primaschema""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "primer_checksum",
                "domain_of": ["PrimerScheme"],
                "slot_uri": "GENEPIO:0100675",
            }
        },
    )
    reference_checksum: Optional[str] = Field(
        default=None,
        description="""Checksum for the reference FASTA file, in format checksum_type:checksum, where checksum_type is lowercase name of checksum generator e.g. primaschema""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "reference_checksum",
                "domain_of": ["PrimerScheme"],
            }
        },
    )

    @field_validator("name")
    def pattern_name(cls, v):
        pattern = re.compile(r"^[\da-z0-9_.-]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid name format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid name format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator("version")
    def pattern_version(cls, v):
        pattern = re.compile(r"^v\d+\.\d+\.\d+(-[a-z0-9]+)?$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid version format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid version format: {v}"
            raise ValueError(err_msg)
        return v


class Vendor(ConfiguredBaseModel):
    """
    Vendor of the primers described in the amplicon scheme or a kit containing these primers
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {
            "class_uri": "GENEPIO:0100674",
            "from_schema": "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme",
        }
    )

    organisation_name: str = Field(
        default=...,
        description="""The name of the vendor""",
        json_schema_extra={
            "linkml_meta": {"alias": "organisation_name", "domain_of": ["Vendor"]}
        },
    )
    home_page: Optional[str] = Field(
        default=None,
        description="""A link to the home page of the vendor""",
        json_schema_extra={
            "linkml_meta": {"alias": "home_page", "domain_of": ["Vendor"]}
        },
    )
    kit_name: Optional[str] = Field(
        default=None,
        description="""Vendor specific kit name for primer kit""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "kit_name",
                "domain_of": ["Vendor"],
                "slot_uri": "GENEPIO:0100693",
            }
        },
    )


class Contributor(ConfiguredBaseModel):
    """
    Person or organisation who contributed to primerscheme development
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {
            "class_uri": "IAO:contributor",
            "from_schema": "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme",
        }
    )

    name: str = Field(
        default=...,
        description="""The name of the person or organisation""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "name",
                "domain_of": ["PrimerScheme", "Contributor", "Mask"],
                "slot_uri": "IAO:0000590",
            }
        },
    )
    orcid: Optional[str] = Field(
        default=None,
        description="""ORCID ('Open Researcher and Contributor IDentifier') of a person""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "orcid",
                "domain_of": ["Contributor"],
                "slot_uri": "IAO:0000708",
            }
        },
    )
    email: Optional[str] = Field(
        default=None,
        description="""Contact email""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "email",
                "domain_of": ["Contributor"],
                "slot_uri": "IAO:0000429",
            }
        },
    )

    @field_validator("email")
    def pattern_email(cls, v):
        pattern = re.compile(r"^[\w\-\.]+@([\w-]+\.)+[\w-]{2,}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid email format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid email format: {v}"
            raise ValueError(err_msg)
        return v


class Mask(ConfiguredBaseModel):
    """
    A region to mask out, with zero-based, half open coordinates
    """

    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta(
        {
            "from_schema": "https://github.com/pha4ge/primer-schemes/schemas/primer-scheme"
        }
    )

    reference: str = Field(
        default=...,
        description="""Name (ID) of the reference sequence""",
        json_schema_extra={
            "linkml_meta": {"alias": "reference", "domain_of": ["Mask"]}
        },
    )
    name: str = Field(
        default=...,
        description="""Name of the region""",
        json_schema_extra={
            "linkml_meta": {
                "alias": "name",
                "domain_of": ["PrimerScheme", "Contributor", "Mask"],
            }
        },
    )
    start: int = Field(
        default=...,
        description="""Start coordinate of the region""",
        ge=1,
        json_schema_extra={"linkml_meta": {"alias": "start", "domain_of": ["Mask"]}},
    )
    end: int = Field(
        default=...,
        description="""End coordination of the region""",
        ge=1,
        json_schema_extra={"linkml_meta": {"alias": "end", "domain_of": ["Mask"]}},
    )


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
PrimerScheme.model_rebuild()
Vendor.model_rebuild()
Contributor.model_rebuild()
Mask.model_rebuild()
