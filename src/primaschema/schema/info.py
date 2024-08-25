from __future__ import annotations
from enum import Enum
import re
from typing import List, Optional
from pydantic.version import VERSION as PYDANTIC_VERSION

if int(PYDANTIC_VERSION[0]) >= 2:
    from pydantic import BaseModel, ConfigDict, Field, field_validator
else:
    from pydantic import BaseModel, Field

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


class SchemeStatus(str, Enum):
    """
    Status of this amplicon primer scheme
    """

    PUBLISHED = "PUBLISHED"
    DEPRECATED = "DEPRECATED"
    DRAFT = "DRAFT"


class PrimerScheme(ConfiguredBaseModel):
    """
    A tiled amplicon PCR primer scheme definition
    """

    schema_version: str = Field(
        ...,
        description="""The version of the schema used to create this scheme definition""",
    )
    name: str = Field(
        ..., description="""The canonical name of the primer scheme (lowercase)"""
    )
    amplicon_size: int = Field(
        ...,
        description="""The length (in base pairs) of an amplicon in the primer scheme""",
        ge=1,
    )
    version: str = Field(...)
    organism: str = Field(
        ...,
        description="""The organism against which this primer scheme is targeted. Lowercase, e.g. sars-cov-2""",
    )
    source_url: Optional[str] = Field(
        None,
        description="""Source URL of primer scheme BED file, if available, e.g. GitHub repository URL""",
    )
    definition_url: Optional[str] = Field(
        None,
        description="""GitHub URL of PHA4GE compatible primer scheme scheme definition""",
    )
    aliases: Optional[List[str]] = Field(
        default_factory=list, description="""Aliases for primer scheme name"""
    )
    license: Optional[str] = Field(
        None, description="""License under which the primer scheme is distributed"""
    )
    status: Optional[SchemeStatus] = Field(
        "PUBLISHED",
        description="""The status of this primer scheme (e.g. published, deprecated)""",
    )
    derived_from: Optional[str] = Field(
        None,
        description="""Canonical name of the primer scheme from which this scheme was derived""",
    )
    developers: List[str] = Field(
        default_factory=list,
        description="""Persons or organisations responsible for developing the primer scheme""",
    )
    citations: Optional[List[str]] = Field(
        default_factory=list,
        description="""URLs of publications describing the scheme (DOIs preferred when available)""",
    )
    notes: Optional[List[str]] = Field(
        default_factory=list, description="""Notes about the amplicon primer scheme"""
    )
    vendors: Optional[List[Vendor]] = Field(
        default_factory=list,
        description="""Vendors where one can purchase the primers described in the amplicon scheme or a kit containing these primers""",
    )
    masks: Optional[List[Mask]] = Field(
        default_factory=list,
        description="""Regions of the reference genome that should be masked out with N""",
    )
    primer_checksum: Optional[str] = Field(
        None,
        description="""Checksum for the primer scheme BED file, in format checksum_type:checksum, where checksum_type is lowercase name of checksum generator e.g. primaschema""",
    )
    reference_checksum: Optional[str] = Field(
        None,
        description="""Checksum for the reference FASTA file, in format checksum_type:checksum, where checksum_type is lowercase name of checksum generator e.g. primaschema""",
    )

    @field_validator("name")
    def pattern_name(cls, v):
        pattern = re.compile(r"^[\da-z0-9_.-]+$")
        if isinstance(v, list):
            for element in v:
                if not pattern.match(element):
                    raise ValueError(f"Invalid name format: {element}")
        elif isinstance(v, str):
            if not pattern.match(v):
                raise ValueError(f"Invalid name format: {v}")
        return v

    @field_validator("version")
    def pattern_version(cls, v):
        pattern = re.compile(r"^[\da-z0-9_.-]+$")
        if isinstance(v, list):
            for element in v:
                if not pattern.match(element):
                    raise ValueError(f"Invalid version format: {element}")
        elif isinstance(v, str):
            if not pattern.match(v):
                raise ValueError(f"Invalid version format: {v}")
        return v


class Vendor(ConfiguredBaseModel):
    """
    Vendor of the primers described in the amplicon scheme or a kit containing these primers
    """

    organisation_name: str = Field(..., description="""The name of the vendor""")
    home_page: Optional[str] = Field(
        None, description="""A link to the home page of the vendor"""
    )
    kit_name: Optional[str] = Field(
        None, description="""Vendor specific kit name for primer kit"""
    )


class Mask(ConfiguredBaseModel):
    """
    A region to mask out, with zero-based, half open coordinates
    """

    reference: str = Field(..., description="""Name (ID) of the reference sequence""")
    name: str = Field(..., description="""Name of the region""")
    start: int = Field(..., description="""Start coordinate of the region""", ge=1)
    end: int = Field(..., description="""End coordination of the region""", ge=1)


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
PrimerScheme.model_rebuild()
Vendor.model_rebuild()
Mask.model_rebuild()
