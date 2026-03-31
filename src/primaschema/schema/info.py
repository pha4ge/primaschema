from __future__ import annotations

import re
import sys
from datetime import (
    date,
    datetime,
    time
)
from decimal import Decimal
from enum import Enum
from typing import (
    Any,
    ClassVar,
    Literal,
    Optional,
    Union
)

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    SerializationInfo,
    SerializerFunctionWrapHandler,
    field_validator,
    model_serializer
)


metamodel_version = "None"
version = "1.0.0-alpha"


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        serialize_by_alias = True,
        validate_by_name = True,
        validate_assignment = True,
        validate_default = True,
        extra = "forbid",
        arbitrary_types_allowed = True,
        use_enum_values = True,
        strict = False,
    )

    @model_serializer(mode='wrap', when_used='unless-none')
    def treat_empty_lists_as_none(
            self, handler: SerializerFunctionWrapHandler,
            info: SerializationInfo) -> dict[str, Any]:
        if info.exclude_none:
            _instance = self.model_copy()
            for field, field_info in type(_instance).model_fields.items():
                if getattr(_instance, field) == [] and not(
                        field_info.is_required()):
                    setattr(_instance, field, None)
        else:
            _instance = self
        return handler(_instance, info)



class LinkMLMeta(RootModel):
    root: dict[str, Any] = {}
    model_config = ConfigDict(frozen=True)

    def __getattr__(self, key:str):
        return getattr(self.root, key)

    def __getitem__(self, key:str):
        return self.root[key]

    def __setitem__(self, key:str, value):
        self.root[key] = value

    def __contains__(self, key:str) -> bool:
        return key in self.root


linkml_meta = None

class SchemeStatus(str, Enum):
    """
    Status of this amplicon primer scheme
    """
    DRAFT = "DRAFT"
    TESTED = "TESTED"
    VALIDATED = "VALIDATED"
    DEPRECATED = "DEPRECATED"
    WITHDRAWN = "WITHDRAWN"


class SchemeLicense(str, Enum):
    """
    License under which the primer scheme is distributed
    """
    CC0_1FULL_STOP0 = "CC0-1.0"
    CC_BY_4FULL_STOP0 = "CC-BY-4.0"
    CC_BY_SA_4FULL_STOP0 = "CC-BY-SA-4.0"
    CC_BY_NC_4FULL_STOP0 = "CC-BY-NC-4.0"
    CC_BY_NC_SA_4FULL_STOP0 = "CC-BY-NC-SA-4.0"
    CC_BY_ND_4FULL_STOP0 = "CC-BY-ND-4.0"
    CC_BY_NC_ND_4FULL_STOP0 = "CC-BY-NC-ND-4.0"


class SchemeTag(str, Enum):
    """
    Tag for this primer scheme
    """
    WASTE_WATER = "WASTE-WATER"
    CLINICAL = "CLINICAL"
    FULL_GENOME = "FULL-GENOME"
    MULTI_TARGET = "MULTI-TARGET"
    PANEL = "PANEL"
    QPCR = "QPCR"



class PrimerScheme(ConfiguredBaseModel):
    """
    A tiled amplicon PCR primer scheme definition
    """
    schema_version: str = Field(default=..., description="""The version of the schema used to create this scheme definition""")
    name: str = Field(default=..., description="""The canonical name of the primer scheme (lowercase)""")
    amplicon_size: int = Field(default=..., description="""The length (in base pairs) of an amplicon in the primer scheme""", ge=1)
    version: str = Field(default=..., description="""The semantic version of the scheme (v.{x}.{y}.{z})""")
    contributors: list[Contributor] = Field(default=..., description="""Individuals, organisations, or institutions that have contributed to the development""", min_length=1)
    target_organisms: list[TargetOrganism] = Field(default=..., description="""The organism against which this primer scheme is targeted.""", min_length=1)
    aliases: Optional[list[str]] = Field(default=[], description="""Aliases for primer scheme name""")
    license: Optional[SchemeLicense] = Field(default='CC-BY-SA-4.0', description="""License under which the primer scheme is distributed""")
    status: SchemeStatus = Field(default=..., description="""The status of this primer scheme (e.g. published, deprecated)""")
    tags: Optional[list[SchemeTag]] = Field(default=[], description="""Tags to describe the primer scheme""")
    derived_from: Optional[str] = Field(default=None, description="""Canonical name of the primer scheme from which this scheme was derived""")
    citations: Optional[list[str]] = Field(default=[], description="""URLs of publications describing the scheme (DOIs preferred when available)""")
    notes: Optional[list[str]] = Field(default=[], description="""Notes about the amplicon primer scheme""")
    vendors: Optional[list[Vendor]] = Field(default=[], description="""Vendors where one can purchase the primers described in the amplicon scheme or a kit containing these primers""")
    algorithm: Optional[Algorithm] = Field(default=None, description="""The algorithm (if any) used to generate this primer scheme""")
    checksums: Optional[Checksums] = Field(default=None, description="""SHA256 checksums for scheme files""")
    ref_selections: Optional[list[RefSelection]] = Field(default=[], description="""Optional reference selections""")
    date_created: Optional[date] = Field(default=None, description="""Date the primer scheme was originally created by its authors""")
    date_added: Optional[date] = Field(default=None, description="""Date the scheme was added to this registry""")

    @field_validator('name')
    def pattern_name(cls, v):
        pattern=re.compile(r"^[\da-z0-9_.-]+$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid name format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid name format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('version')
    def pattern_version(cls, v):
        pattern=re.compile(r"^v\d+\.\d+\.\d+(-[a-z0-9]+)?$")
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
    organisation_name: str = Field(default=..., description="""The name of the vendor""")
    home_page: Optional[str] = Field(default=None, description="""A link to the home page of the vendor""")
    kit_name: Optional[str] = Field(default=None, description="""Vendor specific kit name for primer kit""")


class Contributor(ConfiguredBaseModel):
    """
    Person or organisation who contributed to primer scheme development
    """
    name: str = Field(default=..., description="""The name of the person or organisation""")
    orcid_id: Optional[str] = Field(default=None, description="""ORCID iD ('Open Researcher and Contributor IDentifier') of a person""")
    email: Optional[str] = Field(default=None, description="""Contact email""")

    @field_validator('email')
    def pattern_email(cls, v):
        pattern=re.compile(r"^[a-zA-Z0-9._%+\-]+@[a-zA-Z0-9.\-]+\.[a-zA-Z]{2,}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid email format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid email format: {v}"
            raise ValueError(err_msg)
        return v


class Algorithm(ConfiguredBaseModel):
    """
    Algorithm used to generate the primer scheme
    """
    name: str = Field(default=..., description="""The name of the Algorithm""")
    version: Optional[str] = Field(default=None, description="""The version of the Algorithm""")


class Checksums(ConfiguredBaseModel):
    """
    SHA256 checksums for scheme files
    """
    primer_sha256: Optional[str] = Field(default=None, description="""SHA256 checksum for the primer scheme BED file""")
    reference_sha256: Optional[str] = Field(default=None, description="""SHA256 checksum for the reference FASTA file""")

    @field_validator('primer_sha256')
    def pattern_primer_sha256(cls, v):
        pattern=re.compile(r"^[a-fA-F0-9]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid primer_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid primer_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('reference_sha256')
    def pattern_reference_sha256(cls, v):
        pattern=re.compile(r"^[a-fA-F0-9]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid reference_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid reference_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class RefSelection(ConfiguredBaseModel):
    file_sha256: Optional[str] = Field(default=None, description="""SHA256 checksum for the reference selection file""")
    file_name: Optional[str] = Field(default=None, description="""File name of the reference selection file""")
    chromosome: Optional[str] = Field(default=None, description="""The chromosome in the primerbed file for which this provides reference selection""")

    @field_validator('file_sha256')
    def pattern_file_sha256(cls, v):
        pattern=re.compile(r"^[a-fA-F0-9]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid file_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid file_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class TargetOrganism(ConfiguredBaseModel):
    """
    The organisms targeted by this primer scheme
    """
    common_name: Optional[str] = Field(default=None, description="""The common name of the organism""")
    ncbi_tax_id: Optional[str] = Field(default=None, description="""A stable unique identifier for each taxon (for a species, a family, an order, or any other group) in the NCBI taxonomy database""")


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
PrimerScheme.model_rebuild()
Vendor.model_rebuild()
Contributor.model_rebuild()
Algorithm.model_rebuild()
Checksums.model_rebuild()
RefSelection.model_rebuild()
TargetOrganism.model_rebuild()
