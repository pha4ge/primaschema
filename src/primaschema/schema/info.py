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


linkml_meta = LinkMLMeta({'default_curi_maps': ['semweb_context'],
     'default_prefix': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme/',
     'default_range': 'string',
     'description': 'Data model for tiling primer scheme definitions',
     'id': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme',
     'imports': ['linkml:types'],
     'name': 'primer-scheme',
     'prefixes': {'EDAM': {'prefix_prefix': 'EDAM',
                           'prefix_reference': 'http://edamontology.org/data_'},
                  'GENEPIO': {'prefix_prefix': 'GENEPIO',
                              'prefix_reference': 'http://purl.obolibrary.org/obo/GENEPIO_'},
                  'IAO': {'prefix_prefix': 'IAO',
                          'prefix_reference': 'http://purl.obolibrary.org/obo/IAO_'},
                  'ORCID': {'prefix_prefix': 'ORCID',
                            'prefix_reference': 'http://identifiers.org/orcid/'},
                  'linkml': {'prefix_prefix': 'linkml',
                             'prefix_reference': 'https://w3id.org/linkml/'},
                  'schema': {'prefix_prefix': 'schema',
                             'prefix_reference': 'http://schema.org/'}},
     'source_file': '/Users/kentcg/primaschema/src/primaschema/schema/info.yml'} )

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
    CC_BY_SA_4FULL_STOP0 = "CC-BY-SA-4.0"


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
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme',
         'tree_root': True})

    schema_version: str = Field(default=..., description="""The version of the schema used to create this scheme definition""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    name: str = Field(default=..., description="""The canonical name of the primer scheme (lowercase)""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme', 'Contributor', 'Algorithm'],
         'slot_uri': 'GENEPIO:0001456'} })
    amplicon_size: int = Field(default=..., description="""The length (in base pairs) of an amplicon in the primer scheme""", ge=1, json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'GENEPIO:0001449'} })
    version: str = Field(default=..., description="""The semantic version of the scheme (v.{x}.{y}.{z})""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme', 'Algorithm']} })
    contributors: list[Contributor] = Field(default=..., description="""Individuals, organisations, or institutions that have contributed to the development""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    target_organisms: list[TargetOrganism] = Field(default=..., description="""The organism against which this primer scheme is targeted.""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'EDAM:1869'} })
    aliases: Optional[list[str]] = Field(default=[], description="""Aliases for primer scheme name""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'GENEPIO:0100670'} })
    license: Optional[SchemeLicense] = Field(default='CC-BY-SA-4.0', description="""License under which the primer scheme is distributed""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'ifabsent': 'SchemeLicense(CC-BY-SA-4.0)'} })
    status: SchemeStatus = Field(default=..., description="""The status of this primer scheme (e.g. published, deprecated)""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'GENEPIO:0100681'} })
    tags: Optional[list[SchemeTag]] = Field(default=[], description="""Tags to describe the primerscheme""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    derived_from: Optional[str] = Field(default=None, description="""Canonical name of the primer scheme from which this scheme was derived""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'GENEPIO:0100671'} })
    citations: Optional[list[str]] = Field(default=[], description="""URLs of publications describing the scheme (DOIs preferred when available)""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'IAO:0000301'} })
    notes: Optional[list[str]] = Field(default=[], description="""Notes about the amplicon primer scheme""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'GENEPIO:0100672'} })
    vendors: Optional[list[Vendor]] = Field(default=[], description="""Vendors where one can purchase the primers described in the amplicon scheme or a kit containing these primers""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    algorithm: Optional[Algorithm] = Field(default=None, description="""The algorithm (if any) used to generate this primerscheme""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    primer_checksum: Optional[str] = Field(default=None, description="""Checksum for the primer scheme BED file, in format checksum_type:checksum, where checksum_type is lowercase name of checksum generator e.g. primaschema""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme'], 'slot_uri': 'GENEPIO:0100675'} })
    primer_file_sha256: Optional[str] = Field(default=None, description="""SHA256 checksum for the primer scheme BED file""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    reference_checksum: Optional[str] = Field(default=None, description="""Checksum for the reference FASTA file, in format checksum_type:checksum, where checksum_type is lowercase name of checksum generator e.g. primaschema""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    reference_file_sha256: Optional[str] = Field(default=None, description="""SHA256 checksum for the reference FASTA file""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })
    ref_selections: Optional[list[RefSelection]] = Field(default=[], description="""Optional reference selections""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme']} })

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

    @field_validator('primer_file_sha256')
    def pattern_primer_file_sha256(cls, v):
        pattern=re.compile(r"^[a-fA-F0-9]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid primer_file_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid primer_file_sha256 format: {v}"
            raise ValueError(err_msg)
        return v

    @field_validator('reference_file_sha256')
    def pattern_reference_file_sha256(cls, v):
        pattern=re.compile(r"^[a-fA-F0-9]{64}$")
        if isinstance(v, list):
            for element in v:
                if isinstance(element, str) and not pattern.match(element):
                    err_msg = f"Invalid reference_file_sha256 format: {element}"
                    raise ValueError(err_msg)
        elif isinstance(v, str) and not pattern.match(v):
            err_msg = f"Invalid reference_file_sha256 format: {v}"
            raise ValueError(err_msg)
        return v


class Vendor(ConfiguredBaseModel):
    """
    Vendor of the primers described in the amplicon scheme or a kit containing these primers
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'GENEPIO:0100674',
         'from_schema': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme'})

    organisation_name: str = Field(default=..., description="""The name of the vendor""", json_schema_extra = { "linkml_meta": {'domain_of': ['Vendor']} })
    home_page: Optional[str] = Field(default=None, description="""A link to the home page of the vendor""", json_schema_extra = { "linkml_meta": {'domain_of': ['Vendor']} })
    kit_name: Optional[str] = Field(default=None, description="""Vendor specific kit name for primer kit""", json_schema_extra = { "linkml_meta": {'domain_of': ['Vendor'], 'slot_uri': 'GENEPIO:0100693'} })


class Contributor(ConfiguredBaseModel):
    """
    Person or organisation who contributed to primerscheme development
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'IAO:contributor',
         'from_schema': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme'})

    name: str = Field(default=..., description="""The name of the person or organisation""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme', 'Contributor', 'Algorithm'],
         'slot_uri': 'IAO:0000590'} })
    orcid_id: Optional[str] = Field(default=None, description="""ORCID iD ('Open Researcher and Contributor IDentifier') of a person""", json_schema_extra = { "linkml_meta": {'domain_of': ['Contributor'], 'slot_uri': 'IAO:0000708'} })
    email: Optional[str] = Field(default=None, description="""Contact email""", json_schema_extra = { "linkml_meta": {'domain_of': ['Contributor'], 'slot_uri': 'IAO:0000429'} })

    @field_validator('email')
    def pattern_email(cls, v):
        pattern=re.compile(r"^[\w\-\.]+@([\w-]+\.)+[\w-]{2,}$")
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
    Algorithm used to generate the primerscheme
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'IAO:0000064',
         'from_schema': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme'})

    name: str = Field(default=..., description="""The name of the Algorithm""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme', 'Contributor', 'Algorithm']} })
    version: Optional[str] = Field(default=None, description="""The version of the Algorithm""", json_schema_extra = { "linkml_meta": {'domain_of': ['PrimerScheme', 'Algorithm'], 'slot_uri': 'IAO:0000129'} })


class RefSelection(ConfiguredBaseModel):
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'from_schema': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme'})

    file_sha256: Optional[str] = Field(default=None, description="""SHA256 checksum for the reference selection file""", json_schema_extra = { "linkml_meta": {'domain_of': ['RefSelection']} })
    file_name: Optional[str] = Field(default=None, description="""File name of the reference selection file""", json_schema_extra = { "linkml_meta": {'domain_of': ['RefSelection']} })
    chromosome: Optional[str] = Field(default=None, description="""The chromosome in the primerbed file for which this provides reference selection""", json_schema_extra = { "linkml_meta": {'domain_of': ['RefSelection']} })

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
    The organisms targeted by this primerscheme
    """
    linkml_meta: ClassVar[LinkMLMeta] = LinkMLMeta({'class_uri': 'EDAM:1869',
         'from_schema': 'https://github.com/pha4ge/primer-schemes/schemas/primer-scheme'})

    common_name: Optional[str] = Field(default=None, description="""The common name of the organism""", json_schema_extra = { "linkml_meta": {'domain_of': ['TargetOrganism'], 'slot_uri': 'EDAM:1874'} })
    ncbi_tax_id: Optional[str] = Field(default=None, description="""A stable unique identifier for each taxon (for a species, a family, an order, or any other group) in the NCBI taxonomy database""", json_schema_extra = { "linkml_meta": {'domain_of': ['TargetOrganism'], 'slot_uri': 'EDAM:1179'} })


# Model rebuild
# see https://pydantic-docs.helpmanual.io/usage/models/#rebuilding-a-model
PrimerScheme.model_rebuild()
Vendor.model_rebuild()
Contributor.model_rebuild()
Algorithm.model_rebuild()
RefSelection.model_rebuild()
TargetOrganism.model_rebuild()
