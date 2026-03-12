from typing import Optional

from pydantic import BaseModel, Field, computed_field

from primaschema import METADATA_FILE_NAME, PRIMER_FILE_NAME, REFERENCE_FILE_NAME
from primaschema.schema.info import (
    Checksums,
    ConfiguredBaseModel,
    Contributor,
    PrimerScheme,
    SchemeLicense,
    SchemeStatus,
    SchemeTag,
    TargetOrganism,
)


class IndexPrimerScheme(BaseModel):
    """
    A subset of PrimerScheme for the index, with additional fields.
    """

    # Taken from PrimerScheme
    name: str = Field(
        default=...,
        description="""The canonical name of the primer scheme (lowercase)""",
    )
    amplicon_size: int = Field(
        default=...,
        description="""The length (in base pairs) of an amplicon in the primer scheme""",
        ge=1,
    )
    version: str = Field(
        default=...,
        description="""The semantic version of the scheme (v.{x}.{y}.{z})""",
    )
    contributors: list[Contributor] = Field(
        default=...,
        description="""Individuals, organisations, or institutions that have contributed to the development""",
    )
    target_organisms: list[TargetOrganism] = Field(
        default=...,
        description="""The organism against which this primer scheme is targeted. Lowercase, e.g. sars-cov-2""",
    )
    license: SchemeLicense = Field(
        default=SchemeLicense.CC_BY_SA_4FULL_STOP0,
        description="""License under which the primer scheme is distributed""",
    )
    status: SchemeStatus = Field(
        default=...,
        description="""The status of this primer scheme (e.g. published, deprecated)""",
    )
    tags: Optional[list[SchemeTag]] = Field(
        default=[],
        description="""Tags to describe the primer scheme""",
    )
    derived_from: Optional[str] = Field(
        default=None,
        description="""Canonical name of the primer scheme from which this scheme was derived""",
    )
    checksums: Optional[Checksums] = Field(
        default=None,
        description="""SHA256 checksums for scheme files""",
    )
    base_url: str = Field(default="", exclude=True)

    # New fields
    @property
    def relative_path(self):
        return f"{self.name}/{self.amplicon_size}/{self.version}"

    @computed_field
    def primer_file_url(self) -> str:
        if self.base_url:
            return f"{self.base_url}/{self.relative_path}/{PRIMER_FILE_NAME}"
        return f"{self.relative_path}/{PRIMER_FILE_NAME}"

    @computed_field
    def reference_file_url(self) -> str:
        if self.base_url:
            return f"{self.base_url}/{self.relative_path}/{REFERENCE_FILE_NAME}"
        return f"{self.relative_path}/{REFERENCE_FILE_NAME}"

    @computed_field
    def info_file_url(self) -> str:
        if self.base_url:
            return f"{self.base_url}/{self.relative_path}/{METADATA_FILE_NAME}"
        return f"{self.relative_path}/{METADATA_FILE_NAME}"

    @classmethod
    def from_primer_scheme(
        cls, scheme: PrimerScheme, base_url: str = ""
    ) -> "IndexPrimerScheme":
        """Create an IndexPrimerScheme from a PrimerScheme instance."""
        return cls(
            name=scheme.name,
            amplicon_size=scheme.amplicon_size,
            version=scheme.version,
            contributors=scheme.contributors,
            target_organisms=scheme.target_organisms,
            license=scheme.license or SchemeLicense.CC_BY_SA_4FULL_STOP0,
            status=scheme.status,
            tags=scheme.tags,
            derived_from=scheme.derived_from,
            checksums=scheme.checksums,
            base_url=base_url,
        )


class PrimerSchemeIndex(ConfiguredBaseModel):
    """
    An index of primer schemes.
    """

    primerschemes: dict[str, dict[int, dict[str, IndexPrimerScheme]]] = Field(
        default_factory=dict,
        description="Index of primer schemes structured as {name: {amplicon_size: {version: scheme}}}",
    )

    def add_index_primer_scheme(self, index: IndexPrimerScheme, strict=True):
        # Get or create the substructure
        name_level = self.primerschemes.setdefault(index.name, {})
        amplicon_size_level = name_level.setdefault(index.amplicon_size, {})

        # If the version already exists, check hashes are consistent if strict
        if index.version in amplicon_size_level and strict:
            original = amplicon_size_level[index.version]

            if original.checksums != index.checksums:
                ValueError(
                    f"checksums have changed for {index.relative_path}. Use Strict == False to allow."
                )

            # If file hashes match update entry.
            amplicon_size_level[index.version] = index
        else:
            # new scheme
            amplicon_size_level[index.version] = index

    def remove_index_primer_scheme(self, index: IndexPrimerScheme) -> bool:
        # Early return if not present
        name_level = self.primerschemes.get(index.name)
        if name_level is None:
            return False

        amplicon_size_level = name_level.get(index.amplicon_size)
        if amplicon_size_level is None:
            return False

        if index.version in amplicon_size_level:
            amplicon_size_level.pop(index.version)
            # Prune the tree
            self.prune_index()  # if this is common could do a targeted prune
            return True
        return False

    def prune_index(self):
        """Removes any empty stems"""
        # Iterate over a copy of keys to allow modification during iteration
        for name in list(self.primerschemes.keys()):
            name_level = self.primerschemes[name]

            for amplicon_size in list(name_level.keys()):
                amplicon_size_level = name_level[amplicon_size]

                # If the version dict is empty, remove the amplicon_size key
                if not amplicon_size_level:
                    del name_level[amplicon_size]

            # If the name dict is empty, remove the name key
            if not name_level:
                del self.primerschemes[name]


def create_index(
    primer_schemes: list[PrimerScheme], base_url: str = ""
) -> PrimerSchemeIndex:
    psi = PrimerSchemeIndex()

    for ps in primer_schemes:
        # Convert to index entry
        mps = IndexPrimerScheme.from_primer_scheme(ps, base_url=base_url)
        psi.add_index_primer_scheme(mps)
    return psi


def update_index(
    primer_schemes: list[PrimerScheme],
    index: PrimerSchemeIndex,
    strict: bool = True,
    base_url: str = "",
):
    for ps in primer_schemes:
        # Convert to index entry
        mps = IndexPrimerScheme.from_primer_scheme(ps, base_url=base_url)
        index.add_index_primer_scheme(mps, strict)
