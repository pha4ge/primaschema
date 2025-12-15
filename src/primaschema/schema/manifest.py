from typing import Any, Literal, cast, Optional
from pydantic import Field, create_model, field_validator, computed_field, BaseModel
from primaschema.schema.info import (
    ConfiguredBaseModel,
    PrimerScheme,
    Contributor,
    SchemeLicense,
    SchemeTag,
    SchemeStatus,
    TargetOrganism,
)

from primaschema import PRIMER_FILE_NAME, REFERENCE_FILE_NAME, METADATA_FILE_NAME


class ManifestPrimerScheme(BaseModel):
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
        description="""Tags to describe the primerscheme""",
    )
    derived_from: Optional[str] = Field(
        default=None,
        description="""Canonical name of the primer scheme from which this scheme was derived""",
    )
    primer_file_sha256: Optional[str] = Field(
        default=None,
        description="""SHA256 checksum for the primer scheme BED file""",
    )
    reference_file_sha256: Optional[str] = Field(
        default=None,
        description="""SHA256 checksum for the reference FASTA file""",
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
    ) -> "ManifestPrimerScheme":
        """Create a ManifestPrimerScheme from a PrimerScheme instance."""
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
            primer_file_sha256=scheme.primer_file_sha256,
            reference_file_sha256=scheme.reference_file_sha256,
            base_url=base_url,
        )


class PrimerSchemeIndex(ConfiguredBaseModel):
    """
    An index of primer primerschemes.
    """

    primerschemes: dict[str, dict[int, dict[str, ManifestPrimerScheme]]] = Field(
        default_factory=dict,
        description="Index of primer primerschemes structured as {name: {amplicon_size: {version: scheme}}}",
    )

    def add_manifest_primer_scheme(self, manifest: ManifestPrimerScheme, strict=True):
        # Get or create the substructure
        name_level = self.primerschemes.setdefault(manifest.name, {})
        amplicon_size_level = name_level.setdefault(manifest.amplicon_size, {})

        # If the version already exists, check hashes are consistent if strict
        if manifest.version in amplicon_size_level and strict:
            original_manifest = amplicon_size_level[manifest.version]

            if original_manifest.primer_file_sha256 != manifest.primer_file_sha256:
                ValueError(
                    f"primer_file_sha256 has changed for {manifest.relative_path}. Original ({original_manifest.primer_file_sha256}) -> New ({manifest.primer_file_sha256}). Use Strict == False to allow."
                )
            if (
                original_manifest.reference_file_sha256
                != manifest.reference_file_sha256
            ):
                ValueError(
                    f"primer_file_sha256 has changed for {manifest.relative_path}. Original ({original_manifest.reference_file_sha256}) -> New ({manifest.reference_file_sha256}). Use Strict == False to allow."
                )

            # If files hashes match update manifest.
            amplicon_size_level[manifest.version] = manifest
        else:
            # new scheme
            amplicon_size_level[manifest.version] = manifest

    def remove_manifest_primer_scheme(self, manifest: ManifestPrimerScheme) -> bool:
        # Early return if not present
        name_level = self.primerschemes.get(manifest.name)
        if name_level is None:
            return False

        amplicon_size_level = name_level.get(manifest.amplicon_size)
        if amplicon_size_level is None:
            return False

        if manifest.version in amplicon_size_level:
            amplicon_size_level.pop(manifest.version)
            # Prune the tree
            self.prune_index()  # if this is common could do a targeted prune
            return True
        return False

    def prune_index(self):
        """Removes any stem without a manifest"""
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
    primerprimerschemes: list[PrimerScheme], base_url: str = ""
) -> PrimerSchemeIndex:
    psi = PrimerSchemeIndex()

    for ps in primerprimerschemes:
        # Convert to manifest
        mps = ManifestPrimerScheme.from_primer_scheme(ps, base_url=base_url)
        psi.add_manifest_primer_scheme(mps)
    return psi


def update_index(
    primerprimerschemes: list[PrimerScheme],
    index: PrimerSchemeIndex,
    strict: bool = True,
    base_url: str = "",
):
    for ps in primerprimerschemes:
        # Convert to manifest
        mps = ManifestPrimerScheme.from_primer_scheme(ps, base_url=base_url)
        index.add_manifest_primer_scheme(mps, strict)
