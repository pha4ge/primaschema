from pydantic import (
    BaseModel,
    ConfigDict,
    NonNegativeInt,
    computed_field,
    model_validator,
)
from typing import Dict, List, Literal, Optional, Tuple


class PrimerModel(BaseModel):
    """A primer as represented by a BED record in Primal Scheme v3 format"""

    model_config = ConfigDict(str_strip_whitespace=True)
    chrom: str
    chrom_start: NonNegativeInt
    chrom_end: NonNegativeInt
    name: str
    pool_name: int
    strand: Literal["+", "-"]
    sequence: str

    @computed_field
    @property
    def name_parts(self) -> List[str]:
        parts = self.name.split("_")
        if len(parts) != 4:
            raise ValueError(
                "Name must be in the format '{name}_{amplicon-number}_{LEFT|RIGHT}_{primer-number}'"
            )
        return parts

    @computed_field
    @property
    def number(self) -> int:
        return int(self.name.split("_")[3])


class AmpliconModel(BaseModel):
    """An amplicon as represented by two or more primer records"""

    primers: List[PrimerModel]

    @computed_field
    @property
    def number(self) -> int:
        return int(self.primers[0].name_parts[1])

    @computed_field
    @property
    def min_start(self) -> int:
        return min(p.chrom_start for p in self.primers)

    @computed_field
    @property
    def max_end(self) -> int:
        return max(p.chrom_end for p in self.primers)


class BedModel(BaseModel):
    """A BED file as represented by a collection of amplicons each comprising primer records"""

    amplicons: Dict[str, List[AmpliconModel]]
    reference_lengths: Optional[Dict[str, int]] = None

    @model_validator(mode="after")
    def check_duplicate_primer_names(self) -> "BedModel":
        primer_names = []
        for chrom, amplicons in self.amplicons.items():
            for amplicon in amplicons:
                for primer in amplicon.primers:
                    primer_names.append(primer.name)
        if not len(primer_names) == len(set(primer_names)):
            raise ValueError("Duplicate primer names detected")
        return self

    @model_validator(mode="after")
    def check_primer_bounds(self) -> "BedModel":
        if self.reference_lengths:
            ref_lens = self.reference_lengths
            for chrom, amplicons in self.amplicons.items():
                for amplicon in amplicons:
                    if chrom not in ref_lens:
                        raise ValueError(
                            f"Reference chromosome {chrom} not found in reference.fasta"
                        )
                    if not amplicon.max_end <= ref_lens[chrom]:
                        raise ValueError(
                            f"A primer in amplicon {amplicon.number} is beyond the bounds of {chrom}"
                        )
        return self

    @model_validator(mode="after")
    def check_primer_tiling(self) -> "BedModel":
        chroms_amplicon_boundaries = {}
        for chrom, amplicons in self.amplicons.items():
            chrom = chrom.partition(" ")[0]
            chroms_amplicon_boundaries[chrom] = [
                (a.min_start, a.max_end) for a in amplicons
            ]
        for chrom, amplicon_boundaries in chroms_amplicon_boundaries.items():
            BedModel.check_tiling(amplicon_boundaries)
        return self

    @staticmethod
    def check_overlap(interval1: Tuple[int, int], interval2: Tuple[int, int]) -> bool:
        """Check if two intervals overlap"""
        return max(interval1[0], interval2[0]) <= min(interval1[1], interval2[1])

    @staticmethod
    def check_tiling(intervals: List[Tuple[int, int]]):
        """Verify if each interval overlaps exactly once with the interval before it and after it,
        and does not overlap with any other interval"""
        n = len(intervals)
        if n < 2:
            raise ValueError("Fewer than two amplicons detected")

        # Check first interval
        if not BedModel.check_overlap(intervals[0], intervals[1]):
            raise ValueError("First and second amplicons do not overlap")
        for j in range(2, n):
            if BedModel.check_overlap(intervals[0], intervals[j]):
                raise ValueError(
                    f"First amplicon overlaps with more than one amplicon ({j})"
                )

        # Check last interval
        if not BedModel.check_overlap(intervals[-1], intervals[-2]):
            raise ValueError("Penultimate and last amplicons do not overlap")
        for j in range(n - 2):
            if BedModel.check_overlap(intervals[-1], intervals[j]):
                raise ValueError(
                    f"Last amplicon overlaps with more than one amplicon ({j})"
                )

        # Check other intervals
        for i in range(1, n - 1):
            if not (
                BedModel.check_overlap(intervals[i], intervals[i - 1])
                and BedModel.check_overlap(intervals[i], intervals[i + 1])
            ):
                raise ValueError(f"Amplicons {i-1} and {i} do not overlap")

            # Check interval overlaps only with previous and next
            for j in range(n):
                if j != i - 1 and j != i + 1 and j != i:
                    if BedModel.check_overlap(intervals[i], intervals[j]):
                        raise ValueError(
                            f"Amplicon {i} overlaps with more than two other amplicons ({j})"
                        )
