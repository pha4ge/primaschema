import hashlib
import re
import shutil
import sys

from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal, Optional, Dict, Tuple

import linkml.validator
import yaml

import altair as alt
import pandas as pd

from Bio import SeqIO

from linkml.generators.pydanticgen import PydanticGenerator

from . import (
    logger,
    MANIFEST_HEADER_PATH,
    CACHE_DIR,
    SCHEMA_DIR,
    SCHEMES_ARCHIVE_URL,
)

from .schema import bed, info
from . import util
import json


SCHEME_BED_FIELDS = ["chrom", "chromStart", "chromEnd", "name", "poolName", "strand"]
PRIMER_BED_FIELDS = SCHEME_BED_FIELDS + ["sequence"]
HASHED_BED_FIELDS = [
    "chrom",
    "chromStart",
    "chromEnd",
    "poolName",
    "strand",
    "sequence",
]
POSITION_FIELDS = ["chromStart", "chromEnd"]
MANDATORY_FILES = ("primer.bed", "reference.fasta", "info.yml")


def hash_string(string: str) -> str:
    """Normalise case, sorting, terminal spaces & return prefixed 64b of SHA256 hex"""
    checksum = hashlib.sha256(str(string).strip().upper().encode()).hexdigest()[:16]
    return f"primaschema:{checksum}"


def parse_scheme_bed(bed_path: Path) -> pd.DataFrame:
    """Parse a 6 column scheme.bed bed file"""
    return pd.read_csv(
        bed_path,
        sep="\t",
        names=SCHEME_BED_FIELDS,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=int,
            strand=str,
        ),
    )


def parse_primer_bed(bed_path: Path) -> pd.DataFrame:
    """Parse a 7 column primer.bed bed file"""
    return pd.read_csv(
        bed_path,
        sep="\t",
        comment="#",
        names=PRIMER_BED_FIELDS,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=int,
            strand=str,
            sequence=str,
        ),
    )


def sort_primer_df(df: pd.DataFrame) -> pd.DataFrame:
    df["amplicon_number"] = df["name"].apply(lambda x: int(x.split("_")[1]))
    return df.sort_values(
        [
            "chrom",
            "amplicon_number",
            "chromStart",
            "chromEnd",
            "poolName",
            "strand",
            "sequence",
        ]
    )[[*PRIMER_BED_FIELDS]]


def normalise_primer_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Removes terminal whitespace and normalises case
    Sorts by chromStart, chromEnd, poolName, strand, sequence
    Removes duplicate records, collapsing alts with same coords if backfilled from ref
    """
    df["sequence"] = df["sequence"].str.strip().str.upper()
    return sort_primer_df(df)


def hash_primer_df(df: pd.DataFrame) -> str:
    """
    Returns prefixed SHA256 digest from stringified dataframe
    """
    normalised_df = normalise_primer_df(df)
    string = normalised_df[[*HASHED_BED_FIELDS]].to_csv(index=False, header=False)
    logger.debug(f"hash_primer_df() {string=}")
    return hash_string(string)


def hash_primer_bed(bed_path: Path):
    """Hash a 7 column primer.bed file"""
    df = parse_primer_bed(bed_path)
    return hash_primer_df(df)


def hash_scheme_bed(bed_path: Path, fasta_path: Path) -> str:
    """
    Hash a 6 column scheme.bed file by first converting to 7 column primer.bed
    """
    logger.info("Hashing scheme.bed using reference backfill")
    ref_record = SeqIO.read(fasta_path, "fasta")
    df = parse_scheme_bed(bed_path)
    records = df.to_dict("records")
    for r in records:
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            r["sequence"] = str(ref_record.seq[start_pos:end_pos])
        elif r["strand"] == "-":
            r["sequence"] = str(ref_record.seq[start_pos:end_pos].reverse_complement())
        else:
            raise RuntimeError(f"Invalid strand for BED record {r}")
    bed7_df = pd.DataFrame(records)
    return hash_primer_df(bed7_df)


def convert_primer_bed_to_scheme_bed(bed_path: Path) -> str:
    df = parse_primer_bed(bed_path).drop("sequence", axis=1)
    return df.to_csv(sep="\t", header=False, index=False)


def convert_scheme_bed_to_primer_bed(bed_path: Path, fasta_path: Path) -> str:
    ids_seqs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    df = parse_scheme_bed(bed_path)
    records = df.to_dict("records")
    for r in records:
        chrom = r["chrom"].partition(" ")[0]  # Use chrom name before first space
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            r["sequence"] = str(ids_seqs[chrom].seq[start_pos:end_pos])
        else:
            r["sequence"] = str(
                ids_seqs[chrom].seq[start_pos:end_pos].reverse_complement()
            )
    df = pd.DataFrame(records)
    return df.to_csv(sep="\t", header=False, index=False)


def convert_vwf_to_primer_bed(vwf_path: Path, chrom: str = "chrom") -> str:
    vwf_df = pd.read_csv(vwf_path, sep="\t")
    bed_records = []
    pool_counter = {}

    for r in vwf_df.to_records("dict"):
        amplicon_name = r["Amplicon_name"]
        primer_name = r["Primer_name"]
        orientation = r["Left_or_right"]
        amplicon_number = int(amplicon_name.split("_")[-1])
        pool_name = 1 if amplicon_number % 2 != 0 else 2
        if amplicon_name not in pool_counter:
            pool_counter[amplicon_name] = 1
        else:
            pool_counter[amplicon_name] += 1
        strand = "+" if orientation == "left" else "-"
        sequence = r["Sequence"]
        start_pos = r["Position"]
        bed_record = {}
        bed_record["chrom"] = chrom
        bed_record["chromStart"] = start_pos
        bed_record["chromEnd"] = start_pos + len(sequence)
        bed_record["name"] = primer_name
        bed_record["poolName"] = str(pool_name)
        bed_record["strand"] = strand
        bed_record["sequence"] = sequence
        bed_records.append(bed_record)

    bed_df = pd.DataFrame(bed_records)
    return bed_df.to_csv(sep="\t", header=False, index=False)


def hash_bed(bed_path: Path) -> str:
    bed_type = infer_bed_type(bed_path)
    if bed_type == "primer":
        checksum = hash_primer_bed(bed_path)
    else:  # bed_type == "scheme"
        checksum = hash_scheme_bed(
            bed_path=bed_path, fasta_path=bed_path.parent / "reference.fasta"
        )
    return checksum


def hash_ref(ref_path: Path):
    chroms_seqs = {}
    for record in SeqIO.parse(ref_path, "fasta"):
        chroms_seqs[record.id] = str(record.seq).upper()
    chroms_seqs_sorted = {key: chroms_seqs[key] for key in sorted(chroms_seqs)}
    string = ""
    for chrom, seq in chroms_seqs_sorted.items():
        string += f">{chrom}\n{seq}\n"
    logger.debug(f"hash_ref() {string=}")
    return hash_string(string.strip())


def count_tsv_columns(bed_path: Path) -> int:
    return len(pd.read_csv(bed_path, sep="\t", comment="#").columns)


def parse_yaml(path: Path) -> dict:
    with open(path, "r") as fh:
        return yaml.safe_load(fh)


def validate_scheme_yaml_with_linkml(path: Path) -> None:
    data = parse_yaml(path)
    report = linkml.validator.validate(data, SCHEMA_DIR / "info.yml", "PrimerScheme")
    if report.results:
        msg = ""
        for result in report.results:
            msg += f"{result.message}\n"
        raise ValueError(msg)


def parse_scheme_yaml(path: Path) -> dict:
    """Parse and validate with Pydantic"""
    data = parse_yaml(path)
    return info.PrimerScheme(**data).model_dump()


def validate_bed_and_ref(
    bed_path: Path, ref_path: Optional[Path] = None
) -> bed.BedModel:
    """Check that primer.bed is a tiled PrimalScheme v3 BED"""
    with open(bed_path) as fh:
        bed_contents = fh.readlines()
    primers = []
    amplicons_primers = defaultdict(list)
    for line in bed_contents:
        if line.startswith("#"):
            continue
        r = line.split("\t")
        primer = bed.PrimerModel(
            chrom=r[0].partition(" ")[0],  # We want to ignore anything after a space
            chrom_start=int(r[1]),
            chrom_end=int(r[2]),
            name=r[3],
            pool_name=int(r[4]),
            strand=r[5],
            sequence=r[6],
        )
        chrom = primer.chrom
        amplicon_number = int(primer.name_parts[1])
        primers.append(primer)
        amplicons_primers[(chrom, amplicon_number)].append(primer)

    amplicons = defaultdict(list)
    for (chrom, amplicon_number), primers in amplicons_primers.items():
        amplicons[chrom].append(bed.AmpliconModel(primers=primers))

    # If a ref_path is supplied, populate the reference_lengths field of BedModel
    if ref_path:
        records = SeqIO.parse(ref_path, "fasta")
        reference_lengths = {r.id.partition(" ")[0]: len(r.seq) for r in records}
    else:
        reference_lengths = None
    return bed.BedModel(amplicons=amplicons, reference_lengths=reference_lengths)


def infer_bed_type(bed_path: Path) -> str:
    bed_columns = count_tsv_columns(bed_path)
    if bed_columns == 7:
        bed_type = "primer"
    elif bed_columns == 6:
        bed_type = "scheme"
    else:
        raise RuntimeError(
            "Bed file shoud have either 6 columns (scheme.bed) or 7 column (primer.bed)"
        )
    return bed_type


def validate(
    scheme_dir: Path,
    full: bool = False,
    ignore_checksums: bool = False,
    recursive: bool = False,
    rebuild: bool = False,
) -> None:
    info_schema_linkml_path = SCHEMA_DIR / "info.yml"
    info_schema_pydantic_path = SCHEMA_DIR / "info.py"
    if recursive:
        for path in Path(scheme_dir).rglob("info.yml"):
            if path.is_file() and path.name == "info.yml":
                validate(
                    scheme_dir=path.parent, full=full, ignore_checksums=ignore_checksums
                )
    else:
        logger.info(f"Validating {scheme_dir}")
        yml_path = Path(scheme_dir / "info.yml")
        bed_path = Path(scheme_dir / "primer.bed")
        ref_path = Path(scheme_dir / "reference.fasta")

        if (
            rebuild
            or not info_schema_pydantic_path.exists()
            or info_schema_pydantic_path.stat().st_size == 0
        ):
            logger.info("Building intermediate pydantic model")
            with open(info_schema_pydantic_path, "w") as fh:
                fh.write(PydanticGenerator(info_schema_linkml_path).serialize())

        if full:
            logger.debug("Validating with linkml model")
            validate_scheme_yaml_with_linkml(path=yml_path)
        else:
            logger.debug("Validating with pydantic model")
            parse_scheme_yaml(yml_path)

        logger.debug("Validating primer.bed and reference.fasta")
        validate_bed_and_ref(bed_path=bed_path, ref_path=ref_path)

        scheme = parse_scheme_yaml(yml_path)
        existing_primer_checksum = scheme.get("primer_checksum")
        existing_reference_checksum = scheme.get("reference_checksum")
        primer_checksum = hash_bed(bed_path)
        reference_checksum = hash_ref(ref_path)
        if (
            existing_primer_checksum
            and not primer_checksum == existing_primer_checksum
            and not ignore_checksums
        ):
            raise RuntimeError(
                f"Calculated and documented primer checksums do not match ({primer_checksum} and {existing_primer_checksum})"
            )
        elif not primer_checksum == existing_primer_checksum:
            logger.warning(
                f"Calculated and documented primer checksums do not match ({primer_checksum} and {existing_primer_checksum})"
            )
        if (
            existing_reference_checksum
            and not reference_checksum == existing_reference_checksum
            and not ignore_checksums
        ):
            raise RuntimeError(
                f"Calculated and documented reference checksums do not match ({reference_checksum} and {existing_reference_checksum})"
            )
        elif not reference_checksum == existing_reference_checksum:
            logger.warning(
                f"Calculated and documented reference checksums do not match ({reference_checksum} and {existing_reference_checksum})"
            )
        logger.info(f"Validated {get_scheme_cname(scheme)}")


def format_primer_bed(bed_path: Path) -> str:
    """Sort a primer bed into a maximally compatible format"""
    df = parse_primer_bed(bed_path)
    return sort_primer_df(df).to_csv(sep="\t", header=False, index=False)


def build(
    scheme_dir: Path,
    out_dir: Path = Path("built"),
    plot: bool = False,
    nested: bool = False,  # Create nested output dir structure
    recursive: bool = False,
    ignore_checksums: bool = False,
) -> None:
    """
    Validate and build a primer scheme given a scheme directory path.
    Optionally do so recursively
    """
    if recursive:
        for path in Path(scheme_dir).rglob("info.yml"):
            if path.is_file() and path.name == "info.yml":
                build(
                    scheme_dir=path.parent,
                    out_dir=out_dir,
                    plot=plot,
                    nested=True,
                    ignore_checksums=ignore_checksums,
                )
    else:
        validate(scheme_dir=scheme_dir, ignore_checksums=ignore_checksums)
        scheme = parse_yaml(scheme_dir / "info.yml")
        scheme_cname = get_scheme_cname(scheme)
        if nested:
            out_dir = Path(out_dir) / Path(scheme_cname)
        else:
            out_dir = Path(out_dir) / get_scheme_cname(scheme, sep="_")
        try:
            out_dir.mkdir(parents=True, exist_ok=True)
        except FileExistsError:
            raise FileExistsError(f"Output directory {out_dir} already exists")
        scheme["primer_checksum"] = hash_bed(scheme_dir / "primer.bed")
        scheme["reference_checksum"] = hash_ref(scheme_dir / "reference.fasta")
        with open(out_dir / "info.yml", "w") as scheme_fh:
            logger.debug(f"Writing info.yml to {out_dir}/info.yml")
            yaml.dump(scheme, scheme_fh, sort_keys=False)
        logger.debug(f"Copying primer.bed to {out_dir}/primer.bed")
        with open(out_dir / "primer.bed", "w") as primer_fh:
            primer_fh.write(format_primer_bed(scheme_dir / "primer.bed"))
        logger.debug(f"Copying reference.fasta to {out_dir}/reference.fasta")
        shutil.copy(scheme_dir / "reference.fasta", out_dir)
        logger.debug(f"Writing scheme.bed to {out_dir}/scheme.bed")
        scheme_bed_str = convert_primer_bed_to_scheme_bed(
            bed_path=out_dir / "primer.bed"
        )
        with open(out_dir.resolve() / "scheme.bed", "w") as fh:
            fh.write(scheme_bed_str)

        if plot:
            plot_primers(
                bed_path=out_dir / "primer.bed", out_path=out_dir / "primer.svg"
            )

        logger.info(f"Built {scheme_cname}")


def get_scheme_cname(scheme: dict, sep: Literal["/", "_"] = "/") -> str:
    organism = str(scheme.get("organism", ""))
    name = str(scheme["name"])
    amplicon_size = str(scheme.get("amplicon_size", ""))
    version = str(scheme["version"])
    return sep.join([organism, name, amplicon_size, version])


def build_manifest(root_dir: Path, out_dir: Path = Path()):
    """Build manifest of schemes inside the specified directory"""

    manifest = parse_yaml(MANIFEST_HEADER_PATH)

    manifest_field_exclude = [
        "schema_version",
    ]

    scheme_path = root_dir / "schemes"
    if not scheme_path.exists():
        scheme_path = root_dir

    schemes = []
    organism_set = set([o["organism"] for o in manifest["organisms"]])
    for scheme_info_path in scheme_path.glob("**/info.yml"):
        scheme = parse_yaml(scheme_info_path)
        for field in manifest_field_exclude:
            if field in scheme:
                del scheme[field]
        if scheme["organism"] not in organism_set:
            logger.warning(
                f"Skipping scheme {scheme['name']} with unknown organism {scheme['organism']}",
            )
        schemes.append(scheme)

    manifest["schemes"] = sorted(schemes, key=get_scheme_cname)

    manifest_file_name = "index.json"
    with open(out_dir / manifest_file_name, "w") as fh:
        logger.info(f"Writing {manifest_file_name} to {out_dir}/{manifest_file_name}")
        json.dump(manifest, fh, indent=4)


def diff(bed1_path: Path, bed2_path: Path, only_positions: bool = False):
    """Show symmetric differences between records in two primer.bed files"""
    df1 = parse_primer_bed(bed1_path).assign(origin="bed1")
    df2 = parse_primer_bed(bed2_path).assign(origin="bed2")
    if only_positions:
        column_subset = POSITION_FIELDS
    else:
        column_subset = PRIMER_BED_FIELDS
    return pd.concat([df1, df2]).drop_duplicates(subset=column_subset, keep=False)


def discordant_primers(scheme_dir: Path) -> pd.DataFrame:
    """Show primer records with sequences not matching the reference sequence"""
    primer_bed_path = scheme_dir / "primer.bed"
    scheme_bed_path = scheme_dir / "scheme.bed"
    with TemporaryDirectory() as temp_dir:
        backfilled_bed_path = Path(temp_dir) / "primer.bed"
        if not scheme_bed_path.exists():
            with open(scheme_bed_path, "w") as fh:
                fh.write(convert_primer_bed_to_scheme_bed(primer_bed_path))
            scheme_bed_path = Path(temp_dir) / "scheme.bed"
        with open(backfilled_bed_path, "w") as fh:
            fh.write(
                convert_scheme_bed_to_primer_bed(
                    bed_path=scheme_dir / "scheme.bed",
                    fasta_path=scheme_dir / "reference.fasta",
                )
            )
        return diff(bed1_path=primer_bed_path, bed2_path=backfilled_bed_path)


def amplicon_intervals(bed_path: Path) -> Dict[str, Dict[str, Tuple[int, int]]]:
    """
    find primer positions for all primers in the bed file and compute maximum
    interval between primers of the same name
    """
    primer_name_re = re.compile(r"^(?P<name>.*)_(LEFT|RIGHT)(_.+)?$")
    all_intervals: dict[str, dict[str, Tuple[int, int]]] = {}
    for line in open(bed_path):
        line_parts = line.strip().split("\t")
        if len(line_parts) < 6:
            # skip lines that don't have at least 6 fields
            continue
        chrom, start, end, name, _, strand = line.strip().split("\t")[:6]
        if chrom not in all_intervals:
            all_intervals[chrom] = {}
        intervals = all_intervals[chrom]
        primer_match = primer_name_re.match(name)
        if not primer_match:
            raise ValueError(f"Invalid primer name {name}")
        primer_name = primer_match.group("name")
        if strand == "+":
            start_pos = int(start)
            end_pos = -1
        if strand == "-":
            start_pos = sys.maxsize
            end_pos = int(end)
        prev_start, prev_end = intervals.get(primer_name, (sys.maxsize, -1))
        intervals[primer_name] = (min(prev_start, start_pos), max(prev_end, end_pos))
    return all_intervals


def plot_primers(bed_path: Path, out_path: Path = Path("primer.html")) -> None:
    """
    Plot amplicon and primer positions from a 7 column primer.bed file
    Requires primers to be named {$scheme_id}_{$amplicon_id}_{LEFT|RIGHT}_{1|2|3…}
    Plots one vertical panel per pool per reference chromosome
    Supported out_path extensions: html (interactive), pdf, png, svg
    """
    primer_df = parse_primer_bed(bed_path)
    primer_df["amplicon"] = primer_df["name"].str.split("_").str[1]
    primer_df["poolName"] = primer_df["poolName"].astype(str)

    amp_df = (
        primer_df.groupby(["chrom", "amplicon", "poolName"])
        .agg(min_start=("chromStart", "min"), max_end=("chromEnd", "max"))
        .reset_index()
    )
    amp_df["is_amplicon"] = True
    primer_df["is_amplicon"] = False
    amp_df = amp_df.rename(columns={"min_start": "chromStart", "max_end": "chromEnd"})
    combined_df = pd.concat([primer_df, amp_df], ignore_index=True)

    primer_marks = (
        alt.Chart(combined_df)
        .transform_filter(alt.datum.is_amplicon == False)  # noqa
        .mark_line(size=15)
        .encode(
            x=alt.X("chromStart:Q", title=None),
            x2="chromEnd:Q",
            y=alt.Y("poolName:O", title="pool", scale=alt.Scale(padding=1)),
            color=alt.Color("strand:N").scale(scheme="set2"),
            tooltip=[
                alt.Tooltip("name:N", title="Primer name"),
                alt.Tooltip("chromStart:Q", title="start"),
                alt.Tooltip("chromEnd:Q", title="end"),
            ],
        )
        .properties(
            width=800,
        )
    )
    amplicon_marks = (
        alt.Chart(combined_df)
        .transform_filter(alt.datum.is_amplicon == True)  # noqa
        .mark_rule(strokeWidth=2)
        .encode(
            x=alt.X("chromStart:Q", title=None),
            x2="chromEnd:Q",
            y=alt.Y("poolName:O", title="pool", scale=alt.Scale(padding=1)),
            tooltip=[
                alt.Tooltip("amplicon:N", title="Amplicon name"),
                alt.Tooltip("chromStart:Q", title="Min primer start"),
                alt.Tooltip("chromEnd:Q", title="Max primer end"),
            ],
        )
        .properties(
            width=800,
        )
    )
    combined_chart = (
        alt.layer(primer_marks, amplicon_marks)
        .facet(row=alt.Row("chrom:O", header=alt.Header(labelOrient="top"), title=""))
        .configure_axis(domain=False, ticks=False)
    )

    combined_chart.interactive().save(str(out_path))
    logger.info(f"Plot saved ({out_path})")


def subset(scheme_dir: Path, chrom: str, out_dir: Path = Path("built")) -> None:
    scheme_dir = Path(scheme_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_chroms = set()
    subset_record = None
    for record in SeqIO.parse(scheme_dir / "reference.fasta", "fasta"):
        reference_chroms.add(record.id)
        if record.id == chrom:
            subset_record = record
    if not subset_record:
        raise ValueError(f"Chrom {chrom} not found in reference.fasta")
    else:
        with open(out_dir / "reference.fasta", "w") as fh:
            SeqIO.write(subset_record, fh, "fasta")

    primers_df = parse_primer_bed(scheme_dir / "primer.bed")
    primers_df["chrom"] = primers_df["chrom"].str.partition(" ")[0]
    primer_chroms = set(primers_df["chrom"].unique())
    subset_primers_df = primers_df.query(f"chrom == '{chrom}'")
    subset_primers_df.to_csv(
        out_dir / "primer.bed", sep="\t", index=False, header=False
    )

    if reference_chroms != primer_chroms:
        logger.info(f"Reference chroms: {reference_chroms}")
        logger.info(f"Primer chroms: {primer_chroms}")
    else:
        logger.info(f"Chroms: {reference_chroms}")

    logger.info(
        f"Wrote subset of {len(subset_primers_df)}/{len(primers_df)} primers for {chrom} to {out_dir.resolve()}"
    )


def synchronise() -> None:
    util.download_github_tarball(SCHEMES_ARCHIVE_URL, CACHE_DIR)


# def fetch(string: str) -> None:
#     with open(CACHE_DIR / "index.json", "r") as fh:
#         manifest = json.load(fh)
