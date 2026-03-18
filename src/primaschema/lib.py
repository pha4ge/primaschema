import hashlib
import json
import logging
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Literal, Optional, Tuple

import httpx

import altair as alt
import linkml.validator
import pandas as pd
import yaml
from linkml.generators.pydanticgen import PydanticGenerator

from primaschema import (
    DEFAULT_SCHEMES_URL,
    INDEX_HEADER_PATH,
    SCHEMA_DIR,
    SCHEME_FILES,
    SCHEME_FILES_EXTRA,
)
from primaschema.schema import bed, info
from primaschema.util import read_fasta_records, reverse_complement, write_fasta_records

logger = logging.getLogger(__name__)

SCHEME_BED_FIELDS = ["chrom", "chromStart", "chromEnd", "name", "poolName", "strand"]
PRIMER_BED_FIELDS = SCHEME_BED_FIELDS + ["sequence"]
POSITION_FIELDS = ["chromStart", "chromEnd"]
MANDATORY_FILES = ("primer.bed", "reference.fasta", "info.yml")


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


def convert_primer_bed_to_scheme_bed(bed_path: Path) -> str:
    df = parse_primer_bed(bed_path).drop("sequence", axis=1)
    return df.to_csv(sep="\t", header=False, index=False)


def convert_scheme_bed_to_primer_bed(bed_path: Path, fasta_path: Path) -> str:
    ids_seqs = {record.id: record.sequence for record in read_fasta_records(fasta_path)}
    df = parse_scheme_bed(bed_path)
    records = df.to_dict("records")
    for r in records:
        chrom = r["chrom"].partition(" ")[0]  # Use chrom name before first space
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            r["sequence"] = str(ids_seqs[chrom][start_pos:end_pos])
        else:
            r["sequence"] = reverse_complement(ids_seqs[chrom][start_pos:end_pos])
    df = pd.DataFrame(records)
    return df.to_csv(sep="\t", header=False, index=False)


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
        records = read_fasta_records(ref_path)
        reference_lengths = {r.id: len(r.sequence) for r in records}
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
        logger.info(f"Validated {get_scheme_cname(scheme)}")


def format_primer_bed(bed_path: Path) -> str:
    """Sort a primer bed into a maximally compatible format"""
    df = parse_primer_bed(bed_path)
    return sort_primer_df(df).to_csv(sep="\t", header=False, index=False)


def get_scheme_cname(scheme: dict, sep: Literal["/", "_"] = "/") -> str:
    target_organisms = scheme.get("target_organisms", [])
    if target_organisms:
        first = target_organisms[0]
        organism = (
            first.get("common_name", "") if isinstance(first, dict) else str(first)
        ) or ""
    else:
        organism = str(scheme.get("organism", ""))
    name = str(scheme["name"])
    amplicon_size = str(scheme.get("amplicon_size", ""))
    version = str(scheme["version"])
    parts = [p for p in [organism, name, amplicon_size, version] if p]
    return sep.join(parts)


def build_index(root_dir: Path, out_dir: Path = Path()):
    """Build index of schemes inside the specified directory"""

    index_data = parse_yaml(INDEX_HEADER_PATH)

    field_exclude = [
        "schema_version",
    ]

    scheme_path = root_dir / "schemes"
    if not scheme_path.exists():
        scheme_path = root_dir

    schemes = []
    organism_set = set([o["organism"] for o in index_data["organisms"]])
    for scheme_info_path in scheme_path.glob("**/info.yml"):
        scheme = parse_yaml(scheme_info_path)
        for field in field_exclude:
            if field in scheme:
                del scheme[field]
        if scheme["organism"] not in organism_set:
            logger.warning(
                f"Skipping scheme {scheme['name']} with unknown organism {scheme['organism']}",
            )
        schemes.append(scheme)

    index_data["schemes"] = sorted(schemes, key=get_scheme_cname)

    index_file_name = "index.json"
    with open(out_dir / index_file_name, "w") as fh:
        logger.info(f"Writing {index_file_name} to {out_dir}/{index_file_name}")
        json.dump(index_data, fh, indent=4)


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
    logger.debug(f"Plot saved ({out_path})")


def subset(scheme_dir: Path, chrom: str, out_dir: Path = Path("built")) -> None:
    scheme_dir = Path(scheme_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_chroms = set()
    subset_record = None
    for record in read_fasta_records(scheme_dir / "reference.fasta"):
        reference_chroms.add(record.id)
        if record.id == chrom:
            subset_record = record
    if not subset_record:
        raise ValueError(f"Chrom {chrom} not found in reference.fasta")
    else:
        write_fasta_records(out_dir / "reference.fasta", [subset_record])

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


def _github_tree_url_to_raw(url: str) -> str:
    """Convert a GitHub tree URL to a raw.githubusercontent.com URL."""
    return url.replace("github.com", "raw.githubusercontent.com").replace(
        "/tree/", "/refs/heads/"
    )


def get_scheme(
    scheme_id: str,
    output: Path = Path("."),
    base_url: str = DEFAULT_SCHEMES_URL,
    all_files: bool = False,
) -> Path:
    """Download a primer scheme by ID (name/amplicon_size/version) from a remote repository."""
    parts = scheme_id.strip("/").split("/")
    if len(parts) != 3:
        raise ValueError(
            f"Invalid scheme_id '{scheme_id}': expected format name/amplicon_size/version"
        )
    name, amplicon_size, version = parts
    raw_base = _github_tree_url_to_raw(base_url.rstrip("/"))
    scheme_url = f"{raw_base}/{scheme_id}"
    output_dir = output / name / amplicon_size / version

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir) / name / amplicon_size / version
        tmp_path.mkdir(parents=True)
        with httpx.Client() as client:
            for filename in SCHEME_FILES:
                url = f"{scheme_url}/{filename}"
                response = client.get(url, follow_redirects=True)
                if response.status_code != 200:
                    raise RuntimeError(
                        f"Failed to download {url}: HTTP {response.status_code}"
                    )
                (tmp_path / filename).write_bytes(response.content)
            if all_files:
                for filename in SCHEME_FILES_EXTRA:
                    url = f"{scheme_url}/{filename}"
                    response = client.get(url, follow_redirects=True)
                    if response.status_code != 200:
                        logger.warning(
                            f"Could not download optional file {filename}: HTTP {response.status_code}"
                        )
                        continue
                    file_path = tmp_path / filename
                    file_path.parent.mkdir(parents=True, exist_ok=True)
                    file_path.write_bytes(response.content)
        # Verify checksums from info.json
        info = json.loads((tmp_path / "info.json").read_text())
        checksums = info.get("checksums", {})
        for field, filename in [
            ("primer_sha256", "primer.bed"),
            ("reference_sha256", "reference.fasta"),
        ]:
            expected = checksums.get(field)
            if expected:
                actual = hashlib.sha256((tmp_path / filename).read_bytes()).hexdigest()
                if actual != expected:
                    raise RuntimeError(
                        f"Checksum mismatch for {filename}: "
                        f"expected {expected}, got {actual}"
                    )
        output_dir.mkdir(parents=True, exist_ok=True)
        shutil.copytree(tmp_path, output_dir, dirs_exist_ok=True)
    logger.info(f"Downloaded scheme {scheme_id} to {output_dir}")
    return output_dir
