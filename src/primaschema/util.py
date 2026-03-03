from hashlib import sha256
from pathlib import Path

from Bio import SeqIO
from primalbedtools.bedfiles import BedLine, BedLineParser, sort_bedlines

from primaschema import METADATA_FILE_NAME


def sha256_checksum(filename: Path):
    """
    Docstring for sha256_checksum

    :param filename: Description
    """
    sha256_hasher = sha256()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256_hasher.update(block)
    return f"sha256:{sha256_hasher.hexdigest()}"


def primaschema_bed_hash(
    bedfile_path: Path | None = None, bedlines: list[BedLine] | None = None
) -> str:
    if bedfile_path is None and bedlines is None:
        raise ValueError("Please provide either an path or BedLines to hash")

    if bedlines is None and bedfile_path is not None:
        _h, bedlines = BedLineParser.from_file(bedfile_path)

    assert bedlines is not None
    s_bedlines = sort_bedlines(bedlines)
    hasher = sha256()
    # extract wanted fields
    for sbedline in s_bedlines:
        hasher.update(
            "\t".join(
                [
                    sbedline.chrom,
                    str(sbedline.start),
                    str(sbedline.end),
                    str(sbedline.pool),
                    sbedline.strand,
                    sbedline.sequence,
                ]
            ).encode()
        )
    # return truncated hash
    return f"primaschema:bed:{hasher.hexdigest()[:16]}"


def primaschema_ref_hash(
    ref_path: Path | None = None,
    seq_records: list | None = None,
) -> str:
    if ref_path is None and seq_records is None:
        raise ValueError("Please provide either an path or SeqRecords to hash")

    if seq_records is None:
        seq_records = list(SeqIO.parse(ref_path, "fasta"))

    # Sort by id
    seq_records.sort(key=lambda x: x.id)

    hasher = sha256()
    for record in seq_records:
        hasher.update(record.id.strip().upper().encode())
        hasher.update(str(record.seq.upper()).encode())

    # return truncated hash
    return f"primaschema:ref:{hasher.hexdigest()[:16]}"


def find_all_info_json(primer_schemes_path: Path):
    return list(primer_schemes_path.rglob(f"*/{METADATA_FILE_NAME}"))
