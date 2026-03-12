from hashlib import sha256
from pathlib import Path

import dnaio

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
    return sha256_hasher.hexdigest()


def read_fasta_records(path: Path) -> list[dnaio.SequenceRecord]:
    with dnaio.open(path) as reader:
        return list(reader)


def write_fasta_records(
    path: Path,
    records: list[dnaio.SequenceRecord],
    line_length: int = 60,
) -> None:
    with dnaio.FastaWriter(path, line_length=line_length) as writer:
        for record in records:
            writer.write(record)


def reverse_complement(sequence: str) -> str:
    return dnaio.SequenceRecord("sequence", sequence).reverse_complement().sequence


def find_all_info_json(primer_schemes_path: Path):
    return list(primer_schemes_path.rglob(f"*/{METADATA_FILE_NAME}"))
