import hashlib
from collections import defaultdict
from pathlib import Path

import defopt
import pandas as pd
from Bio import SeqIO


def hash_sequences(sequences: list[str]) -> str:
    sequences = sorted(sequences)
    concatenated_seqs = ",".join(sequences).upper()
    chars_in_concatenated_seqs = set("".join(sequences))
    allowed_chars = set("ACGT")
    if not chars_in_concatenated_seqs <= allowed_chars:
        raise RuntimeError("Found illegal characters in primer sequences")
    hex_digest = hashlib.sha256(concatenated_seqs.encode()).hexdigest()
    return hex_digest


def hash_primer_bed(bed_path: Path):
    """Hash a 7 column {scheme}.primer.bed file"""
    primer_bed_fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "poolName",
        "strand",
        "sequence",
    ]
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=primer_bed_fields,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            strand=str,
            sequence=str,
        ),
    )
    return hash_sequences(df["sequence"].tolist())


def hash_scheme_bed(bed_path: Path, fasta_path: Path) -> str:
    """
    Hash a 6 column {scheme}.scheme.bed file

    Seqs need to first be retrieved and sorted against in order to normalise prior to hashing
    """
    ref_record = SeqIO.read(fasta_path, "fasta")
    scheme_bed_fields = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "poolName",
        "strand",
    ]
    df = pd.read_csv(
        bed_path,
        sep="\t",
        names=scheme_bed_fields,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            strand=str,
        ),
    )
    primer_sequences = []
    for r in df.to_dict("records"):
        start_pos, end_pos = r["chromStart"], r["chromEnd"]
        if r["strand"] == "+":
            primer_sequences.append(str(ref_record.seq[start_pos:end_pos]))
        else:
            primer_sequences.append(
                str(ref_record.seq[start_pos:end_pos].reverse_complement())
            )
    return hash_sequences(primer_sequences)


def hash_ref(ref_path: Path):
    record = SeqIO.read(ref_path, "fasta")
    hex_digest = hashlib.sha256(str(record.seq).upper().encode()).hexdigest()
    return hex_digest


def build():
    """Build PHA4GE-compatible primer scheme bundle"""
    pass


def main():
    defopt.run(
        {
            "hash-bed": hash_bed,
            "hash-ref": hash_ref,
            "build": build,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()