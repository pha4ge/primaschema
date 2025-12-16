import shutil
import subprocess
import tarfile
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import httpx
from Bio import SeqIO
from primalbedtools.bedfiles import BedLine, BedLineParser, sort_bedlines

from primaschema import METADATA_FILE_NAME, logger


def run(cmd, cwd="./"):  # Helper for CLI testing
    return subprocess.run(
        cmd, cwd=cwd, shell=True, check=True, text=True, capture_output=True
    )


def copy_single_child_dir_to_parent(parent_dir_path: Path) -> None:
    parent_dir = Path(parent_dir_path)
    child_dirs = [
        d for d in parent_dir.iterdir() if d.is_dir() and not d.name.startswith(".")
    ]
    if len(child_dirs) != 1:
        raise FileNotFoundError(
            f"Expected one child directory not starting with a dot, but found {len(child_dirs)}."
        )
    child_dir = child_dirs[0]

    for item in child_dir.iterdir():
        destination = parent_dir / item.name
        if item.is_dir():
            shutil.copytree(item, destination, dirs_exist_ok=True)
        else:
            shutil.copy2(item, destination)

    shutil.rmtree(child_dir)


def download_github_tarball(archive_url: str, out_dir: Path) -> None:
    if not archive_url.endswith(".tar.gz"):
        raise ValueError("Archive URL must end with .tar.gz")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    response = httpx.get(archive_url, follow_redirects=True)
    response.raise_for_status()

    shutil.rmtree(out_dir, ignore_errors=True)
    with tarfile.open(fileobj=BytesIO(response.content), mode="r:gz") as tf_fh:
        tf_fh.extractall(out_dir)

    copy_single_child_dir_to_parent(out_dir)

    logger.info(f"Schemes downloaded and extracted to {out_dir}")


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


def find_all_info_json(primerschemes_repo: Path):
    return list(primerschemes_repo.rglob(f"*/{METADATA_FILE_NAME}"))
