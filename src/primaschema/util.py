import subprocess
import shutil
import tarfile

from io import BytesIO
from pathlib import Path

import httpx

from . import logger


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
