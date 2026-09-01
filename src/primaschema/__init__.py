"""Infrastructure for tiled amplicon PCR primer scheme definitions"""

import logging
from importlib.resources import files as _pkg_files
from pathlib import Path

METADATA_FILE_NAME: str = "info.json"
PRIMER_FILE_NAME: str = "primer.bed"
REFERENCE_FILE_NAME: str = "reference.fasta"
README_FILE_NAME: str = "README.md"

INDEX_FILE_NAME: str = "index.json"

# raw.githubusercontent.com paths are /{owner}/{repo}/{ref}/{path} — the ref is not
# optional, so these must name a branch explicitly rather than relying on a default.
DEFAULT_SCHEMES_URL = (
    "https://raw.githubusercontent.com/pha4ge/primer-schemes/main/schemes"
)
DEFAULT_INDEX_URL = (
    "https://raw.githubusercontent.com/pha4ge/primer-schemes/main/index.json"
)

# Locate schema files via importlib.resources — works for both editable and installed packages.
SCHEMA_DIR = Path(str(_pkg_files("primaschema").joinpath("schema")))
INDEX_SCHEMA_PATH = SCHEMA_DIR / INDEX_FILE_NAME
INDEX_HEADER_PATH = SCHEMA_DIR / "index-header.yml"

logging.getLogger("primaschema").addHandler(logging.NullHandler())
