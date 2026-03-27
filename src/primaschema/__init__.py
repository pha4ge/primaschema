"""Infrastructure for tiled amplicon PCR primer scheme definitions"""

import logging
from importlib.resources import files as _pkg_files
from pathlib import Path

METADATA_FILE_NAME: str = "info.json"
PRIMER_FILE_NAME: str = "primer.bed"
REFERENCE_FILE_NAME: str = "reference.fasta"

INDEX_FILE_NAME: str = "index.json"

DEFAULT_SCHEMES_URL = (
    "https://raw.githubusercontent.com/pha4ge/primer-schemes/v1b/schemes"
)
DEFAULT_INDEX_URL = (
    "https://raw.githubusercontent.com/pha4ge/primer-schemes/v1b/index.json"
)

SCHEME_FILES = [METADATA_FILE_NAME, PRIMER_FILE_NAME, REFERENCE_FILE_NAME]
SCHEME_FILES_EXTRA = ["README.md", "work/primer.svg"]

# Locate schema files via importlib.resources — works for both editable and installed packages.
SCHEMA_DIR = Path(str(_pkg_files("primaschema").joinpath("schema")))
INDEX_SCHEMA_PATH = SCHEMA_DIR / INDEX_FILE_NAME
INDEX_HEADER_PATH = SCHEMA_DIR / "index-header.yml"

logging.getLogger("primaschema").addHandler(logging.NullHandler())
