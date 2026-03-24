"""Infrastructure for tiled amplicon PCR primer scheme definitions"""

import logging
import os
from pathlib import Path

METADATA_FILE_NAME: str = "info.json"
PRIMER_FILE_NAME: str = "primer.bed"
REFERENCE_FILE_NAME: str = "reference.fasta"

INDEX_FILE_NAME: str = "index.json"

DEFAULT_SCHEMES_URL = "https://github.com/pha4ge/primer-schemes/tree/v1b/schemes"
DEFAULT_INDEX_URL = "https://github.com/pha4ge/primer-schemes/tree/v1b/index.json"

SCHEME_FILES = [METADATA_FILE_NAME, PRIMER_FILE_NAME, REFERENCE_FILE_NAME]
SCHEME_FILES_EXTRA = ["README.md", "work/primer.svg"]

PKG_DIR = Path(
    os.environ.get(
        "PRIMASCHEMA_ROOT_PATH", Path(__file__).absolute().parent.parent.parent
    )
)
SCHEMA_DIR = PKG_DIR / "src" / "primaschema" / "schema"
INDEX_SCHEMA_PATH = SCHEMA_DIR / INDEX_FILE_NAME
INDEX_HEADER_PATH = SCHEMA_DIR / "index-header.yml"

logging.getLogger("primaschema").addHandler(logging.NullHandler())
