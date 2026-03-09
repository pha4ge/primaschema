"""Infrastructure for tiled amplicon PCR primer scheme definitions"""

import logging
import os
from pathlib import Path

METADATA_FILE_NAME: str = "info.json"
PRIMER_FILE_NAME: str = "primer.bed"
REFERENCE_FILE_NAME: str = "reference.fasta"

PKG_DIR = Path(
    os.environ.get(
        "PRIMASCHEMA_ROOT_PATH", Path(__file__).absolute().parent.parent.parent
    )
)
SCHEMA_DIR = PKG_DIR / "src" / "primaschema" / "schema"
MANIFEST_SCHEMA_PATH = SCHEMA_DIR / "manifest.json"
MANIFEST_HEADER_PATH = SCHEMA_DIR / "manifest-header.yml"

logging.getLogger("primaschema").addHandler(logging.NullHandler())
