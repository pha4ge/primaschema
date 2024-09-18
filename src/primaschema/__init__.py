"""Infrastructure for tiled amplicon PCR primer scheme definitions"""

import logging
import logging.config
import os

from pathlib import Path

from platformdirs import user_data_dir


__version__ = "1.0.0a0"


SCHEMES_ARCHIVE_URL = os.environ.get(
    "PRIMASCHEMA_SCHEMES_ARCHIVE_URL",
    "https://github.com/pha4ge/primer-schemes/archive/refs/heads/main.tar.gz",
)
CACHE_DIR = Path(
    os.environ.get("PRIMASCHEMA_CACHE_DIR") or user_data_dir("primaschema", "PHA4GE")
)

PKG_DIR = Path(
    os.environ.get(
        "PRIMASCHEMA_ROOT_PATH", Path(__file__).absolute().parent.parent.parent
    )
)
SCHEMA_DIR = PKG_DIR / "src" / "primaschema" / "schema"
MANIFEST_SCHEMA_PATH = SCHEMA_DIR / "manifest.json"
MANIFEST_HEADER_PATH = SCHEMA_DIR / "manifest-header.yml"

logging_config = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "default": {
            "format": "%(name)s %(levelname)s: %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "default",
            "level": "INFO",
            "stream": "ext://sys.stderr",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",
    },
    "loggers": {
        "primaschema": {
            "handlers": ["console"],
            "level": "INFO",
            "propagate": False,
        },
    },
}

logging.config.dictConfig(logging_config)
logger = logging.getLogger("primaschema")

logger.debug(f"{PKG_DIR=} {SCHEMA_DIR=}")
