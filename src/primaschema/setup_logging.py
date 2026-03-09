"""Logging helpers for primaschema."""

import logging
import logging.config
from enum import Enum


class LogLevel(str, Enum):
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


_BASE_LOGGING_CONFIG = {
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


def _resolve_level(log_level: LogLevel | str | None) -> int:
    if log_level is None:
        return logging.INFO

    if isinstance(log_level, LogLevel):
        level_name = log_level.value
    else:
        level_name = str(log_level).upper()

    numeric_level = getattr(logging, level_name, None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")
    return numeric_level


def configure_logging(
    log_level: LogLevel | str | None = None,
) -> int:
    """Configure primaschema logging.

    If log_level is provided it takes precedence over verbose.
    """
    numeric_level = _resolve_level(log_level)

    logging.config.dictConfig(_BASE_LOGGING_CONFIG)
    logger = logging.getLogger("primaschema")
    logger.setLevel(numeric_level)
    for handler in logger.handlers:
        handler.setLevel(numeric_level)

    return numeric_level
