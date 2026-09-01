import gzip

import httpx
import pytest

from primaschema.get_scheme import (
    DownloadError,
    _decompress_gzip_bounded,
    _ensure_https,
    _reject_non_https_request,
    _resolve_timeout,
    _resolve_workers,
)


def test_ensure_https_accepts_https():
    _ensure_https("https://example.com/index.json")


def test_ensure_https_rejects_http():
    with pytest.raises(ValueError, match="Only https URLs are allowed"):
        _ensure_https("http://example.com/index.json")


def test_reject_non_https_request_accepts_https():
    _reject_non_https_request(httpx.Request("GET", "https://example.com"))


def test_reject_non_https_request_rejects_http():
    """This is the hook that fires on every redirect hop, not just the first request."""
    with pytest.raises(ValueError, match="Only https URLs are allowed"):
        _reject_non_https_request(httpx.Request("GET", "http://example.com"))


def test_resolve_timeout_defaults_when_none():
    from primaschema.get_scheme import DEFAULT_HTTP_TIMEOUT_SECONDS

    timeout = _resolve_timeout(None)
    assert timeout.connect == DEFAULT_HTTP_TIMEOUT_SECONDS


def test_resolve_timeout_uses_given_value():
    timeout = _resolve_timeout(5.0)
    assert timeout.connect == 5.0


def test_resolve_workers_uses_default_when_none(monkeypatch):
    monkeypatch.setattr("primaschema.get_scheme.MAX_DOWNLOAD_WORKERS", "4")
    assert _resolve_workers(None) == 4


def test_resolve_workers_falls_back_on_invalid_env(monkeypatch):
    monkeypatch.setattr("primaschema.get_scheme.MAX_DOWNLOAD_WORKERS", "not-a-number")
    from primaschema.get_scheme import DEFAULT_MAX_WORKERS

    assert _resolve_workers(None) == DEFAULT_MAX_WORKERS


def test_resolve_workers_clamps_below_minimum():
    assert _resolve_workers(0) == 1


def test_resolve_workers_clamps_above_limit():
    assert _resolve_workers(1000) == 32


def test_decompress_gzip_bounded_succeeds_within_limit():
    payload = gzip.compress(b"hello world")
    assert _decompress_gzip_bounded(payload, max_bytes=1000) == b"hello world"


def test_decompress_gzip_bounded_raises_when_over_limit():
    """A perfectly ordinary small gzip payload, decompressed against a
    deliberately tiny cap"""
    payload = gzip.compress(b"x" * 1000)
    with pytest.raises(DownloadError, match="exceeds max size"):
        _decompress_gzip_bounded(payload, max_bytes=10)
