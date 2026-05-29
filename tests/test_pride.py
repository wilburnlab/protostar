"""Tests for the PRIDE v3 client (pure parsing + an opt-in live check)."""

from __future__ import annotations

import pytest

from protostar.fetch import pride


def test_parse_file_rewrites_ftp_to_https():
    obj = {
        "fileName": "pool.raw",
        "fileSizeBytes": 123,
        "checksum": "abc123",
        "fileCategory": {"value": "RAW"},
        "publicFileLocations": [
            {"name": "Aspera Protocol", "value": "prd_ascp@fasp.ebi.ac.uk:pride/x/pool.raw"},
            {"name": "FTP Protocol", "value": "ftp://ftp.pride.ebi.ac.uk/pride/x/pool.raw"},
        ],
    }
    pf = pride._parse_file(obj)
    assert pf.file_name == "pool.raw"
    assert pf.size_bytes == 123
    assert pf.sha1 == "abc123"
    assert pf.category == "RAW" and pf.is_raw
    assert pf.https_url == "https://ftp.pride.ebi.ac.uk/pride/x/pool.raw"


def test_parse_file_handles_missing_checksum_and_category():
    pf = pride._parse_file({"fileName": "x.zip", "fileSizeBytes": 0, "checksum": ""})
    assert pf.sha1 is None  # empty string → None (Wilhelm2021 case)
    assert pf.category == "OTHER"
    assert pf.https_url is None


@pytest.mark.network
def test_list_files_live_counts():
    files = pride.list_files("PXD004732", page_size=100)
    # >> one page: confirms pagination doesn't stop at the 100-item cap.
    assert len(files) > 100
    cats = {f.category for f in files}
    assert "RAW" in cats and "SEARCH" in cats
    assert all(f.https_url and f.https_url.startswith("https://") for f in files if f.is_raw)
