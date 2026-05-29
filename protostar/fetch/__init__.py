"""Manifest building and raw-file acquisition.

Canonical path: query the ProteomeXchange/PRIDE API for each dataset's ``.raw``
file list and published checksums, download missing files, verify against the
checksum, re-fetch on mismatch (resumable + repairable). An optional
``--seed-from <dir>`` matches files already present on ESS and moves/hardlinks
them in instead of downloading.
"""

from .manifest import (
    DEFAULT_FETCH_CATEGORIES,
    FileStatus,
    Manifest,
    ManifestEntry,
    build_manifest,
    load_manifest,
    local_path,
    manifest_path,
    reconcile,
    summarize,
    write_manifest,
    write_status,
)
from .net import DownloadResult, md5_of, sha1_of, stream_download, verify
from .pride import PrideFile, list_files
from .seed import SeedResult, seed_dataset

__all__ = [
    "DEFAULT_FETCH_CATEGORIES",
    "DownloadResult",
    "FileStatus",
    "Manifest",
    "ManifestEntry",
    "PrideFile",
    "SeedResult",
    "build_manifest",
    "list_files",
    "load_manifest",
    "local_path",
    "manifest_path",
    "md5_of",
    "reconcile",
    "seed_dataset",
    "sha1_of",
    "stream_download",
    "summarize",
    "verify",
    "write_manifest",
    "write_status",
]
