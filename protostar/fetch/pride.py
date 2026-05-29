"""ProteomeXchange / PRIDE Archive (v3) file listing.

Constellation has no PRIDE client. The v3 API
(``/pride/ws/archive/v3/projects/{accession}/files``) returns a paginated
array; per file we keep ``fileName``, ``fileSizeBytes``, ``checksum`` (SHA-1),
the ``fileCategory`` CV-param value (``RAW`` / ``SEARCH`` / ...), and the
``publicFileLocations`` (FTP + Aspera). We retain every category (the manifest
is a complete record) and rewrite each FTP URL to its resumable HTTPS mirror.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass

_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
_USER_AGENT = "protostar/0.1 (+https://github.com/wilburn-lab/protostar)"
_FTP_PREFIX = "ftp://ftp.pride.ebi.ac.uk/"
_HTTPS_PREFIX = "https://ftp.pride.ebi.ac.uk/"
# The v3 API silently caps pageSize at 100 — requesting more still returns 100,
# so a naive "len(batch) < requested" end-of-list test stops after one page.
_MAX_PAGE_SIZE = 100


@dataclass(frozen=True, slots=True)
class PrideFile:
    """One file record from the PRIDE v3 ``files`` listing."""

    file_name: str
    size_bytes: int
    sha1: str | None
    category: str  # RAW / SEARCH / RESULT / OTHER / ...
    https_url: str | None  # resumable mirror; None if no FTP location published

    @property
    def is_raw(self) -> bool:
        return self.category == "RAW"


def _get_json(url: str, *, timeout: int = 120):
    req = urllib.request.Request(
        url, headers={"User-Agent": _USER_AGENT, "Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read())


def _https_from_locations(locations: "list[dict] | None") -> str | None:
    """Pick the FTP location and rewrite it to the resumable HTTPS mirror.

    ``ftp://ftp.pride.ebi.ac.uk/...`` and ``https://ftp.pride.ebi.ac.uk/...``
    serve the same tree; the HTTPS host advertises ``Accept-Ranges: bytes``.
    """
    for loc in locations or []:
        if loc.get("name") == "FTP Protocol":
            value = loc.get("value", "") or ""
            if value.startswith(_FTP_PREFIX):
                return _HTTPS_PREFIX + value[len(_FTP_PREFIX) :]
            return value or None
    return None


def _parse_file(obj: dict) -> PrideFile:
    category = (obj.get("fileCategory") or {}).get("value") or "OTHER"
    return PrideFile(
        file_name=obj["fileName"],
        size_bytes=int(obj.get("fileSizeBytes") or 0),
        sha1=obj.get("checksum") or None,
        category=category,
        https_url=_https_from_locations(obj.get("publicFileLocations")),
    )


def list_files(
    accession: str,
    *,
    categories: "set[str] | None" = None,
    page_size: int = 100,
    timeout: int = 120,
) -> list[PrideFile]:
    """List files for a PXD accession, paginating until the array is exhausted.

    ``categories`` (e.g. ``{"RAW", "SEARCH"}``) filters the result; ``None``
    keeps every file. Sorted by file name for deterministic manifests.
    """
    effective = min(page_size, _MAX_PAGE_SIZE)
    out: list[PrideFile] = []
    page = 0
    while True:
        url = f"{_API}/projects/{accession}/files?pageSize={effective}&page={page}"
        batch = _get_json(url, timeout=timeout)
        if not batch:
            break
        out.extend(_parse_file(o) for o in batch)
        if len(batch) < effective:
            break
        page += 1
    if categories is not None:
        out = [f for f in out if f.category in categories]
    return sorted(out, key=lambda f: f.file_name)


__all__ = ["PrideFile", "list_files"]
