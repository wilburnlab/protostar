"""Streaming, resumable HTTP downloads + checksum verification.

The one genuinely novel primitive protostar owns. Constellation's
``http_get_bytes`` loads whole responses into memory (fine for small catalog
files, not for 800 MB Thermo ``.raw`` files), and its checksum helpers cover
SHA-256 / MD5 only. PRIDE publishes **SHA-1** and Zenodo publishes **MD5**, so
we add streamed digests and a Range-resumable downloader here.

Download integrity (the PRIDE SHA-1 / Zenodo MD5 checked here) is a separate
concern from conversion provenance (the SHA-256 that
``constellation.massspec.io.thermo.convert`` records in each bundle manifest);
the two are never cross-compared.
"""

from __future__ import annotations

import hashlib
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from constellation.core.progress import (
    ProgressCallback,
    emit_done,
    emit_progress,
    emit_start,
)

_CHUNK = 1 << 20  # 1 MiB, matches constellation's checksum helpers
_USER_AGENT = "protostar/0.1 (+https://github.com/wilburn-lab/protostar)"
_STAGE = "download"

VerifyState = Literal["ok", "corrupt", "missing", "partial"]


@dataclass(frozen=True, slots=True)
class DownloadResult:
    """Outcome of one :func:`stream_download`."""

    path: Path
    n_bytes: int
    sha1: str | None
    resumed: bool


def _part_path(dest: Path) -> Path:
    return dest.parent / (dest.name + ".part")


def _hash_of(path: Path, algo: str, *, chunk: int = _CHUNK) -> str:
    h = hashlib.new(algo)
    with Path(path).open("rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def sha1_of(path: "str | Path", *, chunk: int = _CHUNK) -> str:
    """Streamed SHA-1 hex digest (PRIDE publishes SHA-1 checksums)."""
    return _hash_of(Path(path), "sha1", chunk=chunk)


def md5_of(path: "str | Path", *, chunk: int = _CHUNK) -> str:
    """Streamed MD5 hex digest (Zenodo publishes MD5 checksums)."""
    return _hash_of(Path(path), "md5", chunk=chunk)


def verify(
    path: "str | Path",
    *,
    expected_sha1: str | None = None,
    expected_md5: str | None = None,
    expected_size: int | None = None,
) -> VerifyState:
    """Classify a local file against expectations.

    ``"missing"`` if neither the final file nor a ``.part`` exists,
    ``"partial"`` if only a ``.part`` is present, ``"corrupt"`` if the final
    file fails a size or checksum check, else ``"ok"``. Checks are cheap-first
    (size before hash) so passing only ``expected_size`` never reads the file.
    """
    path = Path(path)
    if not path.exists():
        return "partial" if _part_path(path).exists() else "missing"
    if expected_size is not None and path.stat().st_size != expected_size:
        return "corrupt"
    if expected_sha1 is not None and sha1_of(path).lower() != expected_sha1.lower():
        return "corrupt"
    if expected_md5 is not None and md5_of(path).lower() != expected_md5.lower():
        return "corrupt"
    return "ok"


def stream_download(
    url: str,
    dest: "str | Path",
    *,
    expected_size: int | None = None,
    resume: bool = True,
    compute_sha1: bool = True,
    timeout: int = 600,
    progress_cb: ProgressCallback | None = None,
) -> DownloadResult:
    """Stream ``url`` to ``dest`` with HTTP Range resume + atomic finalize.

    Writes to ``<dest>.part`` and ``os.replace``s to ``dest`` only after a
    clean, full-length transfer, so a final-named file is complete by
    construction (a killed download never masquerades as done). When a
    ``.part`` exists and ``resume`` is set, requests ``Range: bytes=<have>-``;
    if the server honours it (HTTP 206) the partial bytes are kept and the
    SHA-1 is seeded from them, otherwise the transfer restarts from zero.

    Returns the byte count and (when ``compute_sha1``) the SHA-1 of the
    finished file, computed in a single streaming pass — the caller verifies
    it against the published checksum.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = _part_path(dest)

    # Already complete (size matches what we expect) — don't re-download.
    if dest.exists() and expected_size is not None and dest.stat().st_size == expected_size:
        sha = sha1_of(dest) if compute_sha1 else None
        return DownloadResult(dest, dest.stat().st_size, sha, resumed=False)

    have = part.stat().st_size if (resume and part.exists()) else 0
    # A .part at/above the expected size is stale — discard before requesting,
    # else the server answers 416 (Range Not Satisfiable).
    if have and expected_size is not None and have >= expected_size:
        part.unlink()
        have = 0

    headers = {"User-Agent": _USER_AGENT}
    if have:
        headers["Range"] = f"bytes={have}-"
    req = urllib.request.Request(url, headers=headers)

    emit_start(progress_cb, _STAGE, total=expected_size or 0, message=dest.name)

    try:
        resp = urllib.request.urlopen(req, timeout=timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 416 and have:
            # Stale partial we couldn't pre-detect (no expected_size) — restart.
            part.unlink(missing_ok=True)
            return stream_download(
                url,
                dest,
                expected_size=expected_size,
                resume=False,
                compute_sha1=compute_sha1,
                timeout=timeout,
                progress_cb=progress_cb,
            )
        raise

    h = hashlib.sha1() if compute_sha1 else None
    with resp:
        status = getattr(resp, "status", None) or resp.getcode()
        resumed = have > 0 and status == 206
        if resumed:
            if h is not None:  # seed the rolling hash from bytes already on disk
                with part.open("rb") as existing:
                    while block := existing.read(_CHUNK):
                        h.update(block)
            mode, completed = "ab", have
        else:
            mode, completed = "wb", 0  # 200 (Range ignored) or nothing to resume
        with part.open(mode) as out:
            while block := resp.read(_CHUNK):
                out.write(block)
                if h is not None:
                    h.update(block)
                completed += len(block)
                emit_progress(progress_cb, _STAGE, completed=completed, total=expected_size or 0)

    n_bytes = part.stat().st_size
    if expected_size is not None and n_bytes != expected_size:
        raise OSError(
            f"download size mismatch for {dest.name}: got {n_bytes} bytes, "
            f"expected {expected_size}; partial left at {part}"
        )
    os.replace(part, dest)
    emit_done(progress_cb, _STAGE, completed=n_bytes, total=n_bytes, message=dest.name)
    return DownloadResult(dest, n_bytes, h.hexdigest() if h is not None else None, resumed)


__all__ = [
    "DownloadResult",
    "VerifyState",
    "md5_of",
    "sha1_of",
    "stream_download",
    "verify",
]
