"""Fetch the published ProteomeTools ``.msp`` spectral libraries from Zenodo.

The 2019 release (https://www.proteometools.org/index.php?id=53) is archived as
Zenodo record ``15705607`` — eight per-mode zips (FTMS/ITMS × HCD/CID × NCE),
each containing ``.msp`` files, with an MD5 published per zip. We download +
verify the zips and extract the ``.msp`` members.

Ingesting the ``.msp`` into a Constellation ``massspec.library.Library`` and
associating them with the raw acquisitions is a **separate, later** task; this
module stops at verified-and-extracted ``.msp`` on disk.
"""

from __future__ import annotations

import json
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path

from constellation.core.progress import ProgressCallback

from ..fetch import net

DEFAULT_RECORD_ID = "15705607"
_API = "https://zenodo.org/api/records"
_USER_AGENT = "protostar/0.1 (+https://github.com/wilburn-lab/protostar)"


@dataclass(frozen=True, slots=True)
class ZenodoFile:
    key: str  # filename, e.g. FTMS_HCD_20_annotated_2019-11-12.zip
    size_bytes: int
    md5: str  # hex digest (the ``md5:`` prefix stripped)
    url: str  # direct-content download URL

    @property
    def mode(self) -> str:
        """Library mode token, e.g. ``FTMS_HCD_20`` (drops ``_annotated_<date>``)."""
        stem = self.key.split("_annotated_")[0]
        return stem or Path(self.key).stem


def list_library_files(
    record_id: str = DEFAULT_RECORD_ID, *, timeout: int = 120
) -> list[ZenodoFile]:
    """List the library zips for a Zenodo record (sorted by mode token)."""
    req = urllib.request.Request(
        f"{_API}/{record_id}",
        headers={"User-Agent": _USER_AGENT, "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        rec = json.loads(resp.read())
    out: list[ZenodoFile] = []
    for f in rec.get("files", []):
        checksum = f.get("checksum", "") or ""
        md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else checksum
        out.append(
            ZenodoFile(
                key=f["key"],
                size_bytes=int(f.get("size") or 0),
                md5=md5,
                url=(f.get("links", {}) or {}).get("self", ""),
            )
        )
    return sorted(out, key=lambda z: z.mode)


def _select(files: list[ZenodoFile], modes: "set[str] | None") -> list[ZenodoFile]:
    if modes is None:
        return files
    return [f for f in files if f.mode in modes]


def fetch_libraries(
    dest_dir: "str | Path",
    *,
    modes: "set[str] | None" = None,
    record_id: str = DEFAULT_RECORD_ID,
    resume: bool = True,
    progress_cb: ProgressCallback | None = None,
) -> list[Path]:
    """Download (resumable) + MD5-verify the selected library zips.

    Returns the local zip paths. ``modes`` filters by mode token
    (e.g. ``{"FTMS_HCD_28", "ITMS_CID_35"}``); ``None`` fetches all eight.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    files = _select(list_library_files(record_id), modes)
    out: list[Path] = []
    for f in files:
        dest = dest_dir / f.key
        net.stream_download(
            f.url,
            dest,
            expected_size=f.size_bytes,
            resume=resume,
            compute_sha1=False,
            progress_cb=progress_cb,
        )
        state = net.verify(dest, expected_md5=f.md5, expected_size=f.size_bytes)
        if state != "ok":
            raise OSError(f"library {f.key} failed verification ({state}; md5 {f.md5})")
        out.append(dest)
    return out


def extract_msp(zip_paths: list[Path], dest_dir: "str | Path") -> list[Path]:
    """Extract ``.msp`` members from each zip into ``<dest_dir>/<mode>/``.

    The mode token is derived from the zip filename; members are written by
    basename (nested zip directories are flattened).
    """
    dest_dir = Path(dest_dir)
    written: list[Path] = []
    for zp in zip_paths:
        zp = Path(zp)
        mode = ZenodoFile(zp.name, 0, "", "").mode
        out_dir = dest_dir / mode
        with zipfile.ZipFile(zp) as zf:
            members = [m for m in zf.namelist() if m.lower().endswith(".msp")]
            for m in members:
                target = out_dir / Path(m).name
                target.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(m) as src, target.open("wb") as dst:
                    while block := src.read(1 << 20):
                        dst.write(block)
                written.append(target)
    return written


__all__ = [
    "DEFAULT_RECORD_ID",
    "ZenodoFile",
    "extract_msp",
    "fetch_libraries",
    "list_library_files",
]
