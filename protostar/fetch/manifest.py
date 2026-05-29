"""Expected-file manifests + filesystem reconciliation.

Two layers:

* **Layer 1 — expected manifest** (committed at ``config/manifests/{ds}.json``):
  the reproducible record of what *should* exist, built once from the PRIDE v3
  API. Captures every category (RAW + SEARCH + the odd metadata file) so a
  fresh user can fetch + verify without re-querying PRIDE.

* **Layer 2 — reconciliation** (:func:`reconcile`): compares a manifest against
  the canonical data tree and classifies each file ``missing`` / ``partial`` /
  ``present`` / ``corrupt``. Cheap (size-only) by default; ``verify=True`` adds
  the streamed SHA-1 check. An optional :func:`write_status` snapshot lands on
  ESS under ``<data_root>/status/`` (gitignored — never committed).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import net
from .pride import PrideFile, list_files

PRIDE_API_VERSION = "v3"

# Category → data-tree subdirectory. Anything unmapped lands under ``other/``.
CATEGORY_DIRS = {"RAW": "raw", "SEARCH": "search"}
_OTHER_DIR = "other"
DEFAULT_FETCH_CATEGORIES = ("RAW", "SEARCH")


# ── Layer 1: the expected manifest ───────────────────────────────────────


@dataclass(frozen=True, slots=True)
class ManifestEntry:
    file_name: str
    size_bytes: int
    sha1: str | None
    https_url: str | None
    file_category: str

    @classmethod
    def from_pride(cls, pf: PrideFile) -> "ManifestEntry":
        return cls(pf.file_name, pf.size_bytes, pf.sha1, pf.https_url, pf.category)

    @classmethod
    def from_dict(cls, d: dict) -> "ManifestEntry":
        return cls(
            file_name=d["file_name"],
            size_bytes=int(d["size_bytes"]),
            sha1=d.get("sha1"),
            https_url=d.get("https_url"),
            file_category=d["file_category"],
        )

    def to_dict(self) -> dict:
        return {
            "file_name": self.file_name,
            "size_bytes": self.size_bytes,
            "sha1": self.sha1,
            "https_url": self.https_url,
            "file_category": self.file_category,
        }


@dataclass(frozen=True, slots=True)
class Manifest:
    dataset: str
    accession: str
    generated_at: str
    pride_api_version: str
    entries: tuple[ManifestEntry, ...]

    def by_category(self, category: str) -> list[ManifestEntry]:
        return [e for e in self.entries if e.file_category == category]

    def fetch_entries(self, categories: "tuple[str, ...] | None") -> list[ManifestEntry]:
        """Entries to fetch — every category when ``categories`` is ``None``."""
        if categories is None:
            return list(self.entries)
        wanted = set(categories)
        return [e for e in self.entries if e.file_category in wanted]

    def to_dict(self) -> dict:
        n_by_cat: dict[str, int] = {}
        for e in self.entries:
            n_by_cat[e.file_category] = n_by_cat.get(e.file_category, 0) + 1
        return {
            "dataset": self.dataset,
            "accession": self.accession,
            "generated_at": self.generated_at,
            "pride_api_version": self.pride_api_version,
            "n_files": len(self.entries),
            "n_by_category": n_by_cat,
            "files": [e.to_dict() for e in self.entries],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Manifest":
        return cls(
            dataset=d["dataset"],
            accession=d["accession"],
            generated_at=d["generated_at"],
            pride_api_version=d.get("pride_api_version", PRIDE_API_VERSION),
            entries=tuple(ManifestEntry.from_dict(f) for f in d["files"]),
        )


def manifest_path(dataset: str, *, manifest_dir: "str | Path") -> Path:
    return Path(manifest_dir) / f"{dataset}.json"


def build_manifest(dataset: str, accession: str, *, page_size: int = 100) -> Manifest:
    """Query PRIDE for every file under ``accession`` → an expected manifest."""
    files = list_files(accession, page_size=page_size)
    return Manifest(
        dataset=dataset,
        accession=accession,
        generated_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        pride_api_version=PRIDE_API_VERSION,
        entries=tuple(ManifestEntry.from_pride(f) for f in files),
    )


def write_manifest(manifest: Manifest, path: "str | Path") -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest.to_dict(), indent=2) + "\n", encoding="utf-8")
    return path


def load_manifest(path: "str | Path") -> Manifest:
    return Manifest.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


# ── routing + Layer 2: reconciliation ────────────────────────────────────


def subdir_for(category: str) -> str:
    return CATEGORY_DIRS.get(category, _OTHER_DIR)


def local_path(entry: ManifestEntry, data_root: "str | Path", dataset: str) -> Path:
    """Canonical on-disk location: ``<data_root>/<subdir>/<dataset>/<file>``."""
    return Path(data_root) / subdir_for(entry.file_category) / dataset / entry.file_name


@dataclass(frozen=True, slots=True)
class FileStatus:
    entry: ManifestEntry
    state: net.VerifyState  # missing / partial / present(ok) / corrupt
    path: Path
    observed_size: int | None
    observed_sha1: str | None


def reconcile(
    manifest: Manifest,
    data_root: "str | Path",
    *,
    categories: "tuple[str, ...] | None" = DEFAULT_FETCH_CATEGORIES,
    verify: bool = False,
) -> list[FileStatus]:
    """Classify each expected file against the data tree.

    Size-only by default (never reads file bodies); ``verify=True`` adds the
    streamed SHA-1 check where a published checksum exists (Wilhelm2021 has
    none — those stay size-only regardless).
    """
    statuses: list[FileStatus] = []
    for e in manifest.fetch_entries(categories):
        p = local_path(e, data_root, dataset=manifest.dataset)
        state = net.verify(
            p,
            expected_size=e.size_bytes,
            expected_sha1=e.sha1 if verify else None,
        )
        size = p.stat().st_size if p.exists() else None
        sha = net.sha1_of(p) if (verify and state == "ok" and e.sha1) else None
        statuses.append(FileStatus(e, state, p, size, sha))
    return statuses


def summarize(statuses: list[FileStatus]) -> dict[str, dict[str, int]]:
    """Per-category counts keyed by state — the body of the dry-run table."""
    table: dict[str, dict[str, int]] = {}
    for s in statuses:
        row = table.setdefault(s.entry.file_category, {})
        row[s.state] = row.get(s.state, 0) + 1
    return table


def status_path(dataset: str, data_root: "str | Path") -> Path:
    return Path(data_root) / "status" / f"{dataset}.json"


def write_status(
    dataset: str,
    data_root: "str | Path",
    statuses: list[FileStatus],
    *,
    extra: "dict | None" = None,
) -> Path:
    """Persist a Layer-2 snapshot (audit only; gitignored under ``data/``)."""
    path = status_path(dataset, data_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {
        "dataset": dataset,
        "checked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "summary": summarize(statuses),
        "files": [
            {
                "file_name": s.entry.file_name,
                "file_category": s.entry.file_category,
                "state": s.state,
                "observed_size": s.observed_size,
                "observed_sha1": s.observed_sha1,
            }
            for s in statuses
        ],
        **(extra or {}),
    }
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return path


__all__ = [
    "CATEGORY_DIRS",
    "DEFAULT_FETCH_CATEGORIES",
    "FileStatus",
    "Manifest",
    "ManifestEntry",
    "build_manifest",
    "load_manifest",
    "local_path",
    "manifest_path",
    "reconcile",
    "subdir_for",
    "summarize",
    "status_path",
    "write_manifest",
    "write_status",
]
