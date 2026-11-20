"""Bundle-v3 access for protostar's data root.

Converted runs live under ``<data_root>/proc/<dataset>/``. Constellation's bundle
format v3 (constellation #104, 2026-10-31) writes each run in one pass as::

    proc/<dataset>/<stem>/            centroid ("peak") bundle: manifest.json,
                                      peaks/part-*.parquet, scan_metadata.parquet,
                                      acquisition_metadata.parquet[, precursor_major/]
    proc/<dataset>/<stem>.profile/    sibling profile bundle (MS1 profile points)

Bundles written before v3 sit in the June 2026 two-pass layout,
``proc/<dataset>/{centroid,profile}/<stem>/``, with a single ``peaks.parquet``.
Constellation refuses them (``OutdatedBundleError``), and there is no upgrader;
they are kept only until the v3 reconversion of their dataset verifies.

This module is the one place protostar decides which directory holds a run and
whether it is usable. Callers never glob ``*/manifest.json`` or open
``peaks.parquet`` themselves. Those habits double-count ``.profile`` siblings and
silently find nothing on v3. An outdated bundle is always reported, never
skipped quietly.

Glue only: reading, versioning and peak access all delegate to
``constellation.massspec.bundle``.
"""

from __future__ import annotations

import subprocess
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import pyarrow as pa
import pyarrow.parquet as pq
from constellation.massspec.bundle import (
    MANIFEST_FILENAME,
    PROFILE_SUFFIX,
    OutdatedBundleError,
    is_profile_bundle,
    read_bundle_manifest,
)

#: Scan metadata file at a bundle's root (named in ``constellation.massspec.bundle``'s
#: module docstring and the manifest's ``outputs``).
SCAN_METADATA_FILENAME = "scan_metadata.parquet"

#: The June 2026 pre-v3 layout's per-mode parents under ``proc/<dataset>/``.
LEGACY_MODES: tuple[str, ...] = ("centroid", "profile")

BundleState = Literal["current", "outdated", "incomplete", "missing"]


def proc_dir(data_root: str | Path, dataset: str) -> Path:
    """``<data_root>/proc/<dataset>``: the parent of every v3 run bundle."""
    return Path(data_root) / "proc" / dataset


def run_bundle_dir(data_root: str | Path, dataset: str, stem: str) -> Path:
    """Where the v3 centroid bundle for ``stem`` lives (whether or not it exists)."""
    return proc_dir(data_root, dataset) / stem


def legacy_dir(data_root: str | Path, dataset: str, mode: str = "centroid") -> Path:
    """The June pre-v3 parent ``proc/<dataset>/<mode>`` (``mode`` in ``LEGACY_MODES``)."""
    if mode not in LEGACY_MODES:
        raise ValueError(f"mode must be one of {LEGACY_MODES}; got {mode!r}")
    return proc_dir(data_root, dataset) / mode


def bundle_state(bundle_dir: str | Path) -> BundleState:
    """Classify a bundle directory.

    - ``missing``: no directory.
    - ``incomplete``: a directory with no manifest (an unfinished conversion).
    - ``outdated``: a manifest Constellation refuses (pre-v3).
    - ``current``: a manifest Constellation opens.
    """
    d = Path(bundle_dir)
    if not d.is_dir():
        return "missing"
    if not (d / MANIFEST_FILENAME).is_file():
        return "incomplete"
    try:
        read_bundle_manifest(d)
    except OutdatedBundleError:
        return "outdated"
    return "current"


def _is_run_candidate(child: Path) -> bool:
    return (
        child.is_dir()
        and child.name not in LEGACY_MODES
        and not child.name.endswith(PROFILE_SUFFIX)
        and not child.name.startswith(".")
    )


def iter_run_bundles(data_root: str | Path, dataset: str) -> Iterator[Path]:
    """Current v3 centroid run bundles for ``dataset``, sorted by stem.

    Skips profile siblings, the legacy ``centroid/`` and ``profile/`` parents, and
    any directory that is not a current bundle. Use :func:`inventory` to see what
    was skipped and why.
    """
    parent = proc_dir(data_root, dataset)
    if not parent.is_dir():
        return
    for child in sorted(parent.iterdir(), key=lambda p: p.name):
        if not _is_run_candidate(child):
            continue
        if bundle_state(child) == "current" and not is_profile_bundle(child):
            yield child


@dataclass
class Inventory:
    """What is on disk for one dataset, by state. Values are sorted stems."""

    current: list[str] = field(default_factory=list)
    outdated: list[str] = field(default_factory=list)
    incomplete: list[str] = field(default_factory=list)
    with_profile: list[str] = field(default_factory=list)
    legacy: dict[str, int] = field(default_factory=dict)  # mode -> bundle count

    def summary(self) -> str:
        legacy = ", ".join(f"{m}={n}" for m, n in sorted(self.legacy.items()) if n) or "none"
        return (
            f"v3 current={len(self.current)} (with profile={len(self.with_profile)}), "
            f"outdated={len(self.outdated)}, incomplete={len(self.incomplete)}; "
            f"legacy pre-v3: {legacy}"
        )


def inventory(data_root: str | Path, dataset: str) -> Inventory:
    """Classify every run directory under ``proc/<dataset>`` plus the legacy trees."""
    inv = Inventory()
    parent = proc_dir(data_root, dataset)
    if not parent.is_dir():
        return inv
    for child in sorted(parent.iterdir(), key=lambda p: p.name):
        if not _is_run_candidate(child):
            continue
        state = bundle_state(child)
        if state == "current":
            inv.current.append(child.name)
            if bundle_state(child.with_name(child.name + PROFILE_SUFFIX)) == "current":
                inv.with_profile.append(child.name)
        elif state == "outdated":
            inv.outdated.append(child.name)
        elif state == "incomplete":
            inv.incomplete.append(child.name)
    for mode in LEGACY_MODES:
        ld = legacy_dir(data_root, dataset, mode)
        inv.legacy[mode] = sum(1 for m in ld.glob(f"*/{MANIFEST_FILENAME}")) if ld.is_dir() else 0
    return inv


def require_run_bundle(data_root: str | Path, dataset: str, stem: str) -> Path | None:
    """The current v3 bundle for ``stem``, or ``None`` if the run was never converted.

    Raises ``OutdatedBundleError`` if a pre-v3 bundle sits at the v3 location.
    Raises ``RuntimeError`` if the conversion is incomplete. Both are things a
    caller must not skip silently. A run that exists only in the legacy
    ``centroid/`` tree is reported as ``None``: it has not been reconverted.
    """
    bundle = run_bundle_dir(data_root, dataset, stem)
    state = bundle_state(bundle)
    if state == "current":
        return bundle
    if state == "missing":
        return None
    if state == "incomplete":
        raise RuntimeError(f"{bundle} has no {MANIFEST_FILENAME}: an unfinished conversion")
    # Re-raise Constellation's own error (it names the fix: convert again).
    read_bundle_manifest(bundle)
    raise AssertionError("unreachable: an outdated bundle opened")  # pragma: no cover


def read_scan_metadata(bundle_dir: str | Path, columns: list[str] | None = None) -> pa.Table:
    """A current bundle's ``scan_metadata.parquet``, optionally column-projected.

    Only the columns present are read when ``columns`` is given, so callers can ask
    for optional trailer fields across reader versions.
    """
    d = Path(bundle_dir)
    read_bundle_manifest(d)  # refuses an outdated bundle
    path = d / SCAN_METADATA_FILENAME
    if columns is None:
        return pq.read_table(path)
    present = set(pq.read_schema(path).names)
    return pq.read_table(path, columns=[c for c in columns if c in present])


def constellation_provenance() -> str:
    """A one-line record of the Constellation build in use, for run logs.

    Gives the package version plus the git commit when Constellation is an
    editable checkout. That is how the lab (and the OSC convert jobs) install it.
    """
    import constellation

    pkg = Path(constellation.__file__).resolve().parent
    version = getattr(constellation, "__version__", "?")
    sha = "unknown"
    try:
        out = subprocess.run(
            ["git", "-C", str(pkg.parent), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        if out.returncode == 0:
            sha = out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return f"constellation {version} @ {sha} ({pkg})"


__all__ = [
    "BundleState",
    "Inventory",
    "LEGACY_MODES",
    "SCAN_METADATA_FILENAME",
    "bundle_state",
    "constellation_provenance",
    "inventory",
    "iter_run_bundles",
    "legacy_dir",
    "proc_dir",
    "read_scan_metadata",
    "require_run_bundle",
    "run_bundle_dir",
]
