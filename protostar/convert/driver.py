"""Drive Constellation's Thermo ``.raw`` → parquet-bundle conversion.

Thin orchestration over ``constellation.massspec.io.thermo.convert_batch``:
protostar owns only file selection (which ``.raw`` to convert, shard/limit for
SLURM arrays) and result reconciliation. The conversion itself — spawn-mode
workers, skip-on-existing-``manifest.json`` resume, source SHA-256, the
``peaks.parquet`` / ``scan_metadata.parquet`` / ``acquisition_metadata.parquet``
bundle — all lives in Constellation. We never reimplement its skip logic.

Output layout keeps centroid and profile side by side so the same ``.raw`` can
be converted in both modes without collision::

    <data_root>/proc/<dataset>/<centroid|profile>/<stem>/{manifest.json, *.parquet}
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from constellation.core.progress import ProgressCallback
from constellation.massspec.io.thermo import BatchResult, convert_batch

DEFAULT_RT_BIN_WIDTH_S = 60.0
DEFAULT_BATCH_SIZE = 64


def mode_name(profile: bool) -> str:
    return "profile" if profile else "centroid"


def bundle_out_parent(data_root: "str | Path", dataset: str, *, profile: bool) -> Path:
    """Parent dir that holds one ``<stem>/`` bundle per ``.raw`` for this mode."""
    return Path(data_root) / "proc" / dataset / mode_name(profile)


def raw_dir(data_root: "str | Path", dataset: str) -> Path:
    return Path(data_root) / "raw" / dataset


def enumerate_raw(data_root: "str | Path", dataset: str) -> list[Path]:
    """Sorted, de-duplicated ``.raw`` files for a dataset (case-insensitive)."""
    rd = raw_dir(data_root, dataset)
    matches: dict[Path, Path] = {}
    for pattern in ("*.raw", "*.RAW"):
        for p in rd.glob(pattern):
            if p.is_file():
                matches.setdefault(p.resolve(), p)
    return sorted(matches.values(), key=lambda p: p.name)


def plan_conversion(
    paths: list[Path],
    *,
    shard: int | None = None,
    n_shards: int | None = None,
    limit: int | None = None,
) -> list[Path]:
    """Deterministically slice ``paths`` for SLURM-array sharding / smoke tests.

    Sharding takes ``sorted(paths)[shard::n_shards]`` (disjoint, balanced,
    stable across runs); ``limit`` then caps the count.
    """
    out = sorted(paths, key=lambda p: p.name)
    if n_shards is not None and n_shards > 1:
        if shard is None or not (0 <= shard < n_shards):
            raise ValueError(f"shard must be in [0,{n_shards}); got {shard}")
        out = out[shard::n_shards]
    if limit is not None:
        out = out[:limit]
    return out


def run_conversion(
    paths: list[Path],
    out_parent: "str | Path",
    *,
    profile: bool = False,
    n_workers: int = 1,
    force: bool = False,
    rt_bin_width_s: float = DEFAULT_RT_BIN_WIDTH_S,
    batch_size: int = DEFAULT_BATCH_SIZE,
    progress_cb: ProgressCallback | None = None,
) -> list[BatchResult]:
    """Convert ``paths`` into ``out_parent`` via ``convert_batch`` (pass-through).

    Returns ``[]`` for an empty input (``convert_batch`` would otherwise raise).
    """
    if not paths:
        return []
    return convert_batch(
        list(paths),
        Path(out_parent),
        n_workers=n_workers,
        force=force,
        profile=profile,
        rt_bin_width_s=rt_bin_width_s,
        batch_size=batch_size,
        progress_cb=progress_cb,
    )


@dataclass
class ConvertSummary:
    n_ok: int = 0
    n_skipped: int = 0
    n_error: int = 0
    errors: list[tuple[str, str]] = field(default_factory=list)  # (input_name, detail)

    @property
    def n_total(self) -> int:
        return self.n_ok + self.n_skipped + self.n_error


def summarize_results(results: list[BatchResult]) -> ConvertSummary:
    s = ConvertSummary()
    for r in results:
        if r.status == "ok":
            s.n_ok += 1
        elif r.status == "skipped":
            s.n_skipped += 1
        else:
            s.n_error += 1
            s.errors.append((r.input_path.name, r.detail or "?"))
    return s


def converted_stems(data_root: "str | Path", dataset: str, *, profile: bool) -> set[str]:
    """Stems with a completed bundle (``manifest.json`` present) for this mode."""
    out_parent = bundle_out_parent(data_root, dataset, profile=profile)
    if not out_parent.is_dir():
        return set()
    return {p.parent.name for p in out_parent.glob("*/manifest.json")}


__all__ = [
    "ConvertSummary",
    "bundle_out_parent",
    "converted_stems",
    "enumerate_raw",
    "mode_name",
    "plan_conversion",
    "raw_dir",
    "run_conversion",
    "summarize_results",
]
