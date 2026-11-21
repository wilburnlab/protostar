"""Drive Constellation's Thermo ``.raw`` → parquet-bundle conversion (bundle v3).

Thin orchestration over ``constellation.massspec.readers.thermo.convert_batch``.
protostar owns:
- file selection: which ``.raw`` files to convert, and shard/limit for SLURM arrays;
- result reconciliation.

The conversion itself — spawn-mode workers, centroiding, the v3 bundle, the
source SHA-256 — lives in Constellation.

**Layout (v3, one pass).** Each ``.raw`` becomes a centroid bundle plus a sibling
profile bundle::

    <data_root>/proc/<dataset>/<stem>/            # peaks/part-*.parquet, scan/acq metadata
    <data_root>/proc/<dataset>/<stem>.profile/    # MS1 profile points + MS2 stored centroids

The June 2026 pre-v3 bundles under ``proc/<dataset>/{centroid,profile}/`` are
never written to here. They coexist untouched until their dataset's v3
reconversion verifies (see ``protostar.bundles``).

**Resume.**
- Constellation skips any output that already has a ``manifest.json``.
- It would also skip an *outdated* (pre-v3) manifest at the v3 location. That is
  constellation CR-A4, in Constellation's
  ``docs/plans/counter-pt-reconciliation.md``.
- So outdated stems are routed to a separate ``force=True`` call. A pre-v3 run is
  never reported as "already converted".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from constellation.core.progress import ProgressCallback
from constellation.massspec.readers.thermo import (
    DEFAULT_CENTROID_ALGO,
    BatchResult,
    convert_batch,
)

from protostar import bundles

DEFAULT_RT_BIN_WIDTH_S = 60.0
DEFAULT_BATCH_SIZE = 64
#: Constellation's default: also write a precursor-major copy of the peaks (#109),
#: which the sweep's precursor-major layout reads. It costs disk; the reconversion
#: pilot measures it (restart plan, Phase 4) and may change this default.
DEFAULT_PRECURSOR_MAJOR = True


def bundle_out_parent(data_root: str | Path, dataset: str) -> Path:
    """The ``output_parent`` handed to ``convert_batch``: ``proc/<dataset>``."""
    return bundles.proc_dir(data_root, dataset)


def raw_dir(data_root: str | Path, dataset: str) -> Path:
    return Path(data_root) / "raw" / dataset


def enumerate_raw(data_root: str | Path, dataset: str) -> list[Path]:
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


@dataclass
class ConversionPlan:
    """``paths`` split by the state of their v3 output bundle."""

    current: list[Path] = field(default_factory=list)  # skip (already v3)
    pending: list[Path] = field(default_factory=list)  # missing or incomplete
    outdated: list[Path] = field(default_factory=list)  # pre-v3 at the v3 location: force

    @property
    def to_convert(self) -> list[Path]:
        return self.pending + self.outdated


def classify(paths: list[Path], out_parent: str | Path) -> ConversionPlan:
    """Split ``paths`` by the state of ``out_parent/<stem>``."""
    plan = ConversionPlan()
    for p in paths:
        state = bundles.bundle_state(Path(out_parent) / p.stem)
        if state == "current":
            plan.current.append(p)
        elif state == "outdated":
            plan.outdated.append(p)
        else:
            plan.pending.append(p)
    return plan


def run_conversion(
    paths: list[Path],
    out_parent: str | Path,
    *,
    keep_profile: bool = True,
    centroid_algo: str = DEFAULT_CENTROID_ALGO,
    precursor_major: bool = DEFAULT_PRECURSOR_MAJOR,
    n_workers: int = 1,
    force: bool = False,
    rt_bin_width_s: float = DEFAULT_RT_BIN_WIDTH_S,
    batch_size: int = DEFAULT_BATCH_SIZE,
    progress_cb: ProgressCallback | None = None,
) -> list[BatchResult]:
    """Convert ``paths`` into ``out_parent`` (v3, centroid + profile in one pass).

    Without ``force``:
    - current bundles are left to Constellation's own skip;
    - outdated (pre-v3) bundles at the v3 location are re-converted with
      ``force=True``.

    Returns ``[]`` for an empty input (``convert_batch`` would otherwise raise).
    """
    if not paths:
        return []
    common = dict(
        n_workers=n_workers,
        rt_bin_width_s=rt_bin_width_s,
        keep_profile=keep_profile,
        centroid_algo=centroid_algo,
        precursor_major=precursor_major,
        batch_size=batch_size,
        progress_cb=progress_cb,
    )
    out = Path(out_parent)
    if force:
        return convert_batch(list(paths), out, force=True, **common)
    plan = classify(list(paths), out)
    results: list[BatchResult] = []
    normal = plan.current + plan.pending
    if normal:
        results += convert_batch(normal, out, force=False, **common)
    if plan.outdated:
        results += convert_batch(plan.outdated, out, force=True, **common)
    return results


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


def converted_stems(data_root: str | Path, dataset: str) -> set[str]:
    """Stems with a **current** (v3) centroid bundle."""
    return {b.name for b in bundles.iter_run_bundles(data_root, dataset)}


__all__ = [
    "ConversionPlan",
    "ConvertSummary",
    "DEFAULT_PRECURSOR_MAJOR",
    "bundle_out_parent",
    "classify",
    "converted_stems",
    "enumerate_raw",
    "plan_conversion",
    "raw_dir",
    "run_conversion",
    "summarize_results",
]
