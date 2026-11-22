#!/usr/bin/env python
"""Pilot P1: what a global MS1 sweep costs (extraction design, Step 2).

For the runs given (default: Zolg2017 pool 1's four injections), sweep MS1 for growing query
sets, from the run's own pool up to every sequence in the search space at z = 1..6. Each sweep
uses Constellation's ``QuerySummaryReducer``, which writes one row per query and keeps no
trace: this is the cost of a matches-only screen. A second, small set of sweeps uses a no-op
reducer that builds the dense (zero-filled) view, to bound what a dense screen would cost.

Measured per sweep: query compile time and size, sweep wall time, peak memory (parent and
workers), matches, and spill bytes. The JSON report is rewritten after every sweep, so a job
that hits its wall time keeps what it finished.

No model or scoring code: Constellation's compile, sweep and reducers, timed.

Usage (OSC, inside an allocation with ``--cpus-per-task`` ≥ the largest ``--workers``)::

    python pipelines/experiments/pilots/p1_screen_cost.py --data-root <data_root> \
        --out <data_root>/pilots/extraction/p1 --workers 4 8 16
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import shutil
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _common  # noqa: E402,F401  (pins libstdc++ before constellation imports)

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from protostar import bundles, pools  # noqa: E402
from protostar import search_space as ss  # noqa: E402

POOL1 = [
    "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2",  # the nb42b reference file
    "01625b_GA1-TUM_first_pool_1_01_01-2xIT_2xHCD-1h-R1",
    "01625b_GA1-TUM_first_pool_1_01_01-3xHCD-1h-R1",
    "01625b_GA1-TUM_first_pool_1_01_01-ETD-1h-R2",
]
ISOTOPES = (-1, 0, 1, 2, 3)


class DenseCost:
    """A no-op reducer that builds the dense MS1 view chunk by chunk and keeps only its size."""

    name = "dense_cost"
    schema = pa.schema([pa.field("cells", pa.int64()), pa.field("acquired", pa.int64())])
    gated = False
    zero_fills = True

    def __init__(self, max_cells: int = 1 << 24):
        self.max_cells = max_cells

    def reduce(self, batch) -> pa.Table:
        cells = acquired = 0
        for chunk in batch.dense_chunks(self.max_cells, levels=[1]):
            d = batch.dense(1, queries=chunk)
            cells += int(np.prod(d.shape))
            acquired += int(d.acquired.sum())
        return pa.table({"cells": [cells], "acquired": [acquired]}, schema=self.schema)


def _peak_rss_gb() -> dict[str, float]:
    self_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    child_kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return {"parent_gb": round(self_kb / 2**20, 2), "largest_worker_gb": round(child_kb / 2**20, 2)}


def _du(path: Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.exists() else 0


def _hashed_fraction(ids: np.ndarray, frac: float) -> np.ndarray:
    """A deterministic pseudo-random fraction of ``ids`` (the same ids in every run)."""
    if frac >= 1.0:
        return ids
    keep = np.array(
        [
            int(hashlib.blake2b(int(i).to_bytes(8, "little"), digest_size=4).hexdigest(), 16)
            % 10_000
            for i in ids
        ]
    ) < int(frac * 10_000)
    return ids[keep]


def _load_space(cache: Path) -> pa.Table:
    if cache.exists():
        return pq.read_table(cache)
    space = ss.build_search_space()
    cache.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(space, cache)
    return space


def _one_sweep(bundle, plan, settings, reducer, *, workers, layout, spill_root):
    from constellation.massspec.sweep.engine import sweep

    spill = Path(tempfile.mkdtemp(prefix="p1spill-", dir=spill_root))
    t0 = time.monotonic()
    try:
        res = sweep(bundle, plan, settings, reducers=[reducer], workers=workers, layout=layout,
                    spill_dir=spill, work_dir=spill if workers > 1 else None)  # fmt: skip
        wall = time.monotonic() - t0
        table = res.table(reducer.name)
        spilled = _du(spill)
    finally:
        shutil.rmtree(spill, ignore_errors=True)
    return wall, table, spilled


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--dataset", default="Zolg2017")
    ap.add_argument("--stems", nargs="+", default=POOL1)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument(
        "--space-cache", type=Path, default=None, help="search-space parquet (built once)"
    )
    ap.add_argument("--fractions", nargs="+", type=float, default=[0.0, 0.1, 0.3, 1.0],
                    help="0 = the run's pool and standards only; else that fraction of all sequences")  # fmt: skip
    ap.add_argument("--layouts", nargs="+", default=["scan_major", "precursor_major"])
    ap.add_argument("--workers", nargs="+", type=int, default=[16])
    ap.add_argument("--full-grid-stem", default=POOL1[0],
                    help="the run that gets every fraction x layout x workers; the others get fraction 1.0 only")  # fmt: skip
    ap.add_argument("--dense-fractions", nargs="*", type=float, default=[0.0, 0.01])
    ap.add_argument("--tol-ppm", type=float, default=20.0)
    args = ap.parse_args()

    from constellation.massspec.sweep import ChannelSettings, compile_queries
    from constellation.massspec.sweep.engine import SweepSettings
    from constellation.massspec.sweep.reduce import QuerySummaryReducer
    from constellation.massspec.sweep.rules import rule_named

    args.out.mkdir(parents=True, exist_ok=True)
    spill_root = Path(os.environ.get("TMPDIR", tempfile.gettempdir()))
    report: dict = {
        "provenance": bundles.constellation_provenance(),
        "args": {k: str(v) for k, v in vars(args).items()},
        "host": os.uname().nodename,
        "cpus": len(os.sched_getaffinity(0)),
        "runs": {},
    }
    t0 = time.monotonic()
    space = _load_space(args.space_cache or args.out / "search_space.parquet")
    index = pools.PoolIndex.from_references(space)
    report["space"] = {"sequences": space.num_rows, "load_s": round(time.monotonic() - t0, 1)}
    all_ids = space.column("sequence_id").to_numpy()
    standards = np.flatnonzero(space.column("is_procal").to_numpy(zero_copy_only=False))
    channels = ChannelSettings(isotope_offsets=ISOTOPES)
    settings = SweepSettings(tolerances={1: (args.tol_ppm, "ppm")}, rule=rule_named("nearest"))

    def save():
        (args.out / "p1_report.json").write_text(json.dumps(report, indent=1, default=str))

    for stem in args.stems:
        bundle = bundles.require_run_bundle(args.data_root, args.dataset, stem)
        if bundle is None:
            report["runs"][stem] = {"error": "no current bundle"}
            save()
            continue
        sm = bundles.read_scan_metadata(
            bundle, columns=["level", "scan_window_lower", "scan_window_upper"]
        )
        mz_range = ss.ms1_mz_range(sm)
        own = np.union1d(index.pool(pools.parse_stem(stem).pool), standards)
        rec = report["runs"].setdefault(stem, {"mz_range": mz_range,
                                               "n_ms1_scans": int(pc.sum(pc.equal(sm.column("level"), 1)).as_py()),
                                               "sweeps": []})  # fmt: skip
        full = stem == args.full_grid_stem
        fractions = args.fractions if full else [1.0]
        for frac in fractions:
            subset = own if frac == 0.0 else np.union1d(_hashed_fraction(all_ids, frac), own)
            targets = ss.targets_for_run(space, mz_range=mz_range, subset=subset)
            tc = time.monotonic()
            plan = compile_queries(targets, levels=[1], channels=channels, acquired_ranges=None)
            compile_s = time.monotonic() - tc
            for layout in args.layouts if full else args.layouts[:1]:
                for w in args.workers if full else [max(args.workers)]:
                    entry = {"fraction": frac, "layout": layout, "workers": w, "reducer": "query_summary",
                             "sequences": int(subset.size), "targets": targets.num_rows,
                             "compile_s": round(compile_s, 1)}  # fmt: skip
                    try:
                        wall, table, spilled = _one_sweep(bundle, plan, settings, QuerySummaryReducer(),
                                                          workers=w, layout=layout, spill_root=spill_root)  # fmt: skip
                        entry.update(wall_s=round(wall, 1), rows=table.num_rows,
                                     matched_ms1=int(pc.sum(table.column("n_matched_ms1")).as_py() or 0),
                                     eligible_ms1=int(pc.sum(table.column("n_eligible_ms1")).as_py() or 0),
                                     spill_bytes=spilled)  # fmt: skip
                    except Exception as exc:  # record and carry on: a pilot measures limits
                        entry["error"] = f"{type(exc).__name__}: {exc}"
                    entry.update(_peak_rss_gb())
                    rec["sweeps"].append(entry)
                    print(json.dumps(entry), flush=True)
                    save()
        if full:
            for frac in args.dense_fractions:
                subset = own if frac == 0.0 else np.union1d(_hashed_fraction(all_ids, frac), own)
                targets = ss.targets_for_run(space, mz_range=mz_range, subset=subset)
                plan = compile_queries(targets, levels=[1], channels=channels, acquired_ranges=None)
                w = max(args.workers)
                entry = {"fraction": frac, "layout": "scan_major", "workers": w, "reducer": "dense_cost",
                         "sequences": int(subset.size), "targets": targets.num_rows}  # fmt: skip
                try:
                    wall, table, spilled = _one_sweep(bundle, plan, settings, DenseCost(), workers=w,
                                                      layout="scan_major", spill_root=spill_root)  # fmt: skip
                    entry.update(wall_s=round(wall, 1), cells=int(pc.sum(table.column("cells")).as_py() or 0),
                                 acquired=int(pc.sum(table.column("acquired")).as_py() or 0), spill_bytes=spilled)  # fmt: skip
                except Exception as exc:
                    entry["error"] = f"{type(exc).__name__}: {exc}"
                entry.update(_peak_rss_gb())
                rec["sweeps"].append(entry)
                print(json.dumps(entry), flush=True)
                save()
    report["total_s"] = round(time.monotonic() - t0, 1)
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
