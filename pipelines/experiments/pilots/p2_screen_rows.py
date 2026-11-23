#!/usr/bin/env python
"""Pilot P2, sweep half: per-scan isotope rows for the screen prototype (extraction design).

P2 asks which per-scan features and thresholds separate in-pool from off-pool signal, so the
screen reducer handed to Constellation can be specified from data. This half sweeps a block of
consecutive runs on one instrument and keeps, for every (query, MS1 scan) where at least two of
M0, M+1 and M+2 matched, the raw isotope intensities and m/z errors (M−1..M+3), plus the scan's
rt, injection time, TIC and floor. Features (KL, cosine, mass coherence), regions and
thresholds are computed offline by ``p2_analyze.py``, so they can change without sweeping again.

Queries per run, all at z = 1..6 inside the run's MS1 range:

- every peptide of every pool in the block (each is in-pool in its pool's runs and off-pool or
  carryover in the others);
- PROCAL and the other standards;
- a hashed 2% of the whole search space, as a general off-pool null.

Outputs, per run, under ``--out``: ``rows/<stem>.parquet`` (the isotope rows),
``channels/<stem>.parquet`` (each query's channel m/z and theoretical shares) and
``labels/<stem>.parquet`` (``protostar.pools.labels_for_run``), plus ``targets.parquet``
(target_id → sequence_id, charge) and a JSON report.

The reducer is a short-lived pilot prototype (ledger promotion entry: the screen reducer); it
emits measurements only and does no scoring.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _common  # noqa: E402,F401  (pins libstdc++ before constellation imports)

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from protostar import bundles, pools  # noqa: E402
from protostar import search_space as ss  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _pilot_reducers import ISOTOPES, IsotopeScanRows  # noqa: E402

CENTER = "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2"


def _hashed(ids: np.ndarray, frac: float) -> np.ndarray:
    import hashlib

    h = np.array(
        [
            int(hashlib.blake2b(int(i).to_bytes(8, "little"), digest_size=4).hexdigest(), 16)
            for i in ids
        ]
    )
    return ids[(h % 10_000) < int(frac * 10_000)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--space-cache", type=Path, default=None)
    ap.add_argument("--center-stem", default=CENTER)
    ap.add_argument("--block-half", type=int, default=10)
    ap.add_argument("--null-fraction", type=float, default=0.02)
    ap.add_argument("--layout", default="precursor_major")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--tol-ppm", type=float, default=20.0)
    ap.add_argument(
        "--stems", nargs="*", default=None, help="restrict to these block runs (smoke tests)"
    )
    args = ap.parse_args()

    from constellation.massspec.sweep import ChannelSettings, compile_queries
    from constellation.massspec.sweep.engine import SweepSettings, sweep
    from constellation.massspec.sweep.rules import rule_named

    for sub in ("rows", "channels", "labels"):
        (args.out / sub).mkdir(parents=True, exist_ok=True)
    report: dict = {"provenance": bundles.constellation_provenance(),
                    "args": {k: str(v) for k, v in vars(args).items()}, "runs": {}}  # fmt: skip

    def save():
        (args.out / "p2_report.json").write_text(json.dumps(report, indent=1, default=str))

    space_path = args.space_cache or args.out / "search_space.parquet"
    space = pq.read_table(space_path) if space_path.exists() else ss.build_search_space()
    index = pools.PoolIndex.from_references(space)
    runs = pools.run_table(args.data_root)
    rows = {r["stem"]: r for r in runs.to_pylist()}
    me = rows[args.center_stem]
    block = sorted(
        (s for s, r in rows.items()
         if r["instrument_serial"] == me["instrument_serial"] and r["instrument_rank"] is not None
         and abs(r["instrument_rank"] - me["instrument_rank"]) <= args.block_half),
        key=lambda s: rows[s]["instrument_rank"],
    )  # fmt: skip
    if args.stems:
        block = [s for s in block if s in set(args.stems)]
    block_pools = sorted({rows[s]["pool"] for s in block})
    seq_ids = np.unique(np.concatenate([index.pool(p) for p in block_pools]))
    standards = np.flatnonzero(space.column("is_procal").to_numpy(zero_copy_only=False))
    qc = np.concatenate(list(index.qc.values())) if index.qc else np.array([], np.int64)
    null = _hashed(space.column("sequence_id").to_numpy(), args.null_fraction)
    query_ids = np.unique(np.concatenate([seq_ids, standards, qc, null]))
    report["block"] = [{"stem": s, "rank": rows[s]["instrument_rank"], "pool": rows[s]["pool"],
                        "injection": rows[s]["injection"]} for s in block]  # fmt: skip
    report["queries"] = {"block_pool_sequences": int(seq_ids.size), "null_sequences": int(null.size),
                         "total_sequences": int(query_ids.size)}  # fmt: skip
    save()
    channels = ChannelSettings(isotope_offsets=ISOTOPES)
    settings = SweepSettings(tolerances={1: (args.tol_ppm, "ppm")}, rule=rule_named("nearest"))
    spill_root = Path(os.environ.get("TMPDIR", tempfile.gettempdir()))
    all_targets = []
    for stem in block:
        t0 = time.monotonic()
        bundle = bundles.require_run_bundle(args.data_root, rows[stem]["dataset"], stem)
        if bundle is None:
            report["runs"][stem] = {"error": "no current bundle"}
            save()
            continue
        sm = bundles.read_scan_metadata(
            bundle, columns=["level", "scan_window_lower", "scan_window_upper"]
        )
        targets = ss.targets_for_run(space, mz_range=ss.ms1_mz_range(sm), subset=query_ids)
        all_targets.append(targets.select(["target_id", "peptide_id", "precursor_charge"]))
        plan = compile_queries(targets, levels=[1], channels=channels, acquired_ranges=None)
        spill = Path(tempfile.mkdtemp(prefix="p2spill-", dir=spill_root))
        try:
            res = sweep(bundle, plan, settings, reducers=[IsotopeScanRows()], workers=args.workers,
                        layout=args.layout, spill_dir=spill, work_dir=spill)  # fmt: skip
            scans_t = res.table("isotope_scans")
            chan_t = res.table("query_channels")
        finally:
            import shutil

            shutil.rmtree(spill, ignore_errors=True)
        pq.write_table(scans_t, args.out / "rows" / f"{stem}.parquet", compression="zstd")
        pq.write_table(chan_t, args.out / "channels" / f"{stem}.parquet", compression="zstd")
        labels = pools.labels_for_run(stem, runs, space, index)
        pq.write_table(labels, args.out / "labels" / f"{stem}.parquet", compression="zstd")
        report["runs"][stem] = {
            "targets": targets.num_rows,
            "rows": scans_t.num_rows,
            "rows_bytes": (args.out / "rows" / f"{stem}.parquet").stat().st_size,
            "labels": pools.label_counts(labels),
            "wall_s": round(time.monotonic() - t0, 1),
        }
        print(stem, json.dumps(report["runs"][stem]), flush=True)
        save()
    if all_targets:
        pq.write_table(pa.concat_tables(all_targets).group_by(["target_id", "peptide_id", "precursor_charge"])
                       .aggregate([]), args.out / "targets.parquet")  # fmt: skip
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
