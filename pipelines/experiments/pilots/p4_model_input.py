#!/usr/bin/env python
"""Pilot P4: what the Tier-2 model-input record costs and covers (extraction design).

Tier 2 stores, for a bounded set of (run, sequence), every eligible MS1 scan, empty ones included,
plus the MS2 scans whose isolation covers the envelope. This pilot extracts that record with
today's API, ``extract_xic(levels=[1, 2], rule="all", drop_unmatched=False)``, for one pool's
peptides (z = 1..6, isotopes M−1..M+3, fragments b/y, or c/z for ETD, up to min(3, z−1)) in:

- the pool's own runs (whole gradient);
- off-pool runs from the same instrument block, with windows at the peptide's in-pool PSM RT
  ± ``--half-window`` s.

MS2 tolerance is per MS level in Constellation, but ion-trap MS2 needs about 0.5 Da and Orbitrap MS2
about 20 ppm. So each run is swept twice, and each MS2 row is kept from the sweep that matches its
scan's analyzer.

Measured per run:
- rows and parquet bytes per level, and per (sequence, charge);
- MS2 scans per target;
- how many MS2 scans have their MS1 scan before and after inside the record (the brackets the
  MS1→MS2 model needs);
- the share of each target's isotope envelope inside its MS2 scans' isolation windows.

Outputs: ``p4_report.json`` and ``traces/<stem>.parquet`` (the combined record, for inspection).
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
import _common  # noqa: E402,F401

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from protostar import bundles, pools  # noqa: E402
from protostar import search_space as ss  # noqa: E402
from protostar.experiments.scans import gate_psms, load_psms  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p3_rt_transfer import best_rt  # noqa: E402

CENTER = "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2"
ISOTOPES = (-1, 0, 1, 2, 3)


def _extract(bundle, targets, ion_types, ms2_tol, workers, spill):
    from constellation.massspec.peptide.ions import IonType
    from constellation.massspec.quant.chromatogram import extract_xic
    from constellation.massspec.sweep import ChannelSettings

    channels = ChannelSettings(
        isotope_offsets=ISOTOPES,
        ion_types=tuple(IonType[t] for t in ion_types),
        max_fragment_charge=3,
        fragment_charge_below_precursor=1,
    )
    t0 = time.monotonic()
    res = extract_xic(bundle, targets, acquisition_id=0, levels=[1, 2],
                      tolerances={1: (20.0, "ppm"), 2: ms2_tol}, channels=channels, rule="all",
                      drop_unmatched=False, workers=workers, work_dir=spill, spill_dir=spill,
                      layout="precursor_major")  # fmt: skip
    path = spill / f"trace-{ms2_tol[0]}{ms2_tol[1]}.parquet"
    res.write("xic_trace", path)
    return pq.read_table(path), time.monotonic() - t0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--space-cache", type=Path, required=True)
    ap.add_argument("--center-stem", default=CENTER)
    ap.add_argument("--block-half", type=int, default=10)
    ap.add_argument("--n-off-pool", type=int, default=8)
    ap.add_argument("--half-window", type=float, default=360.0)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--max-pep", type=float, default=0.01)
    ap.add_argument(
        "--limit-sequences", type=int, default=None, help="smoke tests: the first N sequences"
    )
    ap.add_argument(
        "--limit-in-pool-runs",
        type=int,
        default=None,
        help="smoke tests: the first N own-pool runs",
    )
    args = ap.parse_args()
    (args.out / "traces").mkdir(parents=True, exist_ok=True)
    report: dict = {"provenance": bundles.constellation_provenance(),
                    "args": {k: str(v) for k, v in vars(args).items()}, "runs": {}}  # fmt: skip

    def save():
        (args.out / "p4_report.json").write_text(json.dumps(report, indent=1, default=str))

    space = pq.read_table(args.space_cache)
    index = pools.PoolIndex.from_references(space)
    runs = pools.run_table(args.data_root)
    rows = {r["stem"]: r for r in runs.to_pylist()}
    me = rows[args.center_stem]
    pool = me["pool"]
    in_pool = sorted(s for s, r in rows.items() if r["pool"] == pool and not r["noincl"])
    block = sorted(
        (s for s, r in rows.items()
         if r["instrument_serial"] == me["instrument_serial"] and r["instrument_rank"] is not None
         and r["pool"] != pool and abs(r["instrument_rank"] - me["instrument_rank"]) <= args.block_half),
        key=lambda s: abs(rows[s]["instrument_rank"] - me["instrument_rank"]),
    )[: args.n_off_pool]  # fmt: skip
    ids = np.union1d(
        index.pool(pool), np.flatnonzero(space.column("is_procal").to_numpy(zero_copy_only=False))
    )
    if args.limit_sequences:
        ids = ids[: args.limit_sequences]
    if args.limit_in_pool_runs:
        in_pool = in_pool[: args.limit_in_pool_runs]
    residues = space.column("residues").to_numpy(zero_copy_only=False)
    # each peptide's in-pool RT: median over the pool's runs of its best PSM
    tmp = Path(tempfile.mkdtemp(prefix="p4-", dir=os.environ.get("TMPDIR")))
    rts: dict[str, list[float]] = {}
    for s in in_pool:
        z = pools.search_zip_for_stem(args.data_root, rows[s]["dataset"], s)
        if z is not None:
            for q, rt in best_rt(
                gate_psms(load_psms(z, extract_dir=tmp / s), max_pep=args.max_pep), s
            ).items():
                rts.setdefault(q, []).append(rt)
    ref_rt = {i: float(np.median(rts[residues[i]])) for i in ids if residues[i] in rts}
    report["selection"] = {"pool": pool, "in_pool_runs": in_pool, "off_pool_runs": block,
                           "sequences": int(ids.size), "with_in_pool_rt": len(ref_rt)}  # fmt: skip
    save()

    for stem in in_pool + block:
        bundle = bundles.require_run_bundle(args.data_root, rows[stem]["dataset"], stem)
        if bundle is None:
            continue
        sm = bundles.read_scan_metadata(
            bundle, columns=["scan", "level", "rt", "analyzer", "isolation_lower", "isolation_upper",
                             "scan_window_lower", "scan_window_upper"],
        ).to_pandas()  # fmt: skip
        is_in = stem in in_pool
        sel = ids if is_in else np.array(sorted(ref_rt), dtype=np.int64)
        targets = ss.targets_for_run(
            space, mz_range=ss.ms1_mz_range(pa.Table.from_pandas(sm)), subset=sel
        )
        if not is_in:
            centre = np.array([ref_rt[int(p)] for p in targets.column("peptide_id").to_pylist()])
            targets = targets.set_column(targets.schema.get_field_index("rt_start"), "rt_start",
                                         pa.array(centre - args.half_window))  # fmt: skip
            targets = targets.set_column(targets.schema.get_field_index("rt_end"), "rt_end",
                                         pa.array(centre + args.half_window))  # fmt: skip
        ion_types = ("C", "Z") if rows[stem]["injection"] == "ETD" else ("B", "Y")
        spill = Path(tempfile.mkdtemp(prefix="p4sweep-", dir=os.environ.get("TMPDIR")))
        rec: dict = {
            "in_pool": is_in,
            "injection": rows[stem]["injection"],
            "targets": targets.num_rows,
        }
        try:
            ft, wall_ft = _extract(bundle, targets, ion_types, (20.0, "ppm"), args.workers, spill)
            it, wall_it = _extract(bundle, targets, ion_types, (0.5, "Da"), args.workers, spill)
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {exc}"
            report["runs"][stem] = rec
            save()
            continue
        analyzer = dict(zip(sm["scan"], sm["analyzer"].fillna("")))
        ana = lambda t: pa.array([analyzer.get(s, "") for s in t.column("scan").to_pylist()])  # noqa: E731
        lv1 = ft.filter(pc.equal(ft.column("level"), 1))
        ms2_ft = ft.filter(pc.equal(ft.column("level"), 2))
        ms2_it = it.filter(pc.equal(it.column("level"), 2))
        ms2 = pa.concat_tables([
            ms2_ft.filter(pc.invert(pc.match_substring(ana(ms2_ft), "ITMS"))),
            ms2_it.filter(pc.match_substring(ana(ms2_it), "ITMS")),
        ])  # fmt: skip
        combined = pa.concat_tables([lv1, ms2])
        path = args.out / "traces" / f"{stem}.parquet"
        pq.write_table(combined, path, compression="zstd")
        n_seq_z = pc.count_distinct(combined.column("target_id")).as_py()
        # MS2 coverage: brackets and envelope transmission
        ms1_scans = np.sort(sm.loc[sm["level"] == 1, "scan"].to_numpy())
        lv1_by_t: dict[int, set] = {}
        for t_, s_ in zip(
            lv1.column("target_id").to_pylist(), lv1.column("scan").to_pylist(), strict=True
        ):
            lv1_by_t.setdefault(t_, set()).add(s_)
        ms2_keys = set(
            zip(ms2.column("target_id").to_pylist(), ms2.column("scan").to_pylist(), strict=True)
        )
        bracketed = 0
        for t_, s_ in ms2_keys:
            k = np.searchsorted(ms1_scans, s_)
            before = ms1_scans[k - 1] if k > 0 else None
            after = ms1_scans[k] if k < ms1_scans.size else None
            have = lv1_by_t.get(t_, set())
            bracketed += int(before in have and after in have)
        mz_of = dict(
            zip(
                targets.column("target_id").to_pylist(),
                targets.column("precursor_mz").to_pylist(),
                strict=True,
            )
        )
        z_of = dict(
            zip(
                targets.column("target_id").to_pylist(),
                targets.column("precursor_charge").to_pylist(),
                strict=True,
            )
        )
        iso_lo = dict(zip(sm["scan"], sm["isolation_lower"]))
        iso_hi = dict(zip(sm["scan"], sm["isolation_upper"]))
        frac_in = []
        for t_, s_ in ms2_keys:
            lo, hi = iso_lo.get(s_), iso_hi.get(s_)
            if lo is None or hi is None or not np.isfinite(lo) or not np.isfinite(hi):
                continue
            mz0, z = mz_of[t_], z_of[t_]
            iso_mz = mz0 + np.arange(0, 4) * ss.C13_SPACING / z
            w = np.array(
                [0.5, 0.3, 0.15, 0.05]
            )  # rough envelope for a ~1.5 kDa peptide; sizing only
            frac_in.append(float(w[(iso_mz >= lo) & (iso_mz <= hi)].sum()))
        frac_in = np.array(frac_in)
        rec.update(
            rows_ms1=lv1.num_rows,
            rows_ms2=ms2.num_rows,
            bytes=path.stat().st_size,
            targets_with_rows=int(n_seq_z),
            bytes_per_target=round(path.stat().st_size / max(n_seq_z, 1), 1),
            ms2_scans_per_target=round(len(ms2_keys) / max(n_seq_z, 1), 2),
            ms2_rows_itms=int(ms2_it.num_rows),
            bracket_coverage=round(bracketed / max(len(ms2_keys), 1), 4),
            envelope_in_window_median=None
            if frac_in.size == 0
            else round(float(np.median(frac_in)), 3),
            envelope_in_window_lt_half=None
            if frac_in.size == 0
            else round(float((frac_in < 0.5).mean()), 4),
            wall_s=round(wall_ft + wall_it, 1),
        )
        import shutil

        shutil.rmtree(spill, ignore_errors=True)
        report["runs"][stem] = rec
        print(stem, json.dumps(rec), flush=True)
        save()
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
