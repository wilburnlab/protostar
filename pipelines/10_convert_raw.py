#!/usr/bin/env python
"""Stage 10 — convert ``.raw`` into Constellation parquet bundles (``proc/``).

Drives ``constellation.massspec.io.thermo.convert_batch`` over a dataset's
``.raw`` files. Each file becomes a bundle directory
(``manifest.json`` + ``peaks.parquet`` + ``scan_metadata.parquet`` +
``acquisition_metadata.parquet``); per-scan ``filter_string`` is recorded so
fragmentation mode is recoverable downstream without pre-splitting. Resume is
free — bundles whose ``manifest.json`` already exists are skipped.

Centroid (default) and profile bundles live side by side::

    <data_root>/proc/<dataset>/<centroid|profile>/<stem>/...

Renamed from the scaffold's ``10_build_mzpeak.py`` — the output is just
Constellation's convert bundle (no "mzpeak"/HUPO-PSI standard exists).

Usage::

    python pipelines/10_convert_raw.py --dataset Zolg2017 --limit 1     # smoke test
    python pipelines/10_convert_raw.py --dataset Zolg2017 --limit 2 --profile  # validate profile
    python pipelines/10_convert_raw.py --dataset all --workers 40       # full centroid pass
"""

from __future__ import annotations

import argparse

from _common import (
    add_common_args,
    add_dataset_arg,
    data_root,
    load_config,
    make_progress,
    resolve_datasets,
)

from protostar.convert import driver


def _convert_dataset(args: argparse.Namespace, config: dict, dataset: str) -> int:
    droot = data_root(config, args.data_root)
    mode = driver.mode_name(args.profile)

    raws = driver.enumerate_raw(droot, dataset)
    if not raws:
        print(f"  {dataset}: no .raw files under {driver.raw_dir(droot, dataset)} — fetch first")
        return 0

    planned = driver.plan_conversion(
        raws, shard=args.shard, n_shards=args.n_shards, limit=args.limit
    )
    out_parent = driver.bundle_out_parent(droot, dataset, profile=args.profile)
    already = driver.converted_stems(droot, dataset, profile=args.profile)
    todo = [p for p in planned if args.force or p.stem not in already]

    if args.dry_run:
        print(
            f"  {dataset} [{mode}]: {len(raws)} raw, {len(planned)} in plan, "
            f"{len(already)} already converted, {len(todo)} to convert"
        )
        return 0

    print(f"  {dataset} [{mode}]: converting {len(todo)} of {len(planned)} planned → {out_parent}")
    results = driver.run_conversion(
        planned,
        out_parent,
        profile=args.profile,
        n_workers=args.workers,
        force=args.force,
        rt_bin_width_s=args.rt_bin_width_s,
        batch_size=args.batch_size,
        progress_cb=make_progress(args.quiet),
    )
    summary = driver.summarize_results(results)
    print(
        f"  {dataset} [{mode}]: {summary.n_ok} ok, {summary.n_skipped} skipped, "
        f"{summary.n_error} failed"
    )
    for name, detail in summary.errors:
        print(f"    ERROR {name}: {detail}")
    return 1 if summary.n_error else 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    add_common_args(p)
    add_dataset_arg(p)
    p.add_argument(
        "--profile", action="store_true", help="profile mode (raw FT grid); default is centroid"
    )
    p.add_argument(
        "--workers", type=int, default=1, help="parallel converter workers (spawn-mode; default 1)"
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="re-convert bundles that already exist (this mode only)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="report how many would convert vs skip; change nothing",
    )
    p.add_argument(
        "--rt-bin-width-s",
        type=float,
        default=driver.DEFAULT_RT_BIN_WIDTH_S,
        help="peak-table row-group RT chunking (s)",
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=driver.DEFAULT_BATCH_SIZE,
        help="scans per internal ScanBatch",
    )
    p.add_argument(
        "--limit", type=int, default=None, help="cap the number of files converted (smoke tests)"
    )
    p.add_argument("--shard", type=int, default=0, help="SLURM-array shard index (with --n-shards)")
    p.add_argument("--n-shards", type=int, default=None, help="number of SLURM-array shards")
    return p


def main(argv: "list[str] | None" = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)

    if not args.dry_run:
        # CLR pre-flight in the parent (pure file check; convert_batch loads the
        # CLR only inside spawn workers — never fork a CLR-initialised parent).
        from constellation.massspec.io.thermo import require_thermo

        try:
            require_thermo()
        except ImportError as exc:
            print(f"error: Thermo reader unavailable: {exc}")
            return 3

    rc = 0
    for dataset in resolve_datasets(args.dataset, config):
        rc |= _convert_dataset(args, config, dataset)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
