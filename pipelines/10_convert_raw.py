#!/usr/bin/env python
"""Stage 10 — convert ``.raw`` into Constellation parquet bundles (``proc/``), bundle v3.

Drives ``constellation.massspec.readers.thermo.convert_batch`` over a dataset's
``.raw`` files. Each file becomes, in **one pass**, a centroid bundle plus a
sibling profile bundle::

    <data_root>/proc/<dataset>/<stem>/            # peaks/part-*.parquet + scan/acq metadata
    <data_root>/proc/<dataset>/<stem>.profile/    # MS1 profile points (--no-profile to omit)

Per-scan ``filter_string`` is recorded, so fragmentation mode is recoverable
downstream without pre-splitting files.

Resume is free: a current (v3) bundle is skipped. A pre-v3 bundle at the v3
location is re-converted rather than skipped. The June 2026 pre-v3 trees
(``proc/<dataset>/{centroid,profile}/``) are left untouched; ``--dry-run``
reports them.

Usage::

    python pipelines/10_convert_raw.py --dataset Zolg2017 --dry-run      # inventory + plan
    python pipelines/10_convert_raw.py --dataset Zolg2017 --limit 1      # smoke test
    python pipelines/10_convert_raw.py --dataset all --workers 16        # full pass
    python pipelines/10_convert_raw.py --raw-name <stem> --dataset Zolg2017   # one named file
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

from protostar import bundles
from protostar.convert import driver


def _convert_dataset(args: argparse.Namespace, config: dict, dataset: str) -> int:
    droot = data_root(config, args.data_root)

    raws = driver.enumerate_raw(droot, dataset)
    if args.raw_name:
        wanted = set(args.raw_name)
        raws = [p for p in raws if p.stem in wanted or p.name in wanted]
    if not raws:
        print(f"  {dataset}: no matching .raw under {driver.raw_dir(droot, dataset)}")
        return 0

    planned = driver.plan_conversion(
        raws, shard=args.shard, n_shards=args.n_shards, limit=args.limit
    )
    out_parent = driver.bundle_out_parent(droot, dataset)
    plan = driver.classify(planned, out_parent)
    n_todo = len(planned) if args.force else len(plan.to_convert)

    if args.dry_run:
        inv = bundles.inventory(droot, dataset)
        print(f"  {dataset}: {inv.summary()}")
        print(
            f"  {dataset}: {len(raws)} raw, {len(planned)} in plan → "
            f"{len(plan.current)} current (skip), {len(plan.pending)} pending, "
            f"{len(plan.outdated)} outdated (re-convert); {n_todo} to convert"
        )
        return 0

    print(
        f"  {dataset}: converting {n_todo} of {len(planned)} planned → {out_parent} "
        f"(keep_profile={args.keep_profile}, centroid_algo={args.centroid_algo}, "
        f"precursor_major={args.precursor_major})",
        flush=True,
    )
    results = driver.run_conversion(
        planned,
        out_parent,
        keep_profile=args.keep_profile,
        centroid_algo=args.centroid_algo,
        precursor_major=args.precursor_major,
        n_workers=args.workers,
        force=args.force,
        rt_bin_width_s=args.rt_bin_width_s,
        batch_size=args.batch_size,
        progress_cb=make_progress(args.quiet),
    )
    summary = driver.summarize_results(results)
    print(f"  {dataset}: {summary.n_ok} ok, {summary.n_skipped} skipped, {summary.n_error} failed")
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
        "--no-profile",
        dest="keep_profile",
        action="store_false",
        help="do not write the sibling <stem>.profile/ bundle (default: write it)",
    )
    p.add_argument(
        "--centroid-algo",
        default=driver.DEFAULT_CENTROID_ALGO,
        help=f"centroider for profile scans (default: {driver.DEFAULT_CENTROID_ALGO!r})",
    )
    p.add_argument(
        "--precursor-major",
        action=argparse.BooleanOptionalAction,
        default=driver.DEFAULT_PRECURSOR_MAJOR,
        help="also write the precursor-major peak copy (constellation #109)",
    )
    p.add_argument(
        "--workers", type=int, default=1, help="parallel converter workers (spawn-mode; default 1)"
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="re-convert even current v3 bundles",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="report the on-disk inventory and what would convert; change nothing",
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
        "--raw-name",
        action="append",
        metavar="STEM",
        help="convert only this .raw (stem or file name; repeatable) — pilot runs",
    )
    p.add_argument(
        "--limit", type=int, default=None, help="cap the number of files converted (smoke tests)"
    )
    p.add_argument("--shard", type=int, default=0, help="SLURM-array shard index (with --n-shards)")
    p.add_argument("--n-shards", type=int, default=None, help="number of SLURM-array shards")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)

    if not args.dry_run:
        # CLR pre-flight in the parent (pure file check; convert_batch loads the
        # CLR only inside spawn workers — never fork a CLR-initialised parent).
        from constellation.massspec.readers.thermo import require_thermo

        try:
            require_thermo()
        except ImportError as exc:
            print(f"error: Thermo reader unavailable: {exc}")
            return 3
        print(bundles.constellation_provenance(), flush=True)

    rc = 0
    for dataset in resolve_datasets(args.dataset, config):
        rc |= _convert_dataset(args, config, dataset)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
