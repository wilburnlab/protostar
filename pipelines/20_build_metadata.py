#!/usr/bin/env python
"""Stage 20 — acquisition time table.

Reads acquisition datetime + instrument from each stage-10 convert bundle's
``acquisition_metadata.parquet`` / ``manifest.json``, orders runs chronologically
per instrument, and persists a ``constellation.massspec.acquisitions.Acquisitions``
table per dataset (``proc/<dataset>/acquisitions.parquet``) — the substrate for
carryover and batch-effect analysis (acquisition-order reordering,
instrument-specific effects). Consumes the ``proc/`` tree rather than re-opening
every ``.raw``. Resume is free — a dataset whose ``acquisitions.parquet`` exists
is skipped unless ``--force``.

Usage::

    python pipelines/20_build_metadata.py --dataset Zolg2017 --dry-run
    python pipelines/20_build_metadata.py --dataset Zolg2017
    python pipelines/20_build_metadata.py --dataset all --force
"""

from __future__ import annotations

import argparse

from _common import (
    add_common_args,
    add_dataset_arg,
    data_root,
    load_config,
    resolve_datasets,
)

from protostar.metadata import driver


def _build_dataset(args: argparse.Namespace, config: dict, dataset: str) -> int:
    droot = data_root(config, args.data_root)
    bundles = driver.enumerate_bundles(droot, dataset)
    out_path = driver.acquisitions_out_path(droot, dataset)

    if not bundles:
        print(
            f"  {dataset}: no converted bundles under "
            f"{driver.proc_parent(droot, dataset)} — convert first"
        )
        return 0

    if args.dry_run:
        # Build to report instrument/datetime stats, but write nothing.
        acq = driver.build_acquisitions(droot, dataset)
        s = driver.summarize(acq, len(bundles))
        print(
            f"  {dataset}: {s.n_bundles} bundles → {s.n_acquisitions} acquisitions, "
            f"{s.n_instruments} instrument(s), {s.n_datetime_null} undated → {out_path}"
        )
        return 0

    if driver.acquisitions_present(droot, dataset) and not args.force:
        print(f"  {dataset}: acquisitions.parquet present — skipping (use --force to rebuild)")
        return 0

    acq = driver.build_acquisitions(droot, dataset)
    driver.write_acquisitions(acq, out_path)
    s = driver.summarize(acq, len(bundles))
    print(
        f"  {dataset}: wrote {s.n_acquisitions} acquisitions "
        f"({s.n_instruments} instrument(s), {s.n_datetime_null} undated) → {out_path}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    add_common_args(p)
    add_dataset_arg(p)
    p.add_argument(
        "--force",
        action="store_true",
        help="rebuild acquisitions.parquet even if it already exists",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="report bundle/acquisition/instrument counts; write nothing",
    )
    return p


def main(argv: "list[str] | None" = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    rc = 0
    for dataset in resolve_datasets(args.dataset, config):
        rc |= _build_dataset(args, config, dataset)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
