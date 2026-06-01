#!/usr/bin/env python
"""Stage 05 — curate the synthetic-peptide reference from manuscript supplements.

Parses each dataset's published supplementary Excel workbook (a one-time **manual**
input — download from the articles; redistribution is not ours to make) into the
compact, version-controlled ``PEPTIDE_REFERENCE_TABLE`` parquet committed under
``protostar/peptides/data/<dataset>.parquet``. That parquet — sequence ↔ pool ↔
set, the join to the acquisition ``.raw`` stems — is the durable reference shipped
with the repo; the workbooks are not needed again once it is built.

Unlike the fetch stages, the source is local Excel, not a download: point
``--supplements-dir`` at the directory holding the workbooks (default file names
in ``protostar.peptides.sources.SUPPLEMENT_FILES``). Needs the ``peptides`` extra
(``openpyxl``) for parsing; reading the committed parquet later does not.

Usage::

    python pipelines/05_peptide_reference.py --supplements-dir ~/PT_peptide_info --dataset all
    python pipelines/05_peptide_reference.py --supplements-dir ~/PT_peptide_info --dataset Zolg2017 --dry-run
"""

from __future__ import annotations

import argparse
from pathlib import Path

from _common import add_dataset_arg, load_config, resolve_datasets

from protostar.peptides import reference, sources


def _build_dataset(args: argparse.Namespace, dataset: str) -> int:
    fname = sources.SUPPLEMENT_FILES.get(dataset)
    if fname is None:
        print(f"  {dataset}: no supplement mapping (known: {', '.join(sources.SUPPLEMENT_FILES)})")
        return 1
    xlsx = Path(args.supplements_dir).expanduser() / fname
    if not xlsx.is_file():
        print(f"  {dataset}: supplement not found at {xlsx} — download it (manual input)")
        return 1

    table = sources.parse_supplement(xlsx, dataset)
    n = table.num_rows
    n_pools = len({p for p in table.column("pool").to_pylist() if p})
    import pyarrow.compute as pc

    n_qc = pc.sum(pc.cast(table.column("is_qc"), "int64")).as_py()
    sets = sorted(set(table.column("set").to_pylist()))
    out = reference.reference_path(dataset)

    if args.dry_run:
        print(
            f"  {dataset}: {n} peptides, {n_pools} pools, {n_qc} QC, "
            f"sets={sets} → {out} (dry-run, not written)"
        )
        return 0

    reference.write_reference(table, dataset)
    print(f"  {dataset}: wrote {n} peptides ({n_pools} pools, {n_qc} QC, {len(sets)} sets) → {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--config", type=Path, default=None, help="path to datasets.toml")
    add_dataset_arg(p)
    p.add_argument(
        "--supplements-dir",
        type=Path,
        required=True,
        help="directory holding the manually-downloaded supplementary .xlsx workbooks",
    )
    p.add_argument("--dry-run", action="store_true", help="report counts; write no parquet")
    return p


def main(argv: "list[str] | None" = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    rc = 0
    for dataset in resolve_datasets(args.dataset, config):
        rc |= _build_dataset(args, dataset)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
