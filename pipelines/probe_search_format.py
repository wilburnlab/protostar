#!/usr/bin/env python
"""Diagnostic — confirm the PRIDE SEARCH files are MaxQuant output.

Downloads one SEARCH file per dataset to a scratch dir, inspects the zip's
entries, and reports whether it carries the MaxQuant ``txt/`` signature
(``msms.txt`` / ``evidence.txt`` / ``peptides.txt`` / ``parameters.txt`` / ...).
If confirmed, a reader for these (``constellation.massspec.io.maxquant``) is the
next-PR target tracked in ``docs/constellation_contributions.md`` (#11), needed
to later associate search identifications with the raw acquisitions.

Usage::

    python pipelines/probe_search_format.py
    python pipelines/probe_search_format.py --keep --scratch /tmp/pt_search
"""

from __future__ import annotations

import argparse
import shutil
import tempfile
import zipfile
from pathlib import Path

from _common import dataset_specs, load_config, resolve_datasets

from protostar.fetch import net, pride

# Hallmark members of a MaxQuant ``combined/txt/`` export.
_MAXQUANT_SIGNATURE = {
    "msms.txt",
    "evidence.txt",
    "peptides.txt",
    "parameters.txt",
    "proteinGroups.txt",
    "modificationSpecificPeptides.txt",
    "allPeptides.txt",
    "msmsScans.txt",
}
# Strong markers — any one is conclusive.
_CONCLUSIVE = {"msms.txt", "evidence.txt", "parameters.txt"}


def _probe_dataset(dataset: str, accession: str, scratch: Path) -> bool:
    search = [f for f in pride.list_files(accession, categories={"SEARCH"}) if f.https_url]
    if not search:
        print(f"{dataset} ({accession}): no SEARCH files listed")
        return False
    pick = min(search, key=lambda f: f.size_bytes)
    dest = scratch / dataset / pick.file_name
    print(f"{dataset} ({accession}): probing {pick.file_name} ({pick.size_bytes / 1e6:.1f} MB)")
    net.stream_download(pick.https_url, dest, expected_size=pick.size_bytes, compute_sha1=False)

    if not zipfile.is_zipfile(dest):
        print(f"  NOT a zip — cannot classify ({pick.file_name})")
        return False
    with zipfile.ZipFile(dest) as zf:
        basenames = {Path(m).name for m in zf.namelist()}
    found = sorted(basenames & _MAXQUANT_SIGNATURE)
    is_maxquant = bool(basenames & _CONCLUSIVE)
    sample = sorted(basenames)[:8]
    print(f"  entries (sample): {sample}")
    print(f"  MaxQuant signature files found: {found or 'none'}")
    print(f"  VERDICT: {'MaxQuant' if is_maxquant else 'UNKNOWN — inspect manually'}")
    return is_maxquant


def main(argv: "list[str] | None" = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--config", default=None, help="path to datasets.toml")
    p.add_argument(
        "--scratch",
        type=Path,
        default=None,
        help="download dir (default: a temp dir, removed unless --keep)",
    )
    p.add_argument("--keep", action="store_true", help="keep the downloaded probe files")
    p.add_argument("--dataset", action="append", help="dataset(s) to probe (default: all)")
    args = p.parse_args(argv)

    config = load_config(args.config)
    specs = dataset_specs(config)
    scratch = args.scratch or Path(tempfile.mkdtemp(prefix="protostar_search_probe_"))
    scratch.mkdir(parents=True, exist_ok=True)

    verdicts: dict[str, bool] = {}
    try:
        for dataset in resolve_datasets(args.dataset, config):
            verdicts[dataset] = _probe_dataset(dataset, specs[dataset]["accession"], scratch)
            print()
    finally:
        if not args.keep and args.scratch is None:
            shutil.rmtree(scratch, ignore_errors=True)

    all_mq = all(verdicts.values()) and bool(verdicts)
    print("summary:", {d: ("MaxQuant" if v else "unknown") for d, v in verdicts.items()})
    if all_mq:
        print(
            "→ all SEARCH files are MaxQuant; ledger #11 (massspec.io.maxquant reader) confirmed as the target."
        )
    return 0 if all_mq else 1


if __name__ == "__main__":
    raise SystemExit(main())
