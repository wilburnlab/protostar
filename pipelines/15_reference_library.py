#!/usr/bin/env python
"""Stage 15 — fetch the published ProteomeTools ``.msp`` spectral libraries.

Downloads the eight per-mode library zips from Zenodo record 15705607 (the
archived 2019 ProteomeTools release linked from proteometools.org/id=53),
MD5-verifies each, and optionally extracts the ``.msp`` members into
``<data_root>/libraries/<mode>/``.

Ingesting the ``.msp`` into a Constellation ``massspec.library.Library`` and
associating them with the raw acquisitions is a separate, later task — this
stage stops at verified, extracted ``.msp`` on disk.

Usage::

    python pipelines/15_reference_library.py --dry-run
    python pipelines/15_reference_library.py --print-config
    python pipelines/15_reference_library.py --modes all --extract
"""

from __future__ import annotations

import argparse

from _common import data_root, load_config, make_progress

from protostar.library import zenodo


def _parse_modes(value: str | None) -> "set[str] | None":
    if value is None or value.strip().lower() == "all":
        return None
    return {tok.strip() for tok in value.split(",") if tok.strip()}


def _print_config(files: list[zenodo.ZenodoFile]) -> None:
    print("[libraries]")
    print('source_url    = "https://www.proteometools.org/index.php?id=53"')
    print(f'zenodo_record = "{zenodo.DEFAULT_RECORD_ID}"')
    print("files = [")
    for f in files:
        print(f'    "{f.key}",')
    print("]")


def main(argv: "list[str] | None" = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--config", default=None, help="path to datasets.toml")
    p.add_argument("--data-root", default=None, help="override the data root")
    p.add_argument(
        "--modes",
        default="all",
        help="comma list of mode tokens (e.g. FTMS_HCD_28,ITMS_CID_35) or 'all'",
    )
    p.add_argument("--extract", action="store_true", help="extract .msp members after download")
    p.add_argument(
        "--dry-run", action="store_true", help="list the selected library files; download nothing"
    )
    p.add_argument(
        "--print-config",
        action="store_true",
        help="emit the [libraries] block for datasets.toml and exit",
    )
    p.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")
    args = p.parse_args(argv)

    files = zenodo.list_library_files()
    if args.print_config:
        _print_config(files)
        return 0

    modes = _parse_modes(args.modes)
    selected = [f for f in files if modes is None or f.mode in modes]

    if args.dry_run:
        total_gb = sum(f.size_bytes for f in selected) / 1e9
        print(f"selected {len(selected)} library zip(s), {total_gb:.2f} GB:")
        for f in selected:
            print(f"  {f.mode:12s} {f.key}  ({f.size_bytes / 1e6:.0f} MB, md5 {f.md5[:8]}…)")
        return 0

    config = load_config(args.config)
    lib_dir = data_root(config, args.data_root) / "libraries"
    print(f"fetching {len(selected)} library zip(s) → {lib_dir}")
    zips = zenodo.fetch_libraries(lib_dir, modes=modes, progress_cb=make_progress(args.quiet))
    print(f"  downloaded + verified {len(zips)} zip(s)")

    if args.extract:
        msps = zenodo.extract_msp(zips, lib_dir)
        print(f"  extracted {len(msps)} .msp file(s) under {lib_dir}/<mode>/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
