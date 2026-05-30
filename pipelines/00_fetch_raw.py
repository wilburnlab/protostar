#!/usr/bin/env python
"""Stage 00 — fetch ProteomeTools files (canonical: fresh download + verify).

Three actions over the per-dataset PRIDE file set:

* ``--build-manifest`` — query the ProteomeXchange/PRIDE v3 API for every file
  under each accession and write the committed expected manifest
  (``config/manifests/<dataset>.json``).
* ``--seed-from DIR`` — relocate existing local ``.raw`` copies (the lab's
  Cartographer tree) into the canonical data root instead of downloading
  (default ``--seed-mode move``; RAW only).
* default — download every still-missing file (RAW + SEARCH by default),
  verifying each against its published SHA-1 and re-fetching on mismatch
  (resumable + repairable). ``--dry-run`` reports present/missing/corrupt.

Usage::

    python pipelines/00_fetch_raw.py --build-manifest --dataset all
    python pipelines/00_fetch_raw.py --dataset Zolg2017 --dry-run
    python pipelines/00_fetch_raw.py --dataset Zolg2017 \\
        --seed-from /path/to/existing/ProteomeTools   # or bare --seed-from (uses config)
    python pipelines/00_fetch_raw.py --dataset Zolg2017          # download the rest
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from _common import (
    add_common_args,
    add_dataset_arg,
    data_root,
    dataset_specs,
    load_config,
    make_progress,
    resolve_datasets,
    seed_from_default,
)

from protostar.fetch import manifest, net, seed

_CATEGORY_ALIASES = {"raw": "RAW", "search": "SEARCH", "other": "OTHER"}


def _parse_categories(value: str) -> "tuple[str, ...] | None":
    if value.strip().lower() == "all":
        return None
    out = []
    for token in value.split(","):
        token = token.strip().lower()
        if not token:
            continue
        out.append(_CATEGORY_ALIASES.get(token, token.upper()))
    return tuple(out)


def _download_and_verify(
    entry: manifest.ManifestEntry, dest: Path, *, progress_cb=None, retries: int = 1
) -> tuple[str, str | None]:
    """Download one file, verify size + (published) SHA-1, retry once on mismatch.

    Returns ``("ok", None)`` / ``("corrupt", detail)`` / ``("error", detail)``.
    Downloads are always integrity-checked — this is independent of the
    ``--verify`` flag (which only governs re-hashing *existing* files).
    """
    if not entry.https_url:
        return "error", "no download URL published"
    for attempt in range(retries + 1):
        try:
            res = net.stream_download(
                entry.https_url,
                dest,
                expected_size=entry.size_bytes,
                compute_sha1=bool(entry.sha1),
                progress_cb=progress_cb,
            )
        except Exception as exc:  # noqa: BLE001 — surface as a per-file error
            if attempt < retries:
                continue  # retry; stream_download resumes from the .part prefix
            return "error", f"{type(exc).__name__}: {exc}"
        if entry.sha1 and res.sha1 and res.sha1.lower() != entry.sha1.lower():
            dest.unlink(missing_ok=True)  # discard + retry
            if attempt >= retries:
                return "corrupt", f"sha1 mismatch (got {res.sha1[:12]}…)"
            continue
        return "ok", None
    return "corrupt", "exhausted retries"


def _print_table(dataset: str, statuses: list[manifest.FileStatus]) -> None:
    table = manifest.summarize(statuses)
    total = len(statuses)
    print(f"  {dataset}: {total} expected file(s)")
    for cat in sorted(table):
        row = table[cat]
        cells = "  ".join(f"{state}={n}" for state, n in sorted(row.items()))
        print(f"    {cat:10s} {cells}")


def _fetch_dataset(args: argparse.Namespace, config: dict, dataset: str) -> int:
    droot = data_root(config, args.data_root)
    categories = _parse_categories(args.categories)
    progress = make_progress(args.quiet)

    mpath = manifest.manifest_path(dataset, manifest_dir=args.manifest_dir)
    if not mpath.exists():
        raise SystemExit(f"no manifest for {dataset} at {mpath}; run with --build-manifest first")
    m = manifest.load_manifest(mpath)

    if args.seed_from:
        results = seed.seed_dataset(
            m,
            droot,
            args.seed_from,
            mode=args.seed_mode,
            verify=args.verify,
            dry_run=args.dry_run,
            progress_cb=progress,
        )
        seeded = sum(1 for r in results if r.action == "seeded")
        missing = sum(1 for r in results if r.action == "missing_source")
        bad = sum(1 for r in results if r.action in ("size_mismatch", "corrupt"))
        verb = "would seed" if args.dry_run else "seeded"
        print(
            f"  {dataset}: {verb} {seeded} raw file(s); {missing} not found locally, {bad} mismatched"
        )

    statuses = manifest.reconcile(m, droot, categories=categories, verify=args.verify)
    if args.dry_run:
        _print_table(dataset, statuses)
        return 0

    todo = sorted((s.entry for s in statuses if s.state != "ok"), key=lambda e: e.file_name)
    if args.n_shards and args.n_shards > 1:
        todo = todo[args.shard :: args.n_shards]
    if args.limit is not None:
        todo = todo[: args.limit]

    if not todo:
        print(f"  {dataset}: nothing to fetch (all present)")
        manifest.write_status(dataset, droot, statuses)
        return 0

    print(f"  {dataset}: fetching {len(todo)} file(s) with {args.workers} worker(s)")
    outcomes: dict[str, tuple[str, str | None]] = {}

    def work(e: manifest.ManifestEntry) -> tuple[str, tuple[str, str | None]]:
        dest = manifest.local_path(e, droot, dataset)
        cb = progress if args.workers <= 1 else None
        return e.file_name, _download_and_verify(e, dest, progress_cb=cb)

    if args.workers <= 1:
        for e in todo:
            name, outcome = work(e)
            outcomes[name] = outcome
            print(f"    {name} → {outcome[0]}" + (f" ({outcome[1]})" if outcome[1] else ""))
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(work, e): e for e in todo}
            for fut in as_completed(futs):
                name, outcome = fut.result()
                outcomes[name] = outcome
                print(f"    {name} → {outcome[0]}" + (f" ({outcome[1]})" if outcome[1] else ""))

    n_ok = sum(1 for v in outcomes.values() if v[0] == "ok")
    n_bad = len(outcomes) - n_ok
    print(f"  {dataset}: {n_ok} ok, {n_bad} failed")
    manifest.write_status(
        dataset,
        droot,
        manifest.reconcile(m, droot, categories=categories, verify=False),
    )
    return 1 if n_bad else 0


def cmd_build_manifest(args: argparse.Namespace, config: dict) -> int:
    specs = dataset_specs(config)
    for dataset in resolve_datasets(args.dataset, config):
        accession = specs[dataset]["accession"]
        print(f"building manifest for {dataset} ({accession})…")
        m = manifest.build_manifest(dataset, accession)
        path = manifest.write_manifest(
            m, manifest.manifest_path(dataset, manifest_dir=args.manifest_dir)
        )
        print(f"  wrote {path}  ({len(m.entries)} files; {m.to_dict()['n_by_category']})")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    add_common_args(p)
    add_dataset_arg(p)
    p.add_argument(
        "--build-manifest",
        action="store_true",
        help="query PRIDE and (re)write the committed manifests, then exit",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="reconcile against the data tree and print a table; change nothing",
    )
    p.add_argument(
        "--categories",
        default="raw,search",
        help="comma list of categories to fetch (raw,search,other | all). Default: raw,search",
    )
    p.add_argument(
        "--seed-from",
        nargs="?",
        const="__config__",
        default=None,
        help="relocate existing local .raw copies from this dir instead of downloading; bare flag uses [defaults].seed_from from config",
    )
    p.add_argument(
        "--seed-mode",
        choices=["move", "hardlink", "copy"],
        default="move",
        help="how --seed-from places files (default: move)",
    )
    p.add_argument(
        "--verify",
        action="store_true",
        help="also stream the SHA-1 of present/seeded files (slow; size-only otherwise)",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=4,
        help="parallel download workers (network-bound; default 4)",
    )
    p.add_argument(
        "--limit", type=int, default=None, help="cap the number of files acted on (smoke tests)"
    )
    p.add_argument("--shard", type=int, default=0, help="SLURM-array shard index (with --n-shards)")
    p.add_argument("--n-shards", type=int, default=None, help="number of SLURM-array shards")
    return p


def main(argv: "list[str] | None" = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    # --seed-from is opt-in: None = don't seed; bare flag = config default; else the given dir.
    if args.seed_from == "__config__":
        args.seed_from = seed_from_default(config)
        if not args.seed_from:
            raise SystemExit("--seed-from given with no path and no [defaults].seed_from in config")
    if args.seed_from is not None:
        args.seed_from = Path(args.seed_from)
    if args.build_manifest:
        return cmd_build_manifest(args, config)
    rc = 0
    for dataset in resolve_datasets(args.dataset, config):
        rc |= _fetch_dataset(args, config, dataset)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
