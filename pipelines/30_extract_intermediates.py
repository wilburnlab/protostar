#!/usr/bin/env python
"""Stage 30 — extract MS1 isotope chromatograms for the Counter-motivation figures.

For the seed peptides (PROCAL P1–P40 + the prior notebooks' known-clear examples)
across the selected datasets' **non-ETD** acquisitions, extracts each MS1 survey
scan's M+0/M+1/M+2 isotope intensities + signed m/z error over an RT window around
the peptide's elution, at **every charge 1..4**, with the per-scan IIT/TIC — the
"score-ready" intermediate the figure scripts read (extract once, plot many). A
single extraction feeds all four figures: isotope shot noise, charge-dependent
α(z), m/z-vs-N + Student-t, and the interference chromatogram.

Shardable for SLURM arrays: ``--shard i --n-shards N`` processes the search-zip
list ``[i::N]`` deterministically. Search zips extract to node-local ``$TMPDIR``.

Reads **bundle v3** runs (``proc/<dataset>/<stem>/``) through ``protostar.bundles``.
- The peaks are streamed from the bundle; nothing is loaded whole.
- A run that has not been reconverted is counted and reported.
- A pre-v3 bundle at the v3 location raises an error.
- A failed extraction is counted, and the shard exits non-zero.
No run is silently skipped.

Output (one parquet shard per task; wide, one row per peptide·charge·scan):
- identity: dataset, modified_sequence, sequence, psm_charge, charge, raw_file, scan,
  rt, rt_center, is_procal;
- isotopes: mz_mono, I_0/I_1/I_2, sum_I, mz_ppm_0/1/2;
- scan metadata: iit, iit_s, tic, base_peak_intensity, peak_count, analyzer;
- when present in the bundle: intensity_floor, max_ion_time_ms, iit_at_max,
  agc_target, agc_fill, ft_resolution, space_charge_comp_ppm,
  lock_mass_correction_ppm.

No model code lives here — extraction is Constellation, the figure fits live in
``protostar.experiments``.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# _common first: it pins the conda env's libstdc++ (re-exec) before pyarrow loads.
from _common import add_common_args, data_root, load_config, resolve_datasets

# isort: split
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from protostar import bundles
from protostar.experiments import scans
from protostar.intermediates import ms1_extract

#: The figures are built on Zolg2017 (where the reference notebooks ran).
DEFAULT_DATASET = "Zolg2017"

#: Prior notebooks' known-clear example peptides (bare AA), extracted alongside
#: all of PROCAL. Fig 1 shot noise / Fig 3 m/z / Fig 4 interference, respectively.
SEED_PEPTIDES: tuple[str, ...] = (
    "ALFSSITDSEK",
    "ASDLLSGYYIK",
    "FGTGTYAGGEK",
    "YFGYTSDTFGK",
    "GFLDYESTGAK",
    "LYTGAGYDEVK",
    "HLTGLTFDTYK",
)

#: Provenance columns added per (peptide, charge, scan) wide row, taken from the
#: deduped PSM anchor by target_id.
_PROVENANCE = ("sequence", "modified_sequence", "psm_charge", "rt_center", "is_procal")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    add_common_args(p)
    p.add_argument(
        "--dataset",
        action="append",
        metavar="NAME",
        help=f"dataset to extract (repeatable; 'all' for every configured). Default: {DEFAULT_DATASET}.",
    )
    p.add_argument("--peptides", choices=("procal", "seed", "all"), default="seed")
    p.add_argument(
        "--exclude-injection",
        default="ETD",
        help="skip search zips whose name contains this (ETD is outside the modeled modes)",
    )
    p.add_argument("--rt-window", type=float, default=60.0, help="± seconds around each PSM RT")
    p.add_argument(
        "--all-in-window",
        action="store_true",
        help="return every in-tolerance peak per isotope (chimeric/interference; breaks the wide pivot)",
    )
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="output dir for shards (default: <data_root>/proc/intermediates/ms1)",
    )
    return p


def _enumerate_zips(droot: Path, datasets: list[str], exclude: str) -> list[tuple[str, Path]]:
    out: list[tuple[str, Path]] = []
    for ds in datasets:
        for z in sorted((droot / "search" / ds).glob("*.zip")):
            if exclude and exclude.lower() in z.name.lower():
                continue
            out.append((ds, z))
    return out


def _peptide_filter(psms: pa.Table, scope: str) -> pa.Table:
    if scope == "all":
        return psms
    is_procal = pc.fill_null(psms.column("is_procal"), False)
    if scope == "procal":
        return psms.filter(is_procal)
    in_seed = pc.is_in(psms.column("sequence"), value_set=pa.array(SEED_PEPTIDES, pa.string()))
    return psms.filter(pc.or_(is_procal, in_seed))


def _best_per_peptide(psms: pa.Table) -> pa.Table:
    """One PSM per ``modified_sequence`` (lowest pep, then highest score) — the RT
    anchor for the multi-charge MS1 extraction. Deduping keeps the extraction one
    window per peptide; ``charge_range=(1,4)`` then captures every charge."""
    order = pc.sort_indices(psms, sort_keys=[("pep", "ascending"), ("score", "descending")])
    psms = psms.take(order)
    seen: set[str] = set()
    keep: list[int] = []
    for i, s in enumerate(psms.column("modified_sequence").to_pylist()):
        if s is not None and s not in seen:
            seen.add(s)
            keep.append(i)
    return psms.take(pa.array(keep, pa.int64()))


def _extract_bundle(bundle: Path, sub: pa.Table, ds: str, raw_file: str, args) -> pa.Table | None:
    """Extract one v3 bundle's MS1 isotope chromatograms for the deduped PSM anchors,
    returning a wide table with provenance, or ``None`` when nothing matched."""
    scan_meta = bundles.read_scan_metadata(bundle)
    anchors = _best_per_peptide(sub)  # target_id = row index into `anchors`
    if anchors.num_rows == 0:
        return None
    targets = ms1_extract.ms1_targets_from_psms(anchors)
    trace = ms1_extract.extract_ms1_isotopes(
        bundle, targets, rt_window=args.rt_window, all_in_window=args.all_in_window
    )
    wide = ms1_extract.attach_scan_meta(ms1_extract.trace_to_wide(trace), scan_meta)
    if wide.num_rows == 0:
        return None
    tids = wide.column("target_id")
    prov = {
        "dataset": pa.array([ds] * wide.num_rows, pa.string()),
        "raw_file": pa.array([raw_file] * wide.num_rows, pa.string()),
        "sequence": pc.take(anchors.column("sequence"), tids),
        "modified_sequence": pc.take(anchors.column("modified_sequence"), tids),
        "psm_charge": pc.cast(pc.take(anchors.column("charge"), tids), pa.int8()),
        "rt_center": pc.cast(pc.take(anchors.column("retention_time_s"), tids), pa.float64()),
        "is_procal": pc.fill_null(pc.take(anchors.column("is_procal"), tids), False),
    }
    for name, col in prov.items():
        wide = wide.append_column(name, col)
    return wide


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    droot = data_root(config, args.data_root)
    datasets = resolve_datasets(args.dataset, config) if args.dataset else [DEFAULT_DATASET]
    out_dir = args.out if args.out is not None else droot / "proc" / "intermediates" / "ms1"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"stage30_ms1_{args.shard}"

    zips = _enumerate_zips(droot, datasets, args.exclude_injection)
    mine = zips[args.shard :: args.n_shards]
    print(
        f"shard {args.shard}/{args.n_shards}: {len(mine)}/{len(zips)} acquisitions "
        f"(datasets={','.join(datasets)}, peptides={args.peptides})",
        flush=True,
    )

    parts: list[pa.Table] = []
    not_converted: list[str] = []
    failed: list[tuple[str, str]] = []
    for k, (ds, zpath) in enumerate(mine):
        try:
            psms = scans.load_psms(zpath, extract_dir=tmp / zpath.stem)
            psms = scans.tag_calibrants(scans.gate_psms(psms))
            psms = _peptide_filter(psms, args.peptides)
            psms = psms.filter(pc.is_valid(psms.column("modified_sequence")))
        except Exception as exc:  # noqa: BLE001 — robustness across 1000s of zips
            print(f"  [skip search] {zpath.name}: {exc}", flush=True)
            continue
        if psms.num_rows == 0:
            continue
        for raw_file in sorted(set(psms.column("raw_file").to_pylist())):
            # Raises on a pre-v3 bundle at the v3 location (convert again) or an
            # unfinished conversion; None = never converted to v3.
            bundle = bundles.require_run_bundle(droot, ds, raw_file)
            if bundle is None:
                not_converted.append(f"{ds}/{raw_file}")
                continue
            sub = psms.filter(pc.equal(psms.column("raw_file"), raw_file))
            try:
                wide = _extract_bundle(bundle, sub, ds, raw_file, args)
            except Exception as exc:  # noqa: BLE001 — counted + reported; shard exits 1
                failed.append((f"{ds}/{raw_file}", f"{type(exc).__name__}: {exc}"))
                print(f"  [FAILED bundle] {raw_file}: {type(exc).__name__}: {exc}", flush=True)
                continue
            if wide is not None and wide.num_rows:
                parts.append(wide)
        if (k + 1) % 10 == 0 or k + 1 == len(mine):
            n = sum(p.num_rows for p in parts)
            print(f"  [{k + 1}/{len(mine)}] {n} rows", flush=True)

    table = pa.concat_tables(parts, promote_options="default") if parts else _empty()
    out_file = out_dir / f"shard_{args.shard:04d}.parquet"
    pq.write_table(table, out_file)
    print(f"wrote {out_file} ({table.num_rows} rows)", flush=True)
    if not_converted:
        print(
            f"NOTE: {len(not_converted)} identified runs have no v3 bundle (not reconverted): "
            f"{', '.join(not_converted[:5])}{' …' if len(not_converted) > 5 else ''}",
            flush=True,
        )
    if failed:
        print(f"ERROR: {len(failed)} bundles failed extraction:", flush=True)
        for name, detail in failed:
            print(f"  {name}: {detail}", flush=True)
        return 1
    return 0


def _empty() -> pa.Table:
    """A 0-row table when a shard matched nothing (keeps the shard set rectangular)."""
    return pa.table({"dataset": pa.array([], pa.string())})


if __name__ == "__main__":
    sys.exit(main())
