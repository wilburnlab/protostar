#!/usr/bin/env python
"""Stage 30 (experiment) — comprehensive calibrant MS2 extraction (array-parallel).

Extracts, for the QC calibrant peptides across **all non-ETD acquisitions** of the
selected datasets, each MS2 scan's b/y fragment **proportion vector** + total
intensity, with the per-scan **fragmentation mode** (analyzer + activation +
collision energy) deconvolved from ``scan_metadata`` — the instrument truth, so
the targeted injections' interleaved modes (the 9 modeled HCD/CID modes) separate
cleanly. Writes one parquet shard per task that the noise-scaling / metric
experiments read back (extract once, analyze many).

Shardable for SLURM arrays: ``--shard i --n-shards N`` processes search zips
``[i::N]`` deterministically. Search zips are extracted to node-local ``$TMPDIR``
(not ESS) to avoid inode blow-up.

Output rows: dataset, modified_sequence, charge, mode, raw_file, scan,
total_intensity, iit (ion injection time, ms), ion_proxy (= total_intensity·iit,
∝ ion count N since intensity is per-time), props (list<double>, the
basis-aligned fragment proportions).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _common import add_common_args, data_root, load_config

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from protostar.experiments import consensus_assembly, ms2_extract, scans

DEFAULT_DATASETS = ("Zolg2017", "Gessulat2019", "Wilhelm2021")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    add_common_args(p)
    p.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    p.add_argument(
        "--exclude-injection", default="ETD", help="skip search zips whose name contains this"
    )
    p.add_argument("--peptides", choices=("qc", "procal", "all"), default="qc")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument("--out", type=Path, required=True)
    return p


def _enumerate_zips(droot, datasets, exclude):
    out = []
    for ds in datasets:
        for z in sorted((droot / "search" / ds).glob("*.zip")):
            if exclude and exclude.lower() in z.name.lower():
                continue
            out.append((ds, z))
    return out


#: Per-scan scan_metadata fields carried into the rows (alongside the ``mode``
#: string). ``iit``+``tic`` drive the ion-count axis and assigned fraction;
#: ``rt`` localizes co-elution; the AGC/space-charge readouts (``agc_fill``,
#: ``agc_target``, ``space_charge_comp_ppm``) let the high-abundance
#: detector-nonlinearity hypothesis be tested against the instrument's own
#: telemetry; ``base_peak_intensity``/``peak_count`` index saturation/complexity.
_META_FIELDS = (
    "iit",
    "tic",
    "rt",
    "agc_fill",
    "agc_target",
    "space_charge_comp_ppm",
    "base_peak_intensity",
    "peak_count",
)


def _scan_meta(bundle_dir: Path) -> tuple[dict[int, str], dict[int, dict]]:
    """Per-scan ``mode`` (``ANALYZER_activation_CE``) and a metadata panel
    (``_META_FIELDS``) from scan_metadata. Intensity is per-time, so the ion
    count is ``N ∝ I·iit`` — iit puts spectra on a true ion-count axis."""
    t = pq.read_table(
        bundle_dir / "scan_metadata.parquet",
        columns=["scan", "analyzer", "activation_type", "collision_energy", *_META_FIELDS],
    )
    cols = {c: t.column(c).to_pylist() for c in t.column_names}
    modes: dict[int, str] = {}
    meta: dict[int, dict] = {}
    for i in range(t.num_rows):
        s = int(cols["scan"][i])
        a, act, ce = cols["analyzer"][i], cols["activation_type"][i], cols["collision_energy"][i]
        if a and act and ce is not None:
            modes[s] = f"{a}_{act}_{round(ce)}"
        meta[s] = {f: cols[f][i] for f in _META_FIELDS}
    return modes, meta


def _peptide_filter(psms, scope):
    if scope == "all":
        return psms
    col = "is_procal" if scope == "procal" else "is_qc"
    return psms.filter(pc.fill_null(psms.column(col), False))


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    droot = data_root(config, args.data_root)
    args.out.mkdir(parents=True, exist_ok=True)
    tmp = Path(os.environ.get("TMPDIR", "/tmp")) / f"exp30_extract_{args.shard}"

    zips = _enumerate_zips(droot, args.datasets, args.exclude_injection)
    mine = zips[args.shard :: args.n_shards]
    print(f"shard {args.shard}/{args.n_shards}: {len(mine)}/{len(zips)} acquisitions", flush=True)

    basis_cache: dict[str, object] = {}
    rows: list[dict] = []
    for k, (ds, zpath) in enumerate(mine):
        try:
            psms = scans.load_psms(zpath, extract_dir=tmp / zpath.stem)
            psms = scans.tag_calibrants(scans.gate_psms(psms), dataset=ds)
            psms = _peptide_filter(psms, args.peptides)
            psms = psms.filter(pc.is_valid(psms.column("modified_sequence")))
        except Exception as exc:  # noqa: BLE001 — robustness across 1000s of zips
            print(f"  [skip search] {zpath.name}: {exc}", flush=True)
            continue
        if psms.num_rows == 0:
            continue
        for raw_file in set(psms.column("raw_file").to_pylist()):
            bundle = droot / "proc" / ds / "centroid" / raw_file
            if not (bundle / "peaks.parquet").is_file():
                continue
            sub = psms.filter(pc.equal(psms.column("raw_file"), raw_file))
            try:
                peaks = pq.read_table(bundle / "peaks.parquet")
                modes, meta = _scan_meta(bundle)
                trace = ms2_extract.extract_ms2_fragments(peaks, sub)
                scan_channels = ms2_extract.trace_to_scan_channels(trace)
            except Exception as exc:  # noqa: BLE001
                print(f"  [skip bundle] {raw_file}: {exc}", flush=True)
                continue
            for ms, ch, sc, perr, score, dscore, pep in zip(
                sub.column("modified_sequence").to_pylist(),
                sub.column("charge").to_pylist(),
                sub.column("scan").to_pylist(),
                sub.column("mass_error_ppm").to_pylist(),
                sub.column("score").to_pylist(),
                sub.column("delta_score").to_pylist(),
                sub.column("pep").to_pylist(),
            ):
                channels = scan_channels.get(sc)
                mode = modes.get(sc)
                m = meta.get(sc)
                iit = m["iit"] if m else None
                if not channels or mode is None or iit is None:
                    continue
                if ms not in basis_cache:
                    try:
                        basis_cache[ms] = consensus_assembly.basis_for(ms)
                    except Exception:  # noqa: BLE001
                        basis_cache[ms] = None
                basis = basis_cache[ms]
                if basis is None:
                    continue
                # Map fragments to the basis by exact (ion_type, position, charge);
                # the XIC trace's per-fragment mz_error_ppm is carried through
                # (intensity-weighted) — no second m/z match.
                vec, err_ppm = ms2_extract.channels_to_basis(channels, basis)
                total = float(vec.sum())
                if total <= 0:
                    continue
                rows.append(
                    {
                        "dataset": ds,
                        "modified_sequence": ms,
                        "charge": int(ch),
                        "mode": mode,
                        "raw_file": raw_file,
                        "scan": int(sc),
                        "total_intensity": total,
                        "iit": float(iit),
                        "ion_proxy": total * float(iit),  # ∝ N (ion count); intensity is per-time
                        "props": (vec / total).tolist(),
                        # per-channel intensity-weighted signed m/z error (ppm); NaN unmatched.
                        "mz_err_ppm": err_ppm.tolist(),
                        # precursor MS1 mass error (ppm, MaxQuant) — MS1 co-isolation signature.
                        "precursor_mz_err_ppm": perr,
                        # full MS2 TIC -> assigned fraction = total_intensity / scan_tic.
                        "scan_tic": m["tic"],
                        # retention time (s) — co-elution localization (hypothesis 1).
                        "rt": m["rt"],
                        # AGC / space-charge telemetry — the detector-nonlinearity test (hypothesis 2).
                        "agc_fill": m["agc_fill"],
                        "agc_target": m["agc_target"],
                        "space_charge_comp_ppm": m["space_charge_comp_ppm"],
                        "base_peak_intensity": m["base_peak_intensity"],
                        "peak_count": m["peak_count"],
                        # MaxQuant PSM confidence — does KL track ID quality?
                        "score": score,
                        "delta_score": dscore,
                        "pep": pep,
                    }
                )
        if (k + 1) % 10 == 0 or k + 1 == len(mine):
            print(f"  [{k + 1}/{len(mine)}] {len(rows)} spectra", flush=True)

    out_file = args.out / f"shard_{args.shard:04d}.parquet"
    pq.write_table(pa.Table.from_pylist(rows), out_file)
    print(f"wrote {out_file} ({len(rows)} spectra)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
