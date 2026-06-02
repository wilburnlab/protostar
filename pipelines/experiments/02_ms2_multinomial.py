#!/usr/bin/env python
"""exp02 — MS2 fragmentation is multinomial (the mean–variance law).

For high-replicate calibrant peptides (PROCAL / QC, spiked into every run), this
driver characterizes fragment intensities as ``Multinomial(N, p)`` draws and
checks the shot-noise law: across replicate spectra of one peptide, the
per-channel variance of the observed proportions ``p̂_k`` should track
``p̄_k(1 − p̄_k) / N`` with a single fitted slope ``1/N_eff`` per peptide. That
slope, vs total intensity, is the bridge to the MS1 ion count N (Part II).

Pipeline (all statistics imported from Constellation; this script only
orchestrates):

  searches (one per raw file) → gate PSMs → tag calibrants → keep QC/PROCAL →
  XIC level-2 assigned-scan MS2 (analyzer-split tolerance) → per-replicate
  fragment vectors on each peptide's fixed basis → median consensus +
  per_replicate[R,K] → variance law + multinomial-deviance χ² check.

Outputs (under ``<out>``): ``exp02_channels.parquet`` (per-channel pbar / var /
p(1−p)), ``exp02_groups.parquet`` (per-peptide N_eff / R² / deviance), and
figures 2.2 (mean–variance law) + 2.3 (N_eff vs intensity) under
``results/figures/``.

Runs on a compute node (loads each bundle's peaks table). Example:
    python pipelines/experiments/02_ms2_multinomial.py --dataset Zolg2017 \
        --injection DDA --max-acquisitions 30 --peptides qc
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

# pipelines/ on sys.path so `_common` is importable from this nested script;
# `_common` must be imported before constellation (it pins libstdc++ on HPC).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from _common import add_common_args, data_root, load_config

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import torch

from protostar.experiments import consensus_assembly, ms2_extract, scans

REPO_ROOT = Path(__file__).resolve().parents[2]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    add_common_args(p)
    p.add_argument("--dataset", default="Zolg2017", help="dataset (default Zolg2017)")
    p.add_argument(
        "--injection",
        default="DDA",
        help="injection-type token to select search zips (default DDA)",
    )
    p.add_argument(
        "--peptides",
        choices=("qc", "procal", "all"),
        default="qc",
        help="peptide scope: qc (labelled calibrants, default), procal (40), or all",
    )
    p.add_argument("--max-acquisitions", type=int, default=30)
    p.add_argument(
        "--min-replicates",
        type=int,
        default=8,
        help="min distinct acquisitions a (peptide,charge,mode) must recur in",
    )
    p.add_argument("--min-channels", type=int, default=4, help="min nonzero channels to fit")
    p.add_argument("--out", type=Path, default=None, help="output dir (default <data_root>/exp/02)")
    return p


def _select_searches(droot: Path, dataset: str, injection: str, limit: int) -> list[Path]:
    sdir = droot / "search" / dataset
    zips = sorted(p for p in sdir.glob(f"*{injection}*.zip"))
    return zips[:limit]


def _bundle_peaks(droot: Path, dataset: str, raw_file: str):
    bundle = droot / "proc" / dataset / "centroid" / raw_file / "peaks.parquet"
    return pq.read_table(bundle) if bundle.is_file() else None


def _peptide_filter(psms, scope: str):
    if scope == "all":
        return psms
    col = "is_procal" if scope == "procal" else "is_qc"
    return psms.filter(pc.fill_null(psms.column(col), False))


def _fit_through_origin(x: torch.Tensor, y: torch.Tensor) -> tuple[float, float]:
    """LS slope of y = m·x through the origin, plus R²."""
    denom = (x * x).sum()
    if denom <= 0:
        return float("nan"), float("nan")
    m = (x * y).sum() / denom
    ss_res = ((y - m * x) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    return float(m), r2


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    droot = data_root(config, args.data_root)
    out = args.out or (droot / "exp" / "02")
    out.mkdir(parents=True, exist_ok=True)
    extract_root = out / "_search_extract"

    searches = _select_searches(droot, args.dataset, args.injection, args.max_acquisitions)
    print(f"{len(searches)} {args.injection} acquisitions for {args.dataset}", flush=True)

    # accumulate per (modseq, charge, mode): replicate spectra + the raw files seen.
    spectra_by_group: dict[tuple, list] = defaultdict(list)
    files_by_group: dict[tuple, set] = defaultdict(set)
    n_parse_fail = 0

    for i, zpath in enumerate(searches):
        psms = scans.load_psms(zpath, extract_dir=extract_root / zpath.stem)
        psms = scans.gate_psms(psms)
        psms = scans.tag_calibrants(psms, dataset=args.dataset)
        psms = _peptide_filter(psms, args.peptides)
        psms = psms.filter(pc.is_valid(psms.column("modified_sequence")))
        if psms.num_rows == 0:
            continue
        raw_files = set(psms.column("raw_file").to_pylist())
        for raw_file in raw_files:
            sub = psms.filter(pc.equal(psms.column("raw_file"), raw_file))
            peaks = _bundle_peaks(droot, args.dataset, raw_file)
            if peaks is None:
                print(f"  [skip] no bundle for {raw_file}", flush=True)
                continue
            trace = ms2_extract.extract_ms2_fragments(peaks, sub)
            scan_spectra = ms2_extract.trace_to_scan_spectra(trace)
            modseq = sub.column("modified_sequence").to_pylist()
            charge = sub.column("charge").to_pylist()
            mode = sub.column("mode").to_pylist()
            scan = sub.column("scan").to_pylist()
            for ms, ch, md, sc in zip(modseq, charge, mode, scan):
                spec = scan_spectra.get(sc)
                if spec is None or spec[0].numel() == 0:
                    continue
                key = (ms, int(ch), md)
                spectra_by_group[key].append(spec)
                files_by_group[key].add(raw_file)
        print(
            f"[{i + 1}/{len(searches)}] {zpath.stem}: {len(spectra_by_group)} groups so far",
            flush=True,
        )

    # per-group variance law
    chan_rows: list[dict] = []
    group_rows: list[dict] = []
    for key, spectra in spectra_by_group.items():
        modseq, charge, mode = key
        if len(files_by_group[key]) < args.min_replicates:
            continue
        try:
            cons = consensus_assembly.assemble_consensus(modseq, spectra, aggregate="median")
        except Exception as exc:  # noqa: BLE001 — skip unparseable / odd peptidoforms
            n_parse_fail += 1
            print(f"  [skip] {modseq}/{charge} {mode}: {exc}", flush=True)
            continue
        per_rep = cons.per_replicate  # (R, K)
        totals = per_rep.sum(dim=1, keepdim=True).clamp(min=1e-12)
        props = per_rep / totals  # (R, K)
        pbar = props.mean(dim=0)
        var_k = props.var(dim=0, unbiased=True)
        pq1mp = pbar * (1.0 - pbar)
        nonzero = pbar > 0
        if int(nonzero.sum()) < args.min_channels:
            continue
        slope, r2 = _fit_through_origin(pq1mp[nonzero], var_k[nonzero])
        n_eff = 1.0 / slope if slope and slope == slope and slope > 0 else float("nan")
        dev = cons.deviance_from_bulk
        b = cons.basis
        for k in range(b.K):
            if not bool(nonzero[k]):
                continue
            chan_rows.append(
                {
                    "modified_sequence": modseq,
                    "charge": charge,
                    "mode": mode,
                    "ion_type": int(b.ion_type[k]),
                    "position": int(b.position[k]),
                    "fragment_charge": int(b.charge[k]),
                    "pbar": float(pbar[k]),
                    "var": float(var_k[k]),
                    "p_one_minus_p": float(pq1mp[k]),
                    "n_replicates": cons.n_replicates,
                }
            )
        group_rows.append(
            {
                "modified_sequence": modseq,
                "charge": charge,
                "mode": mode,
                "K": int(nonzero.sum()),
                "n_replicates": cons.n_replicates,
                "n_eff": n_eff,
                "fit_r2": r2,
                "mean_total_intensity": float(totals.mean()),
                "deviance_mean": float(dev.mean()),
                "deviance_over_dof": float(dev.mean()) / max(int(nonzero.sum()) - 1, 1),
            }
        )

    chan_tbl = pa.Table.from_pylist(chan_rows)
    group_tbl = pa.Table.from_pylist(group_rows)
    pq.write_table(chan_tbl, out / "exp02_channels.parquet")
    pq.write_table(group_tbl, out / "exp02_groups.parquet")
    print(
        f"\n{group_tbl.num_rows} qualifying peptide-groups "
        f"({chan_tbl.num_rows} channel rows); {n_parse_fail} skipped.",
        flush=True,
    )
    if group_tbl.num_rows:
        med_r2 = float(pc.approximate_median(group_tbl.column("fit_r2")).as_py() or float("nan"))
        med_dof = float(
            pc.approximate_median(group_tbl.column("deviance_over_dof")).as_py() or float("nan")
        )
        print(f"median fit R² = {med_r2:.3f}; median 2N·KL/dof = {med_dof:.2f} (≈1 if multinomial)")

    _make_figures(chan_tbl, group_tbl, REPO_ROOT / "results" / "figures")
    return 0


def _make_figures(chan_tbl, group_tbl, figdir: Path) -> None:
    if group_tbl.num_rows == 0:
        print("no qualifying groups — skipping figures", flush=True)
        return
    figdir.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Fig 2.2 — mean–variance law, faceted by mode; a few top-replicate peptides.
    modes = sorted(set(chan_tbl.column("mode").to_pylist()))
    fig, axes = plt.subplots(
        1, max(len(modes), 1), figsize=(6 * max(len(modes), 1), 5), squeeze=False
    )
    ch = chan_tbl.to_pylist()
    grp = {(_g["modified_sequence"], _g["charge"], _g["mode"]): _g for _g in group_tbl.to_pylist()}
    for ax, mode in zip(axes[0], modes):
        rows = [r for r in ch if r["mode"] == mode]
        # rank peptides in this mode by replicate count; show the top few
        pep_n = defaultdict(int)
        for r in rows:
            pep_n[(r["modified_sequence"], r["charge"])] = r["n_replicates"]
        top = sorted(pep_n, key=pep_n.get, reverse=True)[:6]
        for pep in top:
            pr = [r for r in rows if (r["modified_sequence"], r["charge"]) == pep]
            x = [r["p_one_minus_p"] for r in pr]
            y = [r["var"] for r in pr]
            ax.scatter(x, y, s=14, alpha=0.6, label=f"{pep[0][:10]}+{pep[1]} (R={pep_n[pep]})")
            g = grp.get((pep[0], pep[1], mode))
            if g and g["n_eff"] == g["n_eff"] and g["n_eff"] > 0:
                xs = sorted(x)
                ax.plot(xs, [xi / g["n_eff"] for xi in xs], lw=1, alpha=0.8)
        ax.set_title(f"{mode}: Var[p̂_k] vs p̄_k(1−p̄_k)")
        ax.set_xlabel("p̄_k (1 − p̄_k)")
        ax.set_ylabel("Var across replicates")
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(figdir / "exp02_fig2_2_mean_variance_law.png", dpi=140)
    print(f"wrote {figdir / 'exp02_fig2_2_mean_variance_law.png'}", flush=True)

    # Fig 2.3 — fitted N_eff vs mean total intensity (the bridge).
    g = group_tbl.to_pylist()
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    for mode in modes:
        gm = [r for r in g if r["mode"] == mode and r["n_eff"] == r["n_eff"] and r["n_eff"] > 0]
        ax2.scatter(
            [r["mean_total_intensity"] for r in gm],
            [r["n_eff"] for r in gm],
            s=20,
            alpha=0.6,
            label=mode,
        )
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("mean total fragment intensity (∝ N)")
    ax2.set_ylabel("fitted N_eff")
    ax2.set_title("N_eff vs intensity — the bridge to MS1 N")
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(figdir / "exp02_fig2_3_neff_vs_intensity.png", dpi=140)
    print(f"wrote {figdir / 'exp02_fig2_3_neff_vs_intensity.png'}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
