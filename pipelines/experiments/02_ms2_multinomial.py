#!/usr/bin/env python
"""exp02 — Is MS2 fragmentation noise multinomial or Gaussian? (shot-noise scaling)

A generative-model comparison via a **nearest-neighbour-in-intensity** test on
calibrant peptides, which sidesteps absolute N (and the gain α) entirely.

For one (peptide, charge, mode), pool every MS2 spectrum across runs, normalize
each to a fragment-proportion vector p̂ with a total intensity I (∝ the ion count
N, up to the unknown gain α). Sort by I and compare **adjacent** spectra — they
are quasi-N-matched (tiny ΔI when thousands of spectra span the range), so their
difference isolates *sampling noise at that intensity*. The squared proportion
difference scales as

    E[‖Δp̂‖²]  ≈  A · I^(−β)  +  B ,   β = 1  (multinomial / Poisson shot noise)
                                        β = 2  (additive Gaussian — interference / readout)

with B a constant floor (real run-to-run drift / residual interference). The
exponent β is **dimensionless** — it needs only relative N (intensity ratios),
never α — so this is the α-free version of the shot-noise test. On clean,
interference-free calibrants we expect β ≈ 1 with small B: the multinomial
signature, and the reason KL (the multinomial MLE) is the right comparator while
cosine/L2 (the Gaussian MLE) over-weights intense-ion shot noise. Cf. Du et al.
2008 (Var = N·p + N·p(1−p); ion-trap Poisson-limited).

All statistics here are **descriptive** (binned medians + a log-log slope, with
β=1/β=2 as reference lines) — no generative/likelihood fit lives in protostar.
The driver saves a per-spectrum parquet (intensity + aligned proportions) and a
nearest-neighbour-pairs parquet for offline re-analysis, plus first-look figures.

Runs on a compute node. Example:
    python pipelines/experiments/02_ms2_multinomial.py --dataset Zolg2017 \
        --injection DDA --max-acquisitions 250 --peptides qc
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

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import torch

from constellation.massspec.spectra.consensus import align_to_basis
from protostar.experiments import consensus_assembly, ms2_extract, scans

REPO_ROOT = Path(__file__).resolve().parents[2]
_EPS = 1e-12


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    add_common_args(p)
    p.add_argument("--dataset", default="Zolg2017")
    p.add_argument("--injection", default="DDA", help="injection token selecting search zips")
    p.add_argument("--peptides", choices=("qc", "procal", "all"), default="qc")
    p.add_argument("--max-acquisitions", type=int, default=250)
    p.add_argument("--min-spectra", type=int, default=40, help="min spectra per group to analyze")
    p.add_argument("--rel-di-max", type=float, default=0.15, help="max ΔI/I for an N-matched pair")
    p.add_argument("--n-bins", type=int, default=18, help="log-intensity bins for the median curve")
    p.add_argument("--out", type=Path, default=None)
    return p


def _select_searches(droot, dataset, injection, limit):
    zips = sorted((droot / "search" / dataset).glob(f"*{injection}*.zip"))
    return zips[:limit]


def _bundle_peaks(droot, dataset, raw_file):
    bundle = droot / "proc" / dataset / "centroid" / raw_file / "peaks.parquet"
    return pq.read_table(bundle) if bundle.is_file() else None


def _peptide_filter(psms, scope):
    if scope == "all":
        return psms
    col = "is_procal" if scope == "procal" else "is_qc"
    return psms.filter(pc.fill_null(psms.column(col), False))


def _aligned_matrix(modseq, spectra):
    """Stack replicate (mz_theoretical, intensity) spectra onto the peptide's
    fixed fragment basis → (R, K) tensor, plus the basis."""
    basis = consensus_assembly.basis_for(modseq)
    rows = [
        align_to_basis(mz, inten, basis, tolerance=20.0, tolerance_unit="ppm")
        for mz, inten in spectra
    ]
    return torch.stack(rows), basis


def _nn_pairs(props: torch.Tensor, totals: torch.Tensor, *, rel_di_max: float):
    """Adjacent-in-intensity pairs (quasi-N-matched): returns (pair_intensity,
    squared proportion distance ‖Δp̂‖²) for pairs with ΔI/I ≤ rel_di_max."""
    order = torch.argsort(totals)
    p, t = props[order], totals[order]
    dp2 = ((p[1:] - p[:-1]) ** 2).sum(dim=1)
    pair_i = 0.5 * (t[1:] + t[:-1])
    di = (t[1:] - t[:-1]).abs()
    keep = di <= rel_di_max * pair_i.clamp(min=_EPS)
    return pair_i[keep], dp2[keep]


def _binned_slope(pair_i: torch.Tensor, dp2: torch.Tensor, *, n_bins: int):
    """Median ‖Δp̂‖² per log-intensity bin, the log-log slope over all bins, and
    the slope over the lower-intensity (shot-noise-dominated) half — plus a floor
    estimate B (median of the top bin). Descriptive only."""
    if pair_i.numel() < 8:
        return None
    li = pair_i.log10()
    edges = torch.linspace(li.min(), li.max(), n_bins + 1)
    xs, ys = [], []
    for b in range(n_bins):
        hi_inc = b == n_bins - 1
        m = (li >= edges[b]) & ((li <= edges[b + 1]) if hi_inc else (li < edges[b + 1]))
        if int(m.sum()) >= 3:
            xs.append(float(pair_i[m].median()))
            ys.append(float(dp2[m].median()))
    if len(xs) < 5:
        return None
    lx, ly = np.log10(np.array(xs)), np.log10(np.array(ys))
    slope_all = float(np.polyfit(lx, ly, 1)[0])
    half = max(3, len(xs) // 2)
    slope_low = float(np.polyfit(lx[:half], ly[:half], 1)[0])  # shot-noise regime
    floor = float(np.median(np.array(ys)[-2:]))  # high-I plateau ~ drift/interference
    return {
        "slope_all": slope_all,
        "slope_low": slope_low,
        "floor": floor,
        "bins": list(zip(xs, ys)),
    }


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    droot = data_root(config, args.data_root)
    out = args.out or (droot / "exp" / "02")
    out.mkdir(parents=True, exist_ok=True)

    searches = _select_searches(droot, args.dataset, args.injection, args.max_acquisitions)
    print(f"{len(searches)} {args.injection} acquisitions for {args.dataset}", flush=True)

    spectra_by_group: dict[tuple, list] = defaultdict(list)
    for i, zpath in enumerate(searches):
        psms = scans.load_psms(zpath, extract_dir=out / "_extract" / zpath.stem)
        psms = scans.tag_calibrants(scans.gate_psms(psms), dataset=args.dataset)
        psms = _peptide_filter(psms, args.peptides)
        psms = psms.filter(pc.is_valid(psms.column("modified_sequence")))
        if psms.num_rows == 0:
            continue
        for raw_file in set(psms.column("raw_file").to_pylist()):
            sub = psms.filter(pc.equal(psms.column("raw_file"), raw_file))
            peaks = _bundle_peaks(droot, args.dataset, raw_file)
            if peaks is None:
                continue
            trace = ms2_extract.extract_ms2_fragments(peaks, sub)
            scan_spectra = ms2_extract.trace_to_scan_spectra(trace)
            for ms, ch, md, sc in zip(
                sub.column("modified_sequence").to_pylist(),
                sub.column("charge").to_pylist(),
                sub.column("mode").to_pylist(),
                sub.column("scan").to_pylist(),
            ):
                spec = scan_spectra.get(sc)
                if spec is not None and spec[0].numel() > 0:
                    spectra_by_group[(ms, int(ch), md)].append(spec)
        if (i + 1) % 25 == 0 or i + 1 == len(searches):
            print(f"[{i + 1}/{len(searches)}] {len(spectra_by_group)} groups", flush=True)

    # per-group: align → proportions → nearest-neighbour scaling
    spec_rows, pair_rows, scal_rows = [], [], []
    for (modseq, charge, mode), spectra in spectra_by_group.items():
        if len(spectra) < args.min_spectra:
            continue
        try:
            mat, _ = _aligned_matrix(modseq, spectra)
        except Exception as exc:  # noqa: BLE001 — skip odd peptidoforms
            print(f"  [skip] {modseq}/{charge} {mode}: {exc}", flush=True)
            continue
        totals = mat.sum(dim=1)
        ok = totals > 0
        mat, totals = mat[ok], totals[ok]
        if totals.numel() < args.min_spectra:
            continue
        props = mat / totals[:, None]
        for j in range(totals.numel()):
            spec_rows.append(
                {
                    "modified_sequence": modseq,
                    "charge": charge,
                    "mode": mode,
                    "total_intensity": float(totals[j]),
                    "props": props[j].tolist(),
                }
            )
        pair_i, dp2 = _nn_pairs(props, totals, rel_di_max=args.rel_di_max)
        for j in range(pair_i.numel()):
            pair_rows.append(
                {
                    "modified_sequence": modseq,
                    "charge": charge,
                    "mode": mode,
                    "pair_intensity": float(pair_i[j]),
                    "dp2": float(dp2[j]),
                }
            )
        sc = _binned_slope(pair_i, dp2, n_bins=args.n_bins)
        if sc is not None:
            scal_rows.append(
                {
                    "modified_sequence": modseq,
                    "charge": charge,
                    "mode": mode,
                    "n_spectra": int(totals.numel()),
                    "n_pairs": int(pair_i.numel()),
                    "slope_all": sc["slope_all"],
                    "slope_low": sc["slope_low"],
                    "floor": sc["floor"],
                }
            )

    pq.write_table(pa.Table.from_pylist(spec_rows), out / "exp02_spectra.parquet")
    pq.write_table(pa.Table.from_pylist(pair_rows), out / "exp02_pairs.parquet")
    scal_tbl = pa.Table.from_pylist(scal_rows)
    pq.write_table(scal_tbl, out / "exp02_scaling.parquet")
    print(
        f"\n{scal_tbl.num_rows} peptide-groups; {len(spec_rows)} spectra; {len(pair_rows)} NN pairs",
        flush=True,
    )
    if scal_tbl.num_rows:
        for m in sorted(set(scal_tbl.column("mode").to_pylist())):
            sm = scal_tbl.filter(pc.equal(scal_tbl.column("mode"), m))
            lo = pc.approximate_median(sm.column("slope_low")).as_py()
            al = pc.approximate_median(sm.column("slope_all")).as_py()
            print(
                f"  {m}: {sm.num_rows} peptides, median shot-noise slope = {lo:+.2f} "
                f"(−1 multinomial / −2 Gaussian); slope(all bins) = {al:+.2f}"
            )
    _figures(scal_tbl, pair_rows, REPO_ROOT / "results" / "figures")
    return 0


def _figures(scal_tbl, pair_rows, figdir: Path) -> None:
    if scal_tbl.num_rows == 0:
        print("no groups — skipping figures", flush=True)
        return
    figdir.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    modes = sorted(set(scal_tbl.column("mode").to_pylist()))
    by_pep = defaultdict(list)
    for r in pair_rows:
        by_pep[(r["modified_sequence"], r["charge"], r["mode"])].append(
            (r["pair_intensity"], r["dp2"])
        )

    # Fig A — the "continuous line": binned ‖Δp̂‖² vs intensity for top peptides, with β=1/β=2 refs.
    fig, axes = plt.subplots(1, len(modes), figsize=(6.5 * len(modes), 5.2), squeeze=False)
    for ax, mode in zip(axes[0], modes):
        peps = sorted([k for k in by_pep if k[2] == mode], key=lambda k: -len(by_pep[k]))[:6]
        for k in peps:
            pts = sorted(by_pep[k])
            xi = np.array([p[0] for p in pts])
            yi = np.array([p[1] for p in pts])
            # bin medians
            edges = np.quantile(np.log10(xi), np.linspace(0, 1, 16))
            bx, byy = [], []
            for b in range(len(edges) - 1):
                m = (np.log10(xi) >= edges[b]) & (np.log10(xi) <= edges[b + 1])
                if m.sum() >= 5:
                    bx.append(np.median(xi[m]))
                    byy.append(np.median(yi[m]))
            if len(bx) >= 4:
                ax.plot(
                    bx, byy, "-o", ms=3, lw=1, alpha=0.7, label=f"{k[0][:9]}+{k[1]} (n={len(pts)})"
                )
        # reference slopes anchored to the panel
        xl = np.array(ax.get_xlim())
        if xl[0] > 0:
            x0 = np.array([xi.min(), xi.max()])
            yl = ax.get_ylim()
            anchor = np.median([p[1] for k in peps for p in by_pep[k]])
            ax.plot(
                x0,
                anchor * (x0 / np.median(x0)) ** -1.0,
                "k-",
                lw=1.4,
                alpha=0.7,
                label="β=1 (multinomial)",
            )
            ax.plot(
                x0,
                anchor * (x0 / np.median(x0)) ** -2.0,
                "k--",
                lw=1.0,
                alpha=0.5,
                label="β=2 (Gaussian)",
            )
            ax.set_ylim(yl)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("pair intensity  I  (∝ N)")
        ax.set_ylabel("‖Δp̂‖²  between adjacent-intensity spectra")
        ax.set_title(f"{mode}: shot-noise scaling")
        ax.legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(figdir / "exp02_fig_shot_noise_scaling.png", dpi=145)
    print(f"wrote {figdir / 'exp02_fig_shot_noise_scaling.png'}", flush=True)

    # Fig B — distribution of per-peptide shot-noise slopes, by mode, vs β=1/β=2.
    fig2, ax2 = plt.subplots(figsize=(7, 5))
    for m in modes:
        sl = scal_tbl.filter(pc.equal(scal_tbl.column("mode"), m)).column("slope_low").to_pylist()
        ax2.hist(
            [s for s in sl if s == s],
            bins=24,
            range=(-3, 0.5),
            alpha=0.5,
            label=f"{m} (n={len(sl)})",
        )
    ax2.axvline(-1, color="k", lw=1.4, label="β=1 multinomial")
    ax2.axvline(-2, color="k", ls="--", lw=1.0, label="β=2 Gaussian")
    ax2.set_xlabel("per-peptide shot-noise scaling slope  (d log‖Δp̂‖² / d log I)")
    ax2.set_ylabel("peptides")
    ax2.set_title("Generative-model signature across calibrant peptides")
    ax2.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(figdir / "exp02_fig_slope_distribution.png", dpi=145)
    print(f"wrote {figdir / 'exp02_fig_slope_distribution.png'}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
