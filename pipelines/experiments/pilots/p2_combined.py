#!/usr/bin/env python
"""Pilot P2, third step: a combined screen score and the output volume it allows.

For each setting whose candidate regions ``p2_analyze.py`` saved (``p2_candidates/<setting>/``):

- **Labels.** Positive regions are in-pool regions containing the unmodified peptide's MaxQuant
  PSM time ±30 s; negative regions are those of clean off-pool sequences.
- **Score.** A logistic score over the standardised region features, fitted by IRLS and
  cross-validated with whole runs held out (5 folds), reporting AUC and weights.
- **Volume.** Off-pool regions kept, and rows per run on the full screen (``--full-targets``), at
  thresholds that keep 99/98/95/90% of true regions.

Pilot analysis (edge code). Writes ``p2_combined.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _common  # noqa: E402,F401

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from protostar import pools  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p3_rt_transfer import best_rt  # noqa: E402

from protostar.experiments.scans import gate_psms, load_psms  # noqa: E402

FEATURES = ["log_flux_sum", "n_qual", "abs_median_e", "median_kl", "median_cos", "max_cos", "max_m1_ratio",
            "mad_e", "median_spread"]  # fmt: skip


def irls(x: np.ndarray, y: np.ndarray, iters: int = 25, lam: float = 1e-3) -> np.ndarray:
    xb = np.c_[np.ones(len(x)), x]
    w = np.zeros(xb.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-xb @ w))
        wt = p * (1 - p) + 1e-9
        h = xb.T @ (xb * wt[:, None]) + lam * np.eye(xb.shape[1])
        w += np.linalg.solve(h, xb.T @ (y - p) - lam * w)
    return w


def auc(s: np.ndarray, y: np.ndarray) -> float:
    r = pd.Series(s).rank().to_numpy()
    n1 = int(y.sum())
    n0 = len(y) - n1
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--rows-dir", required=True, type=Path, help="p2_screen_rows.py --out")
    ap.add_argument("--space-cache", type=Path, required=True)
    ap.add_argument("--full-targets", type=float, default=2_476_934)
    ap.add_argument("--max-pep", type=float, default=0.01)
    args = ap.parse_args()

    space = pq.read_table(args.space_cache, columns=["residues"])
    rmap = {s: i for i, s in enumerate(space.column("residues").to_pylist())}
    runs = {r["stem"]: r for r in pools.run_table(args.data_root).to_pylist()}
    analysis = json.loads((args.rows_dir / "p2_analysis.json").read_text())
    tmp = Path(tempfile.mkdtemp(prefix="p2comb-"))
    truth: dict[str, dict[int, float]] = {}
    out: dict = {}
    for sdir in sorted(p for p in (args.rows_dir / "p2_candidates").iterdir() if p.is_dir()):
        frames = []
        for f in sorted(sdir.glob("*.parquet")):
            stem = f.stem
            if stem not in truth:
                z = pools.search_zip_for_stem(args.data_root, runs[stem]["dataset"], stem)
                t = best_rt(
                    gate_psms(load_psms(z, extract_dir=tmp / stem), max_pep=args.max_pep), stem
                )
                truth[stem] = {rmap[s]: rt for s, rt in t.items() if s in rmap}
            reg = pq.read_table(f).to_pandas()
            reg["psm_rt"] = reg["sid"].map(truth[stem])
            reg["pos"] = (
                (reg["label"] == "in_pool")
                & (reg["rt_left"] - 30 <= reg["psm_rt"])
                & (reg["rt_right"] + 30 >= reg["psm_rt"])
            )
            reg["neg"] = reg["label"] == "off_pool"
            reg["run"] = stem
            frames.append(reg[[*FEATURES, "pos", "neg", "run", "target_id"]])
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        df["median_kl"] = df["median_kl"].clip(upper=5)
        x = np.nan_to_num(df[FEATURES].to_numpy(float), nan=0.0)
        xs = (x - x.mean(0)) / (x.std(0) + 1e-9)
        y = df["pos"].to_numpy()
        keep = y | df["neg"].to_numpy()
        folds = df["run"].astype("category").cat.codes.to_numpy() % 5
        score = np.full(len(df), np.nan)
        for k in range(5):
            tr, te = keep & (folds != k), folds == k
            w = irls(xs[tr], y[tr].astype(float))
            score[te] = np.c_[np.ones(te.sum()), xs[te]] @ w
        w_all = irls(xs[keep], y[keep].astype(float))
        name = sdir.name
        base = analysis["settings"].get(name, {}).get("null_regions_per_query")
        rec = {"n_pos": int(y.sum()), "n_neg": int(df["neg"].sum()), "cv_auc": round(auc(score[keep], y[keep]), 4),
               "weights_standardised": dict(zip(["intercept", *FEATURES], np.round(w_all, 3).tolist(), strict=True)),
               "recall_in_pool": analysis["settings"].get(name, {}).get("recall_in_pool"),
               "null_regions_per_query": base, "volume": []}  # fmt: skip
        null = df["neg"].to_numpy()
        for keep_frac in (0.99, 0.98, 0.95, 0.90):
            thr = float(np.quantile(score[y], 1 - keep_frac))
            kept = float((score[null] >= thr).mean())
            rec["volume"].append({"true_regions_kept": keep_frac, "threshold": round(thr, 3), "off_pool_regions_kept": round(kept, 4),
                                  "rows_per_run": None if base is None else round(kept * base * args.full_targets)})  # fmt: skip
        out[name] = rec
        print(
            name,
            json.dumps({k: rec[k] for k in ("cv_auc", "recall_in_pool", "null_regions_per_query")}),
            flush=True,
        )
    (args.rows_dir / "p2_combined.json").write_text(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
