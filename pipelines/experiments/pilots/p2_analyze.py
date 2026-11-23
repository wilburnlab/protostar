#!/usr/bin/env python
"""Pilot P2, analysis half: screen features, regions and thresholds from the isotope rows.

Reads ``p2_screen_rows.py``'s output and, per run, computes features for each row (a query in one
MS1 scan with ≥ 2 of M0..M+2 matched):

==========  ======================================================================================
feature     definition
==========  ======================================================================================
kl          KL(observed ‖ theoretical) over the matched isotopes among M0..M+2, both renormalised
cos         cosine of the same two vectors
m1_ratio    I(M−1) / max I(M0..M+2): a large value means the target is another species' M+1
e_mean      mean ppm error over the matched isotopes, minus the run's PROCAL median offset
e_spread    max − min ppm error over the matched isotopes
flux        Σ I·iit (iit in s): an ion-count proxy up to α(z)
dev         multinomial deviance G² at N = flux / (10·z) (seconds convention; a rough α)
==========  ======================================================================================

Then, for each setting of the qualifying thresholds (kl_max, e_max, gap, m_min), it merges
qualifying scans into regions along the run's MS1 scan order, keeps the top 5 regions per query,
and joins the run's labels and MaxQuant identifications (PEP ≤ 0.01) to report:

- **recall**: identified in-pool sequences with a region containing their PSM RT (±30 s), and the
  same for PROCAL;
- **separation**: per-feature AUC of true in-pool regions (containing the PSM RT) against regions
  of clean off-pool sequences;
- **volume**: region rows per query and per run, extrapolated to the full 2.7M-target screen;
- **carryover**: the rate of a region near the sequence's own-pool RT, by pools back;
- **shadowed**: regions of shadowed sequences, and the fraction within ±30 s of the in-pool RT of
  the sequence they shadow.

Pilot analysis (edge code; pandas allowed). Writes ``p2_analysis.json`` and
``p2_candidates/<setting>.parquet``.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _common  # noqa: E402,F401

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from protostar import bundles, pools  # noqa: E402
from protostar.experiments.scans import gate_psms, load_psms  # noqa: E402

ALPHA_PER_CHARGE = 10.0  # seconds convention; rough (Hsu et al. 2025; nb42b)
TOP_K = 5
FEATURES = ["n_qual", "median_kl", "max_cos", "median_cos", "max_m1_ratio", "abs_median_e", "mad_e",
            "median_spread", "log_flux_sum", "median_dev"]  # fmt: skip


def row_features(
    rows: pd.DataFrame, chans: pd.DataFrame, charge: pd.Series, offset_ppm: float
) -> pd.DataFrame:
    p = chans.pivot_table(index="target_id", columns="isotope", values="abundance", aggfunc="first")
    p = p.reindex(columns=[0, 1, 2]).reindex(rows["target_id"]).to_numpy(dtype=float)
    i = rows[["i_0", "i_1", "i_2"]].to_numpy(dtype=float)
    e = rows[["e_0", "e_1", "e_2"]].to_numpy(dtype=float)
    m = i > 0
    pm = np.where(m, p, 0.0)
    pm = pm / np.clip(pm.sum(1, keepdims=True), 1e-12, None)
    om = i / np.clip(i.sum(1, keepdims=True), 1e-12, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl = np.where(m, om * np.log(np.where(m, om, 1) / np.where(m, pm, 1)), 0.0).sum(1)
        cos = (om * pm).sum(1) / np.sqrt((om**2).sum(1) * (pm**2).sum(1))
    imax = i.max(1)
    em = np.where(m, e, np.nan)
    out = pd.DataFrame(
        {
            "target_id": rows["target_id"].to_numpy(),
            "scan": rows["scan"].to_numpy(),
            "rt": rows["rt"].to_numpy(),
            "kl": kl,
            "cos": cos,
            "m1_ratio": rows["i_m1"].to_numpy() / np.clip(imax, 1e-12, None),
            "e_mean": np.nanmean(em, 1) - offset_ppm,
            "e_spread": np.nanmax(em, 1) - np.nanmin(em, 1),
            "flux": i.sum(1) * rows["iit_ms"].to_numpy() / 1000.0,
        }
    )
    z = charge.reindex(out["target_id"]).to_numpy(dtype=float)
    n_obs = i * (rows["iit_ms"].to_numpy()[:, None] / 1000.0) / (ALPHA_PER_CHARGE * z[:, None])
    n_tot = n_obs.sum(1, keepdims=True)
    n_exp = n_tot * pm
    with np.errstate(divide="ignore", invalid="ignore"):
        out["dev"] = 2 * np.where(
            m & (n_obs > 0), n_obs * np.log(n_obs / np.where(n_exp > 0, n_exp, 1)), 0.0
        ).sum(1)
    return out


def regions(
    feat: pd.DataFrame, ordinal: dict, kl_max: float, e_max: float, gap: int
) -> pd.DataFrame:
    """Every region of qualifying scans (gaps of at most ``gap`` MS1 scans), with its statistics."""
    q = feat[(feat["kl"] <= kl_max) & (feat["e_mean"].abs() <= e_max)].copy()
    if q.empty:
        return pd.DataFrame()
    q["ord"] = q["scan"].map(ordinal)
    q = q.dropna(subset=["ord"]).sort_values(["target_id", "ord"])
    q["ord"] = q["ord"].astype(int)
    new = (q["target_id"] != q["target_id"].shift()) | (q["ord"] - q["ord"].shift() > gap + 1)
    q["region"] = new.cumsum()
    q["w"] = q["flux"] * np.exp(-q["kl"])
    q["abs_dev_e"] = (q["e_mean"] - q.groupby("region")["e_mean"].transform("median")).abs()
    g = q.groupby("region")
    agg = g.agg(
        target_id=("target_id", "first"),
        n_qual=("scan", "size"),
        ord_left=("ord", "min"),
        ord_right=("ord", "max"),
        rt_left=("rt", "min"),
        rt_right=("rt", "max"),
        median_kl=("kl", "median"),
        max_cos=("cos", "max"),
        median_cos=("cos", "median"),
        max_m1_ratio=("m1_ratio", "max"),
        median_e=("e_mean", "median"),
        mad_e=("abs_dev_e", "median"),
        median_spread=("e_spread", "median"),
        flux_sum=("flux", "sum"),
        median_dev=("dev", "median"),
    )
    agg["apex_rt"] = q.loc[g["w"].idxmax(), ["region", "rt"]].set_index("region")["rt"]
    agg["abs_median_e"] = agg["median_e"].abs()
    agg["log_flux_sum"] = np.log10(agg["flux_sum"].clip(lower=1e-12))
    agg["score"] = agg["n_qual"] * np.exp(-agg["median_kl"]) * agg["log_flux_sum"].clip(lower=0)
    return agg.reset_index(drop=True)


def top_regions(agg: pd.DataFrame, m_min: int) -> pd.DataFrame:
    """Regions with at least ``m_min`` qualifying scans, the best ``TOP_K`` per query."""
    r = agg[agg["n_qual"] >= m_min].sort_values(["target_id", "score"], ascending=[True, False])
    r = r.assign(rank=r.groupby("target_id").cumcount() + 1)
    return r[r["rank"] <= TOP_K].reset_index(drop=True)


def auc(pos: np.ndarray, neg: np.ndarray) -> float | None:
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if pos.size == 0 or neg.size == 0:
        return None
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank().to_numpy()
    return float((ranks[: pos.size].sum() - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--rows-dir", required=True, type=Path, help="p2_screen_rows.py --out")
    ap.add_argument("--space-cache", type=Path, required=True)
    ap.add_argument("--kl-max", nargs="+", type=float, default=[0.1, 0.5, np.inf])
    ap.add_argument("--e-max", nargs="+", type=float, default=[10.0, 20.0])
    ap.add_argument("--gap", nargs="+", type=int, default=[2])
    ap.add_argument("--m-min", nargs="+", type=int, default=[2, 3, 5])
    ap.add_argument(
        "--full-targets", type=float, default=2.70e6, help="targets per run in the full screen"
    )
    ap.add_argument("--max-pep", type=float, default=0.01)
    args = ap.parse_args()

    out_dir = args.rows_dir
    (out_dir / "p2_candidates").mkdir(exist_ok=True)
    space = pq.read_table(args.space_cache)
    seq_res = space.column("residues").to_numpy(zero_copy_only=False)
    targets = pq.read_table(out_dir / "targets.parquet").to_pandas()
    charge = targets.set_index("target_id")["precursor_charge"]
    seq_of_target = targets.set_index("target_id")["peptide_id"]
    runs = pools.run_table(args.data_root)
    rowsmeta = {r["stem"]: r for r in runs.to_pylist()}
    procal_ids = set(
        np.flatnonzero(space.column("is_procal").to_numpy(zero_copy_only=False)).tolist()
    )
    settings = list(product(args.kl_max, args.e_max, args.gap, args.m_min))
    report: dict = {"provenance": bundles.constellation_provenance(), "settings": {}, "runs": {}}
    stems = sorted(p.stem for p in (out_dir / "rows").glob("*.parquet"))
    # identifications per run (truth for recall), and each sequence's own-pool RT (for carryover/shadow)
    tmp = Path(tempfile.mkdtemp(prefix="p2psm-"))
    psm_rt: dict[str, dict[int, float]] = {}
    for stem in stems:
        z = pools.search_zip_for_stem(args.data_root, rowsmeta[stem]["dataset"], stem)
        if z is None:
            continue
        psms = gate_psms(load_psms(z, extract_dir=tmp / stem), max_pep=args.max_pep).to_pandas()
        psms = psms[psms["raw_file"].str.lower() == stem.lower()].sort_values(
            ["pep", "score"], ascending=[True, False]
        )
        best = psms.drop_duplicates("sequence")
        rmap = {s: i for i, s in enumerate(seq_res)}
        psm_rt[stem] = {
            rmap[s]: float(rt)
            for s, rt in zip(best["sequence"], best["retention_time_s"])
            if s in rmap
        }
    own_rt: dict[int, list[float]] = {}
    for stem, d in psm_rt.items():
        for sid, rt in d.items():
            own_rt.setdefault(sid, []).append(rt)
    own_med = {k: float(np.median(v)) for k, v in own_rt.items()}

    per_setting: dict[tuple, dict] = {s: {"pos": [], "neg": [], "recall_hit": 0, "recall_n": 0, "procal_hit": 0,
                                          "procal_n": 0, "rows": 0, "null_rows": 0, "null_queries": 0,
                                          "carry": {}, "shadow_rt": [0, 0]} for s in settings}  # fmt: skip
    for stem in stems:
        rows = pq.read_table(out_dir / "rows" / f"{stem}.parquet").to_pandas()
        chans = (
            pq.read_table(out_dir / "channels" / f"{stem}.parquet")
            .to_pandas()
            .drop_duplicates(["target_id", "isotope"])
        )
        labels = (
            pq.read_table(out_dir / "labels" / f"{stem}.parquet")
            .to_pandas()
            .set_index("sequence_id")
        )
        bundle = bundles.require_run_bundle(args.data_root, rowsmeta[stem]["dataset"], stem)
        sm = bundles.read_scan_metadata(bundle, columns=["scan", "level"]).to_pandas()
        ms1 = np.sort(sm.loc[sm["level"] == 1, "scan"].to_numpy())
        ordinal = {int(s): k for k, s in enumerate(ms1)}
        rows["sid"] = rows["target_id"].map(seq_of_target)
        pro = rows[rows["sid"].isin(procal_ids) & (rows["n_iso"] == 3)]
        offset = float(np.nanmedian(pro["e_0"])) if len(pro) else 0.0
        feat = row_features(rows, chans, charge, offset)
        feat["sid"] = feat["target_id"].map(seq_of_target)
        n_queries = rows["target_id"].nunique()
        lab = feat["sid"].map(labels["label"]).fillna("off_pool")
        run_q = chans["target_id"].nunique()
        null_targets = set(chans["target_id"]) - set(
            targets[targets["peptide_id"].isin(labels.index)]["target_id"]
        )  # queries labelled off_pool in this run
        report["runs"][stem] = {"rows": len(rows), "queries_with_rows": int(n_queries), "queries": int(run_q),
                                "ppm_offset": round(offset, 2), "rows_by_label": lab.value_counts().to_dict()}  # fmt: skip
        truth = psm_rt.get(stem, {})
        truth_df = pd.DataFrame({"sid": list(truth), "psm_rt": list(truth.values())})
        truth_df["label"] = truth_df["sid"].map(labels["label"]).fillna("off_pool")
        truth_df["is_procal"] = truth_df["sid"].isin(procal_ids)
        for kl_max, e_max, gap in product(args.kl_max, args.e_max, args.gap):
            base = regions(feat, ordinal, kl_max, e_max, gap)
            for m_min in args.m_min:
                s = (kl_max, e_max, gap, m_min)
                acc = per_setting[s]
                reg = top_regions(base, m_min) if len(base) else base
                acc["rows"] += len(reg)
                if reg.empty:
                    continue
                reg["sid"] = reg["target_id"].map(seq_of_target)
                reg["label"] = reg["sid"].map(labels["label"]).fillna("off_pool")
                # recall: identified in-pool / PROCAL sequences with a region containing the PSM RT ±30 s
                mm = reg[["sid", "rt_left", "rt_right"]].merge(truth_df, on="sid")
                mm["hit"] = (mm["rt_left"] - 30 <= mm["psm_rt"]) & (
                    mm["rt_right"] + 30 >= mm["psm_rt"]
                )
                hit = mm.groupby("sid")["hit"].any()
                tin = truth_df[truth_df["label"] == "in_pool"]["sid"]
                tpro = truth_df[truth_df["is_procal"]]["sid"]
                acc["recall_n"] += len(tin)
                acc["recall_hit"] += int(hit.reindex(tin).fillna(False).sum())
                acc["procal_n"] += len(tpro)
                acc["procal_hit"] += int(hit.reindex(tpro).fillna(False).sum())
                inpool = reg[reg["label"] == "in_pool"].copy()
                inpool["psm_rt"] = inpool["sid"].map(truth)
                ok = (inpool["rt_left"] - 30 <= inpool["psm_rt"]) & (
                    inpool["rt_right"] + 30 >= inpool["psm_rt"]
                )
                null_reg = reg[reg["target_id"].isin(null_targets)]
                acc["pos"].append(inpool[ok][FEATURES])
                acc["neg"].append(null_reg[FEATURES])
                acc["null_rows"] += len(null_reg)
                acc["null_queries"] += len(null_targets)
                carry = reg[reg["label"] == "carryover"]
                if len(carry):
                    lb = labels.loc[carry["sid"], "pools_back"].to_numpy()
                    near = (carry["sid"].map(own_med) - carry["apex_rt"]).abs().to_numpy() <= 60
                    for b, n in zip(lb, near, strict=True):
                        c = acc["carry"].setdefault(int(b), [0, 0])
                        c[0] += int(n)
                        c[1] += 1
                sh = reg[reg["label"] == "shadowed"]
                if len(sh):
                    src = labels.loc[sh["sid"], "shadow_of"].to_numpy()
                    src_rt = np.array([truth.get(int(x), np.nan) for x in src])
                    coinc = np.abs(src_rt - sh["apex_rt"].to_numpy()) <= 30
                    acc["shadow_rt"][0] += int(np.nansum(coinc))
                    acc["shadow_rt"][1] += len(sh)
                if s == (
                    args.kl_max[-1],
                    args.e_max[-1],
                    args.gap[0],
                    args.m_min[0],
                ):  # loosest: keep for inspection
                    pq.write_table(pa.Table.from_pandas(reg, preserve_index=False),
                                   out_dir / "p2_candidates" / f"{stem}.parquet")  # fmt: skip
        print(stem, json.dumps(report["runs"][stem], default=str), flush=True)

    n_runs = max(len(stems), 1)
    for s, acc in per_setting.items():
        pos = pd.concat(acc["pos"]) if acc["pos"] else pd.DataFrame(columns=FEATURES)
        neg = pd.concat(acc["neg"]) if acc["neg"] else pd.DataFrame(columns=FEATURES)
        sign = {f: (-1 if f in ("median_kl", "max_m1_ratio", "abs_median_e", "mad_e", "median_spread", "median_dev") else 1)
                for f in FEATURES}  # fmt: skip
        aucs = {
            f: auc(sign[f] * pos[f].to_numpy(float), sign[f] * neg[f].to_numpy(float))
            for f in FEATURES
        }
        null_per_query = acc["null_rows"] / max(acc["null_queries"], 1)
        report["settings"]["kl{}_e{}_g{}_m{}".format(*s)] = {
            "recall_in_pool": round(acc["recall_hit"] / max(acc["recall_n"], 1), 4),
            "n_in_pool_ids": acc["recall_n"],
            "recall_procal": round(acc["procal_hit"] / max(acc["procal_n"], 1), 4),
            "regions_per_run": round(acc["rows"] / n_runs, 1),
            "null_regions_per_query": round(null_per_query, 4),
            "projected_rows_per_run_full_screen": int(null_per_query * args.full_targets),
            "auc": {k: (None if v is None else round(v, 3)) for k, v in aucs.items()},
            "n_pos": len(pos),
            "n_neg": len(neg),
            "carryover_near_own_rt_by_pools_back": {
                k: [v[0], v[1]] for k, v in sorted(acc["carry"].items())
            },
            "shadowed_regions_rt_coincident": acc["shadow_rt"],
        }
    (out_dir / "p2_analysis.json").write_text(json.dumps(report, indent=1, default=str))
    best = sorted(
        report["settings"].items(),
        key=lambda kv: (-kv[1]["recall_in_pool"], kv[1]["projected_rows_per_run_full_screen"]),
    )
    for name, r in best[:10]:
        print(
            name,
            r["recall_in_pool"],
            r["recall_procal"],
            r["projected_rows_per_run_full_screen"],
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
