"""Datasets for the extraction-pilot viewer, second round (after the PI's chromatogram review).

Reads the pilot outputs under <data_root>/pilots/extraction/ and writes JSON {rows, meta} files,
each carrying this script, to pilots/extraction/viewer2/. Truth is the unmodified peptide's best
MaxQuant PSM. Example chromatograms are extracted fresh from the bundles (rule nearest,
zero-filled) with per-scan m/z errors; identified neighbours are MaxQuant PSMs whose isotope
peaks fall within 20 ppm of an example's channels.
"""

import json
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

CODE_DIR = Path(os.environ["PILOTS_CHECKOUT"])  # a protostar checkout of the pilots branch
sys.path.insert(0, str(CODE_DIR))
sys.path.insert(0, str(CODE_DIR / "pipelines/experiments/pilots"))
from p3_rt_transfer import best_rt  # noqa: E402

from protostar import bundles, pools  # noqa: E402
from protostar import search_space as ss  # noqa: E402
from protostar.experiments.scans import gate_psms, load_psms  # noqa: E402

D = Path(os.environ["PROTOSTAR_DATA_ROOT"])
E = D / "pilots/extraction"
OUT = E / "viewer2"
OUT.mkdir(exist_ok=True)
CODE = Path(__file__).read_text()
REC = os.environ.get("RECOMMENDED", "kl0.1_e5.0_m10.3_s1_g2_m3")
BASE = "kl0.1_e10.0_m1inf_s0_g2_m3"  # the first round's rules, with the corrected truth
IN_POOL_RUN = "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2"
OFF_RUNS = ["01625b_GB1-TUM_first_pool_2_01_01-DDA-1h-R2", "01717a_BF4-TUM_second_pool_42_01_01-DDA-1h-R1",
            "01625b_GC1-TUM_first_pool_3_01_01-DDA-1h-R2"]  # fmt: skip
FULL_TARGETS = 2_476_934
ISO = (-1, 0, 1, 2, 3)
ISO_COL = {-1: "m_minus_1", 0: "m0", 1: "m1", 2: "m2", 3: "m3"}


def clean(v):
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return None if (math.isnan(v) or math.isinf(v)) else v
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def dump(name, rows, inputs):
    rows = [{k: clean(v) for k, v in r.items()} for r in rows]
    meta = [
        {"label": "viewer_data2.py", "value": CODE, "type": "code", "lang": "python"},
        {"label": "Input", "value": inputs},
    ]
    (OUT / f"{name}.json").write_text(json.dumps({"rows": rows, "meta": meta}, allow_nan=False))
    print(name, len(rows), "rows", (OUT / f"{name}.json").stat().st_size, "bytes", flush=True)


def short_run(stem):
    n = pools.parse_stem(stem)
    return f"{n.pool} {n.injection}"


def describe(setting):
    kl, e, m1, s, g, m = setting.split("_")
    return (
        f"|ppm| ≤ {float(e[1:]):g}, "
        + ("M−1 check, " if m1 != "m1inf" else "")
        + ("valley split, " if s == "s1" else "")
        + f"≥ {int(m[1:])} scans"
    )


# P2 settings, AUC, volume --------------------------------------------------------------------
p2 = json.loads((E / "p2/p2_analysis.json").read_text())
comb = json.loads((E / "p2/p2_combined.json").read_text())
rows = []
for name, v in p2["settings"].items():
    kl, e, m1, s, g, m = name.split("_")
    rows.append(
        {
            "setting": name,
            "ppm_max": float(e[1:]),
            "m1_check": "yes" if m1 != "m1inf" else "no",
            "valley_split": "yes" if s == "s1" else "no",
            "min_scans": int(m[1:]),
            "recall_in_pool": v["recall_in_pool"],
            "recall_procal": v["recall_procal"],
            "null_regions_per_query": v["null_regions_per_query"],
            "rows_per_run_full_screen": round(v["null_regions_per_query"] * FULL_TARGETS),
            "recommended": name == REC,
        }
    )
dump(
    "p2_settings",
    rows,
    "pilots/extraction/p2/p2_analysis.json (job 17122676); targets per full screen from P1",
)

LABELS = {"log_flux_sum": "Flux, log Σ I·iit", "n_qual": "Qualifying scans", "abs_median_e": "|Median ppm error|",
          "median_kl": "Median isotope KL", "median_cos": "Median isotope cosine", "max_cos": "Max isotope cosine",
          "median_spread": "ppm spread across isotopes", "mad_e": "ppm MAD across scans", "max_m1_ratio": "M−1 ratio",
          "median_dev": "Count-scaled deviance"}  # fmt: skip
rows = [
    {"feature": LABELS[k], "key": k, "auc": v, "kind": "single"}
    for k, v in p2["settings"][REC]["auc"].items()
    if v is not None
]
rows.append(
    {
        "feature": "Combined (logistic, run-held-out CV)",
        "key": "combined",
        "auc": comb[REC]["cv_auc"],
        "kind": "combined",
    }
)
rows.sort(key=lambda r: -r["auc"])
dump("p2_auc", rows, f"p2_analysis.json and p2_combined.json, setting {REC}")

rows = []
for name, rec in comb.items():
    for v in rec["volume"]:
        rows.append(
            {
                "setting": name,
                "setting_label": describe(name),
                "recommended": name == REC,
                "true_regions_kept": v["true_regions_kept"],
                "off_pool_regions_kept": v["off_pool_regions_kept"],
                "rows_per_run": v["rows_per_run"],
                "cv_auc": rec["cv_auc"],
                "recall_in_pool": rec["recall_in_pool"],
            }
        )
dump("p2_volume", rows, "pilots/extraction/p2/p2_combined.json")

# truth (unmodified PSMs) for the in-pool run and the P2 block ------------------------------------
space = pq.read_table(E / "search_space.parquet").to_pandas()
rmap = {s: i for i, s in enumerate(space["residues"])}
pid = dict(zip(space["proforma"], space["sequence_id"], strict=True))
runs = {r["stem"]: r for r in pools.run_table(D).to_pylist()}
tmp = Path(tempfile.mkdtemp())
psm_cache: dict[str, pd.DataFrame] = {}


def run_psms(stem):
    if stem not in psm_cache:
        z = pools.search_zip_for_stem(D, runs[stem]["dataset"], stem)
        t = load_psms(z, extract_dir=tmp / stem)
        psm_cache[stem] = t
    return psm_cache[stem]


def truth_of(stem):
    t = best_rt(gate_psms(run_psms(stem), max_pep=0.01), stem)
    return {rmap[s]: rt for s, rt in t.items() if s in rmap}


# P2 feature distributions from the recommended setting ------------------------------------------
frames = []
for f in sorted((E / "p2/p2_candidates" / REC).glob("*.parquet")):
    tr = truth_of(f.stem)
    reg = pq.read_table(f).to_pandas()
    reg["psm_rt"] = reg["sid"].map(tr)
    reg["cls"] = np.where((reg["label"] == "in_pool") & (reg["rt_left"] - 30 <= reg["psm_rt"]) & (reg["rt_right"] + 30 >= reg["psm_rt"]),
                          "in-pool (true)", np.where(reg["label"] == "off_pool", "off-pool", ""))  # fmt: skip
    frames.append(reg[reg["cls"] != ""])
allreg = pd.concat(frames, ignore_index=True)
rows = []
for key in (
    "log_flux_sum",
    "n_qual",
    "abs_median_e",
    "median_kl",
    "max_cos",
    "max_m1_ratio",
    "median_dev",
):
    x = allreg[key].astype(float)
    lo, hi = np.nanquantile(x, 0.005), np.nanquantile(x, 0.995)
    if key == "n_qual":
        lo, hi = 2, min(hi, 60)
    if hi <= lo:
        hi = lo + 1
    edges = np.linspace(lo, hi, 31)
    for cls, g in allreg.groupby("cls"):
        v = np.clip(g[key].astype(float).dropna().to_numpy(), lo, hi)
        h, _ = np.histogram(v, bins=edges)
        share = h / max(h.sum(), 1)
        for i in range(len(h)):
            rows.append({"feature": LABELS[key], "key": key, "class": cls, "bin_lo": float(edges[i]), "bin_hi": float(edges[i + 1]),
                         "share": float(share[i]), "n": int(h[i])})  # fmt: skip
dump(
    "p2_feature_dist",
    rows,
    f"p2/p2_candidates/{REC}/*.parquet + unmodified MaxQuant PSMs (PEP <= 0.01)",
)

# P3 ------------------------------------------------------------------------------------------
REG = {"within_pool": "Within pool", "cross_pool": "Cross pool", "cross_dataset": "Cross dataset"}
res = pq.read_table(E / "p3/p3_residuals.parquet").to_pandas()
edges = np.arange(-400, 401, 10)
rows = []
for regime, g in res.groupby("regime"):
    for method in ("procal", "refined"):
        v = np.clip(g[f"resid_{method}"].dropna().to_numpy(), -400, 400)
        h, _ = np.histogram(v, bins=edges)
        for i in range(len(h)):
            rows.append({"regime": REG[regime], "map": "PROCAL only" if method == "procal" else "Refined", "bin_lo": int(edges[i]),
                         "bin_hi": int(edges[i + 1]), "share": float(h[i] / max(h.sum(), 1)), "n": int(h[i])})  # fmt: skip
dump(
    "p3_hist",
    rows,
    "pilots/extraction/p3/p3_residuals.parquet (job 17122677, unmodified truth); clipped to ±400 s",
)
rows = []
for truth_name, path in (
    ("Unmodified peptide", E / "p3/p3_report.json"),
    ("Any form (first round)", E / "p3/p3_report_v1.json"),
):
    rep = json.loads(path.read_text())
    for k, v in rep["summary"].items():
        parts = k.split("/")
        if len(parts) != 2:
            continue
        rows.append({"truth": truth_name, "regime": REG[parts[0]], "map": "PROCAL only" if parts[1] == "procal" else "Refined",
                     "n": v["n"], "median_s": v["median"], "sd_mad_s": v["sd_mad"], "p95_abs_s": v["p95_abs"],
                     "p99_abs_s": v["p99_abs"], "over_300s": v["frac_abs_gt_300s"]})  # fmt: skip
dump("p3_stats", rows, "p3_report.json (job 17122677) and p3_report_v1.json (job 17104327)")

# Examples ------------------------------------------------------------------------------------
truth = truth_of(IN_POOL_RUN)
first_round = [("ALPTLDPK", 2, "Bright, found"), ("SFPLFNR", 2, "Bright, found"), ("IEIIATLK", 2, "Median, found"),
               ("YTASDLFR", 2, "Median, found"), ("MPEFYNR", 2, "Weak, found"), ("EIVEEAIR", 2, "Weak, found"),
               ("AETEALSK", 2, "Missed by the screen"), ("AFLDEMK", 2, "Missed by the screen"), ("HLTGLTFDTYK", 2, "PROCAL")]  # fmt: skip
index = pools.PoolIndex.from_references(pq.read_table(E / "search_space.parquet"))
pool1 = [int(i) for i in index.pool("TUM_first_pool_1")]
rec_in = pq.read_table(E / f"p2/p2_candidates/{REC}/{IN_POOL_RUN}.parquet").to_pandas()
gp = gate_psms(run_psms(IN_POOL_RUN), max_pep=0.01).to_pandas()
gp = gp[(gp["raw_file"].str.lower() == IN_POOL_RUN.lower())]
# confident in-pool IDs that the refined screen misses
miss = []
for sid in pool1:
    if sid not in truth:
        continue
    t = truth[sid]
    reg = rec_in[rec_in["sid"] == sid]
    if not ((reg["rt_left"] - 30 <= t) & (reg["rt_right"] + 30 >= t)).any():
        pep = gp.loc[gp["sequence"] == space.at[sid, "residues"], "pep"].min()
        miss.append((pep, sid))
miss = [s for _, s in sorted(miss)[:3]]
# co-eluting identified near-isobars (5-20 ppm on M0, |dRT| <= 2 min), at z = 2
C13, P = ss.C13_SPACING, ss.PROTON_MASS
gp_mz = gp[["sequence", "modified_sequence", "charge", "mz", "retention_time_s"]].copy()
iso_pairs = []
for sid in pool1:
    if sid not in truth or sid in miss:
        continue
    m0 = (space.at[sid, "monoisotopic_mass"] + 2 * P) / 2
    near = gp_mz[
        (gp_mz["sequence"] != space.at[sid, "residues"])
        & ((gp_mz["retention_time_s"] - truth[sid]).abs() <= 120)
    ]
    ppm = (near["mz"] - m0) / m0 * 1e6
    near = near[(ppm.abs() >= 5) & (ppm.abs() <= 20) & (near["charge"] == 2)]
    if len(near):
        iso_pairs.append((len(near), sid))
iso_pairs = [s for _, s in sorted(iso_pairs, reverse=True)[:2]]
ex = []
for k, (seq, z, kind) in enumerate(first_round):
    ex.append({"example": f"E{k + 1:02d}", "sid": pid[seq], "charge": z, "kind": kind})
for sid in miss:
    ex.append(
        {
            "example": f"E{len(ex) + 1:02d}",
            "sid": sid,
            "charge": 2,
            "kind": "Confident ID, missed by the refined screen",
        }
    )
for sid in iso_pairs:
    ex.append(
        {
            "example": f"E{len(ex) + 1:02d}",
            "sid": sid,
            "charge": 2,
            "kind": "Co-eluting near-isobar",
        }
    )

# chromatograms, fresh from the bundles ---------------------------------------------------------
from constellation.massspec.quant.chromatogram import extract_xic  # noqa: E402
from constellation.massspec.sweep import ChannelSettings  # noqa: E402

offsets = {stem: r["ppm_offset"] for stem, r in p2["runs"].items()}
space_t = pq.read_table(E / "search_space.parquet")
targets = pd.concat(
    [
        ss.targets_for_run(space_t, charges=(e["charge"],), subset=np.array([e["sid"]])).to_pandas()
        for e in ex
    ],
    ignore_index=True,
)
import pyarrow as pa  # noqa: E402

tgt = pa.Table.from_pandas(targets, preserve_index=False)
ex_rows, pt_rows, reg_rows, nb_rows = [], [], [], []
for e in ex:
    sid = e["sid"]
    e["psm_rt"] = truth.get(sid)
    e["proforma"] = space.at[sid, "proforma"]
    ex_rows.append({"example": e["example"], "kind": e["kind"], "sequence": e["proforma"], "charge": e["charge"],
                    "psm_rt_min": None if e["psm_rt"] is None else round(e["psm_rt"] / 60, 3), "in_pool_run": short_run(IN_POOL_RUN)})  # fmt: skip
for stem in [IN_POOL_RUN, *OFF_RUNS]:
    role = "In-pool" if stem == IN_POOL_RUN else "Off-pool"
    bundle = bundles.require_run_bundle(D, runs[stem]["dataset"], stem)
    res_ = extract_xic(bundle, tgt, acquisition_id=0, levels=[1], tolerances={1: (20.0, "ppm")},
                       channels=ChannelSettings(isotope_offsets=ISO), rule="nearest", drop_unmatched=False, layout="precursor_major")  # fmt: skip
    tr = (
        res_.table("xic_trace")
        .select(["target_id", "rt", "isotope", "intensity", "mz_error_ppm"])
        .to_pandas()
    )
    labels = pq.read_table(E / f"p2/labels/{stem}.parquet").to_pandas().set_index("sequence_id")
    off = offsets.get(stem, 0.0)
    ps_all = run_psms(stem).to_pandas()
    ps_all = ps_all[
        (ps_all["raw_file"].str.lower() == stem.lower())
        & (ps_all["pep"] <= 0.05)
        & ~ps_all["is_decoy"].fillna(False)
    ]
    for e in ex:
        tid = e["sid"] * ss.TARGET_ID_STRIDE + e["charge"]
        if e["psm_rt"] is None:
            continue
        w = tr[(tr["target_id"] == tid) & ((tr["rt"] - e["psm_rt"]).abs() <= 360)]
        lab = labels["label"].get(e["sid"], "off_pool")
        wi = w.pivot_table(index="rt", columns="isotope", values="intensity", aggfunc="max")
        we = w.pivot_table(index="rt", columns="isotope", values="mz_error_ppm", aggfunc="first")
        for rt in wi.index:
            r = {
                "example": e["example"],
                "run": short_run(stem),
                "role": role,
                "label": lab,
                "rt_min": round(rt / 60, 4),
            }
            for k in ISO:
                v = wi.at[rt, k] if k in wi.columns else None
                r[ISO_COL[k]] = (
                    None
                    if v is None or (isinstance(v, float) and math.isnan(v))
                    else round(float(v))
                )
            for k in (0, 1):
                v = we.at[rt, k] if (k in we.columns and rt in we.index) else None
                r[f"ppm_{ISO_COL[k]}"] = (
                    None
                    if v is None or (isinstance(v, float) and math.isnan(v))
                    else round(float(v) - off, 2)
                )
            pt_rows.append(r)
        for setting_name, path in (
            (REC, E / f"p2/p2_candidates/{REC}/{stem}.parquet"),
            ("first round (20 ppm)", E / f"p2/p2_candidates_v1/{stem}.parquet"),
        ):
            if not path.exists():
                continue
            c = pq.read_table(path).to_pandas()
            for _, g in c[c["target_id"] == tid].iterrows():
                reg_rows.append({"example": e["example"], "run": short_run(stem), "setting": "refined" if setting_name == REC else "first round",
                                 "rank": int(g["rank"]), "rt_left_min": round(g["rt_left"] / 60, 3), "rt_right_min": round(g["rt_right"] / 60, 3),
                                 "apex_min": round(g["apex_rt"] / 60, 3), "n_qual": int(g["n_qual"]), "median_ppm": round(float(g["median_e"]), 2),
                                 "median_kl": round(float(g["median_kl"]), 3)})  # fmt: skip
        # identified neighbours on our channels (any modification, PEP <= 0.05)
        ps = ps_all[(ps_all["retention_time_s"] - e["psm_rt"]).abs() <= 360]
        mass = float(space.at[e["sid"], "monoisotopic_mass"])
        chan = {k: (mass + k * C13 + e["charge"] * P) / e["charge"] for k in ISO}
        found = {}
        for r in ps.itertuples():
            m2 = float(r.mz) * int(r.charge) - int(r.charge) * P
            for kk in ISO:
                mzk = (m2 + kk * C13 + int(r.charge) * P) / int(r.charge)
                for k, cm in chan.items():
                    d = (mzk - cm) / cm * 1e6
                    if abs(d) <= 20:
                        key = (r.modified_sequence or r.sequence, int(r.charge), kk, k)
                        f = found.setdefault(
                            key, {"rts": [], "ppm": d, "self": r.modified_sequence == e["proforma"]}
                        )
                        f["rts"].append(float(r.retention_time_s))
        for (seq, zz, kk, k), f in found.items():
            nb_rows.append({"example": e["example"], "run": short_run(stem), "sequence": seq, "charge": zz, "its_isotope": kk,
                            "our_isotope": k, "ppm": round(f["ppm"], 1), "n_psms": len(f["rts"]), "self": f["self"],
                            "rt_min_first": round(min(f["rts"]) / 60, 3), "rt_min_last": round(max(f["rts"]) / 60, 3),
                            "rt_min_median": round(float(np.median(f["rts"])) / 60, 3)})  # fmt: skip
dump(
    "xic_examples",
    ex_rows,
    f"unmodified MaxQuant PSMs in {IN_POOL_RUN}; first-round examples plus refined-screen misses and near-isobars",
)
dump(
    "xic_points",
    pt_rows,
    "extract_xic on the bundles (MS1, rule nearest, zero-filled, 20 ppm); ppm offset-corrected per run (PROCAL)",
)
dump(
    "xic_regions",
    reg_rows,
    f"p2/p2_candidates/{REC} (refined) and p2_candidates_v1 (first round, 20 ppm)",
)
dump(
    "xic_neighbours",
    nb_rows,
    "MaxQuant PSMs (any modification, PEP <= 0.05) with an isotope within 20 ppm of an example's channel",
)
print("done")
