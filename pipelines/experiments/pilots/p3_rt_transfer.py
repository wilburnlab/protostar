#!/usr/bin/env python
"""Pilot P3: how well a peptide's retention time transfers to another run (extraction design).

Tier 2 extracts off-pool and carryover controls at a peptide's *transferred* reference RT, so the
window width depends on how far a transferred RT lands from where the peptide really elutes.
This pilot measures that with MaxQuant identifications (PEP ≤ 0.01) as truth, in three regimes
of increasing difficulty:

- ``within_pool``: a run's own pool peptides, from the pool's other injections;
- ``cross_pool``: peptides synthesized in two pools, from one pool's runs into the other's;
- ``cross_dataset``: a Wilhelm2021 ``noincl`` re-acquisition, from the Zolg2017 runs of its pool.

The reference clock is PROCAL's published 60-minute-gradient RT. Each run gets two maps with
Constellation's ``rt.align``:

- ``procal``: fitted on PROCAL alone;
- ``refined``: PROCAL plus the run's own pool peptides. Evaluated peptides are held out 5-fold,
  so no peptide calibrates its own prediction.

A peptide's reference RT is the median, over the runs it is transferred *from*, of its
PROCAL-mapped RT. The output is the residuals (observed − predicted, s) per peptide, run, regime
and map, plus summary statistics.

Usage::

    python pipelines/experiments/pilots/p3_rt_transfer.py --data-root <data_root> \
        --out <data_root>/pilots/extraction/p3
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import _common  # noqa: E402,F401  (pins libstdc++ before constellation imports)

import numpy as np  # noqa: E402
import pyarrow as pa  # noqa: E402
import pyarrow.compute as pc  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

from protostar import bundles, pools  # noqa: E402
from protostar import search_space as ss  # noqa: E402
from protostar.experiments.scans import gate_psms, load_psms  # noqa: E402
from protostar.peptides.procal import load_procal  # noqa: E402

CENTER = "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2"
_MINUTES = re.compile(r"^\s*([0-9.]+)")


def procal_reference_s() -> dict[str, float]:
    """PROCAL sequence → published 60-minute-gradient RT, in seconds."""
    out = {}
    t = load_procal()
    for seq, rt in zip(
        t.column("sequence").to_pylist(), t.column("rt_60min_gradient").to_pylist(), strict=True
    ):
        m = _MINUTES.match(rt or "")
        if m:
            out[seq] = float(m.group(1)) * 60.0
    return out


def best_rt(psms: pa.Table, stem: str) -> dict[str, float]:
    """Bare sequence → RT (s) of its best PSM (lowest PEP) in ``stem``."""
    raw = pc.utf8_lower(psms.column("raw_file"))
    psms = psms.filter(pc.equal(raw, stem.lower()))
    order = pc.sort_indices(psms, sort_keys=[("pep", "ascending"), ("score", "descending")])
    out: dict[str, float] = {}
    for s, rt in zip(
        psms.column("sequence").take(order).to_pylist(),
        psms.column("retention_time_s").take(order).to_pylist(),
        strict=True,
    ):
        if s is not None and rt is not None and s not in out:
            out[s] = float(rt)
    return out


#: rt.align's refinement ("auto": modal, about 6.5 s a fit; None: about 0.1 s, nearly the same fit).
REFINE: str | None = "auto"


def fit_map(ref: list[float], obs: list[float]):
    from constellation.massspec.rt import align

    if len(ref) < 8:
        return None
    order = np.argsort(ref)
    pairs = pa.table(
        {"reference_rt": np.asarray(ref)[order], "observed_rt": np.asarray(obs)[order]}
    )
    try:
        return align(pairs, refine=REFINE)
    except Exception:  # noqa: BLE001 - a run whose anchors cannot be fitted is reported as such
        return None


def stats(res: np.ndarray) -> dict:
    a = np.abs(res[np.isfinite(res)])
    if a.size == 0:
        return {"n": 0}
    mad_sd = 1.4826 * float(np.median(np.abs(res - np.median(res))))
    return {
        "n": int(a.size),
        "median": round(float(np.median(res)), 2),
        "sd_mad": round(mad_sd, 2),
        "p95_abs": round(float(np.quantile(a, 0.95)), 1),
        "p99_abs": round(float(np.quantile(a, 0.99)), 1),
        "p999_abs": round(float(np.quantile(a, 0.999)), 1),
        "frac_abs_gt_60s": round(float((a > 60).mean()), 4),
        "frac_abs_gt_300s": round(float((a > 300).mean()), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--space-cache", type=Path, default=None)
    ap.add_argument("--center-stem", default=CENTER)
    ap.add_argument("--block-half", type=int, default=10, help="runs either side of the centre run")
    ap.add_argument("--pool-pairs", type=int, default=5)
    ap.add_argument("--min-shared", type=int, default=50)
    ap.add_argument("--noincl-runs", type=int, default=2)
    ap.add_argument("--max-pep", type=float, default=0.01)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--refine", choices=["auto", "none"], default="auto")
    args = ap.parse_args()
    global REFINE
    REFINE = None if args.refine == "none" else args.refine
    args.out.mkdir(parents=True, exist_ok=True)
    t_start = time.monotonic()
    report: dict = {
        "provenance": bundles.constellation_provenance(),
        "args": {k: str(v) for k, v in vars(args).items()},
    }

    def save():
        (args.out / "p3_report.json").write_text(json.dumps(report, indent=1, default=str))

    runs = pools.run_table(args.data_root)
    rows = {r["stem"]: r for r in runs.to_pylist()}
    space_path = args.space_cache or args.out / "search_space.parquet"
    space = pq.read_table(space_path) if space_path.exists() else ss.build_search_space()
    index = pools.PoolIndex.from_references(space)
    residues = space.column("residues").to_numpy(zero_copy_only=False)
    members = {p: set(residues[ids].tolist()) for p, ids in index.members.items()}
    runs_of_pool: dict[str, list[str]] = defaultdict(list)
    for s, r in rows.items():
        runs_of_pool[r["pool"]].append(s)

    # A: the block of consecutive runs around the centre run on its instrument (whole pools)
    me = rows[args.center_stem]
    block = [
        s for s, r in rows.items()
        if r["instrument_serial"] == me["instrument_serial"] and r["instrument_rank"] is not None
        and abs(r["instrument_rank"] - me["instrument_rank"]) <= args.block_half
    ]  # fmt: skip
    block_pools = sorted({rows[s]["pool"] for s in block})
    # B: pool pairs sharing the most peptides
    seq_pools: dict[str, list[str]] = defaultdict(list)
    for p, seqs in members.items():
        for s in seqs:
            seq_pools[s].append(p)
    pair_n = Counter()
    for ps in seq_pools.values():
        if 2 <= len(ps) <= 4:
            for i in range(len(ps)):
                for j in range(i + 1, len(ps)):
                    pair_n[tuple(sorted((ps[i], ps[j])))] += 1
    pairs = [pq_ for pq_, n in pair_n.most_common() if n >= args.min_shared
             and runs_of_pool.get(pq_[0]) and runs_of_pool.get(pq_[1])][: args.pool_pairs]  # fmt: skip
    # C: noincl re-acquisitions and the Zolg2017 runs of their pools
    noincl = sorted(s for s, r in rows.items() if r["noincl"])[: args.noincl_runs]
    report["selection"] = {"block": sorted(block), "block_pools": block_pools,
                           "pool_pairs": [[a, b, pair_n[(a, b)]] for a, b in pairs], "noincl": noincl}  # fmt: skip
    save()

    needed = set()
    for p in block_pools:
        needed |= set(runs_of_pool[p])
    for a, b in pairs:
        needed |= set(runs_of_pool[a]) | set(runs_of_pool[b])
    for s in noincl:
        needed |= {x for x in runs_of_pool[rows[s]["pool"]]}
    # observed RTs per run
    obs: dict[str, dict[str, float]] = {}
    tmp = Path(tempfile.mkdtemp(prefix="p3psm-", dir=os.environ.get("TMPDIR")))
    missing = []
    for s in sorted(needed):
        z = pools.search_zip_for_stem(args.data_root, rows[s]["dataset"], s)
        if z is None:
            missing.append(s)
            continue
        t0 = time.monotonic()
        psms = gate_psms(load_psms(z, extract_dir=tmp / s), max_pep=args.max_pep)
        obs[s] = best_rt(psms, s)
        print(f"psms {s}: {len(obs[s])} sequences, {time.monotonic() - t0:.1f}s", flush=True)
    report["psm_runs"] = {"loaded": len(obs), "missing_zip": missing}
    save()

    procal = procal_reference_s()
    procal_map = {}
    for s, o in obs.items():
        common = [q for q in procal if q in o]
        t0 = time.monotonic()
        procal_map[s] = fit_map([procal[q] for q in common], [o[q] for q in common])
        print(f"procal map {s}: {len(common)} anchors, {time.monotonic() - t0:.1f}s", flush=True)
    report["procal_maps"] = {
        "fitted": sum(m is not None for m in procal_map.values()),
        "runs": len(procal_map),
    }
    save()

    # each run's identified sequences on the PROCAL reference clock, mapped in one call per run
    ref_of: dict[str, dict[str, float]] = {}
    for s_, o in obs.items():
        if procal_map.get(s_) is None or not o:
            continue
        seqs_ = list(o)
        mapped = procal_map[s_].to_reference(np.array([o[q] for q in seqs_])).numpy()
        ref_of[s_] = dict(zip(seqs_, mapped.tolist(), strict=True))

    def ref_rt(seq: str, from_stems: list[str]) -> float | None:
        vals = [ref_of[f][seq] for f in from_stems if f in ref_of and seq in ref_of[f]]
        return float(np.median(vals)) if vals else None

    out_rows = []

    def evaluate(regime: str, run: str, evaluated: dict[str, float]):
        """``evaluated``: sequence → reference RT, for sequences identified in ``run``."""
        if run not in obs or procal_map.get(run) is None or not evaluated:
            return
        o = obs[run]
        seqs = sorted(evaluated)
        ref = np.array([evaluated[q] for q in seqs])
        truth = np.array([o[q] for q in seqs])
        pred_p = procal_map[run].to_observed(ref).numpy()
        # anchors: the run's own pool peptides not being evaluated, referenced from its pool's other runs
        own = [x for x in runs_of_pool[rows[run]["pool"]] if x != run]
        anchors = {
            q: ref_rt(q, own)
            for q in members.get(rows[run]["pool"], ())
            if q in o and q not in evaluated
        }
        anchors = {q: r for q, r in anchors.items() if r is not None}
        base_ref = [procal[q] for q in procal if q in o] + list(anchors.values())
        base_obs = [o[q] for q in procal if q in o] + [o[q] for q in anchors]
        pred_r = np.full(len(seqs), np.nan)
        fold = np.arange(len(seqs)) % args.folds
        t0 = time.monotonic()
        for f in range(args.folds):
            train = fold != f
            m = fit_map(base_ref + ref[train].tolist(), base_obs + truth[train].tolist())
            if m is not None:
                pred_r[fold == f] = m.to_observed(ref[fold == f]).numpy()
        inj = rows[run]["injection"]
        print(
            f"{regime} {run}: {len(seqs)} evaluated, {len(anchors)} anchors, {time.monotonic() - t0:.1f}s",
            flush=True,
        )
        for q, r_, t_, a_, b_ in zip(seqs, ref, truth, pred_p, pred_r, strict=True):
            out_rows.append({"regime": regime, "run": run, "injection": inj, "sequence": q,
                             "reference_rt": r_, "observed_rt": t_,
                             "resid_procal": t_ - a_, "resid_refined": t_ - b_,
                             "n_anchors": len(anchors)})  # fmt: skip

    for run in sorted(block):
        pool = rows[run]["pool"]
        others = [x for x in runs_of_pool[pool] if x != run]
        ev = {q: ref_rt(q, others) for q in members.get(pool, ()) if run in obs and q in obs[run]}
        evaluate("within_pool", run, {q: r for q, r in ev.items() if r is not None})
    for a, b in pairs:
        shared = members[a] & members[b]
        for src, dst in ((a, b), (b, a)):
            for run in runs_of_pool[dst]:
                if run not in obs:
                    continue
                ev = {q: ref_rt(q, runs_of_pool[src]) for q in shared if q in obs[run]}
                evaluate("cross_pool", run, {q: r for q, r in ev.items() if r is not None})
    for run in noincl:
        pool = rows[run]["pool"]
        src = [x for x in runs_of_pool[pool] if not rows[x]["noincl"]]
        if run not in obs:
            continue
        ev = {q: ref_rt(q, src) for q in members.get(pool, ()) if q in obs[run]}
        evaluate("cross_dataset", run, {q: r for q, r in ev.items() if r is not None})

    table = pa.Table.from_pylist(out_rows)
    if table.num_rows:
        pq.write_table(table, args.out / "p3_residuals.parquet")
    summary = {}
    for regime in ("within_pool", "cross_pool", "cross_dataset"):
        sel = [r for r in out_rows if r["regime"] == regime]
        for method in ("procal", "refined"):
            res = np.array([r[f"resid_{method}"] for r in sel], dtype=float)
            summary[f"{regime}/{method}"] = stats(res)
            for inj in pools.INJECTION_TYPES:
                r_inj = np.array(
                    [r[f"resid_{method}"] for r in sel if r["injection"] == inj], dtype=float
                )
                if r_inj.size:
                    summary[f"{regime}/{method}/{inj}"] = stats(r_inj)
    report["summary"] = summary
    report["total_s"] = round(time.monotonic() - t_start, 1)
    save()
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
