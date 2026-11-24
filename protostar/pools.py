"""Pools, run order and per-run peptide labels: the in-pool / off-pool ground truth.

ProteomeTools synthesized peptides in pools of about 1,000 and acquired each pool as four
injections (DDA, 2xIT_2xHCD, 3xHCD, ETD). Which pool a run holds is in its file name; which
peptides a pool holds is in the published supplements (``protostar.peptides``). Together they
say, for any run and any sequence, whether the sequence *should* be there. The global extraction
uses that instead of decoys to judge its detections (decision log D19, once written).

Three pieces:

- ``parse_stem``: pool, injection type, replicate and ``noincl`` from a ``.raw`` stem, for every
  file of the three datasets.
- ``run_table``: every run with its instrument and acquisition time, ranked per instrument
  **across datasets**. A dataset's own ``acquisition_order`` restarts per dataset, so it cannot
  find the run acquired before the first run of another dataset on the same instrument.
- ``labels_for_run``: for one run, the sequences that are not plain off-pool, each with why:

  ======== ==================================================================================
  label    meaning
  ======== ==================================================================================
  standard spiked into the run: PROCAL, the heavy Pierce standards listed in its pool, or a
           supplement QC peptide (``kind`` says which)
  in_pool  synthesized in the run's pool
  carryover in the pool of one of the ``k`` previous distinct pools on the same instrument
  shadowed off-pool, but its neutral mass matches something expected present (``kind``:
           ``composition`` twin, ``isobaric`` within tolerance, Met ``oxidation``, single-residue
           synthesis ``deletion``)
  ======== ==================================================================================

  Every sequence not listed is ``off_pool``. Precedence is the table's order.

Glue only: labels are ground truth for scoring, not a score.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from constellation.core.sequence.proforma import parse_proforma
from constellation.core.sequence.protein import peptide_masses

from protostar.metadata.driver import acquisitions_out_path
from protostar.peptides.reference import POOL_PREFIX_DATASET, load_reference, pool_prefix
from protostar.search_space import DATASETS, normalize_sequence

#: Injection types, as the file names spell them.
INJECTION_TYPES: tuple[str, ...] = ("DDA", "2xIT_2xHCD", "3xHCD", "ETD")

#: ``<run>_<plate>-<pool>_<n>_<f1>_<f2>-<injection>-<gradient>[noincl]-R<replicate>``.
STEM_RE = re.compile(
    r"^(?P<run>[^_]+)_(?P<plate>[^-]+)-"
    r"(?P<pool>(?P<prefix>[A-Za-z][A-Za-z0-9_]*?)_(?P<pool_n>\d+))_(?P<f1>\d+)_(?P<f2>\d+)-"
    r"(?P<injection>[A-Za-z0-9_]+?)-(?P<gradient>\d+h)(?P<noincl>noincl)?-R(?P<replicate>\d+)$"
)

#: Monoisotopic mass of one oxygen atom (Met oxidation), Da.
OXYGEN_MASS: float = 15.9949146196


@dataclass(frozen=True)
class RunName:
    """What a ``.raw`` stem says about its run."""

    stem: str
    run: str  # leading run code, e.g. ``01625b``
    plate: str  # e.g. ``GA1``
    pool: str  # e.g. ``TUM_first_pool_1``
    pool_prefix: str  # e.g. ``TUM_first_pool``
    pool_dataset: str  # the dataset whose supplement lists the pool's peptides
    pool_part: str  # the two fields after the pool number, e.g. ``01_01``
    injection: str  # one of ``INJECTION_TYPES``
    gradient: str  # e.g. ``1h``
    noincl: bool  # a Wilhelm2021 no-inclusion re-acquisition of a Zolg2017 pool
    replicate: int


def parse_stem(stem: str) -> RunName:
    """Parse a ``.raw`` stem (with or without the extension). Raises ``ValueError`` if it does not
    match, or if its pool prefix or injection type is unknown."""
    stem = stem[:-4] if stem.lower().endswith(".raw") else stem
    m = STEM_RE.match(stem)
    if m is None:
        raise ValueError(f"unrecognised run stem {stem!r}")
    prefix = m["prefix"]
    if prefix not in POOL_PREFIX_DATASET:
        raise ValueError(f"unknown pool prefix {prefix!r} in {stem!r}")
    if m["injection"] not in INJECTION_TYPES:
        raise ValueError(f"unknown injection type {m['injection']!r} in {stem!r}")
    return RunName(
        stem=stem,
        run=m["run"],
        plate=m["plate"],
        pool=m["pool"],
        pool_prefix=prefix,
        pool_dataset=POOL_PREFIX_DATASET[prefix],
        pool_part=f"{m['f1']}_{m['f2']}",
        injection=m["injection"],
        gradient=m["gradient"],
        noincl=m["noincl"] is not None,
        replicate=int(m["replicate"]),
    )


def search_zip_for_stem(data_root: str | Path, dataset: str, stem: str) -> Path | None:
    """The MaxQuant search zip of one run, or ``None`` if there is none.

    Zolg2017 and Wilhelm2021 name a search after the pool and injection
    (``TUM_first_pool_1_01_01_DDA-1h-R2-tryptic.zip``; the enzyme suffix varies); Gessulat2019
    names it after the ``.raw`` stem, sometimes in a different letter case."""
    sdir = Path(data_root) / "search" / dataset
    if not sdir.is_dir():
        return None
    name = parse_stem(stem)
    by_lower = {z.name.lower(): z for z in sdir.glob("*.zip")}
    hit = by_lower.get(f"{name.stem}.zip".lower())
    if hit is not None:
        return hit
    prefix = (
        f"{name.pool}_{name.pool_part}_{name.injection}-{name.gradient}"
        f"{'noincl' if name.noincl else ''}-R{name.replicate}-"
    ).lower()
    hits = sorted(z for k, z in by_lower.items() if k.startswith(prefix))
    if len(hits) > 1:
        raise ValueError(f"several search zips match {stem!r}: {[h.name for h in hits]}")
    return hits[0] if hits else None


RUN_TABLE: pa.Schema = pa.schema(
    [
        pa.field("dataset", pa.string(), nullable=False),
        pa.field("stem", pa.string(), nullable=False),
        pa.field("acquisition_id", pa.int64(), nullable=True),
        pa.field("acquisition_datetime", pa.string(), nullable=True),
        pa.field("instrument_serial", pa.string(), nullable=True),
        pa.field("instrument_model", pa.string(), nullable=True),
        # 1-based rank by acquisition time within instrument_serial, across all datasets.
        pa.field("instrument_rank", pa.int32(), nullable=True),
        pa.field("pool", pa.string(), nullable=False),
        pa.field("pool_dataset", pa.string(), nullable=False),
        pa.field("injection", pa.string(), nullable=False),
        pa.field("noincl", pa.bool_(), nullable=False),
        pa.field("replicate", pa.int16(), nullable=False),
    ],
    metadata={b"schema_name": b"ProtostarRunTable", b"schema_version": b"1"},
)


def run_table_from_acquisitions(acquisitions: dict[str, pa.Table]) -> pa.Table:
    """Build the run table from each dataset's acquisition table (stage 20 output)."""
    rows = []
    for ds, t in sorted(acquisitions.items()):
        for r in t.to_pylist():
            name = parse_stem(r["source_file"])
            rows.append(
                {
                    "dataset": ds,
                    "stem": name.stem,
                    "acquisition_id": r.get("acquisition_id"),
                    "acquisition_datetime": r.get("acquisition_datetime"),
                    "instrument_serial": r.get("instrument_serial"),
                    "instrument_model": r.get("instrument_model"),
                    "instrument_rank": None,
                    "pool": name.pool,
                    "pool_dataset": name.pool_dataset,
                    "injection": name.injection,
                    "noincl": name.noincl,
                    "replicate": name.replicate,
                }
            )
    # Rank per instrument by time (stem breaks ties); runs without a time or serial stay unranked.
    timed = [r for r in rows if r["instrument_serial"] and r["acquisition_datetime"]]
    timed.sort(key=lambda r: (r["instrument_serial"], r["acquisition_datetime"], r["stem"]))
    serial, rank = None, 0
    for r in timed:
        rank = rank + 1 if r["instrument_serial"] == serial else 1
        serial = r["instrument_serial"]
        r["instrument_rank"] = rank
    rows.sort(key=lambda r: (r["dataset"], r["stem"]))
    return pa.Table.from_pylist(rows, schema=RUN_TABLE)


def run_table(data_root: str | Path, datasets: tuple[str, ...] = DATASETS) -> pa.Table:
    """The run table from ``<data_root>/proc/<dataset>/acquisitions.parquet`` (stage 20)."""
    return run_table_from_acquisitions(
        {ds: pq.read_table(acquisitions_out_path(data_root, ds)) for ds in datasets}
    )


@dataclass(frozen=True)
class PreviousPool:
    """A pool acquired before a run on the same instrument."""

    pool: str
    pools_back: int  # 1 = the most recent different pool
    runs_back: int  # how many runs back its latest run was
    dt_s: float | None  # seconds between that run's start and this run's start


def previous_pools(runs: pa.Table, stem: str, k: int = 2, max_runs: int = 16) -> list[PreviousPool]:
    """The ``k`` most recent pools, other than the run's own, acquired before ``stem`` on the
    same instrument (looking back at most ``max_runs`` runs)."""
    rows = runs.to_pylist()
    me = next((r for r in rows if r["stem"] == stem), None)
    if me is None:
        raise KeyError(f"{stem!r} is not in the run table")
    if me["instrument_rank"] is None:
        return []
    same = {
        r["instrument_rank"]: r
        for r in rows
        if r["instrument_serial"] == me["instrument_serial"] and r["instrument_rank"] is not None
    }
    t0 = np.datetime64(me["acquisition_datetime"])
    out: list[PreviousPool] = []
    seen = {me["pool"]}
    for back in range(1, max_runs + 1):
        prev = same.get(me["instrument_rank"] - back)
        if prev is None:
            break
        if prev["pool"] in seen:
            continue
        seen.add(prev["pool"])
        dt = float((t0 - np.datetime64(prev["acquisition_datetime"])) / np.timedelta64(1, "s"))
        out.append(PreviousPool(prev["pool"], len(out) + 1, back, dt))
        if len(out) == k:
            break
    return out


class PoolIndex:
    """Pool → ``sequence_id`` array, plus the supplement QC peptides that have no pool."""

    def __init__(self, members: dict[str, np.ndarray], qc: dict[str, np.ndarray] | None = None):
        self.members = members
        self.qc = qc or {}

    @classmethod
    def from_references(cls, space: pa.Table, datasets: tuple[str, ...] = DATASETS) -> PoolIndex:
        lookup = dict(
            zip(
                space.column("proforma").to_pylist(),
                space.column("sequence_id").to_pylist(),
                strict=True,
            )
        )
        members: dict[str, set[int]] = {}
        qc: dict[str, set[int]] = {}
        for ds in datasets:
            t = load_reference(ds)
            for seq, pool, is_qc, qc_type in zip(
                t.column("sequence").to_pylist(),
                t.column("pool").to_pylist(),
                t.column("is_qc").to_pylist(),
                t.column("qc_type").to_pylist(),
                strict=True,
            ):
                sid = lookup[normalize_sequence(seq)]
                if pool is not None:
                    members.setdefault(pool, set()).add(sid)
                elif is_qc:
                    qc.setdefault(qc_type or "qc", set()).add(sid)
        as_arr = lambda d: {k: np.array(sorted(v), dtype=np.int64) for k, v in d.items()}  # noqa: E731
        return cls(as_arr(members), as_arr(qc))

    def pool(self, pool: str) -> np.ndarray:
        if pool not in self.members:
            raise KeyError(
                f"pool {pool!r} has no peptides in the reference (prefix {pool_prefix(pool)!r})"
            )
        return self.members[pool]


RUN_LABEL_TABLE: pa.Schema = pa.schema(
    [
        pa.field("sequence_id", pa.int64(), nullable=False),
        pa.field("label", pa.string(), nullable=False),  # standard / in_pool / carryover / shadowed
        pa.field("kind", pa.string(), nullable=True),  # standard or shadow kind
        pa.field("source_pool", pa.string(), nullable=True),
        pa.field("pools_back", pa.int16(), nullable=True),  # carryover
        pa.field("runs_back", pa.int16(), nullable=True),  # carryover
        pa.field("dt_s", pa.float64(), nullable=True),  # carryover
        pa.field("shadow_of", pa.int64(), nullable=True),  # the expected-present sequence matched
    ],
    metadata={b"schema_name": b"ProtostarRunLabels", b"schema_version": b"1"},
)

SHADOW_KINDS: tuple[str, ...] = ("composition", "isobaric", "oxidation", "deletion")

_TOKEN = re.compile(r"[A-Z](?:\[[^\]]*\])?")


def _deletion_proformas(proforma: str) -> list[str]:
    tokens = _TOKEN.findall(proforma)
    return sorted({"".join(tokens[:i] + tokens[i + 1 :]) for i in range(len(tokens))} - {""})


def _mass_window_hits(sorted_mass: np.ndarray, query: np.ndarray, tol_ppm: float):
    """For each query mass, the index range of ``sorted_mass`` within ``tol_ppm``: (query idx, hit idx)."""
    lo = np.searchsorted(sorted_mass, query * (1 - tol_ppm * 1e-6), side="left")
    hi = np.searchsorted(sorted_mass, query * (1 + tol_ppm * 1e-6), side="right")
    n = hi - lo
    q = np.repeat(np.arange(query.size), n)
    start = np.repeat(lo, n)
    offs = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    return q, start + offs


def labels_for_run(
    stem: str,
    runs: pa.Table,
    space: pa.Table,
    pools: PoolIndex,
    *,
    k_carryover: int = 2,
    tol_ppm: float = 20.0,
    shadow_kinds: tuple[str, ...] = SHADOW_KINDS,
    max_oxidations: int = 2,
) -> pa.Table:
    """The labelled sequences of one run (``RUN_LABEL_TABLE``); every other sequence is off-pool."""
    me = next((r for r in runs.to_pylist() if r["stem"] == stem), None)
    if me is None:
        raise KeyError(f"{stem!r} is not in the run table")
    sid = space.column("sequence_id").to_numpy()
    if not np.array_equal(sid, np.arange(sid.size)):
        raise ValueError("space must be indexed by sequence_id (as build_search_space writes it)")
    mass = space.column("monoisotopic_mass").to_numpy()
    rows: dict[int, dict] = {}

    def put(ids, **fields):
        for i in np.asarray(ids, dtype=np.int64).tolist():
            if i not in rows:
                rows[i] = {"sequence_id": i, **fields}

    # standard
    put(
        np.flatnonzero(space.column("is_procal").to_numpy(zero_copy_only=False)),
        label="standard",
        kind="procal",
    )
    own = pools.pool(me["pool"])
    heavy = space.column("is_heavy").to_numpy(zero_copy_only=False)
    put(own[heavy[own]], label="standard", kind="prtc", source_pool=me["pool"])
    for qc_type, ids in sorted(pools.qc.items()):
        put(ids, label="standard", kind=qc_type)
    # in_pool
    put(own, label="in_pool", source_pool=me["pool"])
    present = np.array(sorted(rows), dtype=np.int64)
    # carryover
    for prev in previous_pools(runs, stem, k=k_carryover):
        put(
            pools.pool(prev.pool),
            label="carryover",
            source_pool=prev.pool,
            pools_back=prev.pools_back,
            runs_back=prev.runs_back,
            dt_s=prev.dt_s,
        )
    # shadowed: off-pool sequences whose neutral mass matches something expected present
    order = np.argsort(mass, kind="stable")
    sorted_mass = mass[order]
    keys = space.column("composition_key").to_numpy(zero_copy_only=False)
    present_keys = {keys[i]: int(i) for i in present}
    species: list[tuple[str, np.ndarray, np.ndarray]] = []  # (kind, query masses, source ids)
    if "isobaric" in shadow_kinds or "composition" in shadow_kinds:
        species.append(("isobaric", mass[present], present))
    if "oxidation" in shadow_kinds:
        n_met = space.column("n_met").to_numpy()[present]
        for k in range(1, max_oxidations + 1):
            sel = n_met >= k
            species.append(("oxidation", mass[present[sel]] + k * OXYGEN_MASS, present[sel]))
    if "deletion" in shadow_kinds:
        proforma = space.column("proforma").to_pylist()
        dels, src = [], []
        for i in present.tolist():
            for d in _deletion_proformas(proforma[i]):
                dels.append(d)
                src.append(i)
        if dels:
            dm = peptide_masses([parse_proforma(d) for d in dels]).numpy()
            species.append(("deletion", dm, np.asarray(src, dtype=np.int64)))
    for kind, qmass, src in species:
        q, hit = _mass_window_hits(sorted_mass, qmass, tol_ppm)
        hit_ids = order[hit]
        for qi, h in zip(q.tolist(), hit_ids.tolist(), strict=True):
            if h in rows:
                continue
            k = kind
            if kind == "isobaric":
                k = "composition" if keys[h] in present_keys else "isobaric"
                if k not in shadow_kinds:
                    continue
            rows[h] = {
                "sequence_id": h,
                "label": "shadowed",
                "kind": k,
                "source_pool": me["pool"],
                "shadow_of": int(src[qi]),
            }
    out = sorted(rows.values(), key=lambda r: r["sequence_id"])
    return pa.Table.from_pylist(out, schema=RUN_LABEL_TABLE)


def label_counts(labels: pa.Table) -> dict[str, int]:
    """``{label or label:kind: count}`` for a quick look at one run's labels."""
    out: dict[str, int] = {}
    keys = pc.binary_join_element_wise(
        labels.column("label"), pc.fill_null(labels.column("kind"), ""), ":"
    )
    for k in keys.to_pylist():
        k = k.rstrip(":")
        out[k] = out.get(k, 0) + 1
    return dict(sorted(out.items()))


__all__ = [
    "INJECTION_TYPES",
    "OXYGEN_MASS",
    "PoolIndex",
    "PreviousPool",
    "RUN_LABEL_TABLE",
    "RUN_TABLE",
    "RunName",
    "SHADOW_KINDS",
    "STEM_RE",
    "label_counts",
    "labels_for_run",
    "parse_stem",
    "previous_pools",
    "run_table",
    "run_table_from_acquisitions",
    "search_zip_for_stem",
]
