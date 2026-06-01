"""Synthetic-peptide reference: the per-dataset sequence ↔ pool ↔ set map.

The ProteomeTools papers each publish a supplementary table listing every
synthetic peptide and the **pool** it was synthesized into; the pool name embeds
the prefix that the acquisition ``.raw`` stems carry (``TUM_first_pool_1``,
``Thermo_SRM_Pool_1``, ``TUM_HLA_1``, …), so this table is the ground-truth join
from peptide → pool → acquisition. It is the known-answer target source for XIC
extraction / validation and downstream identification work.

The published supplements are large Excel workbooks whose redistribution is not
ours to make; this module **extracts only the facts we need** (sequence, pool,
set, and QC/calibrant typing where given) into a compact, version-controlled
parquet committed under ``protostar/peptides/data/<dataset>.parquet``. The
workbooks themselves are a one-time manual input (see ``sources.py``); the
committed parquet is the durable, always-available reference.

Scope (decisions): bare AA ``sequence`` only — **no charge** (charge is a
per-extraction concern: the XIC target builder sweeps a charge range or takes it
from a search), and **no modified-sequence / m/z** (deferred to point of use, so
TMT/HLA chemistry is not baked in here). Retention time is carried only where a
dataset's supplement provides it (Wilhelm ``Average Retention Time`` / iRT).
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

#: Schema version — bump on any additive column change.
PEPTIDE_REFERENCE_SCHEMA_VERSION: int = 1

#: One row per (dataset, sequence, pool). A peptide can recur across pools/sets,
#: so ``peptide_id`` is row-unique within a dataset, not per-sequence.
PEPTIDE_REFERENCE_TABLE: pa.Schema = pa.schema(
    [
        pa.field("peptide_id", pa.int64(), nullable=False),  # row id within dataset
        pa.field("sequence", pa.string(), nullable=False),  # bare AA, as published
        pa.field("pool", pa.string(), nullable=True),  # e.g. TUM_first_pool_1 (null for QC)
        pa.field("pool_prefix", pa.string(), nullable=True),  # e.g. TUM_first_pool (pool minus _N)
        pa.field("set", pa.string(), nullable=False),  # proteotypic/tmt/hla_class_i/qc/…
        pa.field("is_qc", pa.bool_(), nullable=False),  # QC / calibrant peptide
        pa.field("qc_type", pa.string(), nullable=True),  # JPT-RT / JPT-QC / Pierce-RT (QC only)
        pa.field("rt", pa.float64(), nullable=True),  # published avg RT (min), where given
        pa.field("irt", pa.float64(), nullable=True),  # published iRT, where given
    ],
    metadata={
        b"schema_name": b"PeptideReferenceTable",
        b"schema_version": str(PEPTIDE_REFERENCE_SCHEMA_VERSION).encode("utf-8"),
    },
)

#: Committed reference lives beside this module so it ships with the repo.
_DATA_DIR = Path(__file__).resolve().parent / "data"


def reference_path(dataset: str) -> Path:
    """Path to the committed reference parquet for a dataset."""
    return _DATA_DIR / f"{dataset}.parquet"


def write_reference(table: pa.Table, dataset: str) -> Path:
    """Persist the curated reference table for a dataset (cast to schema)."""
    from constellation.core.io.schemas import cast_to_schema

    out = reference_path(dataset)
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(cast_to_schema(table, PEPTIDE_REFERENCE_TABLE), out)
    return out


def load_reference(dataset: str) -> pa.Table:
    """Load the committed synthetic-peptide reference for a dataset.

    Raises ``FileNotFoundError`` with a build hint if the parquet is absent
    (it is produced by ``pipelines/05_peptide_reference.py`` from the manual
    Excel supplements; see ``sources.py``).
    """
    p = reference_path(dataset)
    if not p.is_file():
        raise FileNotFoundError(
            f"no peptide reference for {dataset!r} at {p}; build it with "
            f"`python pipelines/05_peptide_reference.py --dataset {dataset}` "
            f"(needs the manual Excel supplements — see protostar/peptides/sources.py)"
        )
    return pq.read_table(p)


# ── cross-dataset pool resolution ──────────────────────────────────────
#
# A dataset's *acquired* pools are not always native to that dataset's own
# peptide supplement. Wilhelm2021 re-acquired 60 ProteomeTools pools as
# "no-inclusion" runs (``…-1hnoincl-…`` files, pool tokens ``TUM_first_pool_*`` /
# ``TUM_second_pool_*``); those peptides are published in the **Zolg2017**
# supplement, not Wilhelm's HLA/AspN/LysN sheets. So a run→pool→peptide join for
# Wilhelm2021 must source those pools from Zolg2017. The pool-name prefix is a
# globally-unambiguous owner key (verified: each prefix appears in exactly one
# dataset's reference), so we resolve by prefix rather than duplicating rows.

#: pool-name prefix (pool minus the trailing ``_<n>``) → owning dataset.
POOL_PREFIX_DATASET: dict[str, str] = {
    # Zolg2017 (also the source for Wilhelm2021's no-inclusion re-acquisitions)
    "TUM_first_pool": "Zolg2017",
    "TUM_second_pool": "Zolg2017",
    "TUM_third_pool": "Zolg2017",
    "Thermo_SRM_Pool": "Zolg2017",
    # Gessulat2019
    "TUM_isoform": "Gessulat2019",
    "TUM_proteo_TMT": "Gessulat2019",
    "TUM_second_addon": "Gessulat2019",
    "TUM_missing_first": "Gessulat2019",
    # Wilhelm2021 (native)
    "TUM_HLA": "Wilhelm2021",
    "TUM_HLA2": "Wilhelm2021",
    "TUM_aspn": "Wilhelm2021",
    "TUM_lysn": "Wilhelm2021",
}

_POOL_SUFFIX_RE = re.compile(r"_\d+$")


def pool_prefix(pool: str) -> str:
    """Pool prefix (the pool name minus its trailing ``_<n>``)."""
    return _POOL_SUFFIX_RE.sub("", pool)


def resolve_pool_dataset(pool: str) -> str:
    """Return the dataset whose reference holds ``pool``'s peptides.

    Resolves by pool-name prefix (e.g. a Wilhelm2021 ``TUM_first_pool_107``
    no-inclusion run → ``Zolg2017``). Raises ``KeyError`` on an unknown prefix.
    """
    pref = pool_prefix(pool)
    try:
        return POOL_PREFIX_DATASET[pref]
    except KeyError:
        raise KeyError(
            f"unknown pool prefix {pref!r} (from pool {pool!r}); add it to POOL_PREFIX_DATASET"
        ) from None


def load_pool_targets(pools: "Iterable[str]") -> pa.Table:
    """Peptides for the given pools, each sourced from its **owning** dataset.

    The cross-dataset-aware entry point for run→pool→peptide joins: pass the pool
    tokens a set of acquisitions belongs to (incl. Wilhelm2021 no-inclusion pools)
    and get back the union of their peptides — Wilhelm's borrowed Zolg pools
    included. Each row carries its native columns plus a ``source_dataset`` column
    naming the reference it came from. Loads each needed reference once.
    """
    wanted = sorted(set(pools))
    by_ds: dict[str, set[str]] = {}
    for p in wanted:
        by_ds.setdefault(resolve_pool_dataset(p), set()).add(p)

    parts: list[pa.Table] = []
    for ds, ds_pools in sorted(by_ds.items()):
        t = load_reference(ds)
        sel = t.filter(pc.is_in(t.column("pool"), value_set=pa.array(sorted(ds_pools))))
        sel = sel.append_column("source_dataset", pa.array([ds] * sel.num_rows, type=pa.string()))
        parts.append(sel)
    if not parts:
        empty = PEPTIDE_REFERENCE_TABLE.append(
            pa.field("source_dataset", pa.string(), nullable=False)
        )
        return empty.empty_table()
    return pa.concat_tables(parts)


__all__ = [
    "PEPTIDE_REFERENCE_SCHEMA_VERSION",
    "PEPTIDE_REFERENCE_TABLE",
    "POOL_PREFIX_DATASET",
    "load_pool_targets",
    "load_reference",
    "pool_prefix",
    "reference_path",
    "resolve_pool_dataset",
    "write_reference",
]
