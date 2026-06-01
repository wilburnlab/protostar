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

from pathlib import Path

import pyarrow as pa
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


__all__ = [
    "PEPTIDE_REFERENCE_SCHEMA_VERSION",
    "PEPTIDE_REFERENCE_TABLE",
    "load_reference",
    "reference_path",
    "write_reference",
]
