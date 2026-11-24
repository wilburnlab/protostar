"""The ProteomeTools search space: every synthesized sequence, with mass and composition.

The global extraction (``docs/extraction_design.md``, once written) queries **every** sequence
in **every** run, so the search space is the union of the three datasets' peptide references
(``protostar.peptides``) plus PROCAL. Pool membership is *not* stored here; it is per run and
lives in ``protostar.pools``.

Normalisation. The published supplements spell sequences several ways. Each becomes one
ProForma 2.0 string:

- Heavy-labelled standards: ``-K8-`` / ``K*`` → ``K[UNIMOD:259]`` (13C6 15N2), ``-R10-`` /
  ``R*`` → ``R[UNIMOD:267]`` (13C6 15N4). These are the Pierce retention-time standards spiked into
  every Zolg2017 and Gessulat2019 pool.
- Cysteine is carbamidomethylated (``C[UNIMOD:4]``), as in the MaxQuant searches.
- A lowercase residue is a typo in the supplement and is upper-cased.

Masses and compositions come from Constellation (``core.sequence.protein``). Two sequences are
**composition twins** when their elemental compositions and monoisotopic masses agree (``NG`` and
``GGG``, ``I``/``L`` swaps, permutations). In MS1 they are indistinguishable.

Glue only: no scoring or matching. Targets are handed to Constellation's sweep as an
``XIC_TARGET_TABLE``.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
from constellation.core.sequence.proforma import parse_proforma
from constellation.core.sequence.protein import peptide_compositions, peptide_masses
from constellation.massspec.quant.schemas import XIC_TARGET_TABLE

from protostar.peptides.procal import procal_sequences
from protostar.peptides.reference import load_reference

#: Proton mass (Da), for precursor m/z.
PROTON_MASS: float = 1.007276466812
#: Mass difference between the M+1 and M isotopologues (13C − 12C), Da.
C13_SPACING: float = 1.0033548378

#: The three datasets whose references make up the search space.
DATASETS: tuple[str, ...] = ("Zolg2017", "Gessulat2019", "Wilhelm2021")

#: Charges the global extraction considers (PI, 2026-11-21).
DEFAULT_CHARGES: tuple[int, ...] = (1, 2, 3, 4, 5, 6)

#: Fixed and label modifications, as ProForma.
CARBAMIDOMETHYL_C = "C[UNIMOD:4]"
HEAVY_K = "K[UNIMOD:259]"
HEAVY_R = "R[UNIMOD:267]"

#: ``target_id = sequence_id * TARGET_ID_STRIDE + charge``, stable across runs.
TARGET_ID_STRIDE: int = 16

SEARCH_SPACE_SCHEMA_VERSION: int = 1

#: One row per distinct normalised sequence.
SEARCH_SPACE_TABLE: pa.Schema = pa.schema(
    [
        # Row index in the table, sorted by ``proforma``; stable for a given set of references.
        pa.field("sequence_id", pa.int64(), nullable=False),
        pa.field("proforma", pa.string(), nullable=False),
        # Bare residues, upper-case, no modifications.
        pa.field("residues", pa.string(), nullable=False),
        pa.field("monoisotopic_mass", pa.float64(), nullable=False),
        # Equal for composition twins: elemental composition plus mass (labels change mass only).
        pa.field("composition_key", pa.string(), nullable=False),
        pa.field("n_met", pa.int16(), nullable=False),
        pa.field("length", pa.int16(), nullable=False),
        pa.field("is_heavy", pa.bool_(), nullable=False),
        pa.field("is_procal", pa.bool_(), nullable=False),
        # Datasets whose reference lists the sequence (PROCAL alone: empty).
        pa.field("datasets", pa.list_(pa.string()), nullable=False),
    ],
    metadata={
        b"schema_name": b"ProtostarSearchSpace",
        b"schema_version": str(SEARCH_SPACE_SCHEMA_VERSION).encode("utf-8"),
    },
)

_HEAVY_DASHED = re.compile(r"-(K8|R10)-")
_HEAVY_STAR = re.compile(r"([KR])\*")
_RESIDUES = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def normalize_sequence(sequence: str) -> str:
    """The ProForma string for a sequence as a supplement spells it.

    >>> normalize_sequence("DIPVP-K8-P-K8-")
    'DIPVPK[UNIMOD:259]PK[UNIMOD:259]'
    >>> normalize_sequence("GILFVGSGVSGGEEGAR*")
    'GILFVGSGVSGGEEGAR[UNIMOD:267]'
    >>> normalize_sequence("ACK")
    'AC[UNIMOD:4]K'
    """
    s = _HEAVY_DASHED.sub(lambda m: "\x00" if m.group(1) == "K8" else "\x01", sequence)
    s = _HEAVY_STAR.sub(lambda m: "\x00" if m.group(1) == "K" else "\x01", s)
    s = s.upper()
    bare = s.replace("\x00", "K").replace("\x01", "R")
    if not _RESIDUES.match(bare):
        raise ValueError(f"unrecognised residues in sequence {sequence!r}")
    s = s.replace("C", CARBAMIDOMETHYL_C)
    return s.replace("\x00", HEAVY_K).replace("\x01", HEAVY_R)


def residues_of(proforma: str) -> str:
    """The bare residues of a normalised ProForma string."""
    return re.sub(r"\[[^\]]*\]", "", proforma)


def _masses_and_keys(proformas: Sequence[str]) -> tuple[np.ndarray, list[str]]:
    pfs = [parse_proforma(p) for p in proformas]
    masses = peptide_masses(pfs).numpy()
    comps = peptide_compositions(pfs).numpy()
    keys = [f"{c.tobytes().hex()}|{m:.4f}" for c, m in zip(comps, masses, strict=True)]
    return masses, keys


def space_table(owners: dict[str, set[str]], procal: set[str] = frozenset()) -> pa.Table:
    """A ``SEARCH_SPACE_TABLE`` from normalised ProForma → owning datasets (and the PROCAL set)."""
    proformas = sorted(set(owners) | set(procal))
    masses, keys = _masses_and_keys(proformas)
    residues = [residues_of(p) for p in proformas]
    return pa.table(
        {
            "sequence_id": pa.array(np.arange(len(proformas), dtype=np.int64)),
            "proforma": pa.array(proformas, pa.string()),
            "residues": pa.array(residues, pa.string()),
            "monoisotopic_mass": pa.array(masses, pa.float64()),
            "composition_key": pa.array(keys, pa.string()),
            "n_met": pa.array([r.count("M") for r in residues], pa.int16()),
            "length": pa.array([len(r) for r in residues], pa.int16()),
            "is_heavy": pa.array([HEAVY_K in p or HEAVY_R in p for p in proformas]),
            "is_procal": pa.array([p in procal for p in proformas]),
            "datasets": pa.array(
                [sorted(owners.get(p, ())) for p in proformas], pa.list_(pa.string())
            ),
        },
        schema=SEARCH_SPACE_TABLE,
    )


def build_search_space(
    datasets: Iterable[str] = DATASETS, *, include_procal: bool = True
) -> pa.Table:
    """Every distinct normalised sequence of the given datasets' references (plus PROCAL)."""
    owners: dict[str, set[str]] = {}
    for ds in datasets:
        for seq in set(load_reference(ds).column("sequence").to_pylist()):
            owners.setdefault(normalize_sequence(seq), set()).add(ds)
    procal = {normalize_sequence(s) for s in procal_sequences()} if include_procal else set()
    return space_table(owners, procal)


def sequence_ids(space: pa.Table, sequences: Iterable[str]) -> np.ndarray:
    """``sequence_id`` of each supplement-spelled sequence (raises on one not in ``space``)."""
    lookup = dict(
        zip(
            space.column("proforma").to_pylist(),
            space.column("sequence_id").to_pylist(),
            strict=True,
        )
    )
    out = []
    for s in sequences:
        p = normalize_sequence(s)
        if p not in lookup:
            raise KeyError(f"{s!r} ({p}) is not in the search space")
        out.append(lookup[p])
    return np.asarray(out, dtype=np.int64)


def ms1_mz_range(scan_metadata: pa.Table) -> tuple[float, float]:
    """The m/z range every MS1 scan of a run acquired: the intersection of their scan-window
    envelopes (``scan_window_lower`` / ``scan_window_upper``)."""
    ms1 = scan_metadata.filter(pc.equal(scan_metadata.column("level"), 1))
    lo = pc.max(ms1.column("scan_window_lower")).as_py()
    hi = pc.min(ms1.column("scan_window_upper")).as_py()
    if lo is None or hi is None or lo >= hi:
        raise ValueError("no MS1 scan windows recorded; pass mz_range explicitly")
    return float(lo), float(hi)


def precursor_mz(mass: np.ndarray, charge: int) -> np.ndarray:
    """Monoisotopic m/z of neutral mass ``mass`` at ``charge``."""
    return (np.asarray(mass, dtype=np.float64) + charge * PROTON_MASS) / charge


def targets_for_run(
    space: pa.Table,
    *,
    charges: Iterable[int] = DEFAULT_CHARGES,
    mz_range: tuple[float, float] | None = None,
    subset: np.ndarray | None = None,
) -> pa.Table:
    """MS1 targets, one per (sequence, charge), as an ``XIC_TARGET_TABLE`` with no RT gating.

    A (sequence, charge) is kept when its monoisotopic m/z lies inside ``mz_range`` (the run's
    MS1 acquired range, ``ms1_mz_range``); without a range every charge is kept. ``subset``
    restricts to those ``sequence_id`` rows (the space is indexed by ``sequence_id``).
    ``target_id = sequence_id * TARGET_ID_STRIDE + charge``; ``peptide_id`` is the sequence_id.
    """
    if subset is not None:
        space = space.take(pa.array(np.asarray(subset, dtype=np.int64)))
    ids = space.column("sequence_id").to_numpy()
    mass = space.column("monoisotopic_mass").to_numpy()
    proforma = space.column("proforma")
    parts = []
    for z in charges:
        if not 1 <= z < TARGET_ID_STRIDE:
            raise ValueError(f"charge {z} outside 1..{TARGET_ID_STRIDE - 1}")
        mz = precursor_mz(mass, z)
        keep = np.ones(mz.size, dtype=bool)
        if mz_range is not None:
            keep = (mz >= mz_range[0]) & (mz <= mz_range[1])
        idx = np.flatnonzero(keep)
        n = idx.size
        parts.append(
            pa.table(
                {
                    "target_id": pa.array(ids[idx] * TARGET_ID_STRIDE + z, pa.int64()),
                    "modified_sequence": proforma.take(pa.array(idx)),
                    "precursor_charge": pa.array(np.full(n, z, dtype=np.int8)),
                    "precursor_mz": pa.array(mz[idx], pa.float64()),
                    "peptide_id": pa.array(ids[idx], pa.int64()),
                }
            )
        )
    table = pa.concat_tables(parts) if parts else pa.table({"target_id": pa.array([], pa.int64())})
    missing = [f for f in XIC_TARGET_TABLE if f.name not in table.column_names]
    for f in missing:
        table = table.append_column(f.name, pa.nulls(table.num_rows, f.type))
    return table.select(XIC_TARGET_TABLE.names).cast(XIC_TARGET_TABLE)


__all__ = [
    "C13_SPACING",
    "DATASETS",
    "DEFAULT_CHARGES",
    "PROTON_MASS",
    "SEARCH_SPACE_SCHEMA_VERSION",
    "SEARCH_SPACE_TABLE",
    "TARGET_ID_STRIDE",
    "build_search_space",
    "ms1_mz_range",
    "normalize_sequence",
    "precursor_mz",
    "residues_of",
    "sequence_ids",
    "space_table",
    "targets_for_run",
]
