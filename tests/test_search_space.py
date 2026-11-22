"""Tests for the ProteomeTools search space: normalisation, masses, twins, targets."""

from __future__ import annotations

import numpy as np
import pyarrow as pa
import pytest
from constellation.core.sequence.proforma import parse_proforma
from constellation.massspec.peptide.mz import precursor_mz as cn_precursor_mz
from constellation.massspec.quant.schemas import XIC_TARGET_TABLE

from protostar import search_space as ss
from protostar.peptides.reference import load_reference


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("PEPTIDEK", "PEPTIDEK"),
        ("ACK", "AC[UNIMOD:4]K"),
        ("DIPVP-K8-P-K8-", "DIPVPK[UNIMOD:259]PK[UNIMOD:259]"),
        ("GILFVGSGVSGGEEGA-R10-", "GILFVGSGVSGGEEGAR[UNIMOD:267]"),
        ("GILFVGSGVSGGEEGAR*", "GILFVGSGVSGGEEGAR[UNIMOD:267]"),
        ("DLPlIIQRL", "DLPLIIQRL"),
    ],
)
def test_normalize_sequence(raw, expected):
    assert ss.normalize_sequence(raw) == expected


def test_normalize_rejects_unknown_residues():
    with pytest.raises(ValueError):
        ss.normalize_sequence("PEPTIDEB")


def test_both_heavy_notations_are_one_sequence():
    assert ss.normalize_sequence("GILFVGSGVSGGEEGA-R10-") == ss.normalize_sequence(
        "GILFVGSGVSGGEEGAR*"
    )


def _space(proformas, procal=()):
    return ss.space_table({p: {"Zolg2017"} for p in proformas}, set(procal))


def test_space_table_masses_and_composition_twins():
    space = _space(
        ["NG", "GGG", "PEPTIDEK", "PEPTLDEK", "GILFVGSGVSGGEEGAR", "GILFVGSGVSGGEEGAR[UNIMOD:267]"]
    )
    rows = {r["proforma"]: r for r in space.to_pylist()}
    # N = GG in composition, so NG and GGG are twins; I/L swaps are twins
    assert rows["NG"]["composition_key"] == rows["GGG"]["composition_key"]
    assert rows["PEPTIDEK"]["composition_key"] == rows["PEPTLDEK"]["composition_key"]
    # a heavy label changes mass only, and must not make light and heavy twins
    light, heavy = rows["GILFVGSGVSGGEEGAR"], rows["GILFVGSGVSGGEEGAR[UNIMOD:267]"]
    assert light["composition_key"] != heavy["composition_key"]
    assert heavy["monoisotopic_mass"] - light["monoisotopic_mass"] == pytest.approx(
        10.008269, abs=1e-5
    )
    assert heavy["is_heavy"] and not light["is_heavy"]
    assert list(space.column("sequence_id").to_numpy()) == list(range(space.num_rows))


def test_targets_match_constellation_precursor_mz_and_prune_by_range():
    space = _space(["PEPTIDEK", "AC[UNIMOD:4]DEFGHIKLMNPQR"])
    t = ss.targets_for_run(space, charges=(1, 2, 3, 4, 5, 6), mz_range=(350.0, 1500.0))
    assert t.schema.equals(XIC_TARGET_TABLE, check_metadata=False)
    for r in t.to_pylist():
        expected = cn_precursor_mz(
            parse_proforma(r["modified_sequence"]), charge=r["precursor_charge"]
        )
        assert r["precursor_mz"] == pytest.approx(expected, abs=1e-6)
        assert 350.0 <= r["precursor_mz"] <= 1500.0
        assert r["target_id"] == r["peptide_id"] * ss.TARGET_ID_STRIDE + r["precursor_charge"]
    # PEPTIDEK (~927 Da): z=1 (928) and z=2 (464.7) in range, z=3 (310) not
    pep = {r["precursor_charge"] for r in t.to_pylist() if r["modified_sequence"] == "PEPTIDEK"}
    assert pep == {1, 2}


def test_targets_subset():
    space = _space(["PEPTIDEK", "PEPTLDEK", "NG"])
    t = ss.targets_for_run(space, charges=(2,), subset=np.array([2]))
    assert t.column("peptide_id").to_pylist() == [2]


def test_ms1_mz_range_intersects_ms1_windows_only():
    sm = pa.table(
        {
            "level": pa.array([1, 2, 1], pa.int8()),
            "scan_window_lower": pa.array([350.0, 100.0, 360.0]),
            "scan_window_upper": pa.array([1500.0, 2000.0, 1450.0]),
        }
    )
    assert ss.ms1_mz_range(sm) == (360.0, 1450.0)


def test_every_reference_sequence_normalises():
    """Guards the supplements' spellings: every sequence of every dataset must normalise, and the
    15 heavy Pierce standards of every Gessulat2019 pool are also in Zolg2017.

    Zolg2017 has one more heavy form: its QC row spells ``DIPVPKPK*`` (C-terminal Lys heavy) where
    the pool rows spell ``DIPVP-K8-P-K8-`` (both Lys heavy). Both stay in the space as distinct
    masses; the screen shows which is present."""
    heavy: dict[str, set[str]] = {}
    for ds in ss.DATASETS:
        for seq in set(load_reference(ds).column("sequence").to_pylist()):
            p = ss.normalize_sequence(seq)
            if "UNIMOD:259" in p or "UNIMOD:267" in p:
                heavy.setdefault(ds, set()).add(p)
    assert len(heavy["Gessulat2019"]) == 15
    assert heavy["Zolg2017"] - heavy["Gessulat2019"] == {"DIPVPKPK[UNIMOD:259]"}
    assert "Wilhelm2021" not in heavy
