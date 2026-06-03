"""Tests for protostar.experiments.scans — PSM gating, calibrant tagging,
recurring-peptide selection, and replicate grouping over a synthetic
PSM_TABLE."""

from __future__ import annotations

import pyarrow as pa
from constellation.massspec.search.schemas import PSM_TABLE

from protostar.experiments import scans

_DEFAULTS = {
    "psm_id": 0,
    "raw_file": "runA",
    "acquisition_id": None,
    "scan": 1,
    "precursor_scan": None,
    "sequence": "PEPTIDEK",
    "modified_sequence": "PEPTIDEK",
    "peptide_id": None,
    "mod_peptide_id": None,
    "evidence_id": None,
    "proteins": None,
    "charge": 2,
    "mz": 500.0,
    "mass": 998.0,
    "mass_error_ppm": 0.0,
    "retention_time_s": 60.0,
    "fragmentation": "HCD",
    "mass_analyzer": "FTMS",
    "psm_type": "MULTI-MSMS",
    "score": 100.0,
    "delta_score": 10.0,
    "pep": 0.001,
    "is_decoy": False,
    "is_contaminant": False,
    "engine": "maxquant",
}


def _psms(rows: list[dict]) -> pa.Table:
    full = []
    for i, r in enumerate(rows):
        d = dict(_DEFAULTS)
        d["psm_id"] = i
        d.update(r)
        full.append(d)
    cols = {name: [d[name] for d in full] for name in PSM_TABLE.names}
    return pa.table(
        {k: pa.array(v, type=PSM_TABLE.field(k).type) for k, v in cols.items()},
        schema=PSM_TABLE,
    )


def test_add_mode_column():
    t = _psms(
        [
            {"mass_analyzer": "FTMS", "fragmentation": "HCD"},
            {"mass_analyzer": "ITMS", "fragmentation": "CID"},
        ]
    )
    assert scans.add_mode_column(t).column("mode").to_pylist() == ["FTMS_HCD", "ITMS_CID"]


def test_gate_psms_filters_each_reason():
    t = _psms(
        [
            {},  # passes
            {"is_decoy": True},
            {"is_contaminant": True},
            {"psm_type": "MULTI-SECPEP"},
            {"pep": 0.5},
            {"score": -1.0},
            {"pep": None},  # null pep fails the numeric cut
        ]
    )
    g = scans.gate_psms(t, max_pep=0.01, min_score=0.0)
    assert g.num_rows == 1


def test_tag_calibrants_procal():
    t = _psms([{"sequence": "YSAHEEHHYDK"}, {"sequence": "PEPTIDEK"}])  # #1 is PROCAL
    tagged = scans.tag_calibrants(t)
    assert tagged.column("is_procal").to_pylist() == [True, False]


def test_select_recurring_min_acquisitions():
    rows = [
        {"modified_sequence": "AAA", "charge": 2, "raw_file": rf} for rf in ("r1", "r2", "r3")
    ] + [{"modified_sequence": "BBB", "charge": 2, "raw_file": "r1"}]
    sel = scans.select_recurring_peptides(_psms(rows), min_acquisitions=2)
    assert set(sel.column("modified_sequence").to_pylist()) == {"AAA"}


def test_select_recurring_top_n():
    rows = [
        {"modified_sequence": "AAA", "charge": 2, "raw_file": rf} for rf in ("r1", "r2", "r3")
    ] + [{"modified_sequence": "BBB", "charge": 2, "raw_file": rf} for rf in ("r1", "r2")]
    sel = scans.select_recurring_peptides(_psms(rows), top_n=1)
    assert set(sel.column("modified_sequence").to_pylist()) == {"AAA"}


def test_replicate_groups():
    t = scans.add_mode_column(
        _psms(
            [
                {
                    "modified_sequence": "AAA",
                    "charge": 2,
                    "mass_analyzer": "FTMS",
                    "fragmentation": "HCD",
                    "raw_file": "r1",
                },
                {
                    "modified_sequence": "AAA",
                    "charge": 2,
                    "mass_analyzer": "FTMS",
                    "fragmentation": "HCD",
                    "raw_file": "r2",
                },
                {
                    "modified_sequence": "AAA",
                    "charge": 2,
                    "mass_analyzer": "ITMS",
                    "fragmentation": "CID",
                    "raw_file": "r1",
                },
            ]
        )
    )
    groups = {
        (k["modified_sequence"], k["charge"], k["mode"]): sub
        for k, sub in scans.replicate_groups(t)
    }
    assert groups[("AAA", 2, "FTMS_HCD")].num_rows == 2
    assert groups[("AAA", 2, "ITMS_CID")].num_rows == 1
