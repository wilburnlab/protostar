"""Tests for pools, cross-dataset run order and per-run labels."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pytest

from protostar import pools
from protostar import search_space as ss

MANIFESTS = Path(__file__).resolve().parents[1] / "config" / "manifests"


@pytest.mark.parametrize(
    ("stem", "pool", "dataset", "injection", "noincl", "part"),
    [
        (
            "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2",
            "TUM_first_pool_1",
            "Zolg2017",
            "DDA",
            False,
            "01_01",
        ),
        (
            "01640c_BD3-Thermo_SRM_Pool_20_01_01-ETD-1h-R2",
            "Thermo_SRM_Pool_20",
            "Zolg2017",
            "ETD",
            False,
            "01_01",
        ),
        (
            "01974c_BA1-TUM_missing_first_1_01_01-3xHCD-1h-R4",
            "TUM_missing_first_1",
            "Gessulat2019",
            "3xHCD",
            False,
            "01_01",
        ),
        ("02445a_BA1-TUM_HLA_1_01_01-DDA-1h-R1", "TUM_HLA_1", "Wilhelm2021", "DDA", False, "01_01"),
        (
            "03174_GF7-TUM_second_pool_122_01_01-3xHCD-1hnoincl-R1",
            "TUM_second_pool_122",
            "Zolg2017",
            "3xHCD",
            True,
            "01_01",
        ),
    ],
)
def test_parse_stem(stem, pool, dataset, injection, noincl, part):
    name = pools.parse_stem(stem + ".raw")
    assert (name.pool, name.pool_dataset, name.injection, name.noincl, name.pool_part) == (
        pool,
        dataset,
        injection,
        noincl,
        part,
    )


def test_parse_stem_hla2_prefix_is_not_hla():
    assert pools.parse_stem("02800a_BA1-TUM_HLA2_7_01_01-DDA-1h-R1").pool_prefix == "TUM_HLA2"


@pytest.mark.parametrize(
    "bad",
    [
        "not_a_stem",
        "01625b_GA1-TUM_unknown_pool_1_01_01-DDA-1h-R2",  # unknown prefix
        "01625b_GA1-TUM_first_pool_1_01_01-SIM-1h-R2",  # unknown injection
    ],
)
def test_parse_stem_rejects(bad):
    with pytest.raises(ValueError):
        pools.parse_stem(bad)


def test_every_manifest_raw_parses():
    counts = {"DDA": 0, "2xIT_2xHCD": 0, "3xHCD": 0, "ETD": 0}
    n = 0
    for path in sorted(MANIFESTS.glob("*.json")):
        for f in json.loads(path.read_text())["files"]:
            if f["file_category"] != "RAW":
                continue
            counts[pools.parse_stem(f["file_name"]).injection] += 1
            n += 1
    assert n == 4213
    assert all(v > 1000 for v in counts.values())


def _acq(rows):
    return pa.Table.from_pylist(
        [
            {
                "acquisition_id": i,
                "source_file": f"{stem}.raw",
                "acquisition_datetime": t,
                "instrument_serial": serial,
                "instrument_model": "Orbitrap Fusion Lumos",
            }
            for i, (stem, t, serial) in enumerate(rows)
        ]
    )


def _stem(pool, injection="DDA", n="01"):
    return f"0{n}000a_BA1-{pool}_01_01-{injection}-1h-R1"


@pytest.fixture
def runs():
    # One instrument shared by two datasets, interleaved in time, plus a second instrument.
    zolg = _acq(
        [
            (_stem("TUM_first_pool_1", "DDA"), "2016-01-01T00:00:00", "FSN1"),
            (_stem("TUM_first_pool_1", "3xHCD"), "2016-01-01T01:00:00", "FSN1"),
            (_stem("TUM_first_pool_2", "DDA"), "2016-01-01T03:00:00", "FSN1"),
            (_stem("TUM_first_pool_3", "DDA"), "2016-01-01T00:30:00", "FSN2"),
        ]
    )
    gess = _acq([(_stem("TUM_isoform_1", "DDA"), "2016-01-01T02:00:00", "FSN1")])
    return pools.run_table_from_acquisitions({"Zolg2017": zolg, "Gessulat2019": gess})


def test_run_table_ranks_across_datasets(runs):
    by_stem = {r["stem"]: r for r in runs.to_pylist()}
    assert by_stem[_stem("TUM_isoform_1")]["instrument_rank"] == 3
    assert by_stem[_stem("TUM_first_pool_2")]["instrument_rank"] == 4
    assert by_stem[_stem("TUM_first_pool_3")]["instrument_rank"] == 1  # other instrument


def test_previous_pools_skips_own_pool_and_crosses_datasets(runs):
    prev = pools.previous_pools(runs, _stem("TUM_first_pool_2"), k=2)
    assert [(p.pool, p.pools_back, p.runs_back) for p in prev] == [
        ("TUM_isoform_1", 1, 1),
        ("TUM_first_pool_1", 2, 2),
    ]
    assert prev[0].dt_s == pytest.approx(3600.0)
    # the 3xHCD run of pool 1 has only its own pool before it
    assert pools.previous_pools(runs, _stem("TUM_first_pool_1", "3xHCD")) == []


def test_labels_for_run(runs):
    space = ss.space_table(
        {
            "PEPTIDEK": {"Zolg2017"},  # pool 2 (the run)
            "PEPTLDEK": {"Zolg2017"},  # composition twin of PEPTIDEK, off-pool
            "MPEPTIDEK": {"Zolg2017"},  # pool 2, has Met
            "LLLLLLLK": {"Gessulat2019"},  # isoform 1 (carryover, 1 pool back)
            "AAAAAAAK": {"Zolg2017"},  # pool 1 (carryover, 2 pools back)
            "WWWWWWWK": {"Zolg2017"},  # pool 3 (other instrument): plain off-pool
            "GILFVGSGVSGGEEGAR[UNIMOD:267]": {"Zolg2017"},  # heavy standard in every pool
            "PEPTIEK": {"Zolg2017"},  # PEPTIDEK minus D: a deletion product, off-pool
        },
        procal={"YSAHEEHHYDK"},
    )
    pid = {r["proforma"]: r["sequence_id"] for r in space.to_pylist()}
    ids = lambda *ps: np.array(sorted(pid[p] for p in ps), dtype=np.int64)  # noqa: E731
    heavy = "GILFVGSGVSGGEEGAR[UNIMOD:267]"
    index = pools.PoolIndex(
        {
            "TUM_first_pool_1": ids("AAAAAAAK", heavy),
            "TUM_first_pool_2": ids("PEPTIDEK", "MPEPTIDEK", heavy),
            "TUM_first_pool_3": ids("WWWWWWWK", heavy),
            "TUM_isoform_1": ids("LLLLLLLK", heavy),
        }
    )
    labels = pools.labels_for_run(_stem("TUM_first_pool_2"), runs, space, index)
    got = {
        space.column("proforma")[r["sequence_id"]].as_py(): (r["label"], r["kind"])
        for r in labels.to_pylist()
    }
    assert got["YSAHEEHHYDK"] == ("standard", "procal")
    assert got[heavy] == ("standard", "prtc")
    assert got["PEPTIDEK"] == ("in_pool", None)
    assert got["MPEPTIDEK"] == ("in_pool", None)
    assert got["LLLLLLLK"] == ("carryover", None)
    assert got["AAAAAAAK"] == ("carryover", None)
    assert got["PEPTLDEK"] == ("shadowed", "composition")
    assert got["PEPTIEK"] == ("shadowed", "deletion")
    assert "WWWWWWWK" not in got  # off_pool is implied
    carry = {
        space.column("proforma")[r["sequence_id"]].as_py(): r
        for r in labels.to_pylist()
        if r["label"] == "carryover"
    }
    assert carry["LLLLLLLK"]["pools_back"] == 1 and carry["AAAAAAAK"]["pools_back"] == 2
    assert pools.label_counts(labels)["in_pool"] == 2


def test_labels_shadow_kinds_can_be_switched_off(runs):
    space = ss.space_table(
        {"PEPTIDEK": {"Zolg2017"}, "PEPTLDEK": {"Zolg2017"}, "PEPTIEK": {"Zolg2017"}}
    )
    pid = {r["proforma"]: r["sequence_id"] for r in space.to_pylist()}
    index = pools.PoolIndex(
        {
            "TUM_first_pool_2": np.array([pid["PEPTIDEK"]]),
            "TUM_first_pool_1": np.array([], dtype=np.int64),
            "TUM_isoform_1": np.array([], dtype=np.int64),
        }
    )
    labels = pools.labels_for_run(
        _stem("TUM_first_pool_2"), runs, space, index, shadow_kinds=("composition",)
    )
    kinds = {r["kind"] for r in labels.to_pylist() if r["label"] == "shadowed"}
    assert kinds == {"composition"}


def test_search_zip_for_stem(tmp_path):
    z = tmp_path / "search" / "Zolg2017"
    z.mkdir(parents=True)
    (z / "TUM_first_pool_1_01_01_DDA-1h-R2-tryptic.zip").write_bytes(b"")
    (z / "TUM_first_pool_1_01_01_DDA-1h-R1-tryptic.zip").write_bytes(b"")
    (z / "TUM_first_pool_107_01_01_DDA-1hnoincl-R2-semitryptic.zip").write_bytes(b"")
    g = tmp_path / "search" / "Gessulat2019"
    g.mkdir(parents=True)
    (g / "01974C_BA1-TUM_missing_first_1_01_01-DDA-1h-R4.zip").write_bytes(b"")
    hit = pools.search_zip_for_stem(
        tmp_path, "Zolg2017", "01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2"
    )
    assert hit.name == "TUM_first_pool_1_01_01_DDA-1h-R2-tryptic.zip"
    noincl = pools.search_zip_for_stem(
        tmp_path, "Zolg2017", "03100_GF1-TUM_first_pool_107_01_01-DDA-1hnoincl-R2"
    )
    assert noincl.name.endswith("noincl-R2-semitryptic.zip")
    gess = pools.search_zip_for_stem(
        tmp_path, "Gessulat2019", "01974c_BA1-TUM_missing_first_1_01_01-DDA-1h-R4"
    )
    assert gess.name.startswith("01974C")
    assert (
        pools.search_zip_for_stem(
            tmp_path, "Zolg2017", "01625b_GA1-TUM_first_pool_1_01_01-ETD-1h-R2"
        )
        is None
    )
