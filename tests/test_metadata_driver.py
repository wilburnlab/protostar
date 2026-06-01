"""Tests for the stage-20 metadata driver (bundle enumeration + acquisition table)."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from protostar.metadata import driver


def _write_bundle(
    data_root: Path,
    dataset: str,
    stem: str,
    *,
    creation_date: str | None = "06/03/2016 19:50:54",
    instrument_serial: str | None = "FSN20108",
    instrument_model: str | None = "Orbitrap Fusion Lumos",
    mode: str = "centroid",
    with_manifest: bool = True,
) -> Path:
    """Create a minimal convert bundle: manifest.json + acquisition_metadata.parquet."""
    bundle = data_root / "proc" / dataset / mode / stem
    bundle.mkdir(parents=True)
    if with_manifest:
        (bundle / "manifest.json").write_text(
            json.dumps({"source_file": f"/data/raw/{dataset}/{stem}.raw"})
        )
    meta = pa.table(
        {
            "creation_date": pa.array([creation_date], type=pa.string()),
            "instrument_serial": pa.array([instrument_serial], type=pa.string()),
            "instrument_model": pa.array([instrument_model], type=pa.string()),
        }
    )
    pq.write_table(meta, bundle / "acquisition_metadata.parquet")
    return bundle


def test_enumerate_bundles_sorted_and_manifest_gated(tmp_path):
    _write_bundle(tmp_path, "Zolg2017", "pool_C")
    _write_bundle(tmp_path, "Zolg2017", "pool_A")
    _write_bundle(tmp_path, "Zolg2017", "pool_B", with_manifest=False)  # no manifest → skipped
    names = [p.name for p in driver.enumerate_bundles(tmp_path, "Zolg2017")]
    assert names == ["pool_A", "pool_C"]


def test_enumerate_bundles_missing_dir_is_empty(tmp_path):
    assert driver.enumerate_bundles(tmp_path, "Nope") == []


def test_normalize_datetime_mmddyyyy_to_iso():
    assert driver.normalize_datetime("06/03/2016 19:50:54") == "2016-06-03T19:50:54"


def test_normalize_datetime_none_passthrough():
    assert driver.normalize_datetime(None) is None


def test_normalize_datetime_malformed_raises():
    with pytest.raises(ValueError, match="does not match expected format"):
        driver.normalize_datetime("2016-06-03 19:50:54")  # ISO, not MM/DD/YYYY


def test_read_acquisition_record_maps_fields(tmp_path):
    b = _write_bundle(tmp_path, "Zolg2017", "pool_A")
    rec = driver.read_acquisition_record(b, 7)
    assert rec["acquisition_id"] == 7
    assert rec["source_file"] == "pool_A.raw"  # basename of manifest source_file
    assert rec["source_kind"] == "thermo_raw"
    assert rec["acquisition_datetime"] == "2016-06-03T19:50:54"
    assert rec["instrument_serial"] == "FSN20108"
    assert "acquisition_order" not in rec  # filled by constellation, not here


def test_build_acquisitions_stable_ids_and_per_instrument_order(tmp_path):
    # Two instruments; datetimes out of input/stem order to prove chronological ranking.
    _write_bundle(
        tmp_path, "Zolg2017", "run_A2", creation_date="06/02/2016 10:00:00", instrument_serial="A"
    )
    _write_bundle(
        tmp_path, "Zolg2017", "run_A1", creation_date="06/01/2016 10:00:00", instrument_serial="A"
    )
    _write_bundle(
        tmp_path, "Zolg2017", "run_B1", creation_date="06/05/2016 10:00:00", instrument_serial="B"
    )
    acq = driver.build_acquisitions(tmp_path, "Zolg2017")
    rows = {r["source_file"]: r for r in acq.table.to_pylist()}
    # acquisition_id is the stem-sorted index: run_A1=0, run_A2=1, run_B1=2
    assert rows["run_A1.raw"]["acquisition_id"] == 0
    assert rows["run_A2.raw"]["acquisition_id"] == 1
    assert rows["run_B1.raw"]["acquisition_id"] == 2
    # order is chronological within instrument: A1 (earlier) → 1, A2 → 2; B1 → 1
    assert rows["run_A1.raw"]["acquisition_order"] == 1
    assert rows["run_A2.raw"]["acquisition_order"] == 2
    assert rows["run_B1.raw"]["acquisition_order"] == 1
    # ISO datetimes
    assert rows["run_A1.raw"]["acquisition_datetime"] == "2016-06-01T10:00:00"
    # unique ids
    ids = acq.table.column("acquisition_id").to_pylist()
    assert len(set(ids)) == len(ids)


def test_build_acquisitions_null_datetime_tolerated(tmp_path):
    _write_bundle(tmp_path, "Zolg2017", "good", creation_date="06/01/2016 10:00:00")
    _write_bundle(tmp_path, "Zolg2017", "undated", creation_date=None)
    acq = driver.build_acquisitions(tmp_path, "Zolg2017")
    s = driver.summarize(acq, n_bundles=2)
    assert s.n_acquisitions == 2
    assert s.n_datetime_null == 1


def test_write_and_present_round_trip(tmp_path):
    _write_bundle(tmp_path, "Zolg2017", "pool_A")
    assert not driver.acquisitions_present(tmp_path, "Zolg2017")
    acq = driver.build_acquisitions(tmp_path, "Zolg2017")
    out = driver.acquisitions_out_path(tmp_path, "Zolg2017")
    driver.write_acquisitions(acq, out)
    assert driver.acquisitions_present(tmp_path, "Zolg2017")
    reloaded = pq.read_table(out)
    assert reloaded.num_rows == 1
    assert "acquisition_order" in reloaded.column_names
