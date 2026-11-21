"""Tests for the stage-20 metadata driver on bundle-v3 run bundles."""

from __future__ import annotations

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from _bundle_fixtures import ms1_scan, write_legacy_bundle, write_run_bundle

from protostar import bundles
from protostar.metadata import driver

DS = "Zolg2017"


def _bundle(root, stem, *, creation_date="06/03/2016 19:50:54", serial="FSN20108", sha=None):
    return write_run_bundle(
        bundles.proc_dir(root, DS),
        stem,
        [ms1_scan(1, 1.0, [500.0], [1.0])],
        run={"creation_date": creation_date, "instrument_serial": serial},
        sha256=sha,
    )


def test_enumerate_bundles_sorted_current_only(tmp_path):
    _bundle(tmp_path, "pool_C")
    _bundle(tmp_path, "pool_A")
    write_legacy_bundle(bundles.proc_dir(tmp_path, DS), "pool_OLD")  # pre-v3: excluded
    write_legacy_bundle(bundles.legacy_dir(tmp_path, DS), "pool_B")  # June tree: excluded
    names = [p.name for p in driver.enumerate_bundles(tmp_path, DS)]
    assert names == ["pool_A", "pool_C"]  # no .profile siblings either


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
    b = _bundle(tmp_path, "pool_A")
    rec = driver.read_acquisition_record(b)
    assert rec["source_file"] == "pool_A.raw"  # basename of the manifest's cluster path
    assert rec["source_kind"] == "thermo_acquisition"  # Constellation's bundle kind
    assert rec["source_sha256"] == "sha-of-pool_A"
    assert rec["acquisition_datetime"] == "2016-06-03T19:50:54"
    assert rec["instrument_serial"] == "FSN20108"
    assert "acquisition_order" not in rec  # filled by constellation, not here


def test_read_acquisition_record_malformed_date_fails_loud(tmp_path):
    b = _bundle(tmp_path, "pool_A", creation_date="2016-06-03 19:50:54")
    with pytest.raises(ValueError, match="does not match expected format"):
        driver.read_acquisition_record(b)


def test_build_acquisitions_stable_ids_and_per_instrument_order(tmp_path):
    # Two instruments; datetimes out of stem order to prove chronological ranking.
    _bundle(tmp_path, "run_A2", creation_date="06/02/2016 10:00:00", serial="A")
    _bundle(tmp_path, "run_A1", creation_date="06/01/2016 10:00:00", serial="A")
    _bundle(tmp_path, "run_B1", creation_date="06/05/2016 10:00:00", serial="B")
    acq = driver.build_acquisitions(tmp_path, DS)
    rows = {r["source_file"]: r for r in acq.table.to_pylist()}
    # acquisition_id follows stem order: run_A1=0, run_A2=1, run_B1=2
    assert rows["run_A1.raw"]["acquisition_id"] == 0
    assert rows["run_A2.raw"]["acquisition_id"] == 1
    assert rows["run_B1.raw"]["acquisition_id"] == 2
    # order is chronological within instrument
    assert rows["run_A1.raw"]["acquisition_order"] == 1
    assert rows["run_A2.raw"]["acquisition_order"] == 2
    assert rows["run_B1.raw"]["acquisition_order"] == 1
    assert rows["run_A1.raw"]["acquisition_datetime"] == "2016-06-01T10:00:00"


def test_build_acquisitions_null_datetime_tolerated(tmp_path):
    _bundle(tmp_path, "good", creation_date="06/01/2016 10:00:00")
    _bundle(tmp_path, "undated", creation_date=None)
    acq = driver.build_acquisitions(tmp_path, DS)
    s = driver.summarize(acq, n_bundles=2)
    assert s.n_acquisitions == 2
    assert s.n_datetime_null == 1


def test_write_and_present_round_trip(tmp_path):
    _bundle(tmp_path, "pool_A")
    assert not driver.acquisitions_present(tmp_path, DS)
    acq = driver.build_acquisitions(tmp_path, DS)
    out = driver.acquisitions_out_path(tmp_path, DS)
    driver.write_acquisitions(acq, out)
    assert driver.acquisitions_present(tmp_path, DS)
    reloaded = pq.read_table(out)
    assert reloaded.num_rows == 1
    assert "acquisition_order" in reloaded.column_names


def test_june_table_is_not_present(tmp_path):
    """A pre-v3 (June) acquisitions table lacks the bundle-schema stamp, so stage 20
    rebuilds it after reconversion instead of silently keeping it."""
    out = driver.acquisitions_out_path(tmp_path, DS)
    out.parent.mkdir(parents=True)
    pq.write_table(pa.table({"acquisition_id": [0]}), out)
    assert not driver.acquisitions_present(tmp_path, DS)
