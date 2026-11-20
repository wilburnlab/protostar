"""Tests for protostar.bundles — v3 run-bundle discovery and state classification."""

from __future__ import annotations

import pytest
from _bundle_fixtures import ms1_scan, write_legacy_bundle, write_run_bundle
from constellation.massspec.bundle import OutdatedBundleError

from protostar import bundles

DS = "Zolg2017"


def _scans():
    return [
        ms1_scan(1, 10.0, [500.0, 500.5], [1e5, 3e4], iit=20.0, max_ion_time_ms=50.0),
        ms1_scan(2, 11.0, [500.0], [2e5], iit=50.0, max_ion_time_ms=50.0),
    ]


@pytest.fixture
def proc(tmp_path):
    p = bundles.proc_dir(tmp_path, DS)
    write_run_bundle(p, "run_B", _scans())
    write_run_bundle(p, "run_A", _scans(), with_profile=False)
    write_legacy_bundle(p, "run_OLD")  # pre-v3 manifest at the v3 location
    (p / "run_PARTIAL").mkdir()  # unfinished conversion: no manifest
    write_legacy_bundle(bundles.legacy_dir(tmp_path, DS, "centroid"), "run_A")
    write_legacy_bundle(bundles.legacy_dir(tmp_path, DS, "centroid"), "run_B")
    write_legacy_bundle(bundles.legacy_dir(tmp_path, DS, "profile"), "run_A")
    return tmp_path


def test_bundle_state_classifies(proc):
    p = bundles.proc_dir(proc, DS)
    assert bundles.bundle_state(p / "run_A") == "current"
    assert bundles.bundle_state(p / "run_B.profile") == "current"
    assert bundles.bundle_state(p / "run_OLD") == "outdated"
    assert bundles.bundle_state(p / "run_PARTIAL") == "incomplete"
    assert bundles.bundle_state(p / "nope") == "missing"


def test_iter_run_bundles_only_current_centroid_runs(proc):
    names = [b.name for b in bundles.iter_run_bundles(proc, DS)]
    # no .profile siblings, no legacy centroid/ or profile/ parents, no outdated/partial
    assert names == ["run_A", "run_B"]


def test_iter_run_bundles_missing_dataset(tmp_path):
    assert list(bundles.iter_run_bundles(tmp_path, "Nope")) == []


def test_inventory_reports_everything(proc):
    inv = bundles.inventory(proc, DS)
    assert inv.current == ["run_A", "run_B"]
    assert inv.with_profile == ["run_B"]
    assert inv.outdated == ["run_OLD"]
    assert inv.incomplete == ["run_PARTIAL"]
    assert inv.legacy == {"centroid": 2, "profile": 1}
    assert "outdated=1" in inv.summary()


def test_require_run_bundle(proc):
    assert bundles.require_run_bundle(proc, DS, "run_A").name == "run_A"
    assert bundles.require_run_bundle(proc, DS, "never_converted") is None
    with pytest.raises(OutdatedBundleError):
        bundles.require_run_bundle(proc, DS, "run_OLD")
    with pytest.raises(RuntimeError, match="unfinished"):
        bundles.require_run_bundle(proc, DS, "run_PARTIAL")


def test_read_scan_metadata_projects_present_columns(proc):
    b = bundles.run_bundle_dir(proc, DS, "run_A")
    t = bundles.read_scan_metadata(b, columns=["scan", "iit", "not_a_column"])
    assert t.column_names == ["scan", "iit"]
    assert t.column("iit").to_pylist() == [20.0, 50.0]


def test_read_scan_metadata_refuses_outdated(proc):
    with pytest.raises(OutdatedBundleError):
        bundles.read_scan_metadata(bundles.proc_dir(proc, DS) / "run_OLD")


def test_legacy_dir_validates_mode(tmp_path):
    with pytest.raises(ValueError):
        bundles.legacy_dir(tmp_path, DS, "peaks")


def test_constellation_provenance_is_one_line():
    line = bundles.constellation_provenance()
    assert line.startswith("constellation ") and "\n" not in line
