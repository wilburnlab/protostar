"""Tests for manifest round-trip, routing, and reconciliation."""

from __future__ import annotations

from pathlib import Path

from protostar.fetch import manifest
from protostar.fetch.manifest import Manifest, ManifestEntry


def _toy_manifest():
    raw = ManifestEntry("pool_A.raw", 100, "a" * 40, "https://x/pool_A.raw", "RAW")
    srch = ManifestEntry("pool_A.zip", 50, "b" * 40, "https://x/pool_A.zip", "SEARCH")
    other = ManifestEntry("sdrf.tsv", 10, "c" * 40, "https://x/sdrf.tsv", "EXPERIMENTAL DESIGN")
    return Manifest("Zolg2017", "PXD004732", "2026-01-01T00:00:00+00:00", "v3", (raw, srch, other))


def test_round_trip_and_counts(tmp_path):
    m = _toy_manifest()
    assert Manifest.from_dict(m.to_dict()) == m
    p = manifest.write_manifest(m, tmp_path / "Zolg2017.json")
    assert manifest.load_manifest(p) == m
    assert m.to_dict()["n_by_category"] == {"RAW": 1, "SEARCH": 1, "EXPERIMENTAL DESIGN": 1}


def test_routing():
    m = _toy_manifest()
    raw, srch, other = m.entries
    assert manifest.local_path(raw, "/data", "Zolg2017") == Path("/data/raw/Zolg2017/pool_A.raw")
    assert manifest.local_path(srch, "/data", "Zolg2017") == Path(
        "/data/search/Zolg2017/pool_A.zip"
    )
    assert manifest.local_path(other, "/data", "Zolg2017") == Path("/data/other/Zolg2017/sdrf.tsv")


def test_fetch_entries_filter():
    m = _toy_manifest()
    assert len(m.fetch_entries(("RAW", "SEARCH"))) == 2
    assert len(m.fetch_entries(("RAW",))) == 1
    assert len(m.fetch_entries(None)) == 3


def test_reconcile_states(tmp_path):
    m = _toy_manifest()
    data_root = tmp_path / "data"

    statuses = manifest.reconcile(m, data_root)  # default RAW+SEARCH
    assert {s.state for s in statuses} == {"missing"}
    assert len(statuses) == 2  # the EXPERIMENTAL DESIGN entry isn't in default categories

    raw = m.entries[0]
    p = manifest.local_path(raw, data_root, "Zolg2017")
    p.parent.mkdir(parents=True)
    p.write_bytes(b"x" * 100)  # correct size
    by_name = {s.entry.file_name: s for s in manifest.reconcile(m, data_root)}
    assert by_name["pool_A.raw"].state == "ok"
    assert by_name["pool_A.zip"].state == "missing"

    p.write_bytes(b"x" * 99)  # wrong size → corrupt
    assert {s.entry.file_name: s.state for s in manifest.reconcile(m, data_root)}[
        "pool_A.raw"
    ] == "corrupt"


def test_summarize():
    m = _toy_manifest()
    statuses = manifest.reconcile(m, "/nonexistent")
    table = manifest.summarize(statuses)
    assert table == {"RAW": {"missing": 1}, "SEARCH": {"missing": 1}}
