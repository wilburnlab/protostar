"""Tests for the convert driver's file selection + result reconciliation."""

from __future__ import annotations

from pathlib import Path

from constellation.massspec.readers.thermo import BatchResult

from protostar.convert import driver


def test_enumerate_raw_sorted_and_filtered(tmp_path):
    rd = tmp_path / "raw" / "Zolg2017"
    rd.mkdir(parents=True)
    for n in ("pool_C.raw", "pool_A.raw", "pool_B.raw", "notes.txt", "UP.RAW"):
        (rd / n).write_bytes(b"")
    names = [p.name for p in driver.enumerate_raw(tmp_path, "Zolg2017")]
    assert names == ["UP.RAW", "pool_A.raw", "pool_B.raw", "pool_C.raw"]


def test_plan_conversion_shard_and_limit(tmp_path):
    paths = [Path(f"p{i}.raw") for i in range(10)]
    s0 = driver.plan_conversion(paths, shard=0, n_shards=3)
    s1 = driver.plan_conversion(paths, shard=1, n_shards=3)
    s2 = driver.plan_conversion(paths, shard=2, n_shards=3)
    assert sorted(p.name for p in s0 + s1 + s2) == sorted(p.name for p in paths)
    assert set(s0).isdisjoint(s1) and set(s1).isdisjoint(s2)
    assert len(driver.plan_conversion(paths, limit=4)) == 4


def test_run_conversion_empty_is_noop(tmp_path):
    assert driver.run_conversion([], tmp_path) == []


def test_bundle_layout_separates_modes(tmp_path):
    c = driver.bundle_out_parent(tmp_path, "Zolg2017", profile=False)
    p = driver.bundle_out_parent(tmp_path, "Zolg2017", profile=True)
    assert c.name == "centroid" and p.name == "profile"
    assert c.parent == p.parent  # proc/Zolg2017/


def test_summarize_results():
    results = [
        BatchResult(Path("a.raw"), Path("a"), "ok", None),
        BatchResult(Path("b.raw"), Path("b"), "skipped", None),
        BatchResult(Path("c.raw"), Path("c"), "error", "RuntimeError: boom"),
    ]
    s = driver.summarize_results(results)
    assert (s.n_ok, s.n_skipped, s.n_error, s.n_total) == (1, 1, 1, 3)
    assert s.errors == [("c.raw", "RuntimeError: boom")]


def test_converted_stems(tmp_path):
    op = driver.bundle_out_parent(tmp_path, "Zolg2017", profile=False)
    (op / "pool_A").mkdir(parents=True)
    (op / "pool_A" / "manifest.json").write_text("{}")
    (op / "pool_B").mkdir(parents=True)  # no manifest → not counted
    assert driver.converted_stems(tmp_path, "Zolg2017", profile=False) == {"pool_A"}
    assert driver.converted_stems(tmp_path, "Zolg2017", profile=True) == set()
