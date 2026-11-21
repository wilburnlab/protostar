"""Tests for the convert driver: file selection, v3 state routing, result reconciliation,
and a contract test on Constellation's ``convert_batch`` keyword API."""

from __future__ import annotations

import inspect
from pathlib import Path

from _bundle_fixtures import ms1_scan, write_legacy_bundle, write_run_bundle
from constellation.massspec.readers.thermo import BatchResult, convert_batch

from protostar import bundles
from protostar.convert import driver

DS = "Zolg2017"


def test_enumerate_raw_sorted_and_filtered(tmp_path):
    rd = tmp_path / "raw" / DS
    rd.mkdir(parents=True)
    for n in ("pool_C.raw", "pool_A.raw", "pool_B.raw", "notes.txt", "UP.RAW"):
        (rd / n).write_bytes(b"")
    names = [p.name for p in driver.enumerate_raw(tmp_path, DS)]
    assert names == ["UP.RAW", "pool_A.raw", "pool_B.raw", "pool_C.raw"]


def test_plan_conversion_shard_and_limit():
    paths = [Path(f"p{i}.raw") for i in range(10)]
    s0 = driver.plan_conversion(paths, shard=0, n_shards=3)
    s1 = driver.plan_conversion(paths, shard=1, n_shards=3)
    s2 = driver.plan_conversion(paths, shard=2, n_shards=3)
    assert sorted(p.name for p in s0 + s1 + s2) == sorted(p.name for p in paths)
    assert set(s0).isdisjoint(s1) and set(s1).isdisjoint(s2)
    assert len(driver.plan_conversion(paths, limit=4)) == 4


def test_run_conversion_empty_is_noop(tmp_path):
    assert driver.run_conversion([], tmp_path) == []


def test_bundle_out_parent_is_proc_dataset(tmp_path):
    # v3: one pass writes <stem>/ and <stem>.profile/ directly under proc/<dataset>
    assert driver.bundle_out_parent(tmp_path, DS) == tmp_path / "proc" / DS


def _layout(tmp_path):
    out = driver.bundle_out_parent(tmp_path, DS)
    write_run_bundle(out, "done", [ms1_scan(1, 1.0, [500.0], [1.0])])
    write_legacy_bundle(out, "stale")  # pre-v3 manifest at the v3 location
    (out / "partial").mkdir()
    # a June-layout bundle of "todo" must NOT count as converted
    write_legacy_bundle(bundles.legacy_dir(tmp_path, DS, "centroid"), "todo")
    paths = [Path(f"{s}.raw") for s in ("done", "stale", "partial", "todo")]
    return out, paths


def test_classify_routes_by_v3_state(tmp_path):
    out, paths = _layout(tmp_path)
    plan = driver.classify(paths, out)
    assert [p.stem for p in plan.current] == ["done"]
    assert [p.stem for p in plan.outdated] == ["stale"]
    assert sorted(p.stem for p in plan.pending) == ["partial", "todo"]
    assert sorted(p.stem for p in plan.to_convert) == ["partial", "stale", "todo"]


def test_run_conversion_forces_only_outdated(tmp_path, monkeypatch):
    out, paths = _layout(tmp_path)
    calls = []

    def fake_convert_batch(paths, output_parent, **kw):
        calls.append(([p.stem for p in paths], kw))
        return [BatchResult(p, output_parent / p.stem, "ok", None) for p in paths]

    monkeypatch.setattr(driver, "convert_batch", fake_convert_batch)
    results = driver.run_conversion(paths, out, n_workers=4, precursor_major=False)
    assert len(results) == 4
    by_force = {kw["force"]: stems for stems, kw in calls}
    assert sorted(by_force[False]) == ["done", "partial", "todo"]  # Constellation skips "done"
    assert by_force[True] == ["stale"]  # never reported "already converted"
    for _, kw in calls:
        assert kw["keep_profile"] is True and kw["precursor_major"] is False
        assert kw["n_workers"] == 4 and "profile" not in kw


def test_run_conversion_force_is_one_call(tmp_path, monkeypatch):
    out, paths = _layout(tmp_path)
    calls = []
    monkeypatch.setattr(
        driver,
        "convert_batch",
        lambda paths, output_parent, **kw: calls.append(kw["force"]) or [],
    )
    driver.run_conversion(paths, out, force=True)
    assert calls == [True]


def test_convert_batch_keyword_contract():
    """The driver's keyword set must exist on Constellation's ``convert_batch``.
    (The June ``profile=`` call broke silently because nothing exercised it.)"""
    params = inspect.signature(convert_batch).parameters
    for kw in (
        "n_workers",
        "force",
        "rt_bin_width_s",
        "keep_profile",
        "centroid_algo",
        "precursor_major",
        "batch_size",
        "progress_cb",
    ):
        assert kw in params, kw
    assert "profile" not in params


def test_summarize_results():
    results = [
        BatchResult(Path("a.raw"), Path("a"), "ok", None),
        BatchResult(Path("b.raw"), Path("b"), "skipped", None),
        BatchResult(Path("c.raw"), Path("c"), "error", "RuntimeError: boom"),
    ]
    s = driver.summarize_results(results)
    assert (s.n_ok, s.n_skipped, s.n_error, s.n_total) == (1, 1, 1, 3)
    assert s.errors == [("c.raw", "RuntimeError: boom")]


def test_converted_stems_counts_only_current_v3(tmp_path):
    _layout(tmp_path)
    assert driver.converted_stems(tmp_path, DS) == {"done"}
