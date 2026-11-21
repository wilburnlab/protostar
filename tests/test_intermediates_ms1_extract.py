"""Tests for protostar.intermediates.ms1_extract — the pure pivot/join transforms.

The Constellation XIC call (``extract_ms1_isotopes``) needs a real peaks bundle, so
it is exercised on OSC; here we test the deterministic transforms that shape its
output into the score-ready wide intermediate, on a synthetic ``XIC_TRACE_TABLE``.
"""

from __future__ import annotations

import math

import pyarrow as pa
from constellation.massspec.quant.schemas import XIC_TARGET_TABLE

from protostar.intermediates import ms1_extract


def _synthetic_trace() -> pa.Table:
    """A tiny long-form MS1 trace: target 0 at z=2 (scans 100, 101) and z=3
    (scan 100, with an unmatched M+2 → intensity 0, null m/z error)."""
    rows = [
        # (target, charge, scan, rt, isotope, intensity, mz_ppm, mz_theo)
        (0, 2, 100, 300.0, 0, 1000.0, 0.5, 500.25),
        (0, 2, 100, 300.0, 1, 300.0, -0.3, 500.75),
        (0, 2, 100, 300.0, 2, 50.0, 1.2, 501.25),
        (0, 2, 101, 301.0, 0, 1200.0, 0.4, 500.25),
        (0, 2, 101, 301.0, 1, 360.0, -0.2, 500.75),
        (0, 2, 101, 301.0, 2, 60.0, 1.1, 501.25),
        (0, 3, 100, 300.0, 0, 200.0, 0.6, 333.84),
        (0, 3, 100, 300.0, 1, 60.0, -0.1, 334.17),
        (0, 3, 100, 300.0, 2, 0.0, None, 334.51),
    ]
    cols = list(zip(*rows))
    return pa.table(
        {
            "target_id": pa.array(cols[0], pa.int64()),
            "precursor_charge": pa.array(cols[1], pa.int8()),
            "scan": pa.array(cols[2], pa.int32()),
            "rt": pa.array(cols[3], pa.float64()),
            "isotope": pa.array(cols[4], pa.int8()),
            "intensity": pa.array(cols[5], pa.float64()),
            "mz_error_ppm": pa.array(cols[6], pa.float64()),
            "mz_theoretical": pa.array(cols[7], pa.float64()),
        }
    )


def _row_dicts(table: pa.Table) -> dict[tuple[int, int, int], dict]:
    d = table.to_pylist()
    return {(r["target_id"], r["precursor_charge"], r["scan"]): r for r in d}


def test_trace_to_wide_pivots_isotopes():
    wide = ms1_extract.trace_to_wide(_synthetic_trace())
    assert wide.num_rows == 3  # (0,2,100), (0,2,101), (0,3,100)
    for c in ("I_0", "I_1", "I_2", "mz_ppm_0", "mz_ppm_1", "mz_ppm_2", "sum_I", "mz_mono"):
        assert c in wide.column_names

    rows = _row_dicts(wide)
    a = rows[(0, 2, 100)]
    assert (a["I_0"], a["I_1"], a["I_2"]) == (1000.0, 300.0, 50.0)
    assert a["sum_I"] == 1350.0
    assert a["mz_mono"] == 500.25
    assert a["rt"] == 300.0
    assert (a["mz_ppm_0"], a["mz_ppm_1"], a["mz_ppm_2"]) == (0.5, -0.3, 1.2)

    b = rows[(0, 2, 101)]
    assert b["sum_I"] == 1620.0

    # The z=3 row has an unmatched M+2: intensity 0, NaN m/z error.
    c = rows[(0, 3, 100)]
    assert c["I_2"] == 0.0
    assert c["sum_I"] == 260.0
    assert math.isnan(c["mz_ppm_2"])
    assert c["mz_ppm_0"] == 0.6


def test_trace_to_wide_empty_is_typed():
    """An empty trace yields a 0-row table with the full typed schema (so shards
    concatenate cleanly)."""
    empty = _synthetic_trace().slice(0, 0)
    wide = ms1_extract.trace_to_wide(empty)
    assert wide.num_rows == 0
    assert wide.schema.field("I_0").type == pa.float64()
    assert wide.schema.field("scan").type == pa.int32()


def test_attach_scan_meta_joins_iit_and_derives_seconds():
    wide = ms1_extract.trace_to_wide(_synthetic_trace())
    scan_meta = pa.table(
        {
            "scan": pa.array([100, 101], pa.int32()),
            "level": pa.array([1, 1], pa.int8()),
            "iit": pa.array([50.0, 60.0], pa.float32()),  # ms
            "tic": pa.array([1.0e6, 1.1e6], pa.float64()),
            "analyzer": pa.array(["FTMS", "FTMS"], pa.string()),
            "peak_count": pa.array([1200, 1300], pa.uint32()),
        }
    )
    joined = ms1_extract.attach_scan_meta(wide, scan_meta)
    assert "iit_s" in joined.column_names
    assert "tic" in joined.column_names
    by_scan = {r["scan"]: r for r in joined.to_pylist()}
    assert by_scan[100]["iit"] == 50.0
    assert math.isclose(by_scan[100]["iit_s"], 0.05, rel_tol=1e-9)
    assert math.isclose(by_scan[101]["iit_s"], 0.06, rel_tol=1e-9)
    assert by_scan[100]["analyzer"] == "FTMS"


def test_attach_scan_meta_tolerates_missing_fields():
    """Only the scan-meta fields that exist are joined (older bundles)."""
    wide = ms1_extract.trace_to_wide(_synthetic_trace())
    scan_meta = pa.table(
        {"scan": pa.array([100, 101], pa.int32()), "tic": pa.array([1.0, 2.0], pa.float64())}
    )
    joined = ms1_extract.attach_scan_meta(wide, scan_meta)
    assert "tic" in joined.column_names
    assert "iit_s" not in joined.column_names  # no iit column → no derived seconds


def test_ms1_targets_from_psms_shape():
    psms = pa.table(
        {
            "modified_sequence": pa.array(["PEPTIDEK", "ELVISLIVESK"], pa.string()),
            "charge": pa.array([2, 3], pa.int64()),
            "retention_time_s": pa.array([300.0, 600.0], pa.float64()),
            "scan": pa.array([100, 250], pa.int64()),
        }
    )
    targets = ms1_extract.ms1_targets_from_psms(psms)
    # Cast to the canonical target contract (nullable extras filled).
    assert targets.schema.equals(XIC_TARGET_TABLE)
    assert targets.column("target_id").to_pylist() == [0, 1]
    assert targets.column("rt_center").to_pylist() == [300.0, 600.0]
    assert targets.column("precursor_charge").to_pylist() == [2, 3]


def test_attach_scan_meta_flags_iit_at_max():
    wide = ms1_extract.trace_to_wide(_synthetic_trace())
    scan_meta = pa.table(
        {
            "scan": pa.array([100, 101], pa.int32()),
            "iit": pa.array([50.0, 12.0], pa.float32()),
            "max_ion_time_ms": pa.array([50.0, 50.0], pa.float64()),
        }
    )
    by_scan = {r["scan"]: r for r in ms1_extract.attach_scan_meta(wide, scan_meta).to_pylist()}
    assert by_scan[100]["iit_at_max"] is True  # injection hit the time limit
    assert by_scan[101]["iit_at_max"] is False


def test_extract_ms1_isotopes_on_v3_bundle(tmp_path):
    """End to end on the real v3 layout: a planted z=2 envelope in three MS1 scans,
    extracted from the bundle directory (streamed), pivoted and joined to scan meta."""
    from _bundle_fixtures import ms1_scan, write_run_bundle
    from constellation.core.sequence.proforma import parse_proforma
    from constellation.massspec.peptide.envelope import peptide_envelope

    from protostar import bundles

    mz, frac = peptide_envelope(parse_proforma("ALFSSITDSEK"), charge=2, n_peaks=3)
    mz, frac = mz.tolist(), frac.tolist()
    scans = [
        ms1_scan(s, rt, mz, [scale * f for f in frac], iit=iit, max_ion_time_ms=50.0)
        for s, rt, scale, iit in (
            (10, 299.0, 1e5, 20.0),
            (11, 300.0, 3e5, 50.0),
            (12, 301.0, 1e5, 22.0),
        )
    ]
    bundle = write_run_bundle(bundles.proc_dir(tmp_path, "Zolg2017"), "run_A", scans)
    psms = pa.table(
        {
            "modified_sequence": ["ALFSSITDSEK"],
            "charge": pa.array([2], pa.int64()),
            "retention_time_s": [300.0],
            "scan": pa.array([11], pa.int64()),
        }
    )
    targets = ms1_extract.ms1_targets_from_psms(psms)
    trace = ms1_extract.extract_ms1_isotopes(bundle, targets, charge_range=(2, 2), rt_window=5.0)
    wide = ms1_extract.attach_scan_meta(
        ms1_extract.trace_to_wide(trace), bundles.read_scan_metadata(bundle)
    )
    rows = {r["scan"]: r for r in wide.to_pylist()}
    assert sorted(rows) == [10, 11, 12]
    mid = rows[11]
    assert math.isclose(mid["sum_I"], 3e5, rel_tol=1e-6)
    assert math.isclose(mid["I_1"] / mid["I_0"], frac[1] / frac[0], rel_tol=1e-6)
    assert abs(mid["mz_ppm_0"]) < 1e-3  # planted exactly on the theoretical m/z
    assert mid["iit_at_max"] is True and rows[10]["iit_at_max"] is False
    assert "intensity_floor" in wide.column_names
