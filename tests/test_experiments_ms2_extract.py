"""Tests for the pure transforms in protostar.experiments.ms2_extract
(the IO/extraction glue is exercised on real bundles, not here)."""

from __future__ import annotations

import pyarrow as pa
import torch
from constellation.massspec.quant.schemas import XIC_TARGET_TABLE

from protostar.experiments import ms2_extract


def test_analyzer_tolerance():
    assert ms2_extract.analyzer_tolerance("FTMS") == (20.0, "ppm")
    assert ms2_extract.analyzer_tolerance("ITMS") == (0.5, "Da")
    assert ms2_extract.analyzer_tolerance("FTMS_HCD") == (20.0, "ppm")
    assert ms2_extract.analyzer_tolerance("ITMS_CID") == (0.5, "Da")
    assert ms2_extract.analyzer_tolerance("") == (20.0, "ppm")  # default


def test_xic_targets_from_psms():
    psms = pa.table(
        {
            "modified_sequence": ["PEPTIDEK", "ELVISK"],
            "charge": pa.array([2, 3], pa.int8()),
            "mz": [500.0, 333.0],
            "retention_time_s": [60.0, 75.0],
            "scan": pa.array([1001, 2002], pa.int32()),
        }
    )
    t = ms2_extract.xic_targets_from_psms(psms)
    assert t.schema.names == XIC_TARGET_TABLE.names
    assert t.column("target_id").to_pylist() == [0, 1]
    assert t.column("scan").to_pylist() == [1001, 2002]
    assert t.column("precursor_charge").to_pylist() == [2, 3]


def test_trace_to_scan_spectra_groups_by_scan():
    trace = pa.table(
        {
            "scan": pa.array([1, 1, 2], pa.int32()),
            "mz_theoretical": [100.0, 200.0, 150.0],
            "intensity": [10.0, 20.0, 5.0],
        }
    )
    out = ms2_extract.trace_to_scan_spectra(trace)
    assert set(out) == {1, 2}
    mz1, in1 = out[1]
    assert torch.allclose(mz1, torch.tensor([100.0, 200.0], dtype=torch.float64))
    assert torch.allclose(in1, torch.tensor([10.0, 20.0], dtype=torch.float64))
    assert out[2][1].tolist() == [5.0]


def test_trace_to_scan_spectra_min_intensity():
    trace = pa.table(
        {
            "scan": pa.array([1, 2], pa.int32()),
            "mz_theoretical": [100.0, 150.0],
            "intensity": [10.0, 5.0],
        }
    )
    out = ms2_extract.trace_to_scan_spectra(trace, min_intensity=6.0)
    assert set(out) == {1}  # scan 2 dropped (intensity 5 <= 6)
