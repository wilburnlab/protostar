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


def test_trace_to_target_spectra_groups_by_target():
    # targets 0 and 1 share scan 1; grouping must isolate them by target_id,
    # never pool by scan.
    trace = pa.table(
        {
            "target_id": pa.array([0, 0, 1], pa.int64()),
            "scan": pa.array([1, 1, 1], pa.int32()),
            "mz_theoretical": [100.0, 200.0, 150.0],
            "intensity": [10.0, 20.0, 5.0],
        }
    )
    out = ms2_extract.trace_to_target_spectra(trace)
    assert set(out) == {0, 1}
    mz0, in0 = out[0]
    assert torch.allclose(mz0, torch.tensor([100.0, 200.0], dtype=torch.float64))
    assert torch.allclose(in0, torch.tensor([10.0, 20.0], dtype=torch.float64))
    assert out[1][1].tolist() == [5.0]  # NOT merged with target 0 despite same scan


def test_trace_to_target_spectra_min_intensity():
    trace = pa.table(
        {
            "target_id": pa.array([0, 1], pa.int64()),
            "scan": pa.array([1, 2], pa.int32()),
            "mz_theoretical": [100.0, 150.0],
            "intensity": [10.0, 5.0],
        }
    )
    out = ms2_extract.trace_to_target_spectra(trace, min_intensity=6.0)
    assert set(out) == {0}  # target 1 dropped (intensity 5 <= 6)


def test_trace_to_target_channels_isolates_chimeric_scan():
    # Two PSMs (target 0, target 1) in the SAME scan, each with a fragment at the
    # SAME (ion_type, position, charge) channel but a DIFFERENT real m/z. Keying on
    # target_id must keep them apart, so a neighbour's fragment is never projected
    # onto this peptide's basis by (ion_type, position, charge) (PR #14 / Codex).
    trace = pa.table(
        {
            "target_id": pa.array([0, 1], pa.int64()),
            "scan": pa.array([5, 5], pa.int32()),
            "ion_type": pa.array([4, 4], pa.int8()),
            "position": pa.array([2, 2], pa.int32()),
            "fragment_charge": pa.array([1, 1], pa.int32()),
            "intensity": [100.0, 999.0],
            "mz_error_ppm": [1.0, -3.0],
            "isotope": pa.array([0, 0], pa.int8()),
        }
    )
    out = ms2_extract.trace_to_target_channels(trace)
    assert set(out) == {0, 1}
    assert out[0][(4, 2, 1)][0] == 100.0  # target 0 keeps only its own fragment...
    assert out[1][(4, 2, 1)][0] == 999.0  # ...not merged into 1099 by shared scan
