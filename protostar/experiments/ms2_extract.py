"""MS2 fragment-spectrum extraction for the replicate groups.

Wraps Constellation's XIC level-2 ``assigned_scans_only`` extraction (the MS2
retrieval path) into per-acquisition fragment spectra. Two facts shape the API:

  * a replicate group spans **many acquisitions** (one converted bundle per raw
    file), so extraction is per-bundle and the driver accumulates across them;
  * within one injection, FTMS (Orbitrap) and ITMS (ion-trap) MS2 scans
    interleave, so the analyzer-appropriate fragment tolerance (FTMS ~20 ppm,
    ITMS ~0.5 Da — required, per the XIC validation) is applied **per PSM**, by
    splitting on ``mass_analyzer``.

The per-scan fragment vectors are keyed back to the consensus fragment basis by
their **theoretical** m/z (an exact key — no second tolerance), so the result
feeds ``constellation.massspec.spectra.consensus.build_consensus`` directly.
No model code lives here.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence

import pyarrow as pa
import torch
from constellation.core.io.schemas import cast_to_schema
from constellation.massspec.peptide.ions import IonType
from constellation.massspec.quant.chromatogram import extract_xic_scan_major
from constellation.massspec.quant.schemas import XIC_TARGET_TABLE, XIC_TRACE_TABLE

#: Analyzer → (tolerance, unit) for MS2 fragment matching (XIC-validated).
ANALYZER_TOLERANCE: dict[str, tuple[float, str]] = {
    "FTMS": (20.0, "ppm"),
    "ITMS": (0.5, "Da"),
}


def analyzer_tolerance(analyzer_or_mode: str) -> tuple[float, str]:
    """Fragment tolerance for an analyzer (``FTMS``/``ITMS``) or a mode label
    (``FTMS_HCD`` → its leading analyzer). Defaults to 20 ppm."""
    analyzer = (analyzer_or_mode or "").split("_", 1)[0].upper()
    return ANALYZER_TOLERANCE.get(analyzer, (20.0, "ppm"))


def xic_targets_from_psms(psms: pa.Table) -> pa.Table:
    """Build an ``XIC_TARGET_TABLE`` from a PSM subset — one target per row,
    anchored on the measured precursor m/z at its assigned scan. ``target_id``
    is the row index into ``psms`` (so traces map straight back)."""
    n = psms.num_rows
    tbl = pa.table(
        {
            "target_id": pa.array(range(n), pa.int64()),
            "modified_sequence": psms.column("modified_sequence"),
            "precursor_charge": psms.column("charge"),
            "precursor_mz": psms.column("mz"),
            "rt_center": psms.column("retention_time_s"),
            "scan": psms.column("scan"),
        }
    )
    return cast_to_schema(tbl, XIC_TARGET_TABLE)


def extract_ms2_fragments(
    peaks: pa.Table,
    psms: pa.Table,
    *,
    acquisition_id: int = 0,
    ion_types: Sequence[IonType] = (IonType.B, IonType.Y),
    max_fragment_charge: int = 2,
) -> pa.Table:
    """Extract assigned-scan MS2 fragments for one bundle's PSMs.

    Splits PSMs by ``mass_analyzer`` so each gets its analyzer-appropriate
    tolerance, runs XIC level-2 ``assigned_scans_only`` on each subset, and
    concatenates the traces (``XIC_TRACE_TABLE``: target_id, scan, ion_type,
    position, fragment_charge, mz_theoretical, intensity, …). ``target_id``
    indexes ``psms``."""
    targets = xic_targets_from_psms(psms)
    analyzers = psms.column("mass_analyzer").to_pylist()
    parts: list[pa.Table] = []
    for analyzer in sorted({(a or "").upper() for a in analyzers}):
        pos = [i for i, a in enumerate(analyzers) if (a or "").upper() == analyzer]
        if not pos:
            continue
        tol, unit = analyzer_tolerance(analyzer)
        parts.append(
            extract_xic_scan_major(
                peaks,
                targets.take(pa.array(pos, pa.int64())),
                acquisition_id=acquisition_id,
                level=2,
                assigned_scans_only=True,
                ion_types=ion_types,
                max_fragment_charge=max_fragment_charge,
                tolerance=tol,
                tolerance_unit=unit,
            )
        )
    if not parts:
        return XIC_TRACE_TABLE.empty_table()
    return parts[0] if len(parts) == 1 else pa.concat_tables(parts)


def trace_to_scan_spectra(
    trace: pa.Table, *, min_intensity: float = 0.0
) -> dict[int, tuple[torch.Tensor, torch.Tensor]]:
    """Group an XIC level-2 trace into per-scan fragment spectra
    ``{scan: (mz_theoretical, intensity)}`` over matched fragments.

    The theoretical m/z is used as the alignment key (it is exact against the
    consensus basis), so ``build_consensus`` maps each fragment to its channel
    with no second tolerance. A pure transform — unit-testable without data."""
    scans = trace.column("scan").to_pylist()
    mz = trace.column("mz_theoretical").to_pylist()
    inten = trace.column("intensity").to_pylist()
    by_scan: dict[int, tuple[list, list]] = defaultdict(lambda: ([], []))
    for s, m, it in zip(scans, mz, inten):
        if it is not None and it > min_intensity and m is not None:
            by_scan[s][0].append(m)
            by_scan[s][1].append(it)
    return {
        s: (
            torch.tensor(ms, dtype=torch.float64),
            torch.tensor(its, dtype=torch.float64),
        )
        for s, (ms, its) in by_scan.items()
    }


__all__ = [
    "ANALYZER_TOLERANCE",
    "analyzer_tolerance",
    "xic_targets_from_psms",
    "extract_ms2_fragments",
    "trace_to_scan_spectra",
]
