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


def trace_to_target_spectra(
    trace: pa.Table, *, min_intensity: float = 0.0
) -> dict[int, tuple[torch.Tensor, torch.Tensor]]:
    """Group an XIC level-2 trace into per-**target (PSM)** fragment spectra
    ``{target_id: (mz_theoretical, intensity)}`` over matched fragments.

    Keyed on ``target_id`` (not ``scan``) so each PSM's fragments stay isolated:
    a chimeric scan can carry rows for several targets, and a neighbour's
    fragments must not be pooled into this peptide's spectrum. The theoretical
    m/z is the alignment key (exact against the consensus basis), so
    ``build_consensus`` maps each fragment to its channel with no second
    tolerance. A pure transform — unit-testable without data."""
    tid = trace.column("target_id").to_pylist()
    mz = trace.column("mz_theoretical").to_pylist()
    inten = trace.column("intensity").to_pylist()
    by_target: dict[int, tuple[list, list]] = defaultdict(lambda: ([], []))
    for t, m, it in zip(tid, mz, inten):
        if it is not None and it > min_intensity and m is not None:
            by_target[int(t)][0].append(m)
            by_target[int(t)][1].append(it)
    return {
        t: (
            torch.tensor(ms, dtype=torch.float64),
            torch.tensor(its, dtype=torch.float64),
        )
        for t, (ms, its) in by_target.items()
    }


def trace_to_target_channels(
    trace: pa.Table, *, min_intensity: float = 0.0, isotope: int = 0
) -> dict[int, dict[tuple[int, int, int], list[float]]]:
    """Group an XIC level-2 trace by **target (PSM)** and fragment-channel
    identity ``(ion_type, position, fragment_charge)``, accumulating intensity
    and the intensity-weighted m/z error.

    Returns ``{target_id: {(ion_type, position, charge): [Σ intensity, Σ
    intensity·mz_error_ppm]}}``. Keying on ``target_id`` (not ``scan``) keeps
    each PSM isolated: a chimeric scan may carry rows for several targets, and a
    neighbour's fragment must not be projected onto this peptide's basis by
    matching (ion_type, position, charge) at a *different* theoretical m/z. The
    XIC extraction already matched observed peaks to theoretical fragments (with
    the analyzer-appropriate tolerance) and recorded ``mz_error_ppm`` per
    fragment, so the error is carried through — no second m/z match. Only the
    monoisotopic peak (``isotope``) is kept, matching the consensus basis. A
    pure transform."""
    tid = trace.column("target_id").to_pylist()
    it_type = trace.column("ion_type").to_pylist()
    pos = trace.column("position").to_pylist()
    fch = trace.column("fragment_charge").to_pylist()
    inten = trace.column("intensity").to_pylist()
    err = trace.column("mz_error_ppm").to_pylist()
    iso = trace.column("isotope").to_pylist()
    out: dict[int, dict[tuple[int, int, int], list[float]]] = defaultdict(
        lambda: defaultdict(lambda: [0.0, 0.0])
    )
    for i in range(trace.num_rows):
        v = inten[i]
        if v is None or v <= min_intensity or (iso[i] is not None and iso[i] != isotope):
            continue
        acc = out[int(tid[i])][(int(it_type[i]), int(pos[i]), int(fch[i]))]
        acc[0] += v
        e = err[i]
        if e is not None:
            acc[1] += v * e
    return {t: dict(ch) for t, ch in out.items()}


def channels_to_basis(
    channels: dict[tuple[int, int, int], list[float]], basis
) -> tuple[torch.Tensor, torch.Tensor]:
    """Project a target's channel accumulator (from :func:`trace_to_target_channels`)
    onto the fixed K-channel ``basis`` order.

    Returns ``(intensity[K], mz_error_ppm[K])`` — the intensity-weighted signed
    m/z error per channel (``NaN`` where unmatched). Mapping is by exact
    ``(ion_type, position, charge)`` identity, so it is unambiguous (no second
    m/z tolerance) and the intensity vector reproduces the theoretical-m/z
    alignment."""
    k = int(basis.K)
    vec = torch.zeros(k, dtype=torch.float64)
    err = torch.full((k,), float("nan"), dtype=torch.float64)
    if not channels:
        return vec, err
    index = {
        (int(a), int(b), int(c)): j
        for j, (a, b, c) in enumerate(
            zip(basis.ion_type.tolist(), basis.position.tolist(), basis.charge.tolist())
        )
    }
    for key, (sum_i, sum_ie) in channels.items():
        j = index.get(key)
        if j is not None and sum_i > 0:
            vec[j] = sum_i
            err[j] = sum_ie / sum_i
    return vec, err


__all__ = [
    "ANALYZER_TOLERANCE",
    "analyzer_tolerance",
    "xic_targets_from_psms",
    "extract_ms2_fragments",
    "trace_to_target_spectra",
    "trace_to_target_channels",
    "channels_to_basis",
]
