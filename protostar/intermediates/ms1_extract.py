"""MS1 isotope-envelope chromatogram extraction for the Counter-motivation figures.

Wraps Constellation's XIC level-1 extraction (M+0/M+1/M+2 isotope traces built
from ``peptide_envelope``) into per-(peptide, charge, scan) **wide** rows carrying
per-isotope intensity + signed m/z error, plus per-scan IIT/TIC — the
"score-ready" intermediate the figure scripts read. RT-windowed and multi-charge
(``charge_range=(1, 4)``) by default, so a single extraction feeds all four
figures: isotope shot noise, charge-dependent α(z), m/z-vs-N + Student-t, and the
interference chromatogram.

No model code lives here (all of that is Constellation). The wide schema keeps
everything a later Counter-score pass needs — per-isotope intensity, m/z error,
IIT, RT, scan — so scores can be joined on (modified_sequence, charge, raw_file,
scan) without re-extraction. A neighbour PSM cannot leak onto this peptide's
envelope: the M+0/M+1/M+2 m/z come from this peptide's own ProForma + charge.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator, Sequence
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
from constellation.core.io.schemas import cast_to_schema
from constellation.massspec.quant.chromatogram import extract_xic_scan_major
from constellation.massspec.quant.schemas import XIC_TARGET_TABLE

#: Number of isotopes (M+0, M+1, M+2) extracted per precursor.
N_ISOTOPES: int = 3

#: Per-scan SCAN_METADATA_TABLE fields carried onto the wide rows (level-1 survey
#: metadata). Only those present in a bundle are joined, so older readers degrade
#: gracefully.
#: - ``iit`` (ms) puts intensity on a true ion-count axis (N ∝ I·iit).
#: - ``tic`` gives the assigned fraction; ``tic × iit`` is the trap-load proxy for the
#:   space-charge test (decision log D10).
#: - ``intensity_floor`` (SCAN_METADATA v5) is the scan's detection floor (D8).
#: - ``max_ion_time_ms`` / ``agc_*`` flag scans whose injection hit the time limit.
#: - ``space_charge_comp_ppm`` / ``lock_mass_correction_ppm`` are Thermo's own per-scan
#:   m/z corrections.
#: - ``ft_resolution`` is the scan's resolution setting.
#: The rest are saturation and complexity context.
_SCAN_META_FIELDS: tuple[str, ...] = (
    "iit",
    "tic",
    "base_peak_intensity",
    "peak_count",
    "analyzer",
    "intensity_floor",
    "max_ion_time_ms",
    "agc_target",
    "agc_fill",
    "ft_resolution",
    "space_charge_comp_ppm",
    "lock_mass_correction_ppm",
)

#: An injection within this fraction of ``max_ion_time_ms`` counts as having hit the
#: limit (the AGC target was not reached).
_IIT_AT_MAX_FRACTION = 0.999


def ms1_targets_from_psms(psms: pa.Table) -> pa.Table:
    """Build an MS1 ``XIC_TARGET_TABLE`` from a PSM subset — one target per row,
    anchored on the PSM's retention time. ``target_id`` is the row index into
    ``psms`` so wide rows map straight back to the PSM (sequence, raw_file, …).

    MS1 queries derive each isotope's m/z from ``modified_sequence`` + the query
    charge, so only the modseq, an RT anchor, and (for provenance) the PSM charge
    are needed. The RT window is applied at extraction (``extract_ms1_isotopes``),
    not here, so the same target sweeps the full elution."""
    n = psms.num_rows
    tbl = pa.table(
        {
            "target_id": pa.array(range(n), pa.int64()),
            "modified_sequence": psms.column("modified_sequence"),
            "precursor_charge": pc.cast(psms.column("charge"), pa.int8()),
            "rt_center": pc.cast(psms.column("retention_time_s"), pa.float64()),
            "scan": pc.cast(psms.column("scan"), pa.int32()),
        }
    )
    return cast_to_schema(tbl, XIC_TARGET_TABLE)


def extract_ms1_isotopes(
    peaks: pa.Table | Path | str,
    targets: pa.Table,
    *,
    acquisition_id: int = 0,
    n_isotopes: int = N_ISOTOPES,
    charge_range: tuple[int, int] | None = (1, 4),
    rt_window: float | None = 60.0,
    all_in_window: bool = False,
    tolerance: float = 20.0,
) -> pa.Table:
    """Extract MS1 isotope traces (M+0 … M+{n−1}) for the targets over an RT window.

    Thin wrapper over ``extract_xic_scan_major(level=1)``. ``peaks`` is a v3 bundle
    directory, which is streamed, or an in-memory peak table (tests). ``charge_range``
    extracts every charge in ``[lo, hi]`` for each peptide (so α(z) sees charges
    the MS2 trigger missed); ``rt_window`` gates each target to
    ``rt_center ± window`` (full elution, not just the assigned scan);
    ``all_in_window=True`` returns every in-tolerance peak per isotope (chimeric /
    interference work — note the wide pivot assumes one peak per cell, so keep it
    ``False`` for the rectangular intermediate). ``drop_unmatched=False`` keeps
    zero-intensity isotopes so every scan has a rectangular M+0/M+1/M+2 row.

    Returns the long ``XIC_TRACE_TABLE`` (one row per target·scan·charge·isotope)."""
    return extract_xic_scan_major(
        peaks,
        targets,
        acquisition_id=acquisition_id,
        level=1,
        n_isotopes=n_isotopes,
        charge_range=charge_range,
        rt_window=rt_window,
        tolerance=tolerance,
        tolerance_unit="ppm",
        all_in_window=all_in_window,
        drop_unmatched=False,
    )


def trace_to_wide(trace: pa.Table, *, n_isotopes: int = N_ISOTOPES) -> pa.Table:
    """Pivot the long MS1 trace into wide per-(target, charge, scan) rows.

    Columns: ``target_id, precursor_charge, scan, rt, mz_mono, sum_I,
    I_0..I_{n−1}, mz_ppm_0..mz_ppm_{n−1}``. The isotope axis is pivoted into the
    per-isotope intensity / signed-m/z-error columns; unmatched isotopes (kept by
    ``drop_unmatched=False``) contribute intensity 0 and NaN m/z error. A pure
    transform — unit-testable on a synthetic ``XIC_TRACE_TABLE`` with no data."""
    tid = trace.column("target_id").to_pylist()
    chg = trace.column("precursor_charge").to_pylist()
    scn = trace.column("scan").to_pylist()
    rt = trace.column("rt").to_pylist()
    iso = trace.column("isotope").to_pylist()
    inten = trace.column("intensity").to_pylist()
    err = trace.column("mz_error_ppm").to_pylist()
    mzt = trace.column("mz_theoretical").to_pylist()

    cells: dict[tuple[int, int, int], dict] = {}
    for i in range(trace.num_rows):
        k = iso[i]
        if k is None or not (0 <= k < n_isotopes):
            continue
        key = (int(tid[i]), int(chg[i]), int(scn[i]))
        cell = cells.get(key)
        if cell is None:
            cell = {
                "rt": rt[i],
                "I": [0.0] * n_isotopes,
                "err": [float("nan")] * n_isotopes,
                "mz_mono": None,
            }
            cells[key] = cell
        v = inten[i]
        cell["I"][k] = float(v) if v is not None else 0.0
        e = err[i]
        cell["err"][k] = float(e) if e is not None else float("nan")
        if k == 0:
            cell["mz_mono"] = mzt[i]

    keys = sorted(cells)
    out: dict[str, pa.Array] = {
        "target_id": pa.array([k[0] for k in keys], pa.int64()),
        "precursor_charge": pa.array([k[1] for k in keys], pa.int8()),
        "scan": pa.array([k[2] for k in keys], pa.int32()),
        "rt": pa.array([cells[k]["rt"] for k in keys], pa.float64()),
        "mz_mono": pa.array([cells[k]["mz_mono"] for k in keys], pa.float64()),
        "sum_I": pa.array([float(sum(cells[k]["I"])) for k in keys], pa.float64()),
    }
    for j in range(n_isotopes):
        out[f"I_{j}"] = pa.array([cells[k]["I"][j] for k in keys], pa.float64())
    for j in range(n_isotopes):
        out[f"mz_ppm_{j}"] = pa.array([cells[k]["err"][j] for k in keys], pa.float64())
    return pa.table(out)


def attach_scan_meta(wide: pa.Table, scan_meta: pa.Table) -> pa.Table:
    """Left-join per-scan metadata (level-1 survey) onto the wide rows by ``scan``.

    Adds two derived columns:
    - ``iit_s`` = ``iit`` ms / 1000, the variance-model units;
    - ``iit_at_max`` = the injection hit ``max_ion_time_ms``, when that field exists.

    Only the ``_SCAN_META_FIELDS`` present in ``scan_meta`` are joined, so it
    degrades gracefully across reader versions. A pure transform."""
    present = [f for f in _SCAN_META_FIELDS if f in scan_meta.column_names]
    # scan_metadata is one row per scan by construction, so the left join stays
    # one-to-one.
    meta = scan_meta.select(["scan", *present])
    idx = meta.schema.get_field_index("scan")
    meta = meta.set_column(idx, "scan", pc.cast(meta.column("scan"), pa.int32()))
    joined = wide.join(meta, keys=["scan"], join_type="left outer")
    if "iit" in joined.column_names:
        iit = pc.cast(joined.column("iit"), pa.float64())
        joined = joined.append_column("iit_s", pc.divide(iit, 1000.0))
        if "max_ion_time_ms" in joined.column_names:
            limit = pc.multiply(
                pc.cast(joined.column("max_ion_time_ms"), pa.float64()), _IIT_AT_MAX_FRACTION
            )
            joined = joined.append_column("iit_at_max", pc.greater_equal(iit, limit))
    return joined


def load_wide(
    root,
    *,
    dataset: str | None = None,
    peptides: Sequence[str] | None = None,
    columns: Sequence[str] | None = None,
) -> pa.Table:
    """Load the wide MS1 intermediate (a directory of ``shard_*.parquet``), with
    optional pushdown filters on ``dataset`` and bare-AA ``sequence``. The
    consumer-side companion to the stage-30 driver, used by the figure scripts."""
    import pyarrow.dataset as pads

    dset = pads.dataset(str(root), format="parquet")
    filt = None
    if dataset is not None:
        filt = pc.field("dataset") == dataset
    if peptides is not None:
        pf = pc.is_in(pc.field("sequence"), value_set=pa.array(list(peptides), pa.string()))
        filt = pf if filt is None else (filt & pf)
    return dset.to_table(columns=list(columns) if columns else None, filter=filt)


def dominant_charge(table: pa.Table) -> pa.Table:
    """Restrict the wide table to each peptide's identified (dominant) charge —
    i.e. ``precursor_charge == psm_charge``. This drops the phantom
    ``charge_range=(1,4)`` extractions where the peptide does not genuinely
    ionize at that charge (random in-tolerance peaks that otherwise contaminate
    the isotope-ratio and m/z-error statistics). No-op if the columns are absent."""
    if "psm_charge" in table.column_names and "precursor_charge" in table.column_names:
        return table.filter(pc.equal(table.column("precursor_charge"), table.column("psm_charge")))
    return table


def to_arrays(table: pa.Table) -> dict:
    """Materialize the per-isotope wide columns into NumPy arrays the estimators
    consume: ``{I (n×K), mz_ppm (n×K), sum_I, rt, iit_s}`` (NaN where absent)."""
    import numpy as np

    n = table.num_rows
    names = set(table.column_names)

    def col(name: str) -> np.ndarray:
        if name not in names:
            return np.full(n, np.nan)
        return np.asarray(table.column(name).to_numpy(zero_copy_only=False), dtype=float)

    def stack(prefix: str) -> np.ndarray:
        return np.stack([col(f"{prefix}{k}") for k in range(N_ISOTOPES)], axis=1)

    return {
        "I": stack("I_"),
        "mz_ppm": stack("mz_ppm_"),
        "sum_I": col("sum_I"),
        "rt": col("rt"),
        "iit_s": col("iit_s"),
    }


def iter_groups(
    table: pa.Table, keys: Sequence[str] = ("modified_sequence", "precursor_charge")
) -> Iterator[tuple[dict, pa.Table]]:
    """Yield ``(key_dict, sub_table)`` per group — the per-(peptide, charge) units
    the figures summarize. Pure (no I/O)."""
    cols = [table.column(k).to_pylist() for k in keys]
    buckets: dict[tuple, list[int]] = defaultdict(list)
    for i, ktuple in enumerate(zip(*cols)):
        buckets[ktuple].append(i)
    for ktuple, rows in buckets.items():
        yield dict(zip(keys, ktuple)), table.take(pa.array(rows, pa.int64()))


__all__ = [
    "N_ISOTOPES",
    "ms1_targets_from_psms",
    "extract_ms1_isotopes",
    "trace_to_wide",
    "attach_scan_meta",
    "load_wide",
    "iter_groups",
    "to_arrays",
    "dominant_charge",
]
