"""Build the acquisition time-table from converted ``proc/`` bundles.

Thin orchestration over Constellation's ``massspec.acquisitions``: protostar
owns only *which* bundles to read and how to map their per-acquisition metadata
into ``ACQUISITION_TABLE`` records; the table container, schema, and
chronological-ordering logic (``Acquisitions.with_acquisition_order``) live in
Constellation. We never reimplement the ordering.

Input is the stage-10 ``proc/`` tree, **not** the ``.raw`` files: each bundle's
``acquisition_metadata.parquet`` already carries the acquisition datetime +
instrument identity. We read the ``centroid`` bundles (every acquisition has
one; ``profile`` is a subset with identical acquisition metadata).

Output is one table per dataset::

    <data_root>/proc/<dataset>/acquisitions.parquet

with ``acquisition_order`` a 1-based chronological rank within each
``instrument_serial`` (per dataset).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq
from constellation.massspec.acquisitions import Acquisitions

#: Thermo ``creation_date`` wire format (``FileHeader.CreationDate`` .NET
#: ``ToString()``); confirmed uniform MM/DD/YYYY across all three datasets.
_THERMO_DATETIME_FMT = "%m/%d/%Y %H:%M:%S"
#: ``source_kind`` recorded for Thermo ``.raw``-derived acquisitions.
SOURCE_KIND = "thermo_raw"


def proc_parent(data_root: "str | Path", dataset: str, *, mode: str = "centroid") -> Path:
    """The ``proc/<dataset>/<mode>`` dir holding one ``<stem>/`` bundle per ``.raw``."""
    return Path(data_root) / "proc" / dataset / mode


def acquisitions_out_path(data_root: "str | Path", dataset: str) -> Path:
    """Where this dataset's acquisition table is written."""
    return Path(data_root) / "proc" / dataset / "acquisitions.parquet"


def enumerate_bundles(
    data_root: "str | Path", dataset: str, *, mode: str = "centroid"
) -> list[Path]:
    """Sorted bundle dirs (those with a ``manifest.json``) under ``proc/<dataset>/<mode>``.

    Sorted by stem so the derived ``acquisition_id`` (index in this list) is
    stable across runs. Mirrors ``convert.driver.converted_stems`` enumeration.
    """
    parent = proc_parent(data_root, dataset, mode=mode)
    if not parent.is_dir():
        return []
    return sorted((m.parent for m in parent.glob("*/manifest.json")), key=lambda p: p.name)


def normalize_datetime(creation_date: "str | None") -> str | None:
    """Normalize a Thermo ``creation_date`` to ISO-8601, or ``None``.

    ``None`` passes through. A non-conforming string raises ``ValueError``
    (the format is confirmed uniform, so a mismatch is a real signal — fail
    loud rather than silently null it).
    """
    if creation_date is None:
        return None
    try:
        return datetime.strptime(creation_date, _THERMO_DATETIME_FMT).isoformat()
    except ValueError as exc:
        raise ValueError(
            f"creation_date {creation_date!r} does not match expected "
            f"format {_THERMO_DATETIME_FMT!r}"
        ) from exc


def _first(table: object, column: str) -> object | None:
    """First value of ``column`` in a single-row Arrow table, or ``None``."""
    if column not in table.column_names:  # type: ignore[attr-defined]
        return None
    col = table.column(column)  # type: ignore[attr-defined]
    return col[0].as_py() if col.length() else None


def read_acquisition_record(bundle_dir: Path, acquisition_id: int) -> dict[str, object | None]:
    """Map one bundle's ``acquisition_metadata.parquet`` to an ``ACQUISITION_TABLE`` record.

    ``acquisition_order`` is intentionally left out — it is filled canonically
    by ``Acquisitions.with_acquisition_order``.
    """
    meta = pq.read_table(bundle_dir / "acquisition_metadata.parquet")
    # source_file: the .raw stem, taken from the manifest's source_file basename
    # so it matches the original acquisition name (the bundle dir name is the
    # same stem, but the manifest is authoritative).
    manifest = json.loads((bundle_dir / "manifest.json").read_text())
    source_file = Path(str(manifest.get("source_file", bundle_dir.name))).name
    return {
        "acquisition_id": acquisition_id,
        "source_file": source_file,
        "source_kind": SOURCE_KIND,
        "acquisition_datetime": normalize_datetime(_first(meta, "creation_date")),
        "instrument_serial": _first(meta, "instrument_serial"),
        "instrument_model": _first(meta, "instrument_model"),
    }


def build_acquisitions(
    data_root: "str | Path", dataset: str, *, mode: str = "centroid"
) -> Acquisitions:
    """Build the chronologically-ordered ``Acquisitions`` table for a dataset.

    ``acquisition_id`` = index in the stem-sorted bundle list (stable across
    re-runs; per-dataset, not global). ``acquisition_order`` is a 1-based rank
    within each ``instrument_serial``.
    """
    bundles = enumerate_bundles(data_root, dataset, mode=mode)
    records = [read_acquisition_record(b, i) for i, b in enumerate(bundles)]
    return Acquisitions.from_records(records).with_acquisition_order()


def write_acquisitions(acq: Acquisitions, out_path: "str | Path") -> None:
    """Persist the table to ``out_path`` (PyArrow Parquet)."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(acq.table, out)


def acquisitions_present(data_root: "str | Path", dataset: str) -> bool:
    """Whether this dataset's acquisition table already exists (resume/skip)."""
    return acquisitions_out_path(data_root, dataset).is_file()


@dataclass
class MetadataSummary:
    n_bundles: int = 0
    n_acquisitions: int = 0
    n_instruments: int = 0
    n_datetime_null: int = 0


def summarize(acq: Acquisitions, n_bundles: int) -> MetadataSummary:
    serials = acq.table.column("instrument_serial").to_pylist()
    dts = acq.table.column("acquisition_datetime").to_pylist()
    return MetadataSummary(
        n_bundles=n_bundles,
        n_acquisitions=len(acq),
        n_instruments=len({s for s in serials if s is not None}),
        n_datetime_null=sum(1 for d in dts if d is None),
    )


__all__ = [
    "MetadataSummary",
    "SOURCE_KIND",
    "acquisitions_out_path",
    "acquisitions_present",
    "build_acquisitions",
    "enumerate_bundles",
    "normalize_datetime",
    "proc_parent",
    "read_acquisition_record",
    "summarize",
    "write_acquisitions",
]
