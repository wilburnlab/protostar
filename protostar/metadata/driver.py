"""Build the acquisition time-table from converted ``proc/`` bundles (bundle v3).

Thin orchestration over Constellation's ``massspec.acquisitions``. protostar owns
only *which* bundles to read: every current v3 centroid run bundle under
``proc/<dataset>/``, via ``protostar.bundles``.

Constellation owns the rest:
- the record mapping (``acquisition_record_from_bundle``, which reads the
  manifest's source file, SHA-256 and run metadata);
- id assignment (``Acquisitions.register``);
- chronological ordering (``Acquisitions.with_acquisition_order``).

protostar adds one check: a run date that is present but not in Thermo's
``MM/DD/YYYY HH:MM:SS`` form fails loudly. Constellation would quietly null it,
and the format is confirmed uniform across all three datasets, so a mismatch is
a real signal.

Output is one table per dataset::

    <data_root>/proc/<dataset>/acquisitions.parquet

``acquisition_id`` is assigned in stem-sorted order, so it is stable across
re-runs and per dataset. ``acquisition_order`` is a 1-based chronological rank
within each ``instrument_serial``.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq
from constellation.massspec.acquisitions import Acquisitions, acquisition_record_from_bundle
from constellation.massspec.bundle import BUNDLE_SCHEMA_VERSION, read_bundle_manifest

from protostar import bundles

#: Parquet schema-metadata key recording which bundle format the table was built
#: from. The June 2026 tables (pre-v3 bundles) lack it, so they are rebuilt rather
#: than silently reused after reconversion.
BUNDLE_SCHEMA_KEY = b"protostar.bundle_schema_version"

#: Thermo ``creation_date`` wire format (``FileHeader.CreationDate`` .NET
#: ``ToString()``); confirmed uniform MM/DD/YYYY across all three datasets.
_THERMO_DATETIME_FMT = "%m/%d/%Y %H:%M:%S"


def acquisitions_out_path(data_root: str | Path, dataset: str) -> Path:
    """Where this dataset's acquisition table is written."""
    return bundles.proc_dir(data_root, dataset) / "acquisitions.parquet"


def enumerate_bundles(data_root: str | Path, dataset: str) -> list[Path]:
    """Current v3 centroid run bundles, sorted by stem (see ``protostar.bundles``)."""
    return list(bundles.iter_run_bundles(data_root, dataset))


def normalize_datetime(creation_date: str | None) -> str | None:
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


def read_acquisition_record(bundle_dir: Path) -> dict[str, object | None]:
    """One bundle's ``Acquisitions.register`` record.

    This is Constellation's ``acquisition_record_from_bundle``, with
    ``source_file`` reduced to its basename (the manifest may hold a cluster path)
    and the run date validated strictly.
    """
    rec = dict(acquisition_record_from_bundle(bundle_dir))
    rec["source_file"] = Path(str(rec["source_file"])).name
    run = read_bundle_manifest(bundle_dir).get("run_metadata") or {}
    rec["acquisition_datetime"] = normalize_datetime(run.get("creation_date"))
    return rec


def build_acquisitions(data_root: str | Path, dataset: str) -> Acquisitions:
    """Build the chronologically-ordered ``Acquisitions`` table for a dataset."""
    records = [read_acquisition_record(b) for b in enumerate_bundles(data_root, dataset)]
    acq, _ids = Acquisitions.empty().register(records)
    return acq.with_acquisition_order()


def write_acquisitions(acq: Acquisitions, out_path: str | Path) -> None:
    """Persist the table to ``out_path`` (PyArrow Parquet), stamped with the bundle
    schema version it was built from (``BUNDLE_SCHEMA_KEY``)."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    meta = dict(acq.table.schema.metadata or {})
    meta[BUNDLE_SCHEMA_KEY] = str(BUNDLE_SCHEMA_VERSION).encode()
    pq.write_table(acq.table.replace_schema_metadata(meta), out)


def acquisitions_present(data_root: str | Path, dataset: str) -> bool:
    """Whether a **current** acquisition table exists (resume/skip): one built from
    bundles of the current format. A June (pre-v3) table does not count."""
    path = acquisitions_out_path(data_root, dataset)
    if not path.is_file():
        return False
    meta = pq.read_schema(path).metadata or {}
    return meta.get(BUNDLE_SCHEMA_KEY) == str(BUNDLE_SCHEMA_VERSION).encode()


@dataclass
class MetadataSummary:
    n_bundles: int = 0
    n_acquisitions: int = 0
    n_instruments: int = 0
    n_datetime_null: int = 0


def summarize(acq: Acquisitions, n_bundles: int) -> MetadataSummary:
    """Counts for the stage log.

    ``n_acquisitions < n_bundles`` means ``register`` merged runs with the same
    SHA-256 (one file under two names). Check before trusting the table.
    """
    serials = acq.table.column("instrument_serial").to_pylist()
    dts = acq.table.column("acquisition_datetime").to_pylist()
    return MetadataSummary(
        n_bundles=n_bundles,
        n_acquisitions=len(acq),
        n_instruments=len({s for s in serials if s is not None}),
        n_datetime_null=sum(1 for d in dts if d is None),
    )


__all__ = [
    "BUNDLE_SCHEMA_KEY",
    "MetadataSummary",
    "acquisitions_out_path",
    "acquisitions_present",
    "build_acquisitions",
    "enumerate_bundles",
    "normalize_datetime",
    "read_acquisition_record",
    "summarize",
    "write_acquisitions",
]
