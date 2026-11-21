"""Tiny synthetic converted bundles, shaped like Constellation's ``massspec convert``.

Modelled on Constellation's own test helper (``tests/_massspec_fixtures.py``,
``write_bundle`` / ``write_bundle_manifest``), so protostar's stages can be tested
on the real v3 layout without Thermo DLLs. ``write_legacy_bundle`` writes the June
2026 pre-v3 shape: manifest ``schema_version`` 1 and a single ``peaks.parquet``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from constellation.massspec.bundle import BUNDLE_SCHEMA_VERSION, PROFILE_BUNDLE_KIND
from constellation.massspec.schemas import (
    ACQUISITION_METADATA_TABLE,
    MS_PEAK_TABLE,
    SCAN_METADATA_TABLE,
)

DEFAULT_RUN = {
    "instrument_model": "Orbitrap Fusion Lumos",
    "instrument_serial": "FSN20108",
    "creation_date": "06/03/2016 19:50:54",
}


def ms1_scan(scan: int, rt: float, mz: list[float], intensity: list[float], **meta: Any) -> dict:
    """One level-1 scan: peaks plus scan-metadata fields (``iit`` etc. via ``meta``)."""
    return {
        "scan": scan,
        "level": 1,
        "rt": rt,
        "mz": np.asarray(mz, dtype=float),
        "intensity": np.asarray(intensity, dtype=float),
        **meta,
    }


def _peaks_table(scans: list[dict]) -> pa.Table:
    rows: dict[str, list[Any]] = {f.name: [] for f in MS_PEAK_TABLE}
    for s in scans:
        n = int(s["mz"].size)
        for f in MS_PEAK_TABLE:
            if f.name in ("mz", "intensity"):
                rows[f.name].extend(s[f.name].tolist())
            else:
                rows[f.name].extend([s.get(f.name)] * n)
    return pa.table(
        {k: pa.array(v, MS_PEAK_TABLE.field(k).type) for k, v in rows.items()},
        schema=MS_PEAK_TABLE,
    )


def _scan_meta_table(scans: list[dict]) -> pa.Table:
    meta_rows = []
    for s in scans:
        row = {k: v for k, v in s.items() if k not in ("mz", "intensity")}
        row["peak_count"] = int(s["mz"].size)
        if s["mz"].size:
            row.setdefault("tic", float(s["intensity"].sum()))
            row.setdefault("intensity_floor", float(s["intensity"].min()))
        meta_rows.append(row)
    return pa.Table.from_pylist(meta_rows, schema=SCAN_METADATA_TABLE)


def write_run_bundle(
    parent: Path,
    stem: str,
    scans: list[dict],
    *,
    run: dict[str, Any] | None = None,
    sha256: str | None = None,
    with_profile: bool = True,
) -> Path:
    """A v3 centroid run bundle at ``parent/<stem>`` (+ ``<stem>.profile`` sibling)."""
    bundle = parent / stem
    (bundle / "peaks").mkdir(parents=True)
    pq.write_table(_peaks_table(scans), bundle / "peaks" / "part-00000.parquet")
    pq.write_table(_scan_meta_table(scans), bundle / "scan_metadata.parquet")
    run_meta = {**DEFAULT_RUN, **(run or {})}
    acquisition = {
        "source_format": "thermo_raw",
        "source_path": f"/data/{stem}.raw",
        "constellation_version": "0.0.0",
        "schema_version": 1,
        "n_scans": len(scans),
        **run_meta,
    }
    pq.write_table(
        pa.Table.from_pylist([acquisition], schema=ACQUISITION_METADATA_TABLE),
        bundle / "acquisition_metadata.parquet",
    )
    _write_manifest(bundle, kind="thermo_acquisition", stem=stem, run=run_meta, sha256=sha256)
    if with_profile:
        prof = parent / f"{stem}.profile"
        (prof / "profile").mkdir(parents=True)
        pq.write_table(_peaks_table(scans), prof / "profile" / "part-00000.parquet")
        _write_manifest(prof, kind=PROFILE_BUNDLE_KIND, stem=stem, run=run_meta, sha256=sha256)
    return bundle


def _write_manifest(
    bundle: Path, *, kind: str, stem: str, run: dict[str, Any], sha256: str | None
) -> None:
    profile = kind == PROFILE_BUNDLE_KIND
    (bundle / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": BUNDLE_SCHEMA_VERSION,
                "kind": kind,
                "source_file": f"/cluster/data/raw/{stem}.raw",
                "source_sha256": sha256 or f"sha-of-{stem}",
                "parameters": {"rt_bin_width_s": 60.0},
                "outputs": {
                    ("profile" if profile else "peaks"): ("profile" if profile else "peaks"),
                    "scan_metadata": "scan_metadata.parquet",
                    "acquisition_metadata": "acquisition_metadata.parquet",
                },
                "run_metadata": run,
            }
        )
    )


def write_legacy_bundle(parent: Path, stem: str) -> Path:
    """A June 2026 pre-v3 bundle: manifest ``schema_version`` 1 + ``peaks.parquet``."""
    bundle = parent / stem
    bundle.mkdir(parents=True)
    (bundle / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "thermo_acquisition",
                "source_file": f"/cluster/data/raw/{stem}.raw",
            }
        )
    )
    pq.write_table(_peaks_table([ms1_scan(1, 1.0, [500.0], [1.0])]), bundle / "peaks.parquet")
    return bundle
