"""Acquisition metadata curation.

Build the acquisition time table: read acquisition datetime + instrument from
each converted ``proc/`` bundle's ``acquisition_metadata.parquet``, order runs
chronologically per instrument, and persist as a
``constellation.massspec.acquisitions.Acquisitions`` table for carryover /
batch-effect analysis.
"""

from .driver import (
    MetadataSummary,
    acquisitions_out_path,
    acquisitions_present,
    build_acquisitions,
    enumerate_bundles,
    normalize_datetime,
    proc_parent,
    read_acquisition_record,
    summarize,
    write_acquisitions,
)

__all__ = [
    "MetadataSummary",
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
