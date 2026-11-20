"""``.raw`` -> parquet-bundle conversion drivers (the ``proc/`` tree).

Thin orchestration over Constellation's Thermo ``.raw`` reader
(``massspec.readers.thermo.convert_batch``). Each ``.raw`` becomes, in one pass,
a v3 centroid bundle (``proc/<dataset>/<stem>/``: ``manifest.json``,
``peaks/part-*.parquet``, ``scan_metadata.parquet``,
``acquisition_metadata.parquet``) and a sibling profile bundle
(``proc/<dataset>/<stem>.profile/``). Per-scan filter strings are recorded, so
fragmentation mode is recoverable without pre-splitting files. Read bundles
through ``protostar.bundles``.
"""

from .driver import (
    ConversionPlan,
    ConvertSummary,
    bundle_out_parent,
    classify,
    converted_stems,
    enumerate_raw,
    plan_conversion,
    run_conversion,
    summarize_results,
)

__all__ = [
    "ConversionPlan",
    "ConvertSummary",
    "bundle_out_parent",
    "classify",
    "converted_stems",
    "enumerate_raw",
    "plan_conversion",
    "run_conversion",
    "summarize_results",
]
