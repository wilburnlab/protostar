"""``.raw`` -> parquet-bundle conversion drivers (the ``proc/`` tree).

Thin orchestration over Constellation's Thermo ``.raw`` reader
(``massspec.readers.thermo.convert_batch``). Each ``.raw`` becomes a bundle
directory (``manifest.json`` + ``peaks.parquet`` + ``scan_metadata.parquet`` +
``acquisition_metadata.parquet``); per-scan filter strings are recorded so
fragmentation mode is recoverable without pre-splitting files. Bundles are
rebuilt from scratch for downstream consistency.
"""

from .driver import (
    ConvertSummary,
    bundle_out_parent,
    converted_stems,
    enumerate_raw,
    mode_name,
    plan_conversion,
    run_conversion,
    summarize_results,
)

__all__ = [
    "ConvertSummary",
    "bundle_out_parent",
    "converted_stems",
    "enumerate_raw",
    "mode_name",
    "plan_conversion",
    "run_conversion",
    "summarize_results",
]
