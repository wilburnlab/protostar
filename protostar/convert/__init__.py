"""``.raw`` -> mzpeak Parquet + scanmeta conversion drivers.

Thin orchestration over Constellation's Thermo ``.raw`` reader. mzpeak caches are
rebuilt from scratch for downstream consistency; per-scan filter strings are
recorded so fragmentation mode is recoverable without pre-splitting files.
"""
