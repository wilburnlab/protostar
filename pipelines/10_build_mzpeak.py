#!/usr/bin/env python
"""Stage 10 — build mzpeak (.raw -> mzpeak Parquet + scanmeta).

Converts each ``.raw`` to an mzpeak Parquet cache plus a scanmeta sidecar
(IIT / TIC / filter_string per scan) via Constellation's Thermo reader. mzpeak caches
are **rebuilt from scratch** (no reuse of prior caches) for downstream consistency.
Per-scan filter strings are recorded so fragmentation mode is recoverable downstream
without pre-splitting files. Resume-safe (skip existing); designed to run as a SLURM array.

STATUS: stub. Depends on Constellation Thermo `.raw` reader + mzpeak writer
(ledger items #1, #2 in docs/constellation_contributions.md).
"""

from __future__ import annotations


def main() -> None:
    raise NotImplementedError("Stage 10 (build_mzpeak) not yet implemented — scaffold only.")


if __name__ == "__main__":
    main()
