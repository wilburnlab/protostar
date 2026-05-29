#!/usr/bin/env python
"""Stage 20 — acquisition time table.

Reads acquisition datetime + instrument per run, orders runs chronologically per
instrument, and persists a ``constellation.massspec.acquisitions.Acquisitions`` table —
the substrate for carryover and batch-effect analysis (acquisition-order reordering,
instrument-specific effects). The needed fields already sit in each stage-10 convert
bundle's ``acquisition_metadata.parquet`` / ``manifest.json``, so this consumes the
``proc/`` tree rather than re-opening every ``.raw``.

STATUS: stub. Ledger item #1 (Thermo reader / bundle metadata) is landed; this is the next
data-stage to implement (reads stage-10 bundles).
"""

from __future__ import annotations


def main() -> None:
    raise NotImplementedError("Stage 20 (build_metadata) not yet implemented — scaffold only.")


if __name__ == "__main__":
    main()
