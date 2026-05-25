#!/usr/bin/env python
"""Stage 20 — acquisition time table.

Parses acquisition datetime + instrument per ``.raw`` (Thermo header), orders runs
chronologically per instrument, and persists a
``constellation.massspec.acquisitions.Acquisitions`` table — the substrate for carryover
and batch-effect analysis (acquisition-order reordering, instrument-specific effects).

STATUS: stub. Depends on Constellation Thermo `.raw` header access (ledger item #1).
"""

from __future__ import annotations


def main() -> None:
    raise NotImplementedError("Stage 20 (build_metadata) not yet implemented — scaffold only.")


if __name__ == "__main__":
    main()
