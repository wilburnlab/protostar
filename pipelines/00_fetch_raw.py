#!/usr/bin/env python
"""Stage 00 — fetch raw data (canonical: fresh download + hash verify).

Builds an expected-file manifest from ``config/datasets.toml`` by querying the
ProteomeXchange/PRIDE API per accession (PXD004732 / PXD010595 / PXD021013) for each
dataset's ``.raw`` file list and published checksums, then downloads missing files,
verifies each against its checksum, and re-fetches on mismatch (resumable + repairable).

This script reflects *what a fresh user does*. The optional ``--seed-from <dir>`` flag
hash-matches files already on ESS (e.g. the existing local copy) and
hardlinks/moves them in instead of downloading — a local time-saver, not part of the
reproducible path.

Usage (planned)::

    python pipelines/00_fetch_raw.py --dataset Zolg2017 --dry-run
    python pipelines/00_fetch_raw.py --dataset Zolg2017
    python pipelines/00_fetch_raw.py --dataset Zolg2017 \
        --seed-from /path/to/ProteomeTools

STATUS: stub. Implementation tracked in CLAUDE.md "Status" / docs ledger.
"""

from __future__ import annotations


def main() -> None:
    raise NotImplementedError("Stage 00 (fetch_raw) not yet implemented — scaffold only.")


if __name__ == "__main__":
    main()
