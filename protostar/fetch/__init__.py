"""Manifest building and raw-file acquisition.

Canonical path: query the ProteomeXchange/PRIDE API for each dataset's ``.raw``
file list and published checksums, download missing files, verify against the
checksum, re-fetch on mismatch (resumable + repairable). An optional
``--seed-from <dir>`` hash-matches files already present on ESS and
hardlinks/moves them in instead of downloading.
"""
