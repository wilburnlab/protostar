"""Acquisition metadata curation.

Build the acquisition time table: parse acquisition datetime + instrument per
``.raw`` (Thermo header), order runs chronologically per instrument, and persist
as a ``constellation.massspec.acquisitions.Acquisitions`` table for carryover /
batch-effect analysis.
"""
