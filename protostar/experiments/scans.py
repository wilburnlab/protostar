"""PSM-anchored scan selection for the MS2 experiments.

Pure data orchestration over a MaxQuant ``Search``'s ``PSM_TABLE`` — load,
gate, calibrant-tag, select the recurring (high-replicate) peptides, and group
PSMs into replicate sets per (peptide, charge, fragmentation mode). No model
code lives here: all scoring / likelihoods come from ``constellation``.

The grouping key's ``mode`` is the per-PSM ``mass_analyzer`` + ``fragmentation``
(the column filter the project mandates, not a filter-string re-parse) — a
coarse mode (FTMS_HCD, ITMS_CID, …); collision-energy refinement is a later
enrichment from ``scan_metadata``.
"""

from __future__ import annotations

import zipfile
from collections import defaultdict
from collections.abc import Iterator, Sequence
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
from constellation.massspec.io.maxquant import read_maxquant_search

from protostar.peptides import procal_sequences

#: Default replicate-grouping key.
DEFAULT_GROUP_KEY: tuple[str, ...] = ("modified_sequence", "charge", "mode")


def add_mode_column(psms: pa.Table) -> pa.Table:
    """Append a coarse ``mode`` column = ``MASS_ANALYZER_FRAGMENTATION`` (e.g.
    ``FTMS_HCD``), upper-cased, from the per-PSM analyzer + fragmentation."""
    analyzer = pc.utf8_upper(pc.fill_null(psms.column("mass_analyzer"), ""))
    frag = pc.utf8_upper(pc.fill_null(psms.column("fragmentation"), ""))
    mode = pc.binary_join_element_wise(analyzer, frag, "_")
    return psms.append_column("mode", mode)


def load_psms(search_path: str | Path, *, extract_dir: str | Path | None = None) -> pa.Table:
    """Load a MaxQuant search's PSMs (with a ``mode`` column appended).

    ``search_path`` is a directory containing ``msms.txt`` (any of the layouts
    ``read_maxquant_search`` accepts) or a ``.zip`` of one. A zip is extracted to
    ``extract_dir`` (defaults to a sibling ``<stem>/`` next to the zip)."""
    p = Path(search_path)
    if p.is_file() and p.suffix == ".zip":
        dest = Path(extract_dir) if extract_dir is not None else p.with_suffix("")
        dest.mkdir(parents=True, exist_ok=True)
        if not any(dest.rglob("msms.txt")):
            with zipfile.ZipFile(p) as z:
                z.extractall(dest)
        p = dest
    search = read_maxquant_search(p)
    return add_mode_column(search.psms)


def gate_psms(
    psms: pa.Table,
    *,
    max_pep: float = 0.01,
    min_score: float = 0.0,
    drop_decoy: bool = True,
    drop_contaminant: bool = True,
    psm_types: Sequence[str] | None = ("MULTI-MSMS",),
) -> pa.Table:
    """Confidence gate — the cheap, principled guard against trusting every
    assignment. Keeps non-decoy, non-contaminant PSMs of an allowed ``psm_type``
    with ``pep <= max_pep`` and ``score >= min_score``. Null ``pep``/``score``
    rows fail the numeric cuts (dropped). Statistical outlier screening beyond
    this is the deferred consensus-refinement pass, not here."""
    conds: list = []
    if drop_decoy:
        conds.append(pc.invert(psms.column("is_decoy")))
    if drop_contaminant:
        conds.append(pc.invert(psms.column("is_contaminant")))
    if psm_types is not None:
        conds.append(pc.is_in(psms.column("psm_type"), value_set=pa.array(list(psm_types))))
    if max_pep is not None:
        conds.append(pc.less_equal(psms.column("pep"), max_pep))
    if min_score is not None:
        conds.append(pc.greater_equal(psms.column("score"), min_score))
    if not conds:
        return psms
    mask = conds[0]
    for c in conds[1:]:
        mask = pc.and_(mask, c)
    return psms.filter(pc.fill_null(mask, False))


def tag_calibrants(psms: pa.Table, *, dataset: str | None = None) -> pa.Table:
    """Append ``is_procal`` (and, when ``dataset`` given, ``is_qc``) by bare-AA
    sequence membership — annotations, not a filter (the analysis set stays
    data-driven via :func:`select_recurring_peptides`)."""
    seq = psms.column("sequence")
    procal = pa.array(sorted(procal_sequences()), type=pa.string())
    out = psms.append_column("is_procal", pc.is_in(seq, value_set=procal))
    if dataset is not None:
        from protostar.peptides import load_reference

        ref = load_reference(dataset)
        qc_seqs = ref.filter(ref.column("is_qc")).column("sequence")
        out = out.append_column("is_qc", pc.is_in(seq, value_set=qc_seqs.combine_chunks()))
    return out


def select_recurring_peptides(
    psms: pa.Table,
    *,
    key: Sequence[str] = ("modified_sequence", "charge"),
    min_acquisitions: int | None = None,
    top_n: int | None = None,
) -> pa.Table:
    """Filter to the recurring (high-replicate) peptides — those identified
    across many acquisitions, the high-N target set.

    Counts distinct ``raw_file`` per ``key``; keeps peptides with at least
    ``min_acquisitions`` (or the ``top_n`` most-recurring). Returns the PSM
    subset for the surviving peptides."""
    if min_acquisitions is None and top_n is None:
        raise ValueError("pass min_acquisitions or top_n")
    counts = psms.group_by(list(key)).aggregate([("raw_file", "count_distinct")])
    n = counts.column("raw_file_count_distinct")
    if top_n is not None:
        order = pc.sort_indices(counts, sort_keys=[("raw_file_count_distinct", "descending")])
        counts = counts.take(order.slice(0, top_n))
    else:
        counts = counts.filter(pc.greater_equal(n, min_acquisitions))
    keep = counts.select(list(key))
    # inner-join the PSMs to the surviving keys.
    return psms.join(keep, keys=list(key), join_type="inner")


def replicate_groups(
    psms: pa.Table, *, key: Sequence[str] = DEFAULT_GROUP_KEY
) -> Iterator[tuple[dict, pa.Table]]:
    """Yield ``(key_dict, psm_subset)`` for each replicate group — all PSMs
    sharing the grouping key (default: peptide, charge, mode). Groups span
    acquisitions (the replicate dimension)."""
    cols = [psms.column(c).to_pylist() for c in key]
    buckets: dict[tuple, list[int]] = defaultdict(list)
    for i, ktuple in enumerate(zip(*cols)):
        buckets[ktuple].append(i)
    for ktuple, rows in buckets.items():
        yield dict(zip(key, ktuple)), psms.take(pa.array(rows, pa.int64()))


__all__ = [
    "DEFAULT_GROUP_KEY",
    "add_mode_column",
    "load_psms",
    "gate_psms",
    "tag_calibrants",
    "select_recurring_peptides",
    "replicate_groups",
]
