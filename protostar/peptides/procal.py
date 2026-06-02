"""PROCAL calibrant reference — the 40 synthetic standard peptides.

PROCAL (Zolg et al., *Proteomics* 2017; DOI 10.1002/pmic.201700263) is a set of
40 synthetic peptides spiked into **every** ProteomeTools acquisition — the
highest-replicate, widest-dynamic-range target set in the corpus, and so the
high-N proving ground for the MS2 (and later MS1) statistics. PROCAL is a subset
of the JPT calibrants already flagged ``is_qc`` in the per-dataset reference
(``reference.py``); this module makes the exact 40 available as a first-class
annotation (the broader "present in all samples" calibrant set is selected
data-driven from the searches — PROCAL is one labelled subset of it).

The source is supplementary Table S1 of the PROCAL paper (sheet
``PROCAL Sequences``). **openpyxl cannot open that specific workbook** (it ships a
malformed defined-name spec that makes openpyxl drop every sheet), so
``parse_procal_supplement`` reads the sheet directly from the OOXML zip — a small,
dependency-free reader sufficient for this one flat table. As with the other
supplements the workbook is a one-time manual input; the committed
``data/procal.parquet`` is the durable, always-available asset.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.sax.saxutils import unescape

import pyarrow as pa
import pyarrow.parquet as pq

#: Schema version — bump on any additive column change.
PROCAL_SCHEMA_VERSION: int = 1

#: One row per PROCAL peptide (40 total).
PROCAL_TABLE: pa.Schema = pa.schema(
    [
        pa.field("peptide_number", pa.int32(), nullable=False),  # 1..40, as published
        pa.field("sequence", pa.string(), nullable=False),  # bare AA
        pa.field("mz_2plus", pa.float64(), nullable=True),  # (M+2H)2+
        pa.field("mz_3plus", pa.float64(), nullable=True),  # (M+3H)3+
        pa.field("expected_charge", pa.string(), nullable=True),  # e.g. "2;3" or "2"
        pa.field("ssrcalc_hi", pa.float64(), nullable=True),  # SSRCalc HI (pred., 2015)
        pa.field("rt_60min_gradient", pa.string(), nullable=True),  # "9.96+/-0.17" (raw)
    ],
    metadata={
        b"schema_name": b"ProcalTable",
        b"schema_version": str(PROCAL_SCHEMA_VERSION).encode("utf-8"),
    },
)

#: Default supplement filename (the manual Excel input).
PROCAL_SUPPLEMENT_FILE: str = "pmic12732-sup-0002-table-s1.xlsx"
_PROCAL_SHEET_NAME: str = "PROCAL Sequences"

#: Committed parquet lives beside this module so it ships with the repo.
_DATA_DIR = Path(__file__).resolve().parent / "data"


# ── minimal OOXML reader (openpyxl chokes on this workbook) ─────────────


def _col_letter(ref: str) -> str:
    return re.match(r"([A-Z]+)", ref).group(1)


def _read_sheet_grid(xlsx_path: Path, sheet_name: str) -> dict[tuple[int, str], str]:
    """Return ``{(row, col_letter): value}`` for one sheet, resolving shared
    strings. A small, dependency-free reader for flat supplement tables."""
    with zipfile.ZipFile(xlsx_path) as z:
        wb = z.read("xl/workbook.xml").decode("utf-8", "replace")
        rels = z.read("xl/_rels/workbook.xml.rels").decode("utf-8", "replace")
        # Resolve sheet name → r:id → target, order-independent in attributes.
        rid = None
        for sm in re.finditer(r"<sheet\b([^>]*)/?>", wb):
            attrs = sm.group(1)
            nm = re.search(r'name="([^"]+)"', attrs)
            idm = re.search(r'r:id="([^"]+)"', attrs)
            if nm and idm and nm.group(1) == sheet_name:
                rid = idm.group(1)
                break
        if rid is None:
            raise ValueError(f"sheet {sheet_name!r} not found in {xlsx_path}")
        rel_target = {
            m.group("id"): m.group("t")
            for m in re.finditer(
                r'<Relationship\b(?=[^>]*Id="(?P<id>[^"]+)")(?=[^>]*Target="(?P<t>[^"]+)")[^>]*/?>',
                rels,
            )
        }
        target = rel_target.get(rid)
        if target is None:
            raise ValueError(f"no rels target for {rid!r} in {xlsx_path}")
        target = target.lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        sheet_xml = z.read(target).decode("utf-8", "replace")
        try:
            ss_xml = z.read("xl/sharedStrings.xml").decode("utf-8", "replace")
            shared = re.findall(r"<t[^>]*>(.*?)</t>", ss_xml, re.S)
        except KeyError:
            shared = []

    grid: dict[tuple[int, str], str] = {}
    for rowm in re.finditer(r'<row[^>]*r="(\d+)"[^>]*>(.*?)</row>', sheet_xml, re.S):
        rn = int(rowm.group(1))
        for cm in re.finditer(r'<c r="([A-Z]+\d+)"([^>]*)>(.*?)</c>', rowm.group(2), re.S):
            ref, attrs, inner = cm.groups()
            tm = re.search(r't="([^"]+)"', attrs)
            ctype = tm.group(1) if tm else None
            if ctype == "inlineStr":  # <is><t>...</t></is>
                ts = re.findall(r"<t[^>]*>(.*?)</t>", inner, re.S)
                if not ts:
                    continue
                val = "".join(ts)
            else:  # numeric / formula-string ("str") / shared-string ("s")
                vm = re.search(r"<v>(.*?)</v>", inner, re.S)
                if vm is None:
                    continue
                val = shared[int(vm.group(1))] if ctype == "s" else vm.group(1)
            grid[(rn, _col_letter(ref))] = unescape(val)
    return grid


def _to_float(s: str | None) -> float | None:
    if s is None:
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def parse_procal_supplement(xlsx_path: str | Path) -> pa.Table:
    """Parse PROCAL Table S1 (sheet ``PROCAL Sequences``) → ``PROCAL_TABLE``.

    Columns (published layout): A peptide number, B sequence, C (M+2H)2+,
    D (M+3H)3+, E expected charge, F SSRCalc HI, H 60-min-gradient RT. Rows
    where column B is a bare-AA sequence are kept.
    """
    grid = _read_sheet_grid(Path(xlsx_path), _PROCAL_SHEET_NAME)
    rows = sorted({r for r, _ in grid})
    cols = {
        "peptide_number": [],
        "sequence": [],
        "mz_2plus": [],
        "mz_3plus": [],
        "expected_charge": [],
        "ssrcalc_hi": [],
        "rt_60min_gradient": [],
    }
    for r in rows:
        seq = grid.get((r, "B"))
        if not seq or not re.fullmatch(r"[A-Z]+", seq):
            continue  # header / title / blank rows
        num = grid.get((r, "A"))
        cols["peptide_number"].append(
            int(num) if num and num.isdigit() else len(cols["sequence"]) + 1
        )
        cols["sequence"].append(seq)
        cols["mz_2plus"].append(_to_float(grid.get((r, "C"))))
        cols["mz_3plus"].append(_to_float(grid.get((r, "D"))))
        cols["expected_charge"].append(grid.get((r, "E")))
        cols["ssrcalc_hi"].append(_to_float(grid.get((r, "F"))))
        cols["rt_60min_gradient"].append(grid.get((r, "H")))

    return pa.table(
        {k: pa.array(v, type=PROCAL_TABLE.field(k).type) for k, v in cols.items()},
        schema=PROCAL_TABLE,
    )


# ── committed-asset access ──────────────────────────────────────────────


def procal_path() -> Path:
    """Path to the committed PROCAL parquet."""
    return _DATA_DIR / "procal.parquet"


def write_procal(table: pa.Table) -> Path:
    """Persist the PROCAL table (cast to schema)."""
    from constellation.core.io.schemas import cast_to_schema

    out = procal_path()
    out.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(cast_to_schema(table, PROCAL_TABLE), out)
    return out


def load_procal() -> pa.Table:
    """Load the committed PROCAL reference (pyarrow only).

    Raises ``FileNotFoundError`` with a build hint if absent (rebuild from the
    manual Excel supplement via ``parse_procal_supplement`` /
    ``pipelines/05_peptide_reference.py``)."""
    p = procal_path()
    if not p.is_file():
        raise FileNotFoundError(
            f"no PROCAL reference at {p}; rebuild from the supplement "
            f"({PROCAL_SUPPLEMENT_FILE}) via parse_procal_supplement / "
            f"`python pipelines/05_peptide_reference.py --procal <xlsx>`"
        )
    return pq.read_table(p)


def procal_sequences() -> frozenset[str]:
    """The 40 PROCAL peptide sequences (bare AA) — the ``is_procal`` annotation set."""
    return frozenset(load_procal().column("sequence").to_pylist())


__all__ = [
    "PROCAL_SCHEMA_VERSION",
    "PROCAL_TABLE",
    "PROCAL_SUPPLEMENT_FILE",
    "parse_procal_supplement",
    "procal_path",
    "write_procal",
    "load_procal",
    "procal_sequences",
]
