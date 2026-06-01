"""Tests for the synthetic-peptide reference (parse + schema + load).

The real manuscript workbooks are large + un-redistributed, so the parser is
exercised against tiny synthetic xlsx fixtures that mimic each dataset's sheet
layout (Pool name + Sequence columns; Zolg's quoted sheet titles + QC sheet;
Wilhelm's Identifications RT sheet).
"""

from __future__ import annotations

import openpyxl
import pyarrow.compute as pc
import pytest

from protostar.peptides import reference, sources
from protostar.peptides.reference import PEPTIDE_REFERENCE_TABLE


def _write_xlsx(path, sheets: dict[str, list[list]]):
    """Write a workbook from {sheet_title: [header_row, *data_rows]}."""
    wb = openpyxl.Workbook()
    wb.remove(wb.active)  # drop the default sheet
    for title, rows in sheets.items():
        ws = wb.create_sheet(title=title)
        for row in rows:
            ws.append(row)
    wb.save(path)


def _zolg_like(path):
    # Mirrors Zolg: 3 peptide sheets (quoted titles) + a QC sheet (no pool).
    _write_xlsx(
        path,
        {
            '"Proteotypic Set"': [
                ["Pool name", "Set Name", "Sequence", "QC type", "Gene names"],
                ["TUM_first_pool_1", "Proteotypic", "FSPVLGR", "", "AAAS"],
                ["TUM_first_pool_1", "Proteotypic", "AYASLFR", "", "ABCA4"],
                ["TUM_first_pool_2", "Proteotypic", "ELVISLIVESK", "", "TP53"],
            ],
            '"Missing Gene Set"': [
                ["Pool name", "Set Name", "Sequence", "QC type", "Gene names"],
                ["TUM_second_pool_1", "Missing", "AAAKAPR", "", "LMAN1L"],
            ],
            '"SRMAtlas Set"': [
                ["Pool name", "Set Name", "Sequence", "QC type", "Gene names"],
                ["Thermo_SRM_Pool_1", "SRMAtlas", "AAAHHYGAQCDKPNK", "", "NDUFA8"],
            ],
            "Quality Control Peptides": [
                ["Sequence", "Type", "QC type"],
                ["TASGVGGFSTK", "Retention Time Peptide", "JPT-RT"],
                ["LGGNEQVTR", "PRTC Pierce", "Pierce-RT"],
            ],
        },
    )


def _wilhelm_like(path):
    _write_xlsx(
        path,
        {
            "HLA Class I": [
                ["Pool name", "Sequence"],
                ["TUM_HLA_1", "EEIRKTFNI"],
                ["TUM_HLA_1", "ARTILEENI"],
            ],
            "HLA Class II": [["Pool name", "Sequence"], ["TUM_HLA2_1", "QRDHSAIPVINRAQ"]],
            "AspN": [["Pool name", "Sequence"], ["TUM_aspn_1", "DNLPEAG"]],
            "LysN": [["Pool name", "Sequence"], ["TUM_lysn_1", "KSPLPSQ"]],
            "Identifications": [
                ["Sequence", "Average Retention Time", "Average iRT"],
                ["EEIRKTFNI", 29.56, 53.09],
                ["DNLPEAG", "NaN", "NaN"],  # NaN must become null, not crash
            ],
        },
    )


# ── schema ─────────────────────────────────────────────────────────────


def test_schema_has_version_and_core_fields():
    names = set(PEPTIDE_REFERENCE_TABLE.names)
    assert {"peptide_id", "sequence", "pool", "pool_prefix", "set", "is_qc"} <= names
    assert PEPTIDE_REFERENCE_TABLE.metadata[b"schema_version"] == b"1"


# ── parsing ────────────────────────────────────────────────────────────


def test_parse_zolg_like(tmp_path):
    x = tmp_path / "zolg.xlsx"
    _zolg_like(x)
    t = sources.parse_supplement(x, "Zolg2017")
    assert t.num_rows == 7  # 3 + 1 + 1 proteotypic-ish + 2 QC
    # peptide_id is a contiguous 0-based row index
    assert t.column("peptide_id").to_pylist() == list(range(7))
    # set labels from the sheet map
    assert set(t.column("set").to_pylist()) == {"proteotypic", "missing_gene", "srmatlas", "qc"}
    # pool_prefix strips the trailing _N
    row = {c: t.column(c)[0].as_py() for c in t.column_names}
    assert row["pool"] == "TUM_first_pool_1" and row["pool_prefix"] == "TUM_first_pool"


def test_parse_qc_flagging(tmp_path):
    x = tmp_path / "zolg.xlsx"
    _zolg_like(x)
    t = sources.parse_supplement(x, "Zolg2017")
    qc = t.filter(t.column("is_qc"))
    assert qc.num_rows == 2
    assert set(qc.column("qc_type").to_pylist()) == {"JPT-RT", "Pierce-RT"}
    assert qc.column("pool").null_count == 2  # QC sheet has no pool
    # non-QC rows are flagged False with null qc_type
    nonqc = t.filter(pc.invert(t.column("is_qc")))
    assert nonqc.column("qc_type").null_count == nonqc.num_rows


def test_parse_wilhelm_rt_carry_and_nan(tmp_path):
    x = tmp_path / "wilhelm.xlsx"
    _wilhelm_like(x)
    t = sources.parse_supplement(x, "Wilhelm2021")
    by_seq = {r["sequence"]: r for r in t.to_pylist()}
    assert by_seq["EEIRKTFNI"]["rt"] == pytest.approx(29.56)
    assert by_seq["EEIRKTFNI"]["irt"] == pytest.approx(53.09)
    assert by_seq["DNLPEAG"]["rt"] is None  # "NaN" string → null
    assert by_seq["QRDHSAIPVINRAQ"]["rt"] is None  # absent from Identifications → null


def test_parse_missing_sheet_raises(tmp_path):
    x = tmp_path / "bad.xlsx"
    _write_xlsx(x, {"Wrong Sheet": [["Pool name", "Sequence"], ["P_1", "PEPTIDEK"]]})
    with pytest.raises(ValueError, match="not in workbook"):
        sources.parse_supplement(x, "Zolg2017")


# ── round-trip ─────────────────────────────────────────────────────────


def test_write_and_load_round_trip(tmp_path, monkeypatch):
    x = tmp_path / "wilhelm.xlsx"
    _wilhelm_like(x)
    t = sources.parse_supplement(x, "Wilhelm2021")
    # redirect the committed-data dir to tmp so the test writes nothing real
    monkeypatch.setattr(reference, "_DATA_DIR", tmp_path / "data")
    out = reference.write_reference(t, "Wilhelm2021")
    assert out.is_file()
    loaded = reference.load_reference("Wilhelm2021")
    assert loaded.num_rows == t.num_rows
    assert loaded.schema.equals(PEPTIDE_REFERENCE_TABLE)


def test_load_missing_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(reference, "_DATA_DIR", tmp_path / "empty")
    with pytest.raises(FileNotFoundError, match="no peptide reference"):
        reference.load_reference("Zolg2017")
