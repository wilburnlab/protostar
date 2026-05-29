"""Tests for the Zenodo library client (mode parsing, selection, extraction)."""

from __future__ import annotations

import zipfile

import pytest

from protostar.library import zenodo
from protostar.library.zenodo import ZenodoFile


def test_mode_token():
    z = ZenodoFile("FTMS_HCD_28_annotated_2019-11-12.zip", 0, "", "")
    assert z.mode == "FTMS_HCD_28"
    assert ZenodoFile("ITMS_CID_35_annotated_2019-11-13.zip", 0, "", "").mode == "ITMS_CID_35"


def test_select():
    files = [
        ZenodoFile("FTMS_HCD_28_annotated_2019-11-12.zip", 1, "", ""),
        ZenodoFile("ITMS_CID_35_annotated_2019-11-13.zip", 1, "", ""),
    ]
    assert len(zenodo._select(files, None)) == 2
    sel = zenodo._select(files, {"ITMS_CID_35"})
    assert [f.mode for f in sel] == ["ITMS_CID_35"]


def test_extract_msp_flattens_into_mode_dir(tmp_path):
    zp = tmp_path / "FTMS_HCD_28_annotated_2019-11-12.zip"
    with zipfile.ZipFile(zp, "w") as zf:
        zf.writestr("nested/dir/FTMS_HCD_28.msp", "Name: PEP/2\nNum peaks: 0\n")
        zf.writestr("readme.txt", "ignore")
    out = zenodo.extract_msp([zp], tmp_path / "libraries")
    assert len(out) == 1
    assert out[0].name == "FTMS_HCD_28.msp"
    assert out[0].parent.name == "FTMS_HCD_28"  # mode subdir
    assert out[0].read_text().startswith("Name: PEP/2")


@pytest.mark.network
def test_list_library_files_live():
    files = zenodo.list_library_files()
    assert len(files) == 8
    assert all(len(f.md5) == 32 and f.url for f in files)
    assert {f.mode for f in files} >= {"FTMS_HCD_28", "ITMS_CID_35"}
