"""Tests for the streaming resumable downloader + checksum helpers."""

from __future__ import annotations

import hashlib

import pytest

from protostar.fetch import net


def test_hashes_and_verify(tmp_path):
    f = tmp_path / "x.bin"
    payload = b"protostar" * 1000
    f.write_bytes(payload)
    assert net.sha1_of(f) == hashlib.sha1(payload).hexdigest()
    assert net.md5_of(f) == hashlib.md5(payload).hexdigest()

    assert net.verify(f, expected_size=len(payload)) == "ok"
    assert net.verify(f, expected_sha1=hashlib.sha1(payload).hexdigest()) == "ok"
    assert net.verify(f, expected_size=len(payload) + 1) == "corrupt"
    assert net.verify(f, expected_sha1="deadbeef") == "corrupt"
    assert net.verify(tmp_path / "missing.bin") == "missing"

    part = tmp_path / "y.bin.part"
    part.write_bytes(b"half")
    assert net.verify(tmp_path / "y.bin") == "partial"


def test_download_full(http_server, tmp_path):
    payload = bytes(range(256)) * 400
    (http_server.serve_dir / "a.bin").write_bytes(payload)
    dest = tmp_path / "a.bin"
    res = net.stream_download(http_server.url("a.bin"), dest, expected_size=len(payload))
    assert res.resumed is False
    assert dest.read_bytes() == payload
    assert res.sha1 == hashlib.sha1(payload).hexdigest()
    assert not (tmp_path / "a.bin.part").exists()


def test_download_skip_when_complete(http_server, tmp_path):
    payload = b"already-here" * 100
    (http_server.serve_dir / "b.bin").write_bytes(payload)
    dest = tmp_path / "b.bin"
    dest.write_bytes(payload)  # pre-existing, correct size
    res = net.stream_download(http_server.url("b.bin"), dest, expected_size=len(payload))
    assert res.resumed is False and res.n_bytes == len(payload)


def test_download_resume_206(http_server, tmp_path):
    payload = bytes(range(256)) * 1000
    (http_server.serve_dir / "c.bin").write_bytes(payload)
    dest = tmp_path / "c.bin"
    part = tmp_path / "c.bin.part"
    part.write_bytes(payload[: len(payload) // 3])  # simulate a partial download
    res = net.stream_download(http_server.url("c.bin"), dest, expected_size=len(payload))
    assert res.resumed is True
    assert dest.read_bytes() == payload
    assert res.sha1 == hashlib.sha1(payload).hexdigest()


def test_download_stale_oversized_part(http_server, tmp_path):
    payload = b"data" * 500
    (http_server.serve_dir / "d.bin").write_bytes(payload)
    dest = tmp_path / "d.bin"
    (tmp_path / "d.bin.part").write_bytes(b"x" * (len(payload) + 10))  # oversized
    res = net.stream_download(http_server.url("d.bin"), dest, expected_size=len(payload))
    assert res.resumed is False
    assert dest.read_bytes() == payload


def test_download_server_ignores_range(http_server, tmp_path):
    payload = b"ignore-range" * 100
    (http_server.serve_dir / "e.bin").write_bytes(payload)
    http_server.ignore_range = True  # server returns 200 even with a Range request
    dest = tmp_path / "e.bin"
    (tmp_path / "e.bin.part").write_bytes(payload[:50])
    res = net.stream_download(http_server.url("e.bin"), dest, expected_size=len(payload))
    assert res.resumed is False  # fell back to a clean full download
    assert dest.read_bytes() == payload


def test_download_size_mismatch_raises(http_server, tmp_path):
    (http_server.serve_dir / "f.bin").write_bytes(b"short")
    with pytest.raises(OSError, match="size mismatch"):
        net.stream_download(http_server.url("f.bin"), tmp_path / "f.bin", expected_size=999)
