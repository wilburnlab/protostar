"""Shared pytest fixtures.

A tiny Range-capable HTTP server lets the downloader tests exercise full /
resume (206) / unsatisfiable (416) / ignore-Range (200) paths fully offline.
Tests that hit the real PRIDE / Zenodo endpoints are marked ``network`` and
skipped unless ``--run-network`` is passed.
"""

from __future__ import annotations

import http.server
import threading
from dataclasses import dataclass
from pathlib import Path

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-network",
        action="store_true",
        default=False,
        help="run tests that hit live PRIDE / Zenodo endpoints",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "network: hits live external endpoints")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-network"):
        return
    skip = pytest.mark.skip(reason="needs --run-network")
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip)


class _RangeHandler(http.server.BaseHTTPRequestHandler):
    def log_message(self, *args):  # silence
        return

    def do_GET(self):
        root: Path = self.server.serve_dir  # type: ignore[attr-defined]
        target = root / self.path.lstrip("/")
        if not target.is_file():
            self.send_error(404)
            return
        data = target.read_bytes()
        rng = self.headers.get("Range")
        if rng and not self.server.ignore_range:  # type: ignore[attr-defined]
            start = int(rng.split("=", 1)[1].split("-", 1)[0] or 0)
            if start >= len(data):
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{len(data)}")
                self.end_headers()
                return
            chunk = data[start:]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
            self.send_header("Content-Length", str(len(chunk)))
            self.send_header("Accept-Ranges", "bytes")
            self.end_headers()
            self.wfile.write(chunk)
        else:
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Accept-Ranges", "bytes")
            self.end_headers()
            self.wfile.write(data)


@dataclass
class Server:
    base_url: str
    serve_dir: Path
    _httpd: http.server.ThreadingHTTPServer

    def url(self, name: str) -> str:
        return f"{self.base_url}/{name}"

    @property
    def ignore_range(self) -> bool:
        return self._httpd.ignore_range  # type: ignore[attr-defined]

    @ignore_range.setter
    def ignore_range(self, value: bool) -> None:
        self._httpd.ignore_range = value  # type: ignore[attr-defined]


@pytest.fixture
def http_server(tmp_path: Path):
    serve_dir = tmp_path / "srv"
    serve_dir.mkdir()
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _RangeHandler)
    httpd.serve_dir = serve_dir  # type: ignore[attr-defined]
    httpd.ignore_range = False  # type: ignore[attr-defined]
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    host, port = httpd.server_address
    try:
        yield Server(f"http://{host}:{port}", serve_dir, httpd)
    finally:
        httpd.shutdown()
        httpd.server_close()
        thread.join(timeout=5)
