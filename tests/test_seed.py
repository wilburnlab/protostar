"""Tests for --seed-from relocation logic."""

from __future__ import annotations

from protostar.fetch import seed
from protostar.fetch.manifest import Manifest, ManifestEntry, local_path


def _manifest():
    raw = ManifestEntry("pool_A.raw", 100, "a" * 40, "https://x/pool_A.raw", "RAW")
    missing = ManifestEntry("pool_B.raw", 100, "b" * 40, "https://x/pool_B.raw", "RAW")
    return Manifest("Zolg2017", "PXD004732", "2026-01-01T00:00:00+00:00", "v3", (raw, missing))


def _cartographer(root, name, size):
    src = root / "Zolg2017" / "raw" / name
    src.parent.mkdir(parents=True, exist_ok=True)
    src.write_bytes(b"x" * size)
    return src


def test_seed_move(tmp_path):
    m = _manifest()
    data_root, seed_from = tmp_path / "data", tmp_path / "cart"
    src = _cartographer(seed_from, "pool_A.raw", 100)
    results = {
        r.entry.file_name: r for r in seed.seed_dataset(m, data_root, seed_from, mode="move")
    }
    assert results["pool_A.raw"].action == "seeded"
    assert results["pool_B.raw"].action == "missing_source"
    assert local_path(m.entries[0], data_root, "Zolg2017").is_file()
    assert not src.exists()  # move empties the source

    # idempotent re-run: dest now present
    again = {r.entry.file_name: r for r in seed.seed_dataset(m, data_root, seed_from, mode="move")}
    assert again["pool_A.raw"].action == "present"


def test_seed_hardlink_keeps_source(tmp_path):
    m = _manifest()
    data_root, seed_from = tmp_path / "data", tmp_path / "cart"
    src = _cartographer(seed_from, "pool_A.raw", 100)
    seed.seed_dataset(m, data_root, seed_from, mode="hardlink")
    dest = local_path(m.entries[0], data_root, "Zolg2017")
    assert dest.is_file() and src.exists()
    assert dest.stat().st_ino == src.stat().st_ino  # same inode


def test_seed_size_mismatch_not_seeded(tmp_path):
    m = _manifest()
    data_root, seed_from = tmp_path / "data", tmp_path / "cart"
    _cartographer(seed_from, "pool_A.raw", 99)  # wrong size
    results = {
        r.entry.file_name: r for r in seed.seed_dataset(m, data_root, seed_from, mode="move")
    }
    assert results["pool_A.raw"].action == "size_mismatch"
    assert not local_path(m.entries[0], data_root, "Zolg2017").exists()


def test_seed_dry_run_changes_nothing(tmp_path):
    m = _manifest()
    data_root, seed_from = tmp_path / "data", tmp_path / "cart"
    src = _cartographer(seed_from, "pool_A.raw", 100)
    results = {
        r.entry.file_name: r for r in seed.seed_dataset(m, data_root, seed_from, dry_run=True)
    }
    assert results["pool_A.raw"].action == "seeded"  # would-seed
    assert src.exists() and not local_path(m.entries[0], data_root, "Zolg2017").exists()
