"""Shared CLI plumbing for the ``pipelines/`` stage scripts.

Config loading (``config/datasets.toml``), dataset-selection arguments, and a
progress-callback factory. Kept tiny on purpose — the real work lives in the
``protostar`` package; these scripts are just argparse front-ends.

The scripts assume ``protostar`` is importable (``pip install -e .`` per
``environment.yml``); this module additionally puts the repo root on
``sys.path`` so a fresh checkout works before any install.
"""

from __future__ import annotations

import argparse
import os
import sys
import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from constellation.core.progress import (  # noqa: E402 — after sys.path bootstrap
    NullProgress,
    ProgressCallback,
    StreamProgress,
)

DEFAULT_CONFIG = REPO_ROOT / "config" / "datasets.toml"
DEFAULT_MANIFEST_DIR = REPO_ROOT / "config" / "manifests"
DEFAULT_OSC_CONFIG = REPO_ROOT / "config" / "osc.toml"


def load_config(path: "str | Path | None" = None) -> dict:
    path = Path(path) if path is not None else DEFAULT_CONFIG
    with open(path, "rb") as f:
        return tomllib.load(f)


def dataset_specs(config: dict) -> dict[str, dict]:
    """Map dataset name → its config table (accession, n_pools, ...)."""
    return dict(config.get("datasets", {}))


def load_osc_config(path: "str | Path | None" = None) -> dict:
    """Lab-specific OSC settings from ``config/osc.toml`` (gitignored; optional).

    Returns ``{}`` when absent — e.g. a fresh public clone that ships only
    ``config/osc.example.toml``.
    """
    p = Path(path) if path is not None else DEFAULT_OSC_CONFIG
    if not p.exists():
        return {}
    with open(p, "rb") as f:
        return tomllib.load(f)


def data_root(config: dict, override: "str | Path | None" = None) -> Path:
    """Resolve the data root: --data-root > $PROTOSTAR_DATA_ROOT > osc.toml > datasets.toml."""
    if override is not None:
        return Path(override)
    if env := os.environ.get("PROTOSTAR_DATA_ROOT"):
        return Path(env)
    if osc := load_osc_config().get("paths", {}).get("data"):
        return Path(osc)
    return Path(config.get("defaults", {}).get("data_root", "data"))


def seed_from_default(config: dict) -> str | None:
    """--seed-from default: $PROTOSTAR_SEED_FROM > osc.toml [paths].seed > datasets.toml."""
    return (
        os.environ.get("PROTOSTAR_SEED_FROM")
        or load_osc_config().get("paths", {}).get("seed")
        or config.get("defaults", {}).get("seed_from")
    )


def resolve_datasets(values: "list[str] | None", config: dict) -> list[str]:
    """Expand the ``--dataset`` selection (``all`` → every configured dataset)."""
    specs = dataset_specs(config)
    if not values or "all" in values:
        return list(specs)
    unknown = [v for v in values if v not in specs]
    if unknown:
        raise SystemExit(
            f"unknown dataset(s): {', '.join(unknown)}; choose from {', '.join(specs)} or 'all'"
        )
    return list(values)


def add_dataset_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--dataset",
        action="append",
        metavar="NAME",
        help="dataset to act on (repeatable; 'all' for every configured dataset). Default: all.",
    )


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="path to datasets.toml")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="override the data root (default: [defaults].data_root from config)",
    )
    parser.add_argument(
        "--manifest-dir",
        type=Path,
        default=DEFAULT_MANIFEST_DIR,
        help="directory holding the committed per-dataset manifests",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")


def make_progress(quiet: bool) -> ProgressCallback:
    return NullProgress() if quiet else StreamProgress()


__all__ = [
    "DEFAULT_CONFIG",
    "DEFAULT_MANIFEST_DIR",
    "REPO_ROOT",
    "add_common_args",
    "add_dataset_arg",
    "data_root",
    "dataset_specs",
    "load_config",
    "make_progress",
    "resolve_datasets",
    "seed_from_default",
]
