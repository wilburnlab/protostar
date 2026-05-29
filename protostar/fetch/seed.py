"""Relocate existing local ``.raw`` copies into the canonical data tree.

The lab already holds the ProteomeTools ``.raw`` files at
``<seed_from>/<dataset>/raw/`` (the Cartographer tree). Since protostar fully
supersedes that work and the destination is the same ESS filesystem, the
default is a **move** (instant, 0 extra bytes, empties the source). ``hardlink``
keeps the source as a backup; ``copy`` is the cross-filesystem fallback.

Files are matched against the expected manifest by name + exact size (a
truncated copy fails the size check); ``verify=True`` adds the SHA-1 check.
Only RAW is seeded — the Cartographer tree has no SEARCH outputs, so those are
fetched fresh.
"""

from __future__ import annotations

import errno
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal

from constellation.core.progress import ProgressCallback, emit_done, emit_progress, emit_start

from . import net
from .manifest import Manifest, ManifestEntry, local_path

SeedMode = Literal["move", "hardlink", "copy"]
SeedAction = Literal["seeded", "present", "missing_source", "size_mismatch", "corrupt"]
_STAGE = "seed"


@dataclass(frozen=True, slots=True)
class SeedResult:
    entry: ManifestEntry
    action: SeedAction
    source: Path | None
    dest: Path


def _candidate_sources(entry: ManifestEntry, seed_from: Path, dataset: str) -> Iterator[Path]:
    # Canonical Cartographer layout first, then looser fallbacks.
    yield seed_from / dataset / "raw" / entry.file_name
    yield seed_from / dataset / entry.file_name
    yield seed_from / entry.file_name


def find_source(entry: ManifestEntry, seed_from: "str | Path", dataset: str) -> Path | None:
    for cand in _candidate_sources(entry, Path(seed_from), dataset):
        if cand.is_file():
            return cand
    return None


def _place(src: Path, dest: Path, mode: SeedMode) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        if mode == "move":
            os.replace(src, dest)
        elif mode == "hardlink":
            os.link(src, dest)
        elif mode == "copy":
            shutil.copy2(src, dest)
        else:  # pragma: no cover - guarded by argparse choices
            raise ValueError(f"unknown seed mode: {mode!r}")
    except OSError as exc:
        if exc.errno == errno.EXDEV:
            raise OSError(
                f"cannot {mode} across filesystems ({src} -> {dest}); re-run with --seed-mode copy"
            ) from exc
        raise


def seed_dataset(
    manifest: Manifest,
    data_root: "str | Path",
    seed_from: "str | Path",
    *,
    mode: SeedMode = "move",
    categories: tuple[str, ...] = ("RAW",),
    verify: bool = False,
    dry_run: bool = False,
    progress_cb: ProgressCallback | None = None,
) -> list[SeedResult]:
    """Match manifest entries to local copies and relocate them.

    ``dry_run`` classifies without touching the filesystem. Returns one
    :class:`SeedResult` per considered entry.
    """
    entries = manifest.fetch_entries(categories)
    emit_start(progress_cb, _STAGE, total=len(entries), message=f"{manifest.dataset} ({mode})")
    results: list[SeedResult] = []
    for i, e in enumerate(entries, 1):
        dest = local_path(e, data_root, manifest.dataset)
        action: SeedAction
        src: Path | None = find_source(e, seed_from, manifest.dataset)

        if dest.is_file() and dest.stat().st_size == e.size_bytes:
            action, src = "present", None
        elif src is None:
            action = "missing_source"
        elif src.stat().st_size != e.size_bytes:
            action = "size_mismatch"
        elif verify and e.sha1 and net.sha1_of(src).lower() != e.sha1.lower():
            action = "corrupt"
        else:
            action = "seeded"
            if not dry_run:
                if dest.exists():  # wrong-size leftover from an aborted seed
                    dest.unlink()
                _place(src, dest, mode)
        results.append(SeedResult(e, action, src, dest))
        emit_progress(
            progress_cb,
            _STAGE,
            completed=i,
            total=len(entries),
            message=f"{e.file_name} → {action}",
        )

    n_seeded = sum(1 for r in results if r.action == "seeded")
    n_missing = sum(1 for r in results if r.action == "missing_source")
    emit_done(
        progress_cb,
        _STAGE,
        completed=len(results),
        total=len(results),
        message=f"{n_seeded} seeded, {n_missing} not found locally",
    )
    return results


__all__ = ["SeedAction", "SeedMode", "SeedResult", "find_source", "seed_dataset"]
