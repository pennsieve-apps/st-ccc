"""
Finding the inputs.

Linked layout (placement: subfolder): INPUT_DIR/spatial/<sample_id>.h5ad,
INPUT_DIR/labels/labels.csv and INPUT_DIR/metadata/samples.csv. Without links
the same files are found by name anywhere under INPUT_DIR (symlinks followed).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger("st-ccc")


@dataclass
class Inputs:
    h5ad: dict[str, Path]                       # sample_id (file stem) -> path, sorted
    labels_csv: Path | None
    samples_csv: Path | None
    files_seen: list[Path] = field(default_factory=list)


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinks (INPUT_DIR entries are symlinks), sorted."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        dirnames[:] = sorted(d for d in dirnames if not d.startswith("."))
        found.extend(Path(dirpath) / f for f in sorted(filenames) if not f.startswith("."))
    return found


def _find(root: Path, port: str, predicate) -> list[Path]:
    """Files of a port: its own folder when it exists, else anywhere under INPUT_DIR."""
    folder = root / port
    base = folder if folder.is_dir() else root
    return [p for p in walk(base) if predicate(p)]


def _one(root: Path, port: str, name: str) -> Path | None:
    hits = _find(root, port, lambda p: p.name.lower() == name)
    if len(hits) > 1:
        log.warning("%d files named %s; using %s", len(hits), name, hits[0])
    return hits[0] if hits else None


def discover(input_dir: Path) -> Inputs:
    files = _find(input_dir, "spatial", lambda p: p.suffix.lower() == ".h5ad")
    h5ad: dict[str, Path] = {}
    for p in files:
        if p.stem in h5ad:
            log.warning("two .h5ad files for sample %s; using %s", p.stem, h5ad[p.stem])
            continue
        h5ad[p.stem] = p
    return Inputs(dict(sorted(h5ad.items())), _one(input_dir, "labels", "labels.csv"),
                  _one(input_dir, "metadata", "samples.csv"), files)
