# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""SHA-256 manifest of the granules, outputs, and source behind a run."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def source_sha256(repo_dir: str | Path | None = None) -> str:
    root = Path(repo_dir or Path(__file__).resolve().parents[2])
    paths = [root / "pyproject.toml", root / "Makefile"]
    source_types = {".py", ".html", ".css", ".js"}
    paths += sorted(p for p in (root / "src").rglob("*") if p.suffix in source_types)
    paths += sorted(p for p in (root / "scripts").rglob("*.py"))
    h = hashlib.sha256()
    for path in paths:
        h.update(str(path.relative_to(root)).encode())
        h.update(b"\0")
        h.update(bytes.fromhex(sha256_file(path)))
    return h.hexdigest()


def build_provenance(out_dir: str | Path, paths: list[str]) -> dict:
    out = Path(out_dir)
    granules = [
        {"file": Path(p).name, "bytes": Path(p).stat().st_size, "sha256": sha256_file(p)}
        for p in sorted(paths)
    ]
    outputs = [
        {"file": p.name, "sha256": sha256_file(p)}
        for p in sorted(out.glob("*.json"))
        if p.name != "provenance.json"
    ]
    record = {
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_sha256": source_sha256(),
        "granules": granules,
        "outputs": outputs,
    }
    (out / "provenance.json").write_text(json.dumps(record, indent=2))
    return record
