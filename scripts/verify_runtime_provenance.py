#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs" / "runtime_provenance.tsv"
SKIP_PARTS = {"__pycache__", "build", ".git"}
SKIP_NAMES = {".DS_Store"}
SKIP_SUFFIXES = {".pyc", ".pyo", ".so", ".dylib", ".pyd", ".dll", ".o"}
SKIP_RELS = {"lib/distancelib.c", "lib/qcprot.c"}


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_sha256(path: Path) -> str:
    h = hashlib.sha256()
    for item in sorted(p for p in path.rglob("*") if p.is_file()):
        rel_parts = item.relative_to(path).parts
        rel = item.relative_to(path).as_posix()
        if any(part in SKIP_PARTS for part in rel_parts):
            continue
        if rel in SKIP_RELS or item.name in SKIP_NAMES or item.suffix in SKIP_SUFFIXES:
            continue
        size = item.stat().st_size
        h.update(
            rel.encode("utf-8")
            + b"\0"
            + file_sha256(item).encode("ascii")
            + b"\0"
            + str(size).encode("ascii")
            + b"\n"
        )
    return h.hexdigest()


def verify(manifest: Path) -> list[str]:
    failures: list[str] = []
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hash_type = row["hash_type"]
            if hash_type == "not_vendored":
                continue
            path = ROOT / row["local_path"]
            if hash_type == "file_sha256":
                if not path.is_file():
                    failures.append(f"{row['component']}: missing file {row['local_path']}")
                    continue
                got = file_sha256(path)
            elif hash_type == "tree_sha256":
                if not path.is_dir():
                    failures.append(f"{row['component']}: missing directory {row['local_path']}")
                    continue
                got = tree_sha256(path)
            else:
                failures.append(f"{row['component']}: unsupported hash_type {hash_type}")
                continue
            if got != row["hash_value"]:
                failures.append(f"{row['component']}: checksum mismatch for {row['local_path']}")
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify vendored runtime provenance checksums.")
    parser.add_argument("--manifest", default=str(MANIFEST))
    args = parser.parse_args(argv)
    failures = verify(Path(args.manifest))
    if failures:
        for failure in failures:
            print(f"FAIL {failure}", file=sys.stderr)
        return 1
    print("PASS runtime provenance checksums")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
