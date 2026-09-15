#!/usr/bin/env python3
"""Independent four-line FASTQ record counter.

Prints one decimal count and a newline, matching `z-fastq count` on valid
four-line input. Gzip is sniffed from the 1f 8b magic, not the suffix.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def line_payload(line: bytes) -> bytes:
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    return line


def count_records(path: Path) -> tuple[int, int]:
    reads = 0
    bases = 0
    with path.open("rb") as raw:
        magic = raw.read(2)
        raw.seek(0)
        stream: object
        if magic == b"\x1f\x8b":
            ctx = gzip.GzipFile(fileobj=raw)
        else:
            ctx = raw
        close = magic == b"\x1f\x8b"
        try:
            handle = ctx if close else raw
            while True:
                lines = tuple(handle.readline() for _ in range(4))
                if not lines[0]:
                    if any(lines[1:]):
                        raise ValueError(f"partial FASTQ boundary in {path}")
                    break
                if any(not line for line in lines[1:]):
                    raise ValueError(f"truncated FASTQ record {reads} in {path}")
                header, sequence, plus, quality = tuple(map(line_payload, lines))
                if not header.startswith(b"@"):
                    raise ValueError(f"record {reads} has no FASTQ header in {path}")
                if not plus.startswith(b"+"):
                    raise ValueError(f"record {reads} has no FASTQ plus line in {path}")
                if len(sequence) != len(quality):
                    raise ValueError(
                        f"record {reads} has unequal sequence and quality in {path}"
                    )
                reads += 1
                bases += len(sequence)
        finally:
            if close:
                handle.close()
    return reads, bases


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", type=Path, help="plain or gzip FASTQ")
    ap.add_argument(
        "--bases",
        action="store_true",
        help="print records<TAB>bases instead of records only",
    )
    args = ap.parse_args()
    if not args.path.is_file():
        print(f"error: file not found: {args.path}", file=sys.stderr)
        raise SystemExit(1)
    try:
        reads, bases = count_records(args.path)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
    if args.bases:
        sys.stdout.write(f"{reads}\t{bases}\n")
    else:
        sys.stdout.write(f"{reads}\n")


if __name__ == "__main__":
    main()
