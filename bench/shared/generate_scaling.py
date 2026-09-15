#!/usr/bin/env python3
"""Generate synthetic scaling FASTQ into bench/shared/cache/scaling/.

Families:
  size        fixed 100k records, grow read length to target file size
  reads       fixed 150 bp, grow record count
  gzip        gzip of reads_fixed_* at compresslevel 6

Generate-if-missing unless --force. Stamp fingerprints parameters.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "bench/shared/cache/scaling"
STAMP = OUT_DIR / ".stamp"

DNA = "ACGTACGTAC" * 16
QUAL_BYTE = "I"
SIZE_MBS = (1, 5, 10, 25, 50, 100, 250, 500)
READS_FIXED_COUNTS = (100_000, 250_000, 500_000, 1_000_000)
SIZE_READ_COUNT = 100_000
FIXED_READ_LEN = 150
GZIP_LEVEL = 6
HEADER_WIDTH = 17  # "@r" + 15-digit index
SCHEMA = "fastq-scaling.v1"


def record_bytes(seq_len: int) -> int:
    # @r{i:015d}\n + seq\n + +\n + qual\n
    return HEADER_WIDTH + 4 + 2 * seq_len


def seq_for(length: int) -> str:
    if length <= 0:
        raise SystemExit("error: sequence length must be positive")
    repeats, rest = divmod(length, len(DNA))
    return DNA * repeats + DNA[:rest]


def write_fastq(path: Path, count: int, seq_len: int) -> None:
    seq = seq_for(seq_len)
    qual = QUAL_BYTE * seq_len
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for i in range(1, count + 1):
            handle.write(f"@r{i:015d}\n{seq}\n+\n{qual}\n")


def gzip_copy(src: Path, dest: Path) -> None:
    with src.open("rb") as incoming, gzip.open(dest, "wb", compresslevel=GZIP_LEVEL, mtime=0) as outgoing:
        while True:
            chunk = incoming.read(1024 * 1024)
            if not chunk:
                break
            outgoing.write(chunk)


def length_for_size_mb(mb: int) -> int:
    total = mb * 1024 * 1024
    per = max(record_bytes(1), total // SIZE_READ_COUNT)
    seq_len = max(1, (per - HEADER_WIDTH - 4) // 2)
    return seq_len


def expected_paths(mode: str) -> list[Path]:
    paths: list[Path] = []
    if mode in ("all", "size"):
        paths.extend(OUT_DIR / f"size_{mb}mb.fastq" for mb in SIZE_MBS)
    if mode in ("all", "reads", "seq"):
        paths.extend(OUT_DIR / f"reads_fixed_{c}.fastq" for c in READS_FIXED_COUNTS)
    if mode in ("all", "gzip"):
        paths.extend(OUT_DIR / f"reads_fixed_{c}.fastq.gz" for c in READS_FIXED_COUNTS)
    return paths


def stamp_payload() -> str:
    parts = [
        f"schema={SCHEMA}",
        f"size_mbs={','.join(map(str, SIZE_MBS))}",
        f"reads_fixed_counts={','.join(map(str, READS_FIXED_COUNTS))}",
        f"size_read_count={SIZE_READ_COUNT}",
        f"fixed_read_len={FIXED_READ_LEN}",
        f"gzip_level={GZIP_LEVEL}",
        f"header_width={HEADER_WIDTH}",
        f"dna_hash={hashlib.sha256(DNA.encode()).hexdigest()[:16]}",
        f"qual_byte={QUAL_BYTE}",
    ]
    return "\n".join(parts) + "\n"


def stamp_ok() -> bool:
    return STAMP.is_file() and STAMP.read_text(encoding="utf-8") == stamp_payload()


def outputs_ready(mode: str) -> bool:
    return stamp_ok() and all(path.is_file() for path in expected_paths(mode))


def generate(*, mode: str, force: bool) -> None:
    if mode not in ("all", "size", "reads", "seq", "gzip"):
        raise SystemExit(f"error: mode must be all|size|reads|gzip, got {mode!r}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if outputs_ready(mode) and not force:
        print(f"scaling cache ok: {OUT_DIR} ({mode})")
        return

    if mode in ("all", "size"):
        for mb in SIZE_MBS:
            path = OUT_DIR / f"size_{mb}mb.fastq"
            if force or not path.is_file():
                seq_len = length_for_size_mb(mb)
                print(f"  size_{mb}mb.fastq  n={SIZE_READ_COUNT}  len={seq_len}")
                write_fastq(path, SIZE_READ_COUNT, seq_len)

    if mode in ("all", "reads", "seq", "gzip"):
        for count in READS_FIXED_COUNTS:
            path = OUT_DIR / f"reads_fixed_{count}.fastq"
            if force or not path.is_file():
                print(f"  reads_fixed_{count}.fastq  n={count}  len={FIXED_READ_LEN}")
                write_fastq(path, count, FIXED_READ_LEN)

    if mode in ("all", "gzip"):
        for count in READS_FIXED_COUNTS:
            src = OUT_DIR / f"reads_fixed_{count}.fastq"
            dest = OUT_DIR / f"reads_fixed_{count}.fastq.gz"
            if force or not dest.is_file():
                if not src.is_file():
                    write_fastq(src, count, FIXED_READ_LEN)
                print(f"  gzip {dest.name}  level={GZIP_LEVEL}")
                gzip_copy(src, dest)

    STAMP.write_text(stamp_payload(), encoding="utf-8")
    print(f"wrote scaling FASTQ -> {OUT_DIR} ({mode})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--mode",
        choices=("all", "size", "reads", "seq", "gzip"),
        default="all",
        help="which fixture families to ensure (seq is an alias of reads)",
    )
    ap.add_argument("--force", action="store_true", help="rebuild even if present")
    args = ap.parse_args()
    generate(mode=args.mode, force=args.force)


if __name__ == "__main__":
    main()
