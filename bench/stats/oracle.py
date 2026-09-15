#!/usr/bin/env python3
"""Independent four-line FASTQ stats aggregator.

Prints the same human fields as `z-fastq stats` on valid four-line input
(Phred+33). Gzip is sniffed from the 1f 8b magic, not the suffix.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

RATIO_SCALE = 1_000_000
HUMAN_KEYS = (
    "reads",
    "bases",
    "min_length",
    "max_length",
    "mean_length",
    "a",
    "c",
    "g",
    "t",
    "n",
    "other_bases",
    "gc_fraction",
    "quality_sum",
    "mean_quality",
    "q20_bases",
    "q20_fraction",
    "q30_bases",
    "q30_fraction",
)


def line_payload(line: bytes) -> bytes:
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    return line


def ratio_text(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "-"
    rounded = (numerator * RATIO_SCALE + denominator // 2) // denominator
    return f"{rounded // RATIO_SCALE}.{rounded % RATIO_SCALE:06d}"


def collect(path: Path) -> dict[str, int]:
    reads = 0
    bases = 0
    min_length = 0
    max_length = 0
    a = c = g = t = n = other_bases = 0
    quality_sum = 0
    q20_bases = 0
    q30_bases = 0
    with path.open("rb") as raw:
        magic = raw.read(2)
        raw.seek(0)
        close = magic == b"\x1f\x8b"
        handle: object = gzip.GzipFile(fileobj=raw) if close else raw
        try:
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
                length = len(sequence)
                if reads == 0:
                    min_length = length
                    max_length = length
                else:
                    if length < min_length:
                        min_length = length
                    if length > max_length:
                        max_length = length
                reads += 1
                bases += length
                for base, encoded in zip(sequence, quality):
                    if base in (ord("A"), ord("a")):
                        a += 1
                    elif base in (ord("C"), ord("c")):
                        c += 1
                    elif base in (ord("G"), ord("g")):
                        g += 1
                    elif base in (ord("T"), ord("t")):
                        t += 1
                    elif base in (ord("N"), ord("n")):
                        n += 1
                    else:
                        other_bases += 1
                    if not 33 <= encoded <= 126:
                        raise ValueError(
                            f"record {reads - 1} has invalid Phred+33 byte {encoded}"
                        )
                    score = encoded - 33
                    quality_sum += score
                    if score >= 20:
                        q20_bases += 1
                    if score >= 30:
                        q30_bases += 1
        finally:
            if close:
                handle.close()
    return {
        "reads": reads,
        "bases": bases,
        "min_length": min_length,
        "max_length": max_length,
        "a": a,
        "c": c,
        "g": g,
        "t": t,
        "n": n,
        "other_bases": other_bases,
        "quality_sum": quality_sum,
        "q20_bases": q20_bases,
        "q30_bases": q30_bases,
    }


def format_human(values: dict[str, int], input_label: str) -> str:
    reads = values["reads"]
    bases = values["bases"]
    lines = [
        f"input: {input_label}",
        f"reads: {reads}",
        f"bases: {bases}",
        f"min_length: {values['min_length'] if reads else '-'}",
        f"max_length: {values['max_length'] if reads else '-'}",
        f"mean_length: {ratio_text(bases, reads)}",
        f"a: {values['a']}",
        f"c: {values['c']}",
        f"g: {values['g']}",
        f"t: {values['t']}",
        f"n: {values['n']}",
        f"other_bases: {values['other_bases']}",
        f"gc_fraction: {ratio_text(values['g'] + values['c'], values['a'] + values['c'] + values['g'] + values['t'])}",
        f"quality_sum: {values['quality_sum']}",
        f"mean_quality: {ratio_text(values['quality_sum'], bases)}",
        f"q20_bases: {values['q20_bases']}",
        f"q20_fraction: {ratio_text(values['q20_bases'], bases)}",
        f"q30_bases: {values['q30_bases']}",
        f"q30_fraction: {ratio_text(values['q30_bases'], bases)}",
    ]
    return "\n".join(lines) + "\n"


def parse_human(text: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        if key == "input" or not key:
            continue
        fields[key] = value.strip()
    return fields


def parse_tsv_row(text: str) -> dict[str, str]:
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("expected a TSV header and one data row")
    header = lines[0].split("\t")
    row = lines[1].split("\t")
    if len(row) < len(header):
        raise ValueError("TSV data row is shorter than the header")
    return {name: row[idx] for idx, name in enumerate(header)}


def require_equal(label: str, left: str, right: str) -> None:
    if left != right:
        raise ValueError(f"{label}: expected {left!r}, got {right!r}")


def as_int(label: str, text: str) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError(f"{label} is not an integer: {text!r}") from exc


def as_float(label: str, text: str) -> float:
    try:
        return float(text)
    except ValueError as exc:
        raise ValueError(f"{label} is not a number: {text!r}") from exc


def agree_human(left_text: str, right_text: str) -> None:
    left = parse_human(left_text)
    right = parse_human(right_text)
    missing = [key for key in HUMAN_KEYS if key not in left or key not in right]
    if missing:
        raise ValueError("missing fields: " + ", ".join(missing))
    extra = sorted((set(left) | set(right)) - set(HUMAN_KEYS) - {"input"})
    if extra:
        raise ValueError("unexpected fields: " + ", ".join(extra))
    for key in HUMAN_KEYS:
        require_equal(key, left[key], right[key])


def overlap_seqkit(human_text: str, seqkit_text: str) -> None:
    human = parse_human(human_text)
    row = parse_tsv_row(seqkit_text)
    require_equal("reads/num_seqs", human["reads"], row["num_seqs"])
    require_equal("bases/sum_len", human["bases"], row["sum_len"])
    require_equal("min_length/min_len", human["min_length"], row["min_len"])
    require_equal("max_length/max_len", human["max_length"], row["max_len"])
    if "sum_n" in row:
        require_equal("n/sum_n", human["n"], row["sum_n"])
    q20 = 100.0 * as_float("q20_fraction", human["q20_fraction"])
    q30 = 100.0 * as_float("q30_fraction", human["q30_fraction"])
    gc = 100.0 * as_float("gc_fraction", human["gc_fraction"])
    if abs(as_float("Q20(%)", row["Q20(%)"]) - q20) >= 1.0:
        raise ValueError(f"Q20(%): seqkit {row['Q20(%)']!r} vs 100*q20_fraction {q20:.6f}")
    if abs(as_float("Q30(%)", row["Q30(%)"]) - q30) >= 1.0:
        raise ValueError(f"Q30(%): seqkit {row['Q30(%)']!r} vs 100*q30_fraction {q30:.6f}")
    if abs(as_float("GC(%)", row["GC(%)"]) - gc) >= 0.1:
        raise ValueError(f"GC(%): seqkit {row['GC(%)']!r} vs 100*gc_fraction {gc:.6f}")


def overlap_seqfu(human_text: str, seqfu_text: str) -> None:
    human = parse_human(human_text)
    row = parse_tsv_row(seqfu_text)
    require_equal("reads/#Seq", human["reads"], row["#Seq"])
    require_equal("bases/Total bp", human["bases"], row["Total bp"])
    require_equal("min_length/Min", human["min_length"], row["Min"])
    require_equal("max_length/Max", human["max_length"], row["Max"])
    mean = as_float("mean_length", human["mean_length"])
    avg = as_float("Avg", row["Avg"])
    if abs(avg - mean) >= 0.01:
        raise ValueError(f"Avg: seqfu {row['Avg']!r} vs mean_length {human['mean_length']}")
    if "%GC" in row:
        gc = as_float("gc_fraction", human["gc_fraction"])
        seqfu_gc = as_float("%GC", row["%GC"])
        if abs(seqfu_gc - gc) >= 0.01:
            raise ValueError(f"%GC: seqfu {row['%GC']!r} vs gc_fraction {human['gc_fraction']}")


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", type=Path, nargs="?", help="plain or gzip FASTQ")
    ap.add_argument(
        "--agree-human",
        nargs=2,
        metavar=("LEFT", "RIGHT"),
        help="exit 0 if two human stats dumps match ignoring input:",
    )
    ap.add_argument(
        "--overlap-seqkit",
        nargs=2,
        metavar=("HUMAN", "SEQKIT"),
        help="exit 0 if SeqKit stats -a TSV overlaps human length/GC/Q20/Q30",
    )
    ap.add_argument(
        "--overlap-seqfu",
        nargs=2,
        metavar=("HUMAN", "SEQFU"),
        help="exit 0 if SeqFu stats TSV overlaps human length (and GC if present)",
    )
    args = ap.parse_args()
    try:
        if args.agree_human:
            agree_human(read_text(Path(args.agree_human[0])), read_text(Path(args.agree_human[1])))
            return
        if args.overlap_seqkit:
            overlap_seqkit(
                read_text(Path(args.overlap_seqkit[0])),
                read_text(Path(args.overlap_seqkit[1])),
            )
            return
        if args.overlap_seqfu:
            overlap_seqfu(
                read_text(Path(args.overlap_seqfu[0])),
                read_text(Path(args.overlap_seqfu[1])),
            )
            return
        if args.path is None:
            ap.error("path is required unless a compare flag is set")
        if not args.path.is_file():
            print(f"error: file not found: {args.path}", file=sys.stderr)
            raise SystemExit(1)
        sys.stdout.write(format_human(collect(args.path), str(args.path)))
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
