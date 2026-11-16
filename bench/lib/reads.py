#!/usr/bin/env python3
"""Read-level checks of a trimmed pair, for tools with their own trimming rules.

Usage:
  reads.py hashes FASTQ OUT        sorted 64-bit fingerprints of each record (name, sequence, qualities)
  reads.py agree HASHES FASTQ      prints: records of FASTQ found in HASHES, records in HASHES, records in FASTQ
  reads.py valid IN1 IN2 OUT1 OUT2 prints accept, or wrong: and the first problem

A record's name is its header up to the first space. valid accepts a trimmed
pair when its records keep the input's pairs in the input's order, both mates
of a pair together, and each output read is a stretch of its input read with
the same quality bytes. It decides nothing about where a tool cuts.
FASTQ files may be plain or gzip; output is plain text on stdout.
"""

import array
import gzip
import hashlib
import itertools
import sys


def records(path):
    with open(path, "rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    stream = gzip.open(path, "rb") if gz else open(path, "rb", buffering=1 << 20)
    with stream:
        while True:
            header = stream.readline()
            if not header:
                return
            seq = stream.readline().rstrip(b"\r\n")
            stream.readline()
            qual = stream.readline().rstrip(b"\r\n")
            yield header[1:].split(None, 1)[0] if header[1:].strip() else b"", seq, qual


def key(name, seq, qual):
    return int.from_bytes(hashlib.blake2b(name + b"\n" + seq + b"\n" + qual, digest_size=8).digest(), "little")


def hashes(path, out):
    values = array.array("Q", sorted(key(*r) for r in records(path)))
    with open(out, "wb") as target:
        values.tofile(target)


def agree(ref_path, path):
    ref = array.array("Q")
    with open(ref_path, "rb") as source:
        ref.frombytes(source.read())
    got = sorted(key(*r) for r in records(path))
    i = j = found = 0
    while i < len(ref) and j < len(got):
        if ref[i] == got[j]:
            found += 1
            i += 1
            j += 1
        elif ref[i] < got[j]:
            i += 1
        else:
            j += 1
    print(found, len(ref), len(got))


def cut_of(source, kept):
    """Whether KEPT (name, seq, qual) is a stretch of SOURCE with the same qualities."""
    _, s_seq, s_qual = source
    _, k_seq, k_qual = kept
    if len(k_seq) != len(k_qual):
        return False
    start = s_seq.find(k_seq)
    while start >= 0:
        if s_qual[start:start + len(k_qual)] == k_qual:
            return True
        start = s_seq.find(k_seq, start + 1)
    return False


def valid(in1, in2, out1, out2):
    inputs = zip(records(in1), records(in2))
    outputs = itertools.zip_longest(records(out1), records(out2))
    pairs = 0
    for o1, o2 in outputs:
        pairs += 1
        if o1 is None or o2 is None:
            return "wrong: mates differ in number"
        for i1, i2 in inputs:
            if i1[0] == o1[0] or i2[0] == o2[0]:
                break
        else:
            return f"wrong: pair {pairs} is not an input pair in input order"
        if i1[0] != o1[0] or i2[0] != o2[0]:
            return f"wrong: pair {pairs} mates out of step"
        if not cut_of(i1, o1) or not cut_of(i2, o2):
            return f"wrong: pair {pairs} is not a cut of its input reads"
    return "accept"


def main():
    command, args = sys.argv[1], sys.argv[2:]
    if command == "hashes":
        hashes(*args)
    elif command == "agree":
        agree(*args)
    elif command == "valid":
        print(valid(*args))
    else:
        sys.exit(f"reads.py: unknown command {command}")


if __name__ == "__main__":
    main()
