#!/usr/bin/env python3
"""Rough speed of z-fastq's kernel tiers on the machine it runs on.

Usage: smoke.py WORK_DIR [REV...]

Builds this checkout's ReleaseFast binary, a -Dcpu=native build, a build for each CPU level this machine runs
(x86-64-v3 and x86-64-v4 on x86-64 when the release binary accepts ZFASTQ_TIER=avx2 or avx512; neoverse_n1 on
Linux aarch64 with that level's features; generic on macOS aarch64), and the ReleaseFast binary of each REV that
builds here (oldest first, so a slowdown can be pinned to the commit that brought it), each with its own Zig
cache directory, since a local cache shared across worktrees can hand one worktree's binary to another with the
same build script, and lists each binary's hash. It generates fixed inputs of 60 to 330 MB (a block of reads repeated, so each run takes
tenths of a second to seconds), checks that every build and every kernel tier the binary
runs here (ZFASTQ_TIER) writes the same bytes, then times each one five times, interleaved, and prints median
wall times as a Markdown table (also appended to $GITHUB_STEP_SUMMARY). Exits 1 when outputs differ.
Shared CI runners make the numbers rough: only large differences mean anything, and no number is compared
across runs.
"""

import hashlib
import gzip
import os
import platform
import random
import statistics
import subprocess
import sys
import time

TIERS = ("baseline", "avx2", "avx512")
ROUNDS = 5
ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def run(argv, cwd=None):
    subprocess.run(argv, cwd=cwd, check=True)


def build(source, prefix, *flags):
    env = {key: value for key, value in os.environ.items() if key != "ZIG_LOCAL_CACHE_DIR"}
    try:
        subprocess.run(["zig", "build", "-Doptimize=ReleaseFast", "--prefix", prefix, "--cache-dir", prefix + "-cache", *flags],
                       cwd=source, env=env, check=True)
    except subprocess.CalledProcessError:
        return None
    return os.path.join(prefix, "bin", "z-fastq")


def file_digest(path):
    with open(path, "rb") as built:
        return hashlib.sha256(built.read()).hexdigest()[:12]


def cpu_levels(binary):
    machine = platform.machine().lower()
    if machine in ("x86_64", "amd64"):
        levels = [("x86_64_v3", "avx2"), ("x86_64_v4", "avx512")]
        return [level for level, tier in levels if invoke(binary, tier, ["--version"], subprocess.DEVNULL) == 0]
    if machine in ("aarch64", "arm64"):
        if platform.system() == "Darwin":
            return ["generic"]
        with open("/proc/cpuinfo") as info:
            features = set(next((line.split(":", 1)[1].split() for line in info if line.startswith("Features")), []))
        if {"asimddp", "atomics", "lrcpc", "fphp"} <= features:
            return ["neoverse_n1"]
    return []


def cpu_name():
    if platform.system() == "Darwin":
        return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True).stdout.strip()
    lscpu = subprocess.run(["lscpu"], capture_output=True, text=True).stdout
    for line in lscpu.splitlines():
        if line.startswith("Model name:"):
            return line.split(":", 1)[1].strip()
    return platform.machine()


def read_block(count, length, seed, mate=None):
    rng = random.Random(seed)
    records = []
    for index in range(count):
        sequence = "".join(rng.choice("ACGTACGTACGTACGN") for _ in range(length))
        if index % 3 == 0 and length > 60:
            adapter = ADAPTER2 if mate == 2 else ADAPTER1
            sequence = sequence[: length - len(adapter)] + adapter
        quality = "".join(chr(35 + rng.randrange(40)) for _ in range(length))
        name = f"r{index}/{mate}" if mate else f"r{index}"
        records.append(f"@{name}\n{sequence}\n+\n{quality}\n")
    return "".join(records).encode()


def write_repeated(path, block, repeats):
    with open(path, "wb") as out:
        for _ in range(repeats):
            out.write(block)


def make_inputs(work):
    paths = {name: os.path.join(work, name) for name in ("short.fastq", "long.fastq", "short.fastq.gz", "r1.fastq", "r2.fastq")}
    short = read_block(20000, 150, 1)
    write_repeated(paths["short.fastq"], short, 50)
    write_repeated(paths["long.fastq"], read_block(250, 8000, 2), 60)
    with gzip.open(paths["short.fastq.gz"], "wb", compresslevel=1) as packed:
        for _ in range(10):
            packed.write(short)
    write_repeated(paths["r1.fastq"], read_block(20000, 150, 3, mate=1), 15)
    write_repeated(paths["r2.fastq"], read_block(20000, 150, 4, mate=2), 15)
    return paths


def cases(paths, work):
    out1 = os.path.join(work, "out1.fastq")
    out2 = os.path.join(work, "out2.fastq")
    report = os.path.join(work, "report.json")
    return [
        ("count short", ["count", paths["short.fastq"]], []),
        ("stats short", ["stats", "--json", paths["short.fastq"]], []),
        ("stats long", ["stats", "--json", paths["long.fastq"]], []),
        ("check long", ["check", paths["long.fastq"]], []),
        ("check gzip", ["check", paths["short.fastq.gz"]], []),
        ("sample short", ["sample", "--fraction", "0.1", "--seed", "3", paths["short.fastq"]], []),
        ("interleave", ["interleave", paths["r1.fastq"], paths["r2.fastq"]], []),
        ("process", ["process", "--adapter1", ADAPTER1, "--adapter2", ADAPTER2, "--min-length", "15", "--pair-duplication",
                     "--report", report, "--out1", out1, "--out2", out2, paths["r1.fastq"], paths["r2.fastq"]], [out1, out2, report]),
    ]


def tier_environment(tier):
    env = {key: value for key, value in os.environ.items() if key != "ZFASTQ_TIER"}
    if tier:
        env["ZFASTQ_TIER"] = tier
    return env


def invoke(binary, tier, argv, stdout):
    return subprocess.run([binary, *argv], env=tier_environment(tier), stdout=stdout, stderr=subprocess.DEVNULL).returncode


def output_digest(binary, tier, argv, outputs):
    for path in outputs:
        if os.path.exists(path):
            os.remove(path)
    result = subprocess.run([binary, *argv], env=tier_environment(tier), capture_output=True)
    digest = hashlib.sha256(result.stdout)
    for path in outputs:
        if not os.path.exists(path):
            return result.returncode, "missing " + path
        with open(path, "rb") as produced:
            digest.update(produced.read())
    return result.returncode, digest.hexdigest()


def main():
    work = os.path.abspath(sys.argv[1])
    revisions = sys.argv[2:]
    os.makedirs(work, exist_ok=True)
    variants = {}
    binary = build(".", os.path.join(work, "new"))
    if binary is None:
        sys.exit("this checkout does not build")
    variants["release"] = (binary, None)
    for tier in TIERS:
        if invoke(binary, tier, ["--version"], subprocess.DEVNULL) == 0:
            variants[f"release, ZFASTQ_TIER={tier}"] = (binary, tier)
    for level in cpu_levels(binary):
        built = build(".", os.path.join(work, level), f"-Dcpu={level}")
        if built:
            variants[f"-Dcpu={level}"] = (built, None)
    native = build(".", os.path.join(work, "native"), "-Dcpu=native")
    if native:
        variants["-Dcpu=native"] = (native, None)
    for revision in revisions:
        source = os.path.join(work, f"source-{revision}")
        run(["git", "worktree", "add", "--detach", source, revision])
        built = build(source, os.path.join(work, f"build-{revision}"))
        if built:
            variants[revision[:7]] = (built, None)

    paths = make_inputs(work)
    names = list(variants)
    lines = [f"### {platform.system()} {platform.machine()}: {cpu_name()}", "",
             "Median wall time of 5 interleaved runs, in ms; the ratio is to `release`.", "",
             "Binaries (SHA-256): " + ", ".join(f"{name} `{file_digest(binary)}`" for name, (binary, tier) in variants.items() if tier is None), "",
             "| case | " + " | ".join(names) + " |", "|---" * (len(names) + 1) + "|"]
    mismatch = False
    for label, argv, outputs in cases(paths, work):
        digests = {name: output_digest(*variants[name], argv, outputs) for name in names}
        if len(set(digests.values())) != 1 or digests["release"][0] != 0:
            mismatch = True
            print(f"output differs or fails for {label}: {digests}", file=sys.stderr)
        times = {name: [] for name in names}
        for round_index in range(ROUNDS):
            for name in names[round_index % len(names):] + names[: round_index % len(names)]:
                for path in outputs:
                    if os.path.exists(path):
                        os.remove(path)
                start = time.perf_counter()
                invoke(*variants[name], argv, subprocess.DEVNULL)
                times[name].append(time.perf_counter() - start)
        reference = statistics.median(times["release"])
        cells = [f"{statistics.median(times[name]) * 1000:.0f} ({statistics.median(times[name]) / reference:.2f})" for name in names]
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    report = "\n".join(lines) + "\n"
    print(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as out:
            out.write(report + "\n")
    sys.exit(1 if mismatch else 0)


if __name__ == "__main__":
    main()
