<!-- markdownlint-disable MD033 MD041 -->

<h1 align="center">Z-FASTQ</h1>

<p align="center">
  <strong>Fast, bounded-memory FASTQ work in one Zig executable.</strong>
</p>

<p align="center">
  <a href="https://github.com/eneskemalergin/z-fastq/releases"><img src="https://img.shields.io/badge/version-v0.0.18-2563eb?style=flat-square" alt="Version v0.0.18"></a>
  <a href="https://ziglang.org/download/"><img src="https://img.shields.io/badge/Zig-0.16.0-F7A41D?style=flat-square&amp;logo=zig&amp;logoColor=white" alt="Zig 0.16.0"></a>
  <img src="https://img.shields.io/badge/platform-Linux%20x86--64-64748b?style=flat-square" alt="Supported platform: Linux x86-64">
  <img src="https://img.shields.io/badge/license-not%20selected-94a3b8?style=flat-square" alt="License not selected yet">
</p>

---

z-fastq is a single-core command-line toolkit for working with plain and gzip-compressed FASTQ. It is designed around predictable streaming, bounded memory use, explicit validation, and useful behavior on real sequencing data.

## What it does

- Count FASTQ records from files or standard input.
- Calculate read-length, base-composition, GC, and quality statistics.
- Validate structure, sequence alphabets, quality bytes, and paired-read names.
- Sample records by fraction or exact count, including paired and interleaved input.
- Interleave and deinterleave paired FASTQ.
- Read plain FASTQ and gzip input through the same streaming interface.
- Emit machine-readable JSON for selected statistics and validation workflows.

The CLI is the primary product. A small Zig module is also exported for applications that need the reader, writer, validation, statistics, and I/O building blocks directly.

## Why it exists

I am building z-fastq around a simple constraint: bioinformatics tools should remain practical on ordinary hardware. That means a small executable, bounded streaming memory, and a predictable cost per process.

The current CLI is single-threaded by design. That keeps one invocation easy to reason about and leaves multi-file scheduling to the workflow layer.

## A few honest boundaries

The supported target is currently **Linux x86-64**. Native Windows and other targets are not supported yet.

The accelerated build uses vendored ISA-L for gzip and CRC work and requires NASM when built from source. The ISA-L-disabled build uses the Zig implementation and retains a portable CRC fallback:

```bash
zig build -Disa-l=false
```

The release path is intended to be static. More detailed format guarantees, limits, error codes, machine-readable output, and compatibility notes belong in the project documentation rather than this overview.

## Start

Build with [Zig 0.16.0](https://ziglang.org/download/):

```bash
zig build
zig build test
zig build -Dstatic=true -Doptimize=ReleaseSafe
zig build -Dstatic=true -Doptimize=ReleaseFast
```

Try the main workflows:

```bash
./zig-out/bin/z-fastq count reads.fastq.gz
./zig-out/bin/z-fastq stats reads.fastq.gz
./zig-out/bin/z-fastq check --alphabet iupac reads.fastq.gz
./zig-out/bin/z-fastq sample --fraction 0.10 --seed 11 reads.fastq.gz > sample.fastq
./zig-out/bin/z-fastq interleave reads_R1.fastq.gz reads_R2.fastq.gz > interleaved.fastq
```

Use `z-fastq --help` for the complete command and option reference.

## Commands

| Command        | Purpose                                                    |
| -------------- | ---------------------------------------------------------- |
| `count`        | Count successfully parsed records.                         |
| `stats`        | Report lengths, composition, GC, and quality metrics.      |
| `check`        | Validate FASTQ structure, symbols, qualities, and pairing. |
| `sample`       | Select records or read pairs by fraction or exact count.   |
| `interleave`   | Combine R1 and R2 into a validated interleaved stream.     |
| `deinterleave` | Split validated interleaved reads into two output files.   |

All commands report errors with non-zero exit status. Validation and parsing failures identify the affected record and source location where available.

## Documentation

The README stays intentionally short. The detailed reference covers:

- command options and examples;
- FASTQ and gzip behavior;
- paired-read name policies;
- sampling compatibility and reproducibility;
- machine-readable output and exit statuses;
- resource limits and portability;
- the Zig module API.

See the [project Wiki](https://github.com/eneskemalergin/z-fastq/wiki) for the evolving user documentation.

## License

The project license is still being selected. The vendored ISA-L subset retains its [BSD-3-Clause license](vendor/ISA-L/LICENSE.md), which must accompany binary distributions that contain that code.

---

<p align="center"><em>Four lines hold a life,<br>
Each base travels through the night,<br>
Reads emerge as light.</em></p>
