<!-- markdownlint-disable MD033 MD041 -->

<h1 align="center">Z-FASTQ</h1>

<p align="center">
  <strong>Fast, bounded-memory FASTQ work in one Zig executable.</strong>
</p>

<p align="center">
  <a href="https://github.com/eneskemalergin/z-fastq/releases"><img src="https://img.shields.io/badge/version-v0.0.15-2563eb?style=flat-square" alt="Version v0.0.15"></a>
  <a href="https://github.com/eneskemalergin/z-fastq"><img src="https://img.shields.io/badge/status-polishing-eab308?style=flat-square" alt="Status: polishing"></a>
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

The current CLI is single-threaded by design. That keeps one invocation easy to reason about and leaves multi-file scheduling to the workflow layer. My ambition is to make z-fastq one of the fastest general-purpose FASTQ toolkits without requiring server-class hardware. The long-term goal is to process several files at once without hiding a second layer of threads or multiplying memory use inside every worker.

## A few honest boundaries

The supported release target is currently **Linux x86-64**. Native Windows and other targets are not supported yet.

The accelerated build uses vendored ISA-L for gzip and CRC work and requires NASM when built from source. The ISA-L-disabled build uses the Zig implementation and retains a portable CRC fallback:

```bash
zig build -Disa-l=false
```

The release path is intended to be static. More detailed format guarantees, limits, error codes, JSON schemas, and compatibility notes belong in the project documentation rather than this overview.

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

The README stays intentionally short. The detailed reference will cover:

- command options and examples;
- FASTQ and gzip behavior;
- paired-read name policies;
- sampling compatibility and reproducibility;
- JSON schemas and exit statuses;
- resource limits and portability;
- benchmark methods and results;
- the Zig module API.

See the [project Wiki](https://github.com/eneskemalergin/z-fastq/wiki) for the evolving documentation. Development notes and local benchmark tooling are kept in [`plan/`](plan/) and [`tools/`](tools/).

## Performance

Performance work is central to the project. The goal is not only to finish a FASTQ job quickly, but to do it with a small enough process footprint that several jobs can share a consumer machine.

### Early signal, not a final benchmark

The project is still in polishing. Internals, defaults, and hot paths change often enough that final benchmark tables and public methodology are still being developed. I prefer to share the direction without presenting unfinished numbers as a permanent leaderboard.

Early local comparisons suggest that z-fastq is often **about 1.5x to 5x faster** than selected peers, while using **roughly 2x to 5x less peak RSS** in comparable runs. The latest audited static ReleaseFast binary is under 1 MB, at about 744 kB.

Those ranges vary with the command, plain versus gzip input, read length, pairing mode, sampling mode, and comparison tool. Sampling and output-heavy workflows are still being tuned, so these figures are signals of direction rather than a universal ranking. The aim is a fast, small process that can be scheduled across several files without every worker consuming more of the machine than necessary.

The peer set is deliberately mixed. It includes familiar baseline tools such as `seqtk`, `fqtools`, `BBTools`, `SeqFu`, and `FastQValidator`, alongside newer native implementations, including Rust-based tools such as `Needletail`, `Helicase`, `fq`, `Fasten`, `fqkit`, and `Rasusa`. The familiar tools anchor common practice. The newer native tools are closer to z-fastq's low-level design and give a more useful comparison for throughput and memory. I only compare overlapping work, and keep differences in validation rules, RNGs, output behavior, and child-process accounting visible.

For perspective, a few large-input scale probes looked like this:

| Workload         | Input on disk | FASTQ bytes processed | Observed wall time |  Peak RSS |
| ---------------- | ------------: | --------------------: | -----------------: | --------: |
| `count`          |  1.4 GB plain |          1.4 GB plain |             212 ms |   896 KiB |
| `count`          |   793 MB gzip |       6.14 GB decoded |             3.55 s |   896 KiB |
| `stats`          |  6.1 GB plain |          6.1 GB plain |             2.26 s |   896 KiB |
| `stats`          |   793 MB gzip |       6.14 GB decoded |             4.62 s |   896 KiB |
| `interleave`     |   1.0 GB gzip |       3.58 GB decoded |             5.72 s | 1,048 KiB |
| `check --paired` |   1.0 GB gzip |       3.58 GB decoded |             5.50 s | 1,336 KiB |

Peak RSS here means resident memory measured for the z-fastq process itself. These single-run scale probes are useful for showing the shape of the system, not for making a final cross-platform promise.

Final benchmark documentation will include the measurement method, peer scope, correctness checks, and raw reports when the implementation settles.

## License

The project license is still being selected. The vendored ISA-L subset retains its [BSD-3-Clause license](vendor/ISA-L/LICENSE.md), which must accompany binary distributions that contain that code.

---

<p align="center"><em>Four lines hold a life,<br>
Each base travels through the night,<br>
Reads emerge as light.</em></p>
