<!-- markdownlint-disable MD007 MD024 -->

# Changelog

All notable changes to z-fastq will be documented in this file.

## Unreleased

### Changed

- `sample`, `interleave`, and `deinterleave` take complete records straight from the read buffer and write output in whole pages, which roughly halves their time on plain input without raising peak memory beyond a few KiB.
- `FileSink` writes whole pages and keeps a partial last page buffered until the next write or flush.

### Added

- `FileSink.initStreaming` writes at the file handle's current position, for standard output and pipes.

## [0.1.0]

The first public release. z-fastq is a single-threaded command-line toolkit for plain and gzip-compressed FASTQ on Linux x86-64, built around bounded streaming memory, strict validation, reproducible selection, and a small Zig module for applications that need the reader, writer, validation, and statistics building blocks. FASTQ-producing commands write plain four-line records with LF endings; the CLI does not write gzip output.

### Added

- **Commands.**
  - `count` reports record counts for one or more plain or gzip inputs, or standard input.
  - `stats` reports reads, bases, length range and mean, A/C/G/T/N and other base counts, GC fraction, quality sum and mean, and Q20/Q30 counts and fractions. Undefined values print as `-` in text and `null` in JSON.
  - `check` validates four-line structure, sequence symbols under the IUPAC or ACGTN policy, and Phred+33 quality bytes. `--paired` compares two mate files and `--interleaved` compares adjacent records.
  - `sample` selects records by seeded fraction or exact count. Fraction sampling is one pass; exact-count sampling is two passes and requires regular file paths.
  - `interleave` writes validated R1 then R2 records to standard output. `deinterleave` writes validated alternating records to two exclusively created output files.
- **Input and streaming.**
  - Plain FASTQ, RFC 1952 gzip input including concatenated members, and explicit standard input with `-`. gzip integrity is checked while decoding.
  - Every command streams with a 256 KiB read buffer and memory bounded independently of total input size, except the exact-count selector, whose state grows with the requested count.
  - The default logical line limit is 16 MiB; `--max-line-bytes N` overrides it per command.
  - `count`, `stats`, and `check` accept multiple independent inputs, process them in argument order, continue after an input failure, and return the highest exit class.
- **Paired reads.**
  - Lock-step two-file and interleaved validation, sampling, and transforms use the pair as the unit.
  - The default `illumina` policy recognizes terminal `/1` and `/2` markers and leading second-token mate fields. `--pair-names exact` compares the complete first token.
  - Unequal counts, mismatched names, reversed mate annotations, and an odd interleaved record are reported before a transform can emit the affected pair.
- **Diagnostics.**
  - Stable validation codes: `S001` through `S006` for structure, alphabet, and quality, and `P001` and `P002` for pair names and counts.
  - Diagnostics carry a zero-based record index, a zero-based byte offset in the decompressed FASTQ stream, and a one-based line number in the record when known.
  - Exit classes are `0` success, `1` validation failure, `2` usage error, `3` I/O or allocation failure, and `4` configured limit or arithmetic failure.
- **Machine-readable output.** `stats --json` and `check --json` emit the `z-fastq/stats-v1` and `z-fastq/check-v1` documents. Each result appears in argument order, with handled errors inside the document; the process exit class remains separate. The two v1 documents are frozen for the 0.1.0 interface.
- **Zig module.** The `z-fastq` module exports `Record` and `OwnedRecord`, `Reader` and `Writer`, validation types and `validateRecord`, `Stats`, `count_scan`, plain and gzip I/O adapters, and the documented limits.
- **Build.** The supported build is Linux x86-64 with Zig 0.16.0. ReleaseSafe and ReleaseFast builds are static. gzip decoding uses Zipir 0.2.1, pinned by package hash; the build needs no C compiler, NASM, or system zlib.
