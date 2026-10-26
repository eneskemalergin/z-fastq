<!-- markdownlint-disable MD007 MD024 -->

# Changelog

All notable changes to z-fastq will be documented in this file.

## [0.2.0]

A `process` command that validates paired reads into two plain or gzip outputs with a JSON report, faster record handling in every command, AVX2 kernels selected at run time on CPUs that have them, bounded memory for `stats` and paired sampling on long reads, and stricter handling of standard input and output. Output bytes, JSON documents, and exit classes are unchanged except where a fix below says so.

### Added

- **Commands.** `process` validates two mate files pair by pair, as `check --paired` does, and writes every accepted pair to two new files, plain or gzip (`--gzip`). `--max-pairs N` stops after the first N pairs, and `--report` writes a `z-fastq/process-v1` JSON document with per-mate statistics and any failure. The same inputs and options give byte-identical outputs and report.
- **Library.**
  - `FileSink.initStreaming` writes at the file handle's current position, for standard output and pipes.
  - `Stats.ratios` returns each ratio of `StatsResult` as an exact numerator and denominator, so callers can round without floating point.
- **Build.**
  - `-Dstrip` sets whether debug information is stripped; it defaults to on for ReleaseFast, so a profiling build can keep its symbols.
  - `zig build test-<suite>` runs a single test binary, such as `test-reader`, or `test-src` for the in-file tests; `zig build -l` lists them.

### Changed

- **Speed.**
  - `sample`, `interleave`, and `deinterleave` take complete records straight from the read buffer and write output in whole pages, which roughly halves their time on plain input without raising peak memory beyond a few KiB.
  - Interleaved fraction sampling, exact-count sampling in every mode, and `check --paired` and `--interleaved` take records the same way, and the exact-count output pass skips unselected records with less work: 5 to 20% less time on plain input.
  - On CPUs with AVX2, builds whose target lacks AVX2 select 32-byte kernels at run time for newline search past the first 1 KiB of a line, validation of sequences of 1 KiB or more (about a third less time for `check` on PacBio reads, 15 to 20% for `check` and `sample` on ONT reads), validation of lines longer than the read buffer, and `stats` counting on reads of 1 KiB or more (about a quarter less time on long reads).
  - Fraction sampling validates the records it does not select without copying them: 5 to 25% less time for `sample` and single-file `check` on long reads, with no change in memory.
- **Memory.**
  - `stats` streams a record larger than the read buffer instead of copying it, so its peak memory no longer grows with read length (772 KiB on every input shape, from up to 1.8 MiB on 300 kb reads), long reads take 5 to 25% less time, and `stats` allocates nothing per record.
  - Paired and interleaved fraction sampling draw each pair before reading it and validate unselected pairs without keeping them: on 300 kb read pairs, peak memory for `sample --paired` falls from 2.3 MiB to 0.8 MiB when few pairs are selected. `check --interleaved` and `deinterleave` read a mate whose partner already failed only far enough to report it.
- **Library.** `FileSink` writes whole pages and keeps a partial last page buffered until the next write or `flush`.
- **Build.** `-Dkernel-backend=portable` builds no AVX2 object and selects no kernel at run time; z-fastq's own kernels then use only the instructions of the build target, as Zipir's do.

### Fixed

- **Standard input and output.**
  - `count`, `stats`, and `check` reject standard output appended to one of their inputs (`count a.fastq b.fastq >> b.fastq`) with status `2` before reading anything, as `sample` and `interleave` already did; the input is no longer changed.
  - Standard input named twice, as `- /dev/stdin` or `/dev/stdin /dev/stdin`, is rejected when it is a pipe or FIFO instead of being read once and reported as zero records for the second name.
  - A pipe closed by its reader ends every command with status `3` and `error: standard output: I/O error`, as documented, instead of a SIGPIPE status of 141; `sample`, `stats`, `count`, and the JSON forms print that message on any standard-output failure.
- **Inputs and outputs.**
  - A directory given as an input is reported as `input is a directory` before anything is read, instead of `I/O error`; `deinterleave` then creates no output files.
  - `deinterleave` removes the empty `--out1` it created when `--out2` cannot be created, so a rerun is not blocked by it.
  - Exact-count sampling of a FIFO or directory exits with status `2`, as for standard input, instead of `3`.
- **Options and limits.**
  - A value option given twice (`--count 1 --count 2`) fails with status `2` and `may appear only once`, as `--out1` and `--out2` already did, instead of keeping the last value.
  - `sample --interleaved` and `deinterleave` accept every `--max-line-bytes` value the other commands accept, instead of exiting with status `4` before reading for values of 2^62 - 1 and above, sometimes with no message.
  - `stats` reports an input location past the supported limit as `input location exceeds supported limit`, as the other commands do.

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
