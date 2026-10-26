# Zig library

The CLI is the primary z-fastq product. The repository also registers a small Zig module for applications that want the parser and data-path building blocks directly.

## Import the module

The package module is named `z-fastq`. Once the repository is present as a Zig dependency, import it from application code with:

```zig
const zfastq = @import("z-fastq");
```

The repository's `build.zig` registers the module from `src/root.zig`.

## Public building blocks

The root module exports:

- `Record`, `OwnedRecord`, and `toOwned`.
- `Reader`, `RecordOffsets`, `Writer`, and their error types, including `ParseError`.
- `LintCode`, `codeTag`, `Alphabet`, and `ValidationOptions`.
- `SemanticField`, `SemanticError`, and `validateRecord`.
- `Stats`, `StatsResult`, `StatsError`, `QualityError`, and `decodePhred33`.
- `io.ByteSource` and `io.ByteSink`.
- Plain slice, reader, file, writer, and sink adapters.
- A gzip reader-source adapter.
- The allocation-free `count_scan` module.
- `limits` (also reachable as `io.limits`) with the default line limit and buffer sizes.
- `VERSION`, the package version string.

The CLI's option parsing, pair-name orchestration, sampling selectors, JSON serialization, and deinterleave file policy remain executable code. They are not all part of the root library API.

## Borrowed and owned records

`Reader` returns a borrowed `Record`. Its fields point into reader-owned storage and remain valid until the reader advances or is deinitialized. Use `toOwned` when a record must survive that lifetime:

```zig
const record = (try reader.next()) orelse return;
var owned = try zfastq.toOwned(allocator, record);
defer owned.deinit();
```

This distinction is central to the module. Borrowed records keep streaming work small. Owned records are the explicit escape hatch.

A `Record`'s `id` is the first space- or tab-delimited token of `header` and points into the same bytes. `OwnedRecord.deinit` releases all fields of a copy made by `toOwned`.

## Writer behavior

`Writer` emits four-line FASTQ records with LF endings. It rejects unequal sequence and quality lengths, invalid header or field bytes, embedded line breaks, and terminal carriage returns before writing the record. The sink adapter it writes to must outlive it, and output is complete only after `flush`.

## Validation and statistics

`validateRecord` checks only the sequence alphabet and Phred+33 quality bytes, without allocating. Record structure and equal field lengths are the reader's checks. `decodePhred33` rejects bytes outside ASCII 33 (`!`) through 126 (`~`).

`Stats` allocates nothing. `addRecord` validates a record and changes no counter when it rejects it; `takeLastQualityError` returns the details of the last S006 failure. In `StatsResult`, the length minimum, maximum, and mean are null when no read was added, and a ratio is null when its denominator is zero.

## Reader and I/O behavior

The reader consumes a `ByteSource`, exposes record offsets, and enforces the configured line limit. The caller opens and closes paths and streams; the I/O adapters only borrow them, so the FASTQ record types never handle files. The gzip adapter validates gzip framing while reading through its backing source.

`ByteSource` and `ByteSink` are copied by value and point at their adapter, so the adapter must stay at a stable address and outlive every copy. A read fills a prefix of the destination and returns its length; zero for a nonempty destination means EOF, and a length beyond the destination is rejected. A write consumes the whole slice; a sink without a flush callback treats `flush` as a no-op. `FileSink` writes whole pages and keeps a partial last page buffered until the next write or `flush`.

An error ends the reader. Read its diagnostic with `takeLastError` and do not call `next` again: the position and record index may reflect a partly read record. The same holds for `count_scan.Scanner`: after an error from `feed` or `finishEof`, read `takeLastError` and stop; `record_index` and `byte_offset` are not defined after an error.

## Compatibility

The `0.1.0` module targets Zig `0.16.0` and Linux x86-64. Pin the package version when an application depends on the exact exported names or error shapes, and check release notes before upgrading.

