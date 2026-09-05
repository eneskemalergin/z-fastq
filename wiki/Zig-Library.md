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

- `Record` and `OwnedRecord`.
- `Reader`, `RecordOffsets`, `Writer`, and their error types.
- `LintCode`, `codeTag`, `Alphabet`, and `ValidationOptions`.
- `SemanticField`, `SemanticError`, and `validateRecord`.
- `Stats`, `StatsResult`, `StatsError`, `QualityError`, and `decodePhred33`.
- `io.ByteSource` and `io.ByteSink`.
- Plain slice, reader, file, writer, and sink adapters.
- A gzip reader-source adapter.
- The allocation-free `count_scan` module.

The CLI's option parsing, pair-name orchestration, sampling selectors, JSON serialization, and deinterleave file policy remain executable code. They are not all part of the root library API.

## Borrowed and owned records

`Reader` returns a borrowed `Record`. Its fields point into reader-owned storage and remain valid until the reader advances or is deinitialized. Use `toOwned` when a record must survive that lifetime:

```zig
const record = (try reader.next()) orelse return;
var owned = try zfastq.toOwned(allocator, record);
defer owned.deinit();
```

This distinction is central to the module. Borrowed records keep streaming work small. Owned records are the explicit escape hatch.

## Writer behavior

`Writer` emits four-line FASTQ records with LF endings. It rejects unequal sequence and quality lengths, invalid header or field bytes, embedded line breaks, and terminal carriage returns before writing the record.

## Reader and I/O behavior

The reader consumes a `ByteSource`, exposes record offsets, and enforces the configured line limit. The I/O adapters keep path and stream ownership outside the FASTQ record types. The gzip adapter validates gzip framing while reading through its backing source.

## Compatibility

The current module targets Zig `0.16.0` and Linux x86-64. The API is pre-`0.1.0`. Pin the project revision when an application depends on the exact exported names or error shapes.

