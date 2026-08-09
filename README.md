# z-fastq

Fast Zig FASTQ I/O library and CLI.

## Build

Requires Zig 0.16.0 at `./zig-0.16.0/zig` (use the `./zig` wrapper).

```bash
./zig build
./zig build test
./zig build -Doptimize=ReleaseFast
```

## Usage

```bash
./zig-out/bin/z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
```

Each successfully parsed path or explicit `-` for standard input prints one record count on stdout. Standard input may appear once and is never selected implicitly. Non-zero exit status indicates an error; parse failures include record index, line number, and byte offset on stderr.

Exit status 1 reports invalid FASTQ, 2 reports command-line usage, 3 reports I/O or unexpected allocation failure, and 4 reports configured or arithmetic limits. Untrusted command, option, and path bytes are displayed using printable ASCII, doubled backslashes, and uppercase `\xHH` escapes for all other bytes.

## Library

Import the `z-fastq` module from `src/root.zig`. The current surface provides `Reader`, `Writer`, borrowed `Record`, `OwnedRecord`, structural diagnostics, `count_scan`, plain `io` adapters, shared `limits`, and `VERSION`.

Records returned by `Reader.next()` borrow reader storage until the next reader advance or deinitialization. Use `toOwned()` and later `OwnedRecord.deinit()` when a record must outlive that boundary. Byte-source and byte-sink wrappers are copied, but their referenced adapters must outlive the reader or writer.

## License

TBD
