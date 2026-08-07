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
./zig-out/bin/z-fastq count file.fastq [file2.fastq ...]
```

Each input file prints one record count on stdout (one decimal line per file). Non-zero exit status indicates an error; parse failures include record index, line number, and byte offset on stderr.

## Library

Import the `z-fastq` module from `src/root.zig`: `Reader`, `Writer`, `Record`, `count_scan`, and `io` adapters for plain files and slices.

## License

TBD
