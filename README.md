# z-fastq

Fast Zig FASTQ I/O library and CLI.

## Build

Requires Zig 0.16.0 at `./zig-0.16.0/zig` (use the `./zig` wrapper).

```bash
./zig build
./zig build test
./zig build -Dstatic=true -Doptimize=ReleaseFast
```

Linux x86-64 builds use the vendored ISA-L 2.32.1 stateful inflate and CRC path by default. Building this path from source requires NASM 2.14.01 or newer. Use `-Disa-l=false` to exclude ISA-L, NASM, and its C linkage from the build; other targets select the Zig path automatically. The Zig path runtime-selects PCLMUL CRC-32 on supported x86-64 processors and retains a portable fallback. Both paths stream through bounded storage and validate the same project-owned gzip framing, CRC-32, ISIZE, and concatenated-member behavior.

### Dependency and portability boundary

The project as a whole is not dependency-free. The accelerated build compiles vendored BSD-licensed ISA-L C and x86-64 assembly and needs NASM when built from source. The static release includes ISA-L and its C runtime support in the executable, so users do not install either separately. The ISA-L-disabled path needs no external compression library, C runtime linkage, or assembler beyond the Zig toolchain, but its PCLMUL CRC schedule retains MIT-licensed `crc32fast` provenance. No equivalent ISA-L acceleration is currently integrated or verified for AArch64, RISC-V, Windows, or macOS.

## Usage

```bash
./zig-out/bin/z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
./zig-out/bin/z-fastq stats [--json] [--max-line-bytes N] <path|-> [<path|-> ...]
./zig-out/bin/z-fastq check [--json] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|-> [<path|-> ...]
```

`count` prints one record count for each successfully parsed plain or gzip path, or explicit `-` for standard input. Input bytes select gzip independently of the path suffix. Standard input may appear once and is never selected implicitly. Non-zero exit status indicates an error; parse failures include record index, line number, and byte offset on stderr.

`stats` reports aggregate read lengths, case-insensitive base composition, GC fraction, and Phred+33 mean, Q20, and Q30 metrics. Undefined values print as `-`. Quality bytes outside ASCII 33 through 126 fail with S006 rather than entering the result.

`check` validates FASTQ structure, sequence symbols, and Phred+33 byte range. Its default alphabet accepts upper- and lowercase IUPAC nucleotide symbols; `--alphabet acgtn` selects the narrower A, C, G, T, and N policy. Successful validation is silent. The first failure in each input reports an S001 through S006 code with its zero-based record index and byte offset plus its one-based record line.

`stats --json` and `check --json` emit one provisional versioned JSON document with results in input order. The schema identifiers are `z-fastq/stats-v1` and `z-fastq/check-v1`. The top-level fields are `schema`, `tool`, `byte_strings`, and `results`. Each result contains `input` and an `ok` or `error` status. Successful stats results add the human-mode aggregate fields, successful check results add no payload, and failed results add an error code, message, and nullable record location. Counters remain JSON integers in the unsigned 64-bit range from 0 through 18,446,744,073,709,551,615, while undefined numeric values are null. Consumers backed only by IEEE-754 doubles cannot preserve every integer above 9,007,199,254,740,991 and need arbitrary-precision number handling when such values are possible. Handled input and validation failures stay inside the result array, while invalid command usage remains on stderr. Input labels use the reversible `escaped-bytes-v1` representation rather than assuming UTF-8 paths.

Exit status 1 reports invalid FASTQ, 2 reports command-line usage, 3 reports I/O or unexpected allocation failure, and 4 reports configured or arithmetic limits. Untrusted command, option, and path bytes are displayed using printable ASCII, doubled backslashes, and uppercase `\xHH` escapes for all other bytes.

## Library

Import the `z-fastq` module from `src/root.zig`. The current surface provides `Reader`, `Writer`, borrowed `Record`, `OwnedRecord`, structural and semantic diagnostics, allocation-free record validation with `validateRecord()`, `count_scan`, the checked allocation-free `Stats` accumulator, plain and gzip `io` adapters, shared `limits`, and `VERSION`.

Records returned by `Reader.next()` borrow reader storage until the next reader advance or deinitialization. Use `toOwned()` and later `OwnedRecord.deinit()` when a record must outlive that boundary. Byte-source and byte-sink wrappers are copied, but their referenced adapters must outlive the reader or writer. `io.gzip.ReaderSource` validates RFC 1952 headers, DEFLATE payloads, trailers, and concatenated members while borrowing a `std.Io.Reader` with at least ten buffer bytes.

## License

The project license is not yet selected. The vendored ISA-L subset retains its [BSD-3-Clause license](vendor/ISA-L/LICENSE.md). The native PCLMUL CRC schedule is adapted from `crc32fast` 1.5.0 and retains its [MIT license](vendor/CRC32FAST-LICENSE-MIT). Both notices must accompany binary distributions that contain the corresponding code.
