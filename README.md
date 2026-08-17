# z-fastq

Fast Zig FASTQ I/O library and CLI.

## Build

Requires Zig 0.16.0 at `./zig-0.16.0/zig` (use the `./zig` wrapper).

```bash
./zig build
./zig build test
./zig build -Dstatic=true -Doptimize=ReleaseFast
```

Debug is the development build. ReleaseFast is always single-threaded, static, and stripped; `-Dstatic=true` states that required release configuration explicitly.

The supported Linux x86-64 build uses the vendored ISA-L 2.32.1 stateful inflate and CRC path by default. Building this path from source requires NASM 2.14.01 or newer. Use `-Disa-l=false` to exclude ISA-L, NASM, and its C linkage from the build. The Zig path runtime-selects PCLMUL CRC-32 on supported x86-64 processors and retains a portable fallback. Both paths stream through bounded storage and validate the same project-owned gzip framing, CRC-32, ISIZE, and concatenated-member behavior.

### Dependency and portability boundary

The project as a whole is not dependency-free. The accelerated build compiles vendored BSD-licensed ISA-L C and x86-64 assembly and needs NASM when built from source. The static release includes ISA-L and its C runtime support in the executable, so users do not install either separately. The ISA-L-disabled path needs no external compression library, C runtime linkage, or assembler beyond the Zig toolchain, but its PCLMUL CRC schedule retains MIT-licensed `crc32fast` provenance. No other build target is currently supported.

## Usage

```bash
./zig-out/bin/z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
./zig-out/bin/z-fastq stats [--json] [--max-line-bytes N] <path|-> [<path|-> ...]
./zig-out/bin/z-fastq check [--json] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|-> [<path|-> ...]
./zig-out/bin/z-fastq check --paired [--json] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <R1|-> <R2|->
./zig-out/bin/z-fastq check --interleaved [--json] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|->
./zig-out/bin/z-fastq sample --fraction P [--seed S] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|->
./zig-out/bin/z-fastq sample --count K [--seed S] [--alphabet iupac|acgtn] [--max-line-bytes N] <path>
./zig-out/bin/z-fastq interleave [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <R1|-> <R2|->
```

`count` prints one record count for each successfully parsed plain or gzip path, or explicit `-` for standard input. Input bytes select gzip independently of the path suffix. Standard input may appear once and is never selected implicitly. Non-zero exit status indicates an error; parse failures include record index, line number, and byte offset on stderr.

`stats` reports aggregate read lengths, case-insensitive base composition, GC fraction, and Phred+33 mean, Q20, and Q30 metrics. Undefined values print as `-`. Quality bytes outside ASCII 33 through 126 fail with S006 rather than entering the result.

`check` validates FASTQ structure, sequence symbols, and Phred+33 byte range. Its default alphabet accepts upper- and lowercase IUPAC nucleotide symbols; `--alphabet acgtn` selects the narrower A, C, G, T, and N policy. Successful validation is silent. The first failure in each input reports an S001 through S006 code with its zero-based record index and byte offset plus its one-based record line.

`check --paired` validates two inputs in lock step, and `check --interleaved` validates consecutive records as R1 and R2. The default Illumina name policy recognizes modern second-token mate fields, legacy first-token `/1` and `/2`, and terminal mate suffixes on the second token. Equal unmarked first tokens also pass. `--pair-names exact` instead requires complete first-token equality. P001 reports name or mate-direction disagreement; P002 reports unequal two-file counts or an odd interleaved record. Paired mode accepts standard input on at most one side.

`sample --fraction` selects complete records independently with probability P and writes plain FASTQ with LF endings. P uses the exact decimal grammar documented by `--help`, and the seed is an unsigned 64-bit decimal value that defaults to 11. Selected indexes and output bytes match `seqtk sample` on the screened subset with LF endings, nonempty records, bare plus lines, printable ASCII headers, a signed 64-bit seed, and a P that seqtk also interprets as a fraction. Every record is validated before selection, so fraction 0 cannot hide invalid input. Values that parse to fraction 0 or fraction 1 avoid random-number work.

`sample --count` selects exactly `min(K, N)` records from a regular plain or gzip file containing N valid records. It uses a bounded index reservoir compatible with `seqtk sample -2`, then reopens and validates the file before writing selected records in input order. K is an unsigned 64-bit decimal value. Standard input is rejected because exact mode requires two passes. The command compares file identity, size, modification time, and the two record counts as a best-effort change guard. Concurrent modification that preserves those signals is unsupported. Output is not transactional and may contain a valid prefix if the second pass fails.

`interleave` reads R1 and R2 in lock step, validates both records and their paired names, then writes R1 followed by R2 as plain LF-terminated FASTQ. Plain and gzip inputs may be mixed, and standard input may appear on at most one side. A validation or input failure can leave complete earlier pairs on stdout, while an output failure may leave an arbitrary byte prefix because stdout is not transactional.

`stats --json` and `check --json` emit one provisional versioned JSON document with results in input order. The schema identifiers are `z-fastq/stats-v1` and `z-fastq/check-v1`. The top-level fields are `schema`, `tool`, `byte_strings`, and `results`. Ordinary results contain `input`; two-file paired results contain the ordered `inputs` array. Successful stats results add the human-mode aggregate fields, successful check results add no payload, and failed results add an error. Single-record S errors retain nullable record locations. P001 contains both record positions, 128-byte bounded first-token and normalized-ID prefixes, original lengths, truncation flags, and recognized mate-marker arrays. P002 contains the remaining side and last known record indexes. Counters remain JSON integers in the unsigned 64-bit range from 0 through 18,446,744,073,709,551,615, while undefined numeric values are null. Consumers backed only by IEEE-754 doubles cannot preserve every integer above 9,007,199,254,740,991 and need arbitrary-precision number handling when such values are possible. Handled input and validation failures stay inside the result array, while invalid command usage remains on stderr. Input labels and identifier prefixes use the reversible `escaped-bytes-v1` representation rather than assuming UTF-8.

Exit status 1 reports invalid FASTQ, 2 reports command-line usage, 3 reports I/O or unexpected allocation failure, and 4 reports configured or arithmetic limits. Untrusted command, option, and path bytes are displayed using printable ASCII, doubled backslashes, and uppercase `\xHH` escapes for all other bytes.

## Library

Import the `z-fastq` module from `src/root.zig`. The current surface provides `Reader`, `Writer`, borrowed `Record`, `OwnedRecord`, structural and semantic diagnostics, allocation-free record validation with `validateRecord()`, `count_scan`, the checked allocation-free `Stats` accumulator, plain and gzip `io` adapters, shared `limits`, and `VERSION`.

Records returned by `Reader.next()` borrow reader storage until the next reader advance or deinitialization. Use `toOwned()` and later `OwnedRecord.deinit()` when a record must outlive that boundary. Byte-source and byte-sink wrappers are copied, but their referenced adapters must outlive the reader or writer. `io.gzip.ReaderSource` validates RFC 1952 headers, DEFLATE payloads, trailers, and concatenated members while borrowing a `std.Io.Reader` with at least ten buffer bytes.

## License

The project license is not yet selected. The vendored ISA-L subset retains its [BSD-3-Clause license](vendor/ISA-L/LICENSE.md). The native PCLMUL CRC schedule is adapted from `crc32fast` 1.5.0 and retains its [MIT license](vendor/CRC32FAST-LICENSE-MIT). Both notices must accompany binary distributions that contain the corresponding code.
