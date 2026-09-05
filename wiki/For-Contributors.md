# For contributors

This page is for contributors who update public documentation alongside the source. Keep user-facing instructions on the owning page, and use this map to find the code and tests behind a documented behavior.

## Source map

```text
src/main.zig       CLI parsing, commands, output, and diagnostics
src/root.zig       public Zig exports and version
src/fastq.zig      records, parser, writer, and validation
src/count_scan.zig specialized count scanner
src/stats.zig      aggregate statistics
src/sample.zig     CLI-only sampling grammars and selectors
src/pair.zig       CLI-only paired-name parsing and matching
src/io.zig         byte sources, sinks, limits, and gzip dispatch
src/inflate.zig    native streaming DEFLATE engine
src/crc32.zig      portable and runtime-selected CRC implementation
tests/              parser, CLI, pair, sample, and output contracts
wiki/               public user documentation
```

## Before changing a documented behavior

1. Find the owning source path and its callers.
2. Read the focused tests that establish the current behavior.
3. Classify the claim as observed, inferred, proposed, or unknown.
4. Update the owning public page only when the user-facing contract changes.
5. Re-read the page and the diff before reporting completion.

The page should explain why a user needs the behavior, not just restate a function name.

## Verification

The normal project gates are:

```bash
zig fmt --check src tests build.zig
zig build test --summary all
zig build test -Dstatic=true -Doptimize=ReleaseSafe --summary all
zig build test -Dstatic=true -Doptimize=ReleaseFast --summary all
```

Do not use a timing result as proof of correctness. Performance evidence and behavior checks answer different questions.

## Documentation ownership

- The command pages own user-facing syntax and workflow examples.
- [Input and output](Input-and-Output) owns stream and format behavior.
- [Automation](Automation) owns exit classes and machine output guidance.
- [Limits and supported formats](Limits-and-Supported-Formats) owns current boundaries.
- [Paired reads](Paired-Reads) owns name normalization and pair errors.
- [Zig library](Zig-Library) owns the public module overview.

Link to the owning page instead of duplicating a detailed contract in every recipe or FAQ answer.
