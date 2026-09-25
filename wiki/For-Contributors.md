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

Run the required checks from the repository root with Zig 0.16.0 on Linux x86-64. The default build uses ISA-L and needs NASM; `-Disa-l=false` uses native gzip without NASM. Both engines must pass Debug, then static ReleaseSafe, then static ReleaseFast:

```bash
zig fmt --check src tests build.zig
zig build test --summary all
zig build test -Disa-l=false --summary all
zig build test -Dstatic=true -Doptimize=ReleaseSafe --summary all
zig build test -Disa-l=false -Dstatic=true -Doptimize=ReleaseSafe --summary all
zig build test -Dstatic=true -Doptimize=ReleaseFast --summary all
zig build test -Disa-l=false -Dstatic=true -Doptimize=ReleaseFast --summary all
```

Each test command installs its selected CLI before running the command tests. Run variants sequentially in one checkout because the CLI tests share `zig-out/bin/z-fastq`; parallel variants need separate checkouts. The checks use inline and tracked synthetic fixtures without corpus downloads or peer tools. CLI tests need `pidfd_open`, `/proc`, FIFOs, and hard and symbolic links; a missing facility is a setup failure.

Optional CPU instructions and timing or RSS thresholds are not test prerequisites. `-Disa-l=false` selects a gzip engine; native CRC can still use runtime PCLMUL. Record skipped test names and reasons:

- `[edge] - [single fraction sample]: fraction zero avoids allocation in the default backend` skips with native gzip because it checks an ISA-L-only allocation property. Zero-selection CLI tests still run with both engines.
- `[property] - [CRC-32 PCLMUL]: matches portable folding boundaries` skips if the compiler cannot emit its PCLMUL code or the running CPU lacks PCLMUL. Portable CRC and dispatch tests still run. This skip can occur with either installed engine because the fastq unit tests always use native gzip.

Investigate other skips. A skipped accelerated test does not prove that path ran, and a generic release build on one CPU does not prove fallback on older CPUs. Do not use a timing result as proof of correct output.

## Documentation ownership

- The command pages own user-facing syntax and workflow examples.
- [Input and output](Input-and-Output) owns stream and format behavior.
- [Automation](Automation) owns exit classes and machine output guidance.
- [Limits and supported formats](Limits-and-Supported-Formats) owns current boundaries.
- [Paired reads](Paired-Reads) owns name normalization and pair errors.
- [Zig library](Zig-Library) owns the public module overview.

Link to the owning page instead of duplicating a detailed contract in every recipe or FAQ answer.
