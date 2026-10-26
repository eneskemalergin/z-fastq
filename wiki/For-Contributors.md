# For contributors

This page is for contributors who update public documentation alongside the source. Keep each user-facing instruction on the one page that documents it (listed below), and use this map to find the code and tests behind a documented behavior.

## Source map

```text
src/main.zig       CLI parsing, commands, output, and diagnostics
src/root.zig       public Zig exports and version
src/fastq.zig      records, parser, writer, and validation
src/count_scan.zig specialized count scanner
src/stats.zig      aggregate statistics
src/sample.zig     CLI-only sampling grammars and selectors
src/pair.zig       CLI-only paired-name parsing and matching
src/io.zig         byte sources, sinks, limits, and the Zipir gzip adapter
tests/              parser, CLI, pair, sample, and output behavior tests
wiki/               public user documentation
```

## Before changing a documented behavior

1. Find the source file that implements the behavior and its callers.
2. Read the focused tests that establish the current behavior.
3. Classify the claim as observed, inferred, proposed, or unknown.
4. Update the page that documents the behavior only when the user-facing behavior changes.
5. Re-read the page and the diff before reporting completion.

The page should explain why a user needs the behavior, not just restate a function name.

## Verification

Run the required checks from the repository root with Zig 0.16.0 on Linux x86-64. gzip decoding comes from the pinned Zipir package, so no NASM is needed. Debug runs with both kernel backends, then static ReleaseSafe, then static ReleaseFast:

```bash
zig fmt --check src tests build.zig
zig build test --summary all
zig build test -Dkernel-backend=portable --summary all
zig build test -Dstatic=true -Doptimize=ReleaseSafe --summary all
zig build test -Dstatic=true -Doptimize=ReleaseFast --summary all
```

Each test command installs its selected CLI before running the command tests. Run variants sequentially in one checkout because the CLI tests share `zig-out/bin/z-fastq`; parallel variants need separate checkouts. The checks use inline and tracked synthetic fixtures without corpus downloads or peer tools. CLI tests need `pidfd_open`, `/proc`, FIFOs, and hard and symbolic links; a missing facility is a setup failure.

Optional CPU instructions and timing or RSS thresholds are not test prerequisites. The suite has no expected skips: every test in each of the four runs above must pass. Investigate any skip. Zipir's own CI tests its portable kernels and a baseline-CPU build; the portable Debug run here checks that z-fastq's output does not depend on the kernel choice. Debug builds target the host CPU, so on a host with AVX2 they use AVX2 inline; the release runs target generic x86-64 and are the ones that test the AVX2 kernels selected at run time against their scalar references. Do not use a timing result as proof of correct output.

## Which page documents what

- The command pages document user-facing syntax and workflow examples.
- [Input and output](Input-and-Output) documents stream and format behavior.
- [Automation](Automation) documents exit classes and machine output guidance.
- [Limits and supported formats](Limits-and-Supported-Formats) documents current limits.
- [Paired reads](Paired-Reads) documents name normalization and pair errors.
- [Zig library](Zig-Library) documents the public module.

Link to that page instead of repeating a detailed rule in every recipe or FAQ answer.
