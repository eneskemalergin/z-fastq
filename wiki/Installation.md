# Installation

## Supported build

The current build contract is:

- Linux x86-64.
- Zig `0.16.0`.
- NASM `2.14.01` or newer for the default ISA-L build.

Release builds are static. The default Linux x86-64 release path uses the vendored ISA-L implementation for gzip and CRC work. You can build without ISA-L when NASM is unavailable or when you want to exercise the Zig implementation.

## Build a release binary

From the repository root:

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast
./zig-out/bin/z-fastq --version
```

The current version prints:

```text
z-fastq 0.0.15
```

For a safety-checked static build, use `ReleaseSafe`:

```bash
zig build -Dstatic=true -Doptimize=ReleaseSafe
```

Use `Debug` while changing the project:

```bash
zig build
```

## Build without ISA-L

The fallback path uses the Zig DEFLATE and CRC implementations:

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast -Disa-l=false
```

This does not change the supported target. It changes the gzip and CRC implementation used by the build.

## Run the tests

```bash
zig build test --summary all
```

The test step covers the parser, writer, count scanner, statistics, validation, sampling, paired operations, and installed CLI behavior.

## Use the binary from the build tree

The project does not install a system-wide command. Run the binary from the build output:

```bash
./zig-out/bin/z-fastq stats reads.fastq
```

If you want `z-fastq` on your `PATH`, copy or link `zig-out/bin/z-fastq` into a directory you manage. The project does not provide a package-manager install path yet.
