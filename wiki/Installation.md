# Installation

## Supported build

The supported build is:

- Linux x86-64.
- Zig `0.16.0`.
- Network access for the first build, which downloads the pinned Zipir package.

Release builds are static. gzip input is decoded by [Zipir](https://github.com/eneskemalergin/zipir), a Zig package that `build.zig.zon` pins to its `v0.2.0` release archive and package hash. The build needs no NASM, C compiler, or system zlib.

## Build a release binary

From the repository root:

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast
./zig-out/bin/z-fastq --version
```

The current version prints:

```text
z-fastq 0.0.18
```

For a safety-checked static build, use `ReleaseSafe`:

```bash
zig build -Dstatic=true -Doptimize=ReleaseSafe
```

Use `Debug` while changing the project:

```bash
zig build
```

## Portable checksum kernels

The release binary is compiled for generic x86-64. By default Zipir picks its checksum kernels at run time for the CPU the binary runs on: gzip CRC-32 uses PCLMUL and SSE4.1 when available. To build only the portable kernels:

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast -Dkernel-backend=portable
```

The default is `-Dkernel-backend=dispatch`. This does not change the supported target or the output; the test suite passes with both settings.

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
