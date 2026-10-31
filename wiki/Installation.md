# Installation

Each [release](https://github.com/eneskemalergin/z-fastq/releases) has a ready-made static binary for Linux on x86-64 and on aarch64; anything else is built from source with Zig.

## Download a release

```bash
curl -fLO https://github.com/eneskemalergin/z-fastq/releases/download/v0.2.0/z-fastq-0.2.0-linux-x86_64.tar.gz
curl -fLO https://github.com/eneskemalergin/z-fastq/releases/download/v0.2.0/SHA256SUMS
sha256sum --check --ignore-missing SHA256SUMS
tar -xzf z-fastq-0.2.0-linux-x86_64.tar.gz
./z-fastq-0.2.0-linux-x86_64/z-fastq --version
```

On an ARM machine, use `linux-aarch64` in place of `linux-x86_64`. The archive holds the binary and its license. It needs no libraries, so copy it into a directory on your `PATH`. With the [GitHub CLI](https://cli.github.com/), `gh attestation verify z-fastq-0.2.0-linux-x86_64.tar.gz --repo eneskemalergin/z-fastq` also confirms the archive was built by this repository's release workflow.

## Build from source

The build needs nothing beyond Zig itself: no C compiler, no NASM, and no system zlib, and the release build below is fully static.

**You need**

- Linux on x86-64 or aarch64.
- [Zig](https://ziglang.org/download/) 0.16.0.
- Network access for the first build, which downloads z-fastq's one package dependency, [zipir](https://github.com/eneskemalergin/zipir) 0.3.0. Later builds reuse it.

From the repository root:

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast
./zig-out/bin/z-fastq --version
```

The second line prints `z-fastq 0.2.0`.

- `-Doptimize=ReleaseSafe` keeps Zig's safety checks in the release build.
- A plain `zig build` gives a debug build, for working on the code.

The binary stays in `zig-out/bin/`; there is no package-manager install yet. Copy or link it into a directory on your `PATH` to call it as `z-fastq`.

## Faster code on newer CPUs

The x86-64 release binary runs on any x86-64 CPU and switches to faster code where the CPU allows: PCLMUL and SSE4.1 for the gzip checksum, and AVX2 for line search, validation, and base and quality counting on records of 1 KiB or more. The output is identical either way. The aarch64 release binary runs on any ARMv8-A CPU and selects nothing at run time.

<details>
<summary>Build without run-time selection</summary>

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast -Dkernel-backend=portable
```

This build uses only the instructions of its target CPU. The default is `-Dkernel-backend=dispatch`, and the tests pass with both.

</details>

## Check the build

```bash
zig build test --summary all
```

## Related pages

- [Getting started](Getting-Started): a first run, command by command.
- [Troubleshooting](Troubleshooting): when the build or the binary does not behave.
- [For contributors](For-Contributors): the full set of test runs.
