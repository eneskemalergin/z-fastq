# Installation

Each [release](https://github.com/eneskemalergin/z-fastq/releases) has a ready-made binary for Linux and for macOS, each on x86-64 and on aarch64; anything else is built from source with Zig. The x86-64 binaries need a CPU with the x86-64-v3 set (AVX2, BMI2, and FMA: Intel since 2013, AMD since 2015); on an older CPU the Linux binary stops with a message, and z-fastq is built from source instead.

## Download a release

```bash
curl -fLO https://github.com/eneskemalergin/z-fastq/releases/download/v0.2.1/z-fastq-0.2.1-linux-x86_64.tar.gz
curl -fLO https://github.com/eneskemalergin/z-fastq/releases/download/v0.2.1/SHA256SUMS
sha256sum --check --ignore-missing SHA256SUMS
tar -xzf z-fastq-0.2.1-linux-x86_64.tar.gz
./z-fastq-0.2.1-linux-x86_64/z-fastq --version
```

On an ARM machine, use `linux-aarch64` in place of `linux-x86_64`; on a Mac, use `macos-aarch64` (Apple Silicon) or `macos-x86_64` (Intel), and check the archive with `grep macos-aarch64 SHA256SUMS | shasum -a 256 --check`. The archive holds the binary and its license. The Linux binary is static and the macOS binary needs only the system library every Mac has, so copy it into a directory on your `PATH`. With the [GitHub CLI](https://cli.github.com/), `gh attestation verify z-fastq-0.2.1-linux-x86_64.tar.gz --repo eneskemalergin/z-fastq` also confirms the archive was built by this repository's release workflow.

On Apple Silicon, use `macos-aarch64`: the `macos-x86_64` binary under Rosetta needs macOS 15 or later, the first Rosetta that runs AVX2 code.

The macOS binary is not notarized. An archive fetched with `curl`, as above, runs as it is; macOS blocks a binary unpacked from an archive downloaded in a browser until its quarantine attribute is removed with `xattr -d com.apple.quarantine z-fastq`.

## Build from source

The build needs nothing beyond Zig itself: no C compiler, no NASM, and no system zlib. The release build below is fully static on Linux; on macOS it links only the system library, `libSystem`.

**You need**

- Linux on x86-64 or aarch64, or macOS 13 or later on Apple Silicon or Intel.
- [Zig](https://ziglang.org/download/) 0.16.0.
- Network access for the first build, which downloads z-fastq's one package dependency, [zipir](https://github.com/eneskemalergin/zipir) 0.3.0. Later builds reuse it.

From the repository root:

```bash
zig build -Doptimize=ReleaseFast
./zig-out/bin/z-fastq --version
```

The second line prints `z-fastq 0.2.1`.

- `-Doptimize=ReleaseSafe` keeps Zig's safety checks in the release build.
- `-Dstatic=true` confirms a static Linux release build; it stops with a message on macOS, where no build is static.
- `-Dcpu` names the CPU a release build targets, such as `-Dcpu=x86_64_v3`; the binary then runs only on CPUs with that CPU's features. Without it, a Linux release build targets generic x86-64 or ARMv8-A, and a macOS release build targets Zig's macOS defaults, `core2` (Intel) or `apple_m1` (Apple Silicon).
- A plain `zig build` gives a debug build, for working on the code.

The binary stays in `zig-out/bin/`; there is no package-manager install yet. Copy or link it into a directory on your `PATH` to call it as `z-fastq`.

## Faster code on newer CPUs

The x86-64 release binaries, for Linux and for macOS, are built for the x86-64-v3 set, so every command runs AVX2 code. Where the CPU also has the x86-64-v4 set (AVX-512 F, BW, CD, DQ, and VL) and the operating system enables it, line search, validation, and base and quality counting on records of 1 KiB or more switch to AVX-512 code at run time. A build from source without `-Dcpu` targets generic x86-64 and runs on any x86-64 CPU, switching to AVX2 or AVX-512 code at run time for records of 1 KiB or more only, and PCLMUL and SSE4.1 for the gzip checksum. The output is identical either way. The aarch64 release binaries select nothing at run time: the Linux one runs on any ARMv8-A CPU, the macOS one on any Apple Silicon Mac.

To run slower kernels than the CPU allows, for example to compare speeds or to rule a kernel set out, set `ZFASTQ_TIER=avx2`, or `ZFASTQ_TIER=baseline` in a build for generic x86-64 (the release binaries' lowest tier is `avx2`); `ZFASTQ_TIER=avx512` names the default on a CPU with AVX-512. A tier the binary cannot run on the CPU, or any other value, stops the command with status 2 and names the tiers it can run. Unset, z-fastq picks the fastest tier.

<details>
<summary>Build without run-time selection</summary>

```bash
zig build -Doptimize=ReleaseFast -Dkernel-backend=portable
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
