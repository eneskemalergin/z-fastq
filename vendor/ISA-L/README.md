# ISA-L source boundary

This directory contains the stateful x86-64 inflate and CRC subset of ISA-L 2.32.1 at commit `7c3479e`. The upstream C, assembly, and header files are unchanged from that release. `igzip/inflate_helpers.c` is a plainly marked extraction of three functions from upstream `igzip.c`; it prevents unrelated compression code from entering this decompression-only boundary.

The project compiles the C source directly with Zig and assembles 45 x86-64 entry files with NASM 2.14.01 or newer. z-fastq owns gzip headers, trailers, CRC-32 comparison, ISIZE comparison, concatenated members, input limits, and diagnostics. ISA-L receives only raw DEFLATE payload bytes.

The current build compiles six C translation units and assembles 45 x86-64 entry files. ISA-L's shared CRC dispatcher requires the implementations for its other checksum APIs even though z-fastq selects only gzip CRC. The directory contains 70 files including the license, this record, and source include closure. It omits compression kernels, command-line programs, erasure coding, RAID code, tests, benchmarks, and non-x86 implementations.

This vendor boundary is integrated only for Linux x86-64 glibc and musl. It is an external source dependency and requires NASM when built from source. Static musl release binaries need no separately installed ISA-L or dynamic C runtime, but they still contain ISA-L code and must retain its license notice. Upstream implementations for other architectures are not part of this subset, and their upstream existence does not establish z-fastq support.

The BSD-3-Clause license is in [`LICENSE.md`](LICENSE.md). Its copyright, conditions, and disclaimer must accompany binary distributions in the documentation or other provided materials.

## Update procedure

1. Select a tagged upstream release and record its tag and commit.
1. Review its release notes, license, security status, stateful inflate API, dispatch code, assembler requirement, and supported targets.
1. Recompute the C and NASM include closure from the source list in `build.zig`.
1. Import only that source closure without editing upstream contents.
1. Reconcile the source list, NASM flags, supported targets, state layout, and return-code handling in `build.zig` and `src/io.zig`.
1. Run the Debug, ReleaseFast, native fallback, parity, malformed-input, binary-size, memory, no-thread, and benchmark gates before updating the pin.
