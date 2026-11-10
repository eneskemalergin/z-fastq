//! Process entry of x86-64 Linux builds whose target has the x86-64-v3 set, such as the release binaries. It is
//! compiled for any x86-64 CPU and runs before any other code: it checks that the CPU has the set the build
//! needs (x86-64-v3, or x86-64-v4 for builds that use AVX-512), stops with status 2 and a message when it does
//! not, and otherwise continues at the standard entry with the stack as the kernel left it.

const std = @import("std");
const simd = @import("simd.zig");
const options = @import("build_options");

const MESSAGE = "error: this z-fastq build needs an x86-64 CPU with " ++
    (if (options.require_x86_64_v4) "the x86-64-v4 set (AVX-512)" else "the x86-64-v3 set (AVX2, BMI2, FMA)") ++
    "; build z-fastq from source for this CPU\n";

export fn zfq_x86_64_entry() callconv(.naked) noreturn {
    asm volatile (
        \\ movq %%rsp, %%r12
        \\ andq $-16, %%rsp
        \\ callq %[checkCpu:P]
        \\ movq %%r12, %%rsp
        \\ jmp _start
        :
        : [checkCpu] "X" (&checkCpu),
    );
}

fn checkCpu() callconv(.c) void {
    const supported = if (options.require_x86_64_v4) simd.cpuHasX86V4() else simd.cpuHasX86V3();
    if (supported) return;
    _ = std.os.linux.write(2, MESSAGE, MESSAGE.len);
    std.os.linux.exit_group(2);
}
