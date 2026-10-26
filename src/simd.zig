//! Byte-kernel tiers, their run-time selection, and the newline search kernels.
//!
//! Every byte kernel has a scalar reference and one generic vector body over a
//! comptime lane count. Inline code uses the widest tier the compile target
//! guarantees. When the target lacks AVX2, the same bodies are also compiled
//! into a separate AVX2 object, and long spans call it on CPUs that have AVX2.

const std = @import("std");
const builtin = @import("builtin");
const build_options = @import("build_options");

pub const Tier = enum { baseline, avx2, avx512 };

/// Vector width of each tier in bytes. AVX-512 is not selected yet: it needs
/// CPU detection and its own object, tested on a host that has it. Kernels are
/// tested at its width already.
pub fn tierWidth(comptime tier: Tier) comptime_int {
    return switch (tier) {
        .baseline => 16,
        .avx2 => 32,
        .avx512 => 64,
    };
}

/// The widest tier the compile target guarantees.
pub const INLINE_TIER: Tier = if (std.Target.x86.featureSetHas(builtin.cpu.features, .avx2)) .avx2 else .baseline;
pub const INLINE_LANES = tierWidth(INLINE_TIER);
const BASELINE_LANES = tierWidth(.baseline);

/// Structural block: one 64-bit mask holds the LF positions of 64 bytes.
pub const BLOCK_BYTES = 64;

// A call into the separate AVX2 object cannot be inlined, so shorter spans stay inline.
const AVX2_MIN_BYTES = 1024;

const SELECT_AT_RUN_TIME = build_options.avx2_object;

/// True when a call that will process `len` bytes should take the AVX2 object:
/// the build links it, the span repays the call, and the running CPU has AVX2.
/// A search that cannot know its length first (`findNewline`, the complete
/// record scanner) applies this to a bound or a hint instead.
pub inline fn useAvx2(len: usize) bool {
    return SELECT_AT_RUN_TIME and len >= AVX2_MIN_BYTES and cpuHasAvx2();
}

// --- Kernel table ---

// Each kernel's object entry point is declared once, as the exported function;
// the caller's declaration takes its type from it.

/// Exports `function` from the object under the kernel name `name`. Call from a
/// `comptime` block; outside the object it does nothing.
pub fn exportKernel(comptime name: []const u8, comptime function: anytype) void {
    if (build_options.simd_object) @export(function, .{ .name = symbol(INLINE_TIER, name) });
}

/// The AVX2 object's copy of the kernel exported as `name` from `function`.
/// Reach it only behind `useAvx2`, which is false in builds without the object.
pub fn avx2Kernel(comptime name: []const u8, comptime function: anytype) *const @TypeOf(function.*) {
    return @extern(*const @TypeOf(function.*), .{ .name = symbol(.avx2, name) });
}

fn symbol(comptime tier: Tier, comptime name: []const u8) []const u8 {
    return "zfq_" ++ @tagName(tier) ++ "_" ++ name;
}

comptime {
    if (build_options.simd_object) {
        // The object's root: the record scanner's and statistics' kernels are exported from their files.
        _ = @import("fastq.zig");
        _ = @import("stats.zig");
    }
    exportKernel("find_newline", &findNewlineExport);
    exportKernel("collect_newlines", &collectNewlinesExport);
    exportKernel("contains_newline", &containsNewlineExport);
}

// --- Newline search ---

/// Searches start inline; the bytes after the first `AVX2_MIN_BYTES` go to the
/// object, so a short line never pays for the call.
pub inline fn findNewline(bytes: []const u8) ?usize {
    if (!SELECT_AT_RUN_TIME) return @call(.always_inline, findNewlineFor, .{ INLINE_LANES, bytes });
    const head = @min(bytes.len, AVX2_MIN_BYTES);
    if (@call(.always_inline, findNewlineFor, .{ BASELINE_LANES, bytes[0..head] })) |index| return index;
    if (head == bytes.len) return null;
    if (!useAvx2(bytes.len)) return if (findNewlineFor(BASELINE_LANES, bytes[head..])) |index| head + index else null;
    const rest = avx2Kernel("find_newline", &findNewlineExport)(bytes.ptr + head, bytes.len - head);
    return if (rest == bytes.len - head) null else head + rest;
}

pub inline fn findNewlines(bytes: []const u8, ends: *[4]usize) bool {
    if (!SELECT_AT_RUN_TIME) return @call(.always_inline, collectNewlines, .{ INLINE_LANES, bytes, 0, ends, 0 }) == ends.len;
    const head = @min(bytes.len, AVX2_MIN_BYTES);
    const count = @call(.always_inline, collectNewlines, .{ BASELINE_LANES, bytes[0..head], 0, ends, 0 });
    if (count == ends.len or head == bytes.len) return count == ends.len;
    if (!useAvx2(bytes.len)) return collectNewlines(BASELINE_LANES, bytes[head..], head, ends, count) == ends.len;
    return avx2Kernel("collect_newlines", &collectNewlinesExport)(bytes.ptr + head, bytes.len - head, head, ends, count) == ends.len;
}

pub inline fn containsNewline(first: []const u8, second: ?[]const u8) bool {
    std.debug.assert(second == null or second.?.len == first.len);
    const second_ptr: ?[*]const u8 = if (second) |s| s.ptr else null;
    if (useAvx2(first.len)) {
        return avx2Kernel("contains_newline", &containsNewlineExport)(first.ptr, second_ptr, first.len);
    }
    return @call(.always_inline, containsNewlineFor, .{ INLINE_LANES, first.ptr, second_ptr, first.len });
}

fn findNewlineExport(bytes: [*]const u8, len: usize) callconv(.c) usize {
    return findNewlineFor(INLINE_LANES, bytes[0..len]) orelse len;
}

fn collectNewlinesExport(bytes: [*]const u8, len: usize, base: usize, ends: *[4]usize, count: usize) callconv(.c) usize {
    return collectNewlines(INLINE_LANES, bytes[0..len], base, ends, count);
}

fn containsNewlineExport(first: [*]const u8, second: ?[*]const u8, len: usize) callconv(.c) bool {
    return containsNewlineFor(INLINE_LANES, first, second, len);
}

// --- CPU detection ---

var avx2_state: std.atomic.Value(u8) = .init(0);

fn cpuHasAvx2() bool {
    const state = avx2_state.load(.monotonic);
    if (state != 0) return state == 2;
    const detected = detectAvx2();
    avx2_state.store(if (detected) 2 else 1, .monotonic);
    return detected;
}

fn detectAvx2() bool {
    if (comptime builtin.cpu.arch != .x86_64) return false;
    if (cpuid(0).eax < 7) return false;
    const leaf1 = cpuid(1);
    const has_osxsave = leaf1.ecx & (1 << 27) != 0;
    const has_avx = leaf1.ecx & (1 << 28) != 0;
    // AVX code faults unless the OS saves XMM and YMM state (XCR0 bits 1 and 2).
    if (!has_osxsave or !has_avx or xgetbv0() & 0x6 != 0x6) return false;
    return cpuid(7).ebx & (1 << 5) != 0;
}

const Cpuid = struct { eax: u32, ebx: u32, ecx: u32, edx: u32 };

fn cpuid(leaf: u32) Cpuid {
    var eax: u32 = undefined;
    var ebx: u32 = undefined;
    var ecx: u32 = undefined;
    var edx: u32 = undefined;
    asm volatile ("cpuid"
        : [_] "={eax}" (eax),
          [_] "={ebx}" (ebx),
          [_] "={ecx}" (ecx),
          [_] "={edx}" (edx),
        : [_] "{eax}" (leaf),
          [_] "{ecx}" (@as(u32, 0)),
    );
    return .{ .eax = eax, .ebx = ebx, .ecx = ecx, .edx = edx };
}

fn xgetbv0() u32 {
    return asm volatile (
        \\xorl %%ecx, %%ecx
        \\xgetbv
        : [_] "={eax}" (-> u32),
        :
        : .{ .ecx = true, .edx = true });
}

// --- Kernels ---

/// Bit `i` is set when `block[i]` is LF.
pub inline fn newlineMask(comptime width: usize, block: *const [width]u8) @Int(.unsigned, width) {
    const Bytes = @Vector(width, u8);
    const bytes: Bytes = block.*;
    return @bitCast(bytes == @as(Bytes, @splat('\n')));
}

inline fn blockMask(comptime width: usize, block: *const [BLOCK_BYTES]u8) u64 {
    var mask: u64 = 0;
    inline for (0..BLOCK_BYTES / width) |i| {
        mask |= @as(u64, newlineMask(width, block[i * width ..][0..width])) << @intCast(i * width);
    }
    return mask;
}

// The last 16 bytes are re-read, so the bits of bytes before `pos` are shifted out.
inline fn tailMask(bytes: []const u8, pos: usize) u16 {
    const start = bytes.len - BASELINE_LANES;
    return newlineMask(BASELINE_LANES, bytes[start..][0..BASELINE_LANES]) >> @intCast(pos - start);
}

inline fn anyNewline(comptime lanes: usize, comptime len: usize, first: [*]const u8, second: ?[*]const u8) bool {
    const Bytes = @Vector(lanes, u8);
    const lf: Bytes = @splat('\n');
    var any: @Vector(lanes, bool) = @splat(false);
    inline for (0..len / lanes) |i| {
        const a: Bytes = first[i * lanes ..][0..lanes].*;
        any |= a == lf;
        if (second) |s| {
            const b: Bytes = s[i * lanes ..][0..lanes].*;
            any |= b == lf;
        }
    }
    return @reduce(.Or, any);
}

// A short read's line never holds 192 LF-free bytes; long lines then skip
// LF-free stretches 128 bytes at a time with one mask test.
const SKIP_AFTER_EMPTY_BYTES = 192;

inline fn skipLineBody(comptime lanes: usize, bytes: []const u8, start: usize) usize {
    var pos = start;
    while (bytes.len - pos >= 2 * BLOCK_BYTES and !anyNewline(lanes, 2 * BLOCK_BYTES, bytes.ptr + pos, null)) pos += 2 * BLOCK_BYTES;
    return pos;
}

pub fn findNewlineFor(comptime lanes: usize, bytes: []const u8) ?usize {
    var pos: usize = 0;
    var empty_bytes: usize = 0;
    while (bytes.len - pos >= BLOCK_BYTES) {
        const mask = blockMask(lanes, bytes[pos..][0..BLOCK_BYTES]);
        if (mask != 0) return pos + @ctz(mask);
        pos += BLOCK_BYTES;
        empty_bytes += BLOCK_BYTES;
        if (empty_bytes >= SKIP_AFTER_EMPTY_BYTES) pos = skipLineBody(lanes, bytes, pos);
    }
    while (bytes.len - pos >= BASELINE_LANES) : (pos += BASELINE_LANES) {
        const mask = newlineMask(BASELINE_LANES, bytes[pos..][0..BASELINE_LANES]);
        if (mask != 0) return pos + @ctz(mask);
    }
    if (pos == bytes.len) return null;
    if (bytes.len >= BASELINE_LANES) {
        const mask = tailMask(bytes, pos);
        return if (mask != 0) pos + @ctz(mask) else null;
    }
    for (bytes[pos..], pos..) |byte, index| {
        if (byte == '\n') return index;
    }
    return null;
}

pub fn collectNewlines(comptime lanes: usize, bytes: []const u8, base: usize, ends: *[4]usize, count: usize) usize {
    var n = count;
    var pos: usize = 0;
    var empty_bytes: usize = 0;
    while (bytes.len - pos >= BLOCK_BYTES) {
        var mask = blockMask(lanes, bytes[pos..][0..BLOCK_BYTES]);
        if (mask == 0) {
            pos += BLOCK_BYTES;
            empty_bytes += BLOCK_BYTES;
            if (empty_bytes >= SKIP_AFTER_EMPTY_BYTES) pos = skipLineBody(lanes, bytes, pos);
            continue;
        }
        empty_bytes = 0;
        while (mask != 0) : (mask &= mask - 1) {
            ends[n] = base + pos + @ctz(mask);
            n += 1;
            if (n == ends.len) return n;
        }
        pos += BLOCK_BYTES;
    }
    while (bytes.len - pos >= BASELINE_LANES) : (pos += BASELINE_LANES) {
        var mask = newlineMask(BASELINE_LANES, bytes[pos..][0..BASELINE_LANES]);
        while (mask != 0) : (mask &= mask - 1) {
            ends[n] = base + pos + @ctz(mask);
            n += 1;
            if (n == ends.len) return n;
        }
    }
    if (pos < bytes.len and bytes.len >= BASELINE_LANES) {
        var mask = tailMask(bytes, pos);
        while (mask != 0) : (mask &= mask - 1) {
            ends[n] = base + pos + @ctz(mask);
            n += 1;
            if (n == ends.len) return n;
        }
        return n;
    }
    for (bytes[pos..], pos..) |byte, index| {
        if (byte != '\n') continue;
        ends[n] = base + index;
        n += 1;
        if (n == ends.len) return n;
    }
    return n;
}

pub fn findNewlinesFor(comptime lanes: usize, bytes: []const u8, ends: *[4]usize) bool {
    return collectNewlines(lanes, bytes, 0, ends, 0) == ends.len;
}

pub fn containsNewlineFor(comptime lanes: usize, first: [*]const u8, second: ?[*]const u8, len: usize) bool {
    var pos: usize = 0;
    while (len - pos >= BLOCK_BYTES) : (pos += BLOCK_BYTES) {
        if (anyNewline(lanes, BLOCK_BYTES, first + pos, if (second) |s| s + pos else null)) return true;
    }
    while (len - pos >= BASELINE_LANES) : (pos += BASELINE_LANES) {
        if (anyNewline(BASELINE_LANES, BASELINE_LANES, first + pos, if (second) |s| s + pos else null)) return true;
    }
    if (pos == len) return false;
    if (len >= BASELINE_LANES) {
        const start = len - BASELINE_LANES;
        return anyNewline(BASELINE_LANES, BASELINE_LANES, first + start, if (second) |s| s + start else null);
    }
    for (pos..len) |index| {
        if (first[index] == '\n') return true;
        if (second) |s| if (s[index] == '\n') return true;
    }
    return false;
}

test "[property] - [newline kernels]: every width and the selected path match a scalar scan" {
    var prng: std.Random.DefaultPrng = .init(std.testing.random_seed);
    const random = prng.random();
    var buffer: [700]u8 = undefined;
    for (0..2000) |_| {
        const len = random.uintLessThan(usize, buffer.len + 1);
        const density = random.uintLessThan(u8, 64);
        for (buffer[0..len]) |*byte| {
            byte.* = if (density != 0 and random.uintLessThan(u8, 64) < density) '\n' else 'A' + random.uintLessThan(u8, 26);
        }
        const bytes = buffer[0..len];
        var want_ends: [4]usize = undefined;
        var want_count: usize = 0;
        for (bytes, 0..) |byte, index| {
            if (byte == '\n' and want_count < 4) {
                want_ends[want_count] = index;
                want_count += 1;
            }
        }
        const want_first: ?usize = if (want_count > 0) want_ends[0] else null;
        inline for (.{ BASELINE_LANES, tierWidth(.avx2), tierWidth(.avx512), 0 }) |lanes| {
            try std.testing.expectEqual(want_first, if (lanes == 0) findNewline(bytes) else findNewlineFor(lanes, bytes));
            var ends: [4]usize = undefined;
            const four = if (lanes == 0) findNewlines(bytes, &ends) else findNewlinesFor(lanes, bytes, &ends);
            try std.testing.expectEqual(want_count == 4, four);
            if (four) try std.testing.expectEqualSlices(usize, &want_ends, &ends);
        }
    }
}

test "[property] - [newline kernels]: masks preserve newline positions at every width" {
    inline for (.{ 1, BASELINE_LANES, tierWidth(.avx2), BLOCK_BYTES }) |lanes| {
        const Mask = @Int(.unsigned, lanes);
        var block: [lanes]u8 = @splat('A');
        try std.testing.expectEqual(@as(Mask, 0), newlineMask(lanes, &block));
        for (0..lanes) |lane| {
            block[lane] = '\n';
            try std.testing.expectEqual(@as(Mask, 1) << @intCast(lane), newlineMask(lanes, &block));
            block[lane] = 'A';
        }

        var expected: Mask = 0;
        for ([_]usize{ 0, 1, 15, 16, 31, 32, 47, 48, 62, 63 }) |lane| {
            if (lane >= lanes) continue;
            block[lane] = '\n';
            expected |= @as(Mask, 1) << @intCast(lane);
        }
        try std.testing.expectEqual(expected, newlineMask(lanes, &block));
    }
}

test "[property] - [newline kernels]: two-slice search is exact at every length and alignment" {
    var first_storage: [5 * BLOCK_BYTES + tierWidth(.avx512)]u8 = undefined;
    var second_storage: [5 * BLOCK_BYTES + tierWidth(.avx512)]u8 = undefined;
    for ([_]usize{ 0, 1, BASELINE_LANES - 1, tierWidth(.avx2) - 1, tierWidth(.avx512) - 1 }) |alignment| {
        for (0..4 * BLOCK_BYTES + 3) |len| {
            // LF bytes just outside the slices must never be reported.
            @memset(&first_storage, '\n');
            @memset(&second_storage, '\n');
            const first = first_storage[alignment..][0..len];
            const second = second_storage[alignment..][0..len];
            @memset(first, 0xff);
            @memset(second, '\r');
            for (0..len + 1) |pos| {
                inline for (.{ first, second }, .{ 0xff, '\r' }) |slice, original| {
                    if (pos < len) slice[pos] = '\n';
                    const first_has_lf = std.mem.findScalar(u8, first, '\n') != null;
                    const second_has_lf = std.mem.findScalar(u8, second, '\n') != null;
                    inline for (.{ BASELINE_LANES, tierWidth(.avx2), tierWidth(.avx512) }) |lanes| {
                        try std.testing.expectEqual(first_has_lf, containsNewlineFor(lanes, first.ptr, null, len));
                        try std.testing.expectEqual(first_has_lf or second_has_lf, containsNewlineFor(lanes, first.ptr, second.ptr, len));
                    }
                    try std.testing.expectEqual(first_has_lf, containsNewline(first, null));
                    try std.testing.expectEqual(first_has_lf or second_has_lf, containsNewline(first, second));
                    if (pos < len) slice[pos] = original;
                }
            }
        }
    }
}

test "[edge] - [newline kernels]: four-LF search is exact at vector and block boundaries" {
    const cases = [_][4]usize{
        .{ 14, 15, 30, 31 },
        .{ 15, 16, 31, 32 },
        .{ 30, 31, 62, 63 },
        .{ 31, 32, 63, 64 },
        .{ 63, 64, 127, 128 },
        .{ 126, 127, 128, 129 },
        .{ 300, 301, 500, 501 },
    };
    for (cases) |expected| {
        var data: [520]u8 = @splat('x');
        for (expected) |index| data[index] = '\n';
        data[510] = '\n';
        var actual: [4]usize = undefined;
        inline for (.{ BASELINE_LANES, tierWidth(.avx2), tierWidth(.avx512) }) |lanes| {
            try std.testing.expect(findNewlinesFor(lanes, &data, &actual));
            try std.testing.expectEqual(expected, actual);
            try std.testing.expect(!findNewlinesFor(lanes, data[0..expected[3]], &actual));
        }
        try std.testing.expect(findNewlines(&data, &actual));
        try std.testing.expectEqual(expected, actual);
        try std.testing.expect(!findNewlines(data[0..expected[3]], &actual));
    }
}

test "[property] - [newline kernels]: the run-time selected path is exact around its length gate" {
    // The span, one sentinel LF, and an LF-free span of the same length.
    var storage: [7 * AVX2_MIN_BYTES]u8 = undefined;
    const lengths = [_]usize{ AVX2_MIN_BYTES - 1, AVX2_MIN_BYTES, AVX2_MIN_BYTES + 1, AVX2_MIN_BYTES + 77, 2 * AVX2_MIN_BYTES, 3 * AVX2_MIN_BYTES + 5 };
    const positions = [_]usize{ 0, 63, AVX2_MIN_BYTES - 1, AVX2_MIN_BYTES, AVX2_MIN_BYTES + 1, AVX2_MIN_BYTES + 31, 2 * AVX2_MIN_BYTES - 1, 3 * AVX2_MIN_BYTES };
    for ([_]usize{ 0, 1, tierWidth(.avx2) - 1 }) |alignment| {
        for (lengths) |len| {
            const bytes = storage[alignment..][0..len];
            for (positions) |first_lf| {
                @memset(&storage, 'A');
                // LF bytes just past the span must never be reported.
                storage[alignment + len] = '\n';
                var want_ends: [4]usize = undefined;
                var want_count: usize = 0;
                for ([_]usize{ first_lf, first_lf + 1, first_lf + 100, first_lf + 700 }) |index| {
                    if (index >= len) break;
                    bytes[index] = '\n';
                    want_ends[want_count] = index;
                    want_count += 1;
                }
                try std.testing.expectEqual(if (want_count > 0) want_ends[0] else null, findNewline(bytes));
                var ends: [4]usize = undefined;
                try std.testing.expectEqual(want_count == 4, findNewlines(bytes, &ends));
                if (want_count == 4) try std.testing.expectEqualSlices(usize, &want_ends, &ends);

                const quiet = storage[alignment + len + 1 ..][0..len];
                @memset(quiet, 'C');
                try std.testing.expectEqual(want_count > 0, containsNewline(bytes, null));
                try std.testing.expectEqual(want_count > 0, containsNewline(quiet, bytes));
                try std.testing.expectEqual(want_count > 0, containsNewline(bytes, quiet));
                try std.testing.expect(!containsNewline(quiet, null));
            }
        }
    }
}

test "[unit] - [cpu]: detection agrees with the AVX2 guarantee of the build target" {
    if (comptime builtin.cpu.arch != .x86_64) return error.SkipZigTest;
    if (comptime std.Target.x86.featureSetHas(builtin.cpu.features, .avx2)) try std.testing.expect(cpuHasAvx2());
    try std.testing.expectEqual(cpuHasAvx2(), cpuHasAvx2());
}
