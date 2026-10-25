//! Newline search over 64-byte blocks, with 32-byte AVX2 vectors when the running CPU has them and 16-byte vectors otherwise.

const std = @import("std");
const builtin = @import("builtin");
const build_options = @import("build_options");

const BLOCK_BYTES = 64;
const SSE_LANES = 16;
const AVX2_LANES = 32;
// A call into the separate AVX2 object cannot be inlined, so shorter spans stay on inline 16-byte code.
const AVX2_MIN_BYTES = 1024;

const SELECT_AT_RUN_TIME = build_options.kernel_backend == .dispatch and build_options.simd_avx2 == .object;
const INLINE_LANES = if (build_options.kernel_backend == .dispatch and build_options.simd_avx2 == .direct) AVX2_LANES else SSE_LANES;

pub inline fn findNewline(bytes: []const u8) ?usize {
    if (!SELECT_AT_RUN_TIME) return @call(.always_inline, findNewlineFor, .{ INLINE_LANES, bytes });
    const head = @min(bytes.len, AVX2_MIN_BYTES);
    if (@call(.always_inline, findNewlineFor, .{ SSE_LANES, bytes[0..head] })) |index| return index;
    if (head == bytes.len) return null;
    if (!cpuHasAvx2()) return if (findNewlineFor(SSE_LANES, bytes[head..])) |index| head + index else null;
    const rest = zfq_avx2_find_newline(bytes.ptr + head, bytes.len - head);
    return if (rest == bytes.len - head) null else head + rest;
}

pub inline fn findNewlines(bytes: []const u8, ends: *[4]usize) bool {
    if (!SELECT_AT_RUN_TIME) return @call(.always_inline, collectNewlines, .{ INLINE_LANES, bytes, 0, ends, 0 }) == ends.len;
    const head = @min(bytes.len, AVX2_MIN_BYTES);
    const count = @call(.always_inline, collectNewlines, .{ SSE_LANES, bytes[0..head], 0, ends, 0 });
    if (count == ends.len or head == bytes.len) return count == ends.len;
    if (!cpuHasAvx2()) return collectNewlines(SSE_LANES, bytes[head..], head, ends, count) == ends.len;
    return zfq_avx2_collect_newlines(bytes.ptr + head, bytes.len - head, head, ends, count) == ends.len;
}

pub inline fn containsNewline(first: []const u8, second: ?[]const u8) bool {
    std.debug.assert(second == null or second.?.len == first.len);
    const second_ptr: ?[*]const u8 = if (second) |s| s.ptr else null;
    if (SELECT_AT_RUN_TIME and first.len >= AVX2_MIN_BYTES and cpuHasAvx2()) {
        return zfq_avx2_contains_newline(first.ptr, second_ptr, first.len);
    }
    return @call(.always_inline, containsNewlineFor, .{ INLINE_LANES, first.ptr, second_ptr, first.len });
}

// --- Run-time selection ---

extern fn zfq_avx2_find_newline(bytes: [*]const u8, len: usize) callconv(.c) usize;
extern fn zfq_avx2_collect_newlines(bytes: [*]const u8, len: usize, base: usize, ends: *[4]usize, count: usize) callconv(.c) usize;
extern fn zfq_avx2_contains_newline(first: [*]const u8, second: ?[*]const u8, len: usize) callconv(.c) bool;

/// True when a span of `len` bytes should take the separate AVX2 object: the
/// running CPU has AVX2 and the span is long enough to repay the call.
pub inline fn useAvx2(len: usize) bool {
    return SELECT_AT_RUN_TIME and len >= AVX2_MIN_BYTES and cpuHasAvx2();
}

comptime {
    if (build_options.simd_object) {
        // The record scanner's and statistics' AVX2 copies are compiled into the same object.
        _ = @import("fastq.zig");
        _ = @import("stats.zig");
        @export(&avx2FindNewline, .{ .name = "zfq_avx2_find_newline" });
        @export(&avx2CollectNewlines, .{ .name = "zfq_avx2_collect_newlines" });
        @export(&avx2ContainsNewline, .{ .name = "zfq_avx2_contains_newline" });
    }
}

fn avx2FindNewline(bytes: [*]const u8, len: usize) callconv(.c) usize {
    return findNewlineFor(AVX2_LANES, bytes[0..len]) orelse len;
}

fn avx2CollectNewlines(bytes: [*]const u8, len: usize, base: usize, ends: *[4]usize, count: usize) callconv(.c) usize {
    return collectNewlines(AVX2_LANES, bytes[0..len], base, ends, count);
}

fn avx2ContainsNewline(first: [*]const u8, second: ?[*]const u8, len: usize) callconv(.c) bool {
    return containsNewlineFor(AVX2_LANES, first, second, len);
}

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

inline fn blockMask(comptime lanes: usize, block: *const [BLOCK_BYTES]u8) u64 {
    const Bytes = @Vector(lanes, u8);
    var mask: u64 = 0;
    inline for (0..BLOCK_BYTES / lanes) |i| {
        const bytes: Bytes = block[i * lanes ..][0..lanes].*;
        const lane_mask: @Int(.unsigned, lanes) = @bitCast(bytes == @as(Bytes, @splat('\n')));
        mask |= @as(u64, lane_mask) << @intCast(i * lanes);
    }
    return mask;
}

inline fn chunkMask(chunk: *const [SSE_LANES]u8) u16 {
    const Bytes = @Vector(SSE_LANES, u8);
    const bytes: Bytes = chunk.*;
    return @bitCast(bytes == @as(Bytes, @splat('\n')));
}

// The last 16 bytes are re-read, so the bits of bytes before `pos` are shifted out.
inline fn tailMask(bytes: []const u8, pos: usize) u16 {
    const start = bytes.len - SSE_LANES;
    return chunkMask(bytes[start..][0..SSE_LANES]) >> @intCast(pos - start);
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
    while (bytes.len - pos >= SSE_LANES) : (pos += SSE_LANES) {
        const mask = chunkMask(bytes[pos..][0..SSE_LANES]);
        if (mask != 0) return pos + @ctz(mask);
    }
    if (pos == bytes.len) return null;
    if (bytes.len >= SSE_LANES) {
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
    while (bytes.len - pos >= SSE_LANES) : (pos += SSE_LANES) {
        var mask = chunkMask(bytes[pos..][0..SSE_LANES]);
        while (mask != 0) : (mask &= mask - 1) {
            ends[n] = base + pos + @ctz(mask);
            n += 1;
            if (n == ends.len) return n;
        }
    }
    if (pos < bytes.len and bytes.len >= SSE_LANES) {
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
    while (len - pos >= SSE_LANES) : (pos += SSE_LANES) {
        if (anyNewline(SSE_LANES, SSE_LANES, first + pos, if (second) |s| s + pos else null)) return true;
    }
    if (pos == len) return false;
    if (len >= SSE_LANES) {
        const start = len - SSE_LANES;
        return anyNewline(SSE_LANES, SSE_LANES, first + start, if (second) |s| s + start else null);
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
        inline for (.{ SSE_LANES, AVX2_LANES, 0 }) |lanes| {
            try std.testing.expectEqual(want_first, if (lanes == 0) findNewline(bytes) else findNewlineFor(lanes, bytes));
            var ends: [4]usize = undefined;
            const four = if (lanes == 0) findNewlines(bytes, &ends) else findNewlinesFor(lanes, bytes, &ends);
            try std.testing.expectEqual(want_count == 4, four);
            if (four) try std.testing.expectEqualSlices(usize, &want_ends, &ends);
        }
    }
}

test "[property] - [newline kernels]: two-slice search is exact at every length and alignment" {
    var first_storage: [5 * BLOCK_BYTES + 2]u8 = undefined;
    var second_storage: [5 * BLOCK_BYTES + 2]u8 = undefined;
    for ([_]usize{ 0, 1, SSE_LANES - 1, AVX2_LANES - 1 }) |alignment| {
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
                    inline for (.{ SSE_LANES, AVX2_LANES }) |lanes| {
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
        inline for (.{ SSE_LANES, AVX2_LANES }) |lanes| {
            try std.testing.expect(findNewlinesFor(lanes, &data, &actual));
            try std.testing.expectEqual(expected, actual);
            try std.testing.expect(!findNewlinesFor(lanes, data[0..expected[3]], &actual));
        }
        try std.testing.expect(findNewlines(&data, &actual));
        try std.testing.expectEqual(expected, actual);
        try std.testing.expect(!findNewlines(data[0..expected[3]], &actual));
    }
}

test "[unit] - [cpu]: detection agrees with the AVX2 guarantee of the build target" {
    if (comptime builtin.cpu.arch != .x86_64) return error.SkipZigTest;
    if (comptime std.Target.x86.featureSetHas(builtin.cpu.features, .avx2)) try std.testing.expect(cpuHasAvx2());
    try std.testing.expectEqual(cpuHasAvx2(), cpuHasAvx2());
}
