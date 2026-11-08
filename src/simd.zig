//! Byte-kernel tiers, their run-time selection, and the newline search kernels.
//!
//! Every byte kernel has a scalar reference and one generic vector body over a
//! comptime lane count. Inline code uses the widest tier the compile target
//! guarantees. When an x86-64 target lacks AVX2 or AVX-512, the same bodies are
//! also compiled into separate AVX2 and AVX-512 objects, and long spans call the
//! best one the CPU runs. On aarch64 the baseline tier's 16 lanes are one NEON
//! register. Kernels whose work per call outweighs the call at any length, such
//! as fastp's pair sum in preprocess.zig, use an object whenever one applies.
//! `selectTier` lowers the run-time choice to any tier from the compiled one to
//! the CPU's best; the CLI sets it once from `ZFASTQ_TIER` before any kernel runs.

const std = @import("std");
const builtin = @import("builtin");
const build_options = @import("build_options");

pub const Tier = enum { baseline, avx2, avx512 };

pub fn tierWidth(comptime tier: Tier) comptime_int {
    return switch (tier) {
        .baseline => 16,
        .avx2 => 32,
        .avx512 => 64,
    };
}

pub const INLINE_TIER: Tier = if (builtin.cpu.arch != .x86_64)
    .baseline
else if (targetHas(&.{ .avx512f, .avx512bw, .avx512vl, .avx512dq }))
    .avx512
else if (targetHas(&.{.avx2}))
    .avx2
else
    .baseline;
pub const INLINE_LANES = tierWidth(INLINE_TIER);
const BASELINE_LANES = tierWidth(.baseline);

pub const BLOCK_BYTES = 64;

const OBJECT_MIN_BYTES = 1024;

const AVX2_OBJECT = build_options.avx2_object;
const AVX512_OBJECT = build_options.avx512_object;
const SELECT_AT_RUN_TIME = AVX2_OBJECT or AVX512_OBJECT;

fn targetHas(comptime features: []const std.Target.x86.Feature) bool {
    for (features) |feature| {
        if (!std.Target.x86.featureSetHas(builtin.cpu.features, feature)) return false;
    }
    return true;
}

pub inline fn spanKernel(comptime name: []const u8, comptime function: anytype, len: usize) ?*const @TypeOf(function.*) {
    if (!SELECT_AT_RUN_TIME or len < OBJECT_MIN_BYTES) return null;
    return objectKernel(name, function);
}

pub inline fn callKernel(comptime name: []const u8, comptime function: anytype) ?*const @TypeOf(function.*) {
    if (!SELECT_AT_RUN_TIME) return null;
    return objectKernel(name, function);
}

inline fn objectKernel(comptime name: []const u8, comptime function: anytype) ?*const @TypeOf(function.*) {
    return switch (runTimeTier()) {
        .avx512 => if (AVX512_OBJECT) objectSymbol(.avx512, name, function) else null,
        .avx2 => if (AVX2_OBJECT) objectSymbol(.avx2, name, function) else null,
        .baseline => null,
    };
}

fn objectSymbol(comptime tier: Tier, comptime name: []const u8, comptime function: anytype) *const @TypeOf(function.*) {
    return @extern(*const @TypeOf(function.*), .{ .name = symbol(tier, name) });
}

pub const TierRange = struct { lowest: Tier, highest: Tier };

pub fn tierRange() TierRange {
    return .{ .lowest = INLINE_TIER, .highest = bestTier() };
}

pub fn selectTier(tier: Tier) error{TierUnavailable}!void {
    const range = tierRange();
    if (@intFromEnum(tier) < @intFromEnum(range.lowest) or @intFromEnum(tier) > @intFromEnum(range.highest)) {
        return error.TierUnavailable;
    }
    run_time_tier.store(tierState(tier), .monotonic);
}

pub fn exportKernel(comptime name: []const u8, comptime function: anytype) void {
    if (build_options.simd_object) @export(function, .{ .name = symbol(INLINE_TIER, name) });
}

fn symbol(comptime tier: Tier, comptime name: []const u8) []const u8 {
    return "zfq_" ++ @tagName(tier) ++ "_" ++ name;
}

comptime {
    if (build_options.simd_object) {
        _ = @import("fastq.zig");
        _ = @import("stats.zig");
        _ = @import("preprocess.zig");
    }
    exportKernel("find_newline", &findNewlineExport);
    exportKernel("collect_newlines", &collectNewlinesExport);
    exportKernel("contains_newline", &containsNewlineExport);
}

pub inline fn findNewline(bytes: []const u8) ?usize {
    if (!SELECT_AT_RUN_TIME) return @call(.always_inline, findNewlineFor, .{ INLINE_LANES, bytes });
    const head = @min(bytes.len, OBJECT_MIN_BYTES);
    if (@call(.always_inline, findNewlineFor, .{ INLINE_LANES, bytes[0..head] })) |index| return index;
    if (head == bytes.len) return null;
    const kernel = spanKernel("find_newline", &findNewlineExport, bytes.len) orelse
        return if (findNewlineFor(INLINE_LANES, bytes[head..])) |index| head + index else null;
    const rest = kernel(bytes.ptr + head, bytes.len - head);
    return if (rest == bytes.len - head) null else head + rest;
}

pub inline fn findNewlines(bytes: []const u8, ends: *[4]usize) bool {
    if (!SELECT_AT_RUN_TIME) return @call(.always_inline, collectNewlines, .{ INLINE_LANES, bytes, 0, ends, 0 }) == ends.len;
    const head = @min(bytes.len, OBJECT_MIN_BYTES);
    const count = @call(.always_inline, collectNewlines, .{ INLINE_LANES, bytes[0..head], 0, ends, 0 });
    if (count == ends.len or head == bytes.len) return count == ends.len;
    const kernel = spanKernel("collect_newlines", &collectNewlinesExport, bytes.len) orelse
        return collectNewlines(INLINE_LANES, bytes[head..], head, ends, count) == ends.len;
    return kernel(bytes.ptr + head, bytes.len - head, head, ends, count) == ends.len;
}

pub inline fn containsNewline(first: []const u8, second: ?[]const u8) bool {
    std.debug.assert(second == null or second.?.len == first.len);
    const second_ptr: ?[*]const u8 = if (second) |s| s.ptr else null;
    if (spanKernel("contains_newline", &containsNewlineExport, first.len)) |kernel| {
        return kernel(first.ptr, second_ptr, first.len);
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

var run_time_tier: std.atomic.Value(u8) = .init(0);

fn tierState(tier: Tier) u8 {
    return @as(u8, @intFromEnum(tier)) + 1;
}

fn runTimeTier() Tier {
    const state = run_time_tier.load(.monotonic);
    if (state != 0) return @enumFromInt(state - 1);
    const best = bestTier();
    run_time_tier.store(tierState(best), .monotonic);
    return best;
}

fn bestTier() Tier {
    if (AVX512_OBJECT and detectAvx512()) return .avx512;
    if (AVX2_OBJECT and detectAvx2()) return .avx2;
    return INLINE_TIER;
}

fn detectAvx2() bool {
    if (comptime builtin.cpu.arch != .x86_64) return false;
    if (cpuid(0).eax < 7) return false;
    const leaf1 = cpuid(1);
    const has_osxsave = leaf1.ecx & (1 << 27) != 0;
    const has_avx = leaf1.ecx & (1 << 28) != 0;
    if (!has_osxsave or !has_avx or xgetbv0() & 0x6 != 0x6) return false;
    return cpuid(7).ebx & (1 << 5) != 0;
}

fn detectAvx512() bool {
    if (comptime builtin.cpu.arch != .x86_64) return false;
    if (!detectAvx2()) return false;
    const leaf1_ecx = bits(&.{ 9, 12, 13, 19, 20, 22, 23, 26, 29 });
    const leaf7_ebx = bits(&.{ 3, 8, 16, 17, 28, 30, 31 });
    if (cpuid(1).ecx & leaf1_ecx != leaf1_ecx) return false;
    if (cpuid(7).ebx & leaf7_ebx != leaf7_ebx) return false;
    if (cpuid(0x80000000).eax < 0x80000001) return false;
    const extended_ecx = bits(&.{ 0, 5 });
    if (cpuid(0x80000001).ecx & extended_ecx != extended_ecx) return false;
    return xgetbv0() & 0xe6 == 0xe6;
}

fn bits(comptime positions: []const u5) u32 {
    var mask: u32 = 0;
    for (positions) |position| mask |= @as(u32, 1) << position;
    return mask;
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

pub inline fn movemask(comptime width: usize, matches: @Vector(width, bool)) @Int(.unsigned, width) {
    if (comptime builtin.cpu.arch != .aarch64 or (width != 32 and width != 64)) return @bitCast(matches);
    return neonMovemask(width, matches);
}

const NeonBytes = @Vector(16, u8);

extern fn @"llvm.aarch64.neon.addp.v16i8"(NeonBytes, NeonBytes) NeonBytes;

inline fn neonMovemask(comptime width: usize, matches: @Vector(width, bool)) @Int(.unsigned, width) {
    const weights = std.simd.repeat(width, NeonBytes{ 1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128 });
    const weighted: @Vector(width, u8) = @select(u8, matches, weights, @as(@Vector(width, u8), @splat(0)));
    const addp = @"llvm.aarch64.neon.addp.v16i8";
    const sums = if (width == 64) wide: {
        const quarters = addp(addp(neonPart(width, weighted, 0), neonPart(width, weighted, 1)), addp(neonPart(width, weighted, 2), neonPart(width, weighted, 3)));
        break :wide addp(quarters, quarters);
    } else narrow: {
        const pairs = addp(neonPart(width, weighted, 0), neonPart(width, weighted, 1));
        const quads = addp(pairs, pairs);
        break :narrow addp(quads, quads);
    };
    const lanes: @Vector(2, u64) = @bitCast(sums);
    return @truncate(lanes[0]);
}

inline fn neonPart(comptime width: usize, weighted: @Vector(width, u8), comptime index: usize) NeonBytes {
    const lanes = comptime std.simd.iota(i32, 16) + @as(@Vector(16, i32), @splat(index * 16));
    return @shuffle(u8, weighted, undefined, lanes);
}

pub inline fn newlineMask(comptime width: usize, block: *const [width]u8) @Int(.unsigned, width) {
    const Bytes = @Vector(width, u8);
    const bytes: Bytes = block.*;
    return movemask(width, bytes == @as(Bytes, @splat('\n')));
}

inline fn blockMask(comptime width: usize, block: *const [BLOCK_BYTES]u8) u64 {
    if (comptime builtin.cpu.arch == .aarch64) return newlineMask(BLOCK_BYTES, block);
    var mask: u64 = 0;
    inline for (0..BLOCK_BYTES / width) |i| {
        mask |= @as(u64, newlineMask(width, block[i * width ..][0..width])) << @intCast(i * width);
    }
    return mask;
}

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

fn expectNewlineSearch(bytes: []const u8, second: []const u8) !void {
    var want_ends: [4]usize = undefined;
    var want_count: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte != '\n' or want_count == 4) continue;
        want_ends[want_count] = index;
        want_count += 1;
    }
    const want_first: ?usize = if (want_count > 0) want_ends[0] else null;
    const want_any = want_first != null or std.mem.findScalar(u8, second, '\n') != null;
    var ends: [4]usize = undefined;
    inline for (.{ BASELINE_LANES, tierWidth(.avx2), tierWidth(.avx512) }) |lanes| {
        try std.testing.expectEqual(want_first, findNewlineFor(lanes, bytes));
        try std.testing.expectEqual(want_count == 4, findNewlinesFor(lanes, bytes, &ends));
        if (want_count == 4) try std.testing.expectEqualSlices(usize, &want_ends, &ends);
        try std.testing.expectEqual(want_first != null, containsNewlineFor(lanes, bytes.ptr, null, bytes.len));
        try std.testing.expectEqual(want_any, containsNewlineFor(lanes, bytes.ptr, second.ptr, bytes.len));
    }
    try std.testing.expectEqual(want_first, findNewline(bytes));
    try std.testing.expectEqual(want_count == 4, findNewlines(bytes, &ends));
    if (want_count == 4) try std.testing.expectEqualSlices(usize, &want_ends, &ends);
    try std.testing.expectEqual(want_any, containsNewline(bytes, second));
}

test "[property] - [newline kernels]: every search matches a scalar scan at vector, block, and skip boundaries" {
    const lengths = [_]usize{ 0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 191, 192, 193, 320, 449 };
    var storage: [449 + 2 * BLOCK_BYTES]u8 = undefined;
    var second_storage: [449 + 2 * BLOCK_BYTES]u8 = undefined;
    for ([_]usize{ 0, 1, BASELINE_LANES - 1 }) |alignment| {
        for (lengths) |len| {
            @memset(&storage, '\n');
            @memset(&second_storage, '\n');
            const bytes = storage[alignment..][0..len];
            const second = second_storage[alignment..][0..len];
            @memset(second, '\r');
            @memset(bytes, 'A');
            try expectNewlineSearch(bytes, second);
            if (len == 0) continue;
            for ([_]usize{ 0, len - 1, len / 2 }) |pos| {
                bytes[pos] = '\n';
                try expectNewlineSearch(bytes, second);
                bytes[pos] = 'A';
            }
            second[len - 1] = '\n';
            try expectNewlineSearch(bytes, second);
            second[len - 1] = '\r';
            if (len < 4) continue;
            const run: usize = @min(len - 4, BLOCK_BYTES - 2);
            for (run..run + 4) |pos| bytes[pos] = '\n';
            try expectNewlineSearch(bytes, second);
        }
    }
}

test "[property] - [masks]: movemask sets bit i exactly when lane i matches, at every width" {
    var state: u32 = 5;
    inline for (.{ 16, 32, 64 }) |width| {
        const Mask = @Int(.unsigned, width);
        for (0..200) |_| {
            var expected: Mask = 0;
            var lanes: [width]bool = undefined;
            for (&lanes, 0..) |*lane, index| {
                state = state *% 1103515245 +% 12345;
                lane.* = (state >> 16) % 3 == 0;
                if (lane.*) expected |= @as(Mask, 1) << @intCast(index);
            }
            const matches: @Vector(width, bool) = lanes;
            try std.testing.expectEqual(expected, movemask(width, matches));
        }
        try std.testing.expectEqual(@as(Mask, 0), movemask(width, @as(@Vector(width, bool), @splat(false))));
        try std.testing.expectEqual(~@as(Mask, 0), movemask(width, @as(@Vector(width, bool), @splat(true))));
    }
}

test "[edge] - [newline kernels]: masks set the bit of each LF at every width" {
    inline for (.{ 1, BASELINE_LANES, tierWidth(.avx2), BLOCK_BYTES }) |lanes| {
        const Mask = @Int(.unsigned, lanes);
        var block: [lanes]u8 = @splat('A');
        try std.testing.expectEqual(@as(Mask, 0), newlineMask(lanes, &block));
        var expected: Mask = 0;
        for ([_]usize{ 0, 1, 15, 16, 31, 32, 62, 63 }) |lane| {
            if (lane >= lanes) continue;
            block[lane] = '\n';
            expected |= @as(Mask, 1) << @intCast(lane);
        }
        try std.testing.expectEqual(expected, newlineMask(lanes, &block));
    }
}

test "[property] - [newline kernels]: the run-time selected path is exact around its length gate in every available tier" {
    const saved = run_time_tier.load(.monotonic);
    defer run_time_tier.store(saved, .monotonic);
    var storage: [5 * OBJECT_MIN_BYTES]u8 = undefined;
    const lengths = [_]usize{ OBJECT_MIN_BYTES - 1, OBJECT_MIN_BYTES, OBJECT_MIN_BYTES + 1, 2 * OBJECT_MIN_BYTES + 5 };
    const positions = [_]usize{ 0, OBJECT_MIN_BYTES - 1, OBJECT_MIN_BYTES, OBJECT_MIN_BYTES + 31, 2 * OBJECT_MIN_BYTES };
    const range = tierRange();
    for (std.enums.values(Tier)) |tier| {
        if (@intFromEnum(tier) < @intFromEnum(range.lowest) or @intFromEnum(tier) > @intFromEnum(range.highest)) continue;
        try selectTier(tier);
        for ([_]usize{ 0, 1 }) |alignment| {
            for (lengths) |len| {
                const bytes = storage[alignment..][0..len];
                for (positions) |first_lf| {
                    @memset(&storage, 'A');
                    storage[alignment + len] = '\n';
                    for ([_]usize{ first_lf, first_lf + 1, first_lf + 100, first_lf + 700 }) |index| {
                        if (index < len) bytes[index] = '\n';
                    }
                    const quiet = storage[alignment + len + 1 ..][0..len];
                    @memset(quiet, 'C');
                    try expectNewlineSearch(bytes, quiet);
                    try expectNewlineSearch(quiet, bytes);
                }
            }
        }
    }
}

test "[unit] - [tier]: selection accepts exactly the tiers from the compiled one to the CPU's best" {
    const saved = run_time_tier.load(.monotonic);
    defer run_time_tier.store(saved, .monotonic);
    const range = tierRange();
    try std.testing.expectEqual(INLINE_TIER, range.lowest);
    const best: Tier = if (AVX512_OBJECT and detectAvx512()) .avx512 else if (AVX2_OBJECT and detectAvx2()) .avx2 else INLINE_TIER;
    try std.testing.expectEqual(best, range.highest);
    for (std.enums.values(Tier)) |tier| {
        const before = runTimeTier();
        if (@intFromEnum(tier) < @intFromEnum(range.lowest) or @intFromEnum(tier) > @intFromEnum(range.highest)) {
            try std.testing.expectError(error.TierUnavailable, selectTier(tier));
            try std.testing.expectEqual(before, runTimeTier());
            continue;
        }
        try selectTier(tier);
        try std.testing.expectEqual(tier, runTimeTier());
        const object = tier != INLINE_TIER;
        try std.testing.expectEqual(object, spanKernel("find_newline", &findNewlineExport, OBJECT_MIN_BYTES) != null);
        try std.testing.expectEqual(object, callKernel("find_newline", &findNewlineExport) != null);
        try std.testing.expect(spanKernel("find_newline", &findNewlineExport, OBJECT_MIN_BYTES - 1) == null);
    }
}

test "[unit] - [cpu]: detection agrees with the AVX2 and AVX-512 guarantees of the build target" {
    if (comptime builtin.cpu.arch != .x86_64) return error.SkipZigTest;
    if (comptime targetHas(&.{.avx2})) try std.testing.expect(detectAvx2());
    if (comptime targetHas(&.{ .avx512f, .avx512bw, .avx512cd, .avx512dq, .avx512vl, .bmi, .bmi2, .fma, .f16c, .lzcnt, .movbe })) {
        try std.testing.expect(detectAvx512());
    }
    if (detectAvx512()) try std.testing.expect(detectAvx2());
}
