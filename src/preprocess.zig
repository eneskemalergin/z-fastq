//! Paired read trimming with fastp 1.1.0's rules: the 3' quality window, poly-G tails, and
//! adapters found by mate overlap or by adapter sequence. Every trim keeps a prefix of the read,
//! so a result is one kept length per mate; sequence and quality are cut together by the caller.
//! The overlap and adapter searches score a block of consecutive offsets per pass over the
//! pattern, reading the read through a fixed stack window padded with zero bytes, which no read
//! or adapter base equals. The block kernel uses the inline vector width only, with no AVX2
//! object for builds whose target lacks AVX2.

const std = @import("std");

const simd = @import("simd.zig");

pub const MIN_ADAPTER_BASES = 4;
pub const MAX_ADAPTER_BASES = 128;
pub const MAX_TAIL_WINDOW = 1000;
pub const MAX_TAIL_QUALITY = 30;

const PHRED_OFFSET = 33;
const OVERLAP_REQUIRE = 30;
const OVERLAP_DIFF_LIMIT = 5;
const OVERLAP_DIFF_DIVISOR = 5;
const OVERLAP_COMPARE_LIMIT = 50;
const ADAPTER_MATCH_REQUIRE = 4;
const ADAPTER_LEAD = 4;
const BASES_PER_MISMATCH = 8;
const POLY_G_MAX_MISMATCHES = 5;
const LANES = simd.INLINE_LANES;
const MAX_LANES = simd.tierWidth(.avx512);
const WINDOW_OFFSETS = 1024;
const WINDOW_SPAN = MAX_ADAPTER_BASES + MAX_LANES;
const WINDOW_BYTES = ADAPTER_LEAD + WINDOW_OFFSETS + WINDOW_SPAN;
const THRESHOLD_CLAMP = 1024;
const PRUNE_COLUMNS = 16;

const COMPLEMENT: [256]u8 = table: {
    var table = [_]u8{'N'} ** 256;
    for ("ACGTacgt", "TGCATGCA") |base, complement| table[base] = complement;
    break :table table;
};

pub const Tail = struct {
    window: u16,
    quality: u8,
};

pub const Options = struct {
    adapters: ?[2][]const u8 = null,
    poly_g: ?u32 = null,
    tail: ?Tail = null,

    pub fn trims(self: Options) bool {
        return self.adapters != null or self.poly_g != null or self.tail != null;
    }
};

pub const Step = enum { quality_tail, poly_g, overlap, adapter };

pub const STEP_COUNT = @typeInfo(Step).@"enum".fields.len;

pub const Overlap = struct {
    offset: isize,
    length: usize,
};

pub const Result = struct {
    kept: [2]usize,
    removed: [2][STEP_COUNT]usize = @splat(@splat(0)),
    overlap: ?Overlap = null,
    adapter_cut: [2]bool = .{ false, false },

    fn cut(self: *Result, mate: usize, step: Step, kept: usize) void {
        std.debug.assert(kept <= self.kept[mate]);
        self.removed[mate][@intFromEnum(step)] += self.kept[mate] - kept;
        self.kept[mate] = kept;
    }
};

pub const StepCounts = struct {
    reads: u64 = 0,
    bases: u64 = 0,
};

pub const MateCounts = struct {
    steps: [STEP_COUNT]StepCounts = @splat(.{}),

    pub fn add(self: *MateCounts, removed: [STEP_COUNT]usize) error{Overflow}!void {
        for (&self.steps, removed) |*counts, bases| {
            if (bases == 0) continue;
            counts.reads = try std.math.add(u64, counts.reads, 1);
            counts.bases = try std.math.add(u64, counts.bases, bases);
        }
    }
};

pub const Trimmer = struct {
    options: Options,
    adapter_splats: [2][MAX_ADAPTER_BASES]@Vector(LANES, u8) = undefined,

    pub fn init(options: Options) Trimmer {
        var trimmer: Trimmer = .{ .options = options };
        if (options.adapters) |adapters| {
            for (&trimmer.adapter_splats, adapters) |*splats, adapter| splatPattern(LANES, adapter, splats);
        }
        return trimmer;
    }

    pub fn trimPair(self: *const Trimmer, sequences: [2][]const u8, qualities: [2][]const u8) Result {
        const options = self.options;
        var result: Result = .{ .kept = .{ sequences[0].len, sequences[1].len } };
        if (options.tail) |tail| {
            for (0..2) |mate| result.cut(mate, .quality_tail, tailKept(sequences[mate], qualities[mate], tail));
            if (result.kept[0] == 0 or result.kept[1] == 0) return result;
        }
        if (options.poly_g) |min_run| {
            for (0..2) |mate| result.cut(mate, .poly_g, polyGKept(sequences[mate][0..result.kept[mate]], min_run));
        }
        const adapters = options.adapters orelse return result;
        const reads = [2][]const u8{ sequences[0][0..result.kept[0]], sequences[1][0..result.kept[1]] };
        var windows: [2]Window = .{ .{ .source = reads[0], .complement = false }, .{ .source = reads[1], .complement = false } };
        result.overlap = findOverlap(&windows[0], reads[1]);
        if (result.overlap) |overlap| if (overlap.offset < 0) {
            for (0..2) |mate| {
                result.cut(mate, .overlap, @min(result.kept[mate], overlap.length));
                result.adapter_cut[mate] = true;
            }
            return result;
        };
        for (0..2) |mate| {
            const splats = self.adapter_splats[mate][0..adapters[mate].len];
            if (adapterKept(&windows[mate], adapters[mate], splats)) |kept| {
                result.cut(mate, .adapter, kept);
                result.adapter_cut[mate] = true;
            }
        }
        return result;
    }
};

const Window = struct {
    source: []const u8,
    complement: bool,
    filled: bool = false,
    start: isize = 0,
    bytes: [WINDOW_BYTES + WINDOW_SPAN]u8 = undefined,

    fn at(self: *Window, index: isize) [*]const u8 {
        std.debug.assert(index >= -ADAPTER_LEAD);
        if (!self.filled or index < self.start or index - self.start > WINDOW_OFFSETS) self.fill(index);
        return self.bytes[@intCast(index - self.start + ADAPTER_LEAD)..].ptr;
    }

    fn fill(self: *Window, start: isize) void {
        self.filled = true;
        self.start = start;
        const first = start - ADAPTER_LEAD;
        const zeros: usize = @intCast(@max(-first, 0));
        self.bytes[0 .. 2 * ADAPTER_LEAD].* = @splat(0);
        const begin: usize = @intCast(@max(first, 0));
        const count = @min(self.source.len -| begin, WINDOW_BYTES - zeros);
        const target = self.bytes[zeros..][0..count];
        if (self.complement) {
            const last = self.source.len -% 1 -% begin;
            for (target, 0..) |*byte, index| byte.* = COMPLEMENT[self.source[last - index]];
        } else {
            @memcpy(target, self.source[begin..][0..count]);
        }
        self.bytes[zeros + count ..][0..WINDOW_SPAN].* = @splat(0);
    }
};

fn tailKept(sequence: []const u8, quality: []const u8, tail: Tail) usize {
    const length = sequence.len;
    const window: usize = tail.window;
    if (length <= window) return 0;
    const threshold = window * (PHRED_OFFSET + @as(usize, tail.quality));
    var sum: usize = 0;
    for (quality[length - window + 1 ..]) |byte| sum += byte;
    var end = length - 1;
    while (end >= window) : (end -= 1) {
        sum += quality[end - window + 1];
        if (end < length - 1) sum -= quality[end + 1];
        if (sum >= threshold) break;
    }
    if (end < length - 1) end = end + 1 - window;
    var kept = end + 1;
    while (kept > 0 and sequence[kept - 1] == 'N') kept -= 1;
    return kept;
}

fn polyGKept(sequence: []const u8, min_run: u32) usize {
    const length = sequence.len;
    var mismatches: usize = 0;
    var first_g = length;
    var scanned: usize = 0;
    while (scanned < length) : (scanned += 1) {
        const index = length - 1 - scanned;
        if (sequence[index] == 'G') first_g = index else mismatches += 1;
        const allowed = (scanned + 1) / BASES_PER_MISMATCH;
        if (mismatches > POLY_G_MAX_MISMATCHES or (mismatches > allowed and scanned + 1 >= min_run)) break;
    }
    return if (scanned >= min_run) first_g else length;
}

fn findOverlap(window1: *Window, read2: []const u8) ?Overlap {
    const read1 = window1.source;
    var head2: [OVERLAP_COMPARE_LIMIT]u8 = undefined;
    const head2_len = @min(read2.len, OVERLAP_COMPARE_LIMIT);
    for (head2[0..head2_len], 0..) |*byte, index| byte.* = COMPLEMENT[read2[read2.len - 1 - index]];
    var splats: [OVERLAP_COMPARE_LIMIT]@Vector(LANES, u8) = undefined;
    if (read1.len > OVERLAP_REQUIRE) {
        splatPattern(LANES, head2[0..head2_len], &splats);
        const rule: OverlapRule = .{ .text_len = read1.len, .other_len = read2.len };
        const end: isize = @intCast(read1.len - OVERLAP_REQUIRE);
        if (firstAccepted(LANES, rule, window1, splats[0..head2_len], 0, end)) |found| {
            const offset: usize = @intCast(found);
            return .{ .offset = found, .length = @min(read1.len - offset, read2.len) };
        }
    }
    if (read2.len > OVERLAP_REQUIRE) {
        const head1 = read1[0..@min(read1.len, OVERLAP_COMPARE_LIMIT)];
        splatPattern(LANES, head1, &splats);
        var complement2: Window = .{ .source = read2, .complement = true };
        const rule: OverlapRule = .{ .text_len = read2.len, .other_len = read1.len };
        const end: isize = @intCast(read2.len - OVERLAP_REQUIRE);
        if (firstAccepted(LANES, rule, &complement2, splats[0..head1.len], 0, end)) |found| {
            const shift: usize = @intCast(found);
            return .{ .offset = -found, .length = @min(read1.len, read2.len - shift) };
        }
    }
    return null;
}

fn adapterKept(window: *Window, adapter: []const u8, splats: []const @Vector(LANES, u8)) ?usize {
    const read = window.source;
    const read_len: isize = @intCast(read.len);
    const first: isize = if (adapter.len >= 16) -4 else if (adapter.len >= 12) -3 else if (adapter.len >= 8) -2 else 0;
    const rule: AdapterRule = .{ .read_len = read.len, .adapter_len = adapter.len };
    if (firstAccepted(LANES, rule, window, splats, first, read_len - ADAPTER_MATCH_REQUIRE)) |found| {
        return @intCast(@max(found, 0));
    }
    if (read.len > ADAPTER_MATCH_REQUIRE + 1) {
        const longest = @min(read.len - 1, adapter.len);
        const length = longestGapMatch(read, adapter, longest);
        if (length != 0) return if (length == longest) 0 else read.len - 1 - length;
    }
    if (read.len > ADAPTER_MATCH_REQUIRE) {
        const longest = @min(read.len, adapter.len - 1);
        const length = longestGapMatch(adapter, read, longest);
        if (length != 0) return if (length == longest) 0 else read.len - length;
    }
    return null;
}

const OverlapRule = struct {
    text_len: usize,
    other_len: usize,

    fn need(self: OverlapRule, comptime lanes: comptime_int, offset: isize) @Vector(lanes, u8) {
        const last_remaining = @as(isize, @intCast(self.text_len)) - offset - (lanes - 1);
        if (last_remaining >= @min(self.other_len, OVERLAP_COMPARE_LIMIT)) {
            return @splat(@intCast(overlapNeed(self.other_len)));
        }
        const Wide = @Vector(lanes, i16);
        const remaining: i16 = @intCast(std.math.clamp(@as(isize, @intCast(self.text_len)) - offset, 0, THRESHOLD_CLAMP));
        const other: i16 = @intCast(@min(self.other_len, THRESHOLD_CLAMP));
        const overlap = @max(@min(@as(Wide, @splat(remaining)) - std.simd.iota(i16, lanes), @as(Wide, @splat(other))), @as(Wide, @splat(0)));
        const compared = @min(overlap, @as(Wide, @splat(OVERLAP_COMPARE_LIMIT)));
        const limit = @min(@as(Wide, @splat(OVERLAP_DIFF_LIMIT)), @divTrunc(overlap, @as(Wide, @splat(OVERLAP_DIFF_DIVISOR))));
        return @intCast(compared - limit);
    }

    fn needAt(self: OverlapRule, offset: isize) usize {
        return overlapNeed(@min(self.text_len - @as(usize, @intCast(offset)), self.other_len));
    }
};

fn overlapNeed(overlap: usize) usize {
    return @min(overlap, OVERLAP_COMPARE_LIMIT) - @min(OVERLAP_DIFF_LIMIT, overlap / OVERLAP_DIFF_DIVISOR);
}

const AdapterRule = struct {
    read_len: usize,
    adapter_len: usize,

    fn need(self: AdapterRule, comptime lanes: comptime_int, position: isize) @Vector(lanes, u8) {
        const last_remaining = @as(isize, @intCast(self.read_len)) - position - (lanes - 1);
        if (position >= 0 and last_remaining >= self.adapter_len) {
            return @splat(@intCast(self.adapter_len - self.adapter_len / BASES_PER_MISMATCH));
        }
        const Wide = @Vector(lanes, i16);
        const iota = std.simd.iota(i16, lanes);
        const zero: Wide = @splat(0);
        const remaining: i16 = @intCast(std.math.clamp(@as(isize, @intCast(self.read_len)) - position, 0, THRESHOLD_CLAMP));
        const start: i16 = @intCast(@min(position, THRESHOLD_CLAMP));
        const compared = @max(@min(@as(Wide, @splat(remaining)) - iota, @as(Wide, @splat(@intCast(self.adapter_len)))), zero);
        const skipped = @max(-(@as(Wide, @splat(start)) + iota), zero);
        const allowed = @divTrunc(compared, @as(Wide, @splat(BASES_PER_MISMATCH)));
        return @intCast(@max(compared - skipped - allowed, zero));
    }

    fn needAt(self: AdapterRule, position: isize) usize {
        const compared: usize = @intCast(@min(@as(isize, @intCast(self.read_len)) - position, @as(isize, @intCast(self.adapter_len))));
        const skipped: usize = @intCast(@max(-position, 0));
        return compared - skipped - compared / BASES_PER_MISMATCH;
    }
};

fn splatPattern(comptime lanes: comptime_int, pattern: []const u8, splats: []@Vector(lanes, u8)) void {
    for (pattern, splats[0..pattern.len]) |base, *splat| splat.* = @splat(base);
}

fn firstAccepted(
    comptime lanes: comptime_int,
    rule: anytype,
    window: *Window,
    splats: []const @Vector(lanes, u8),
    first: isize,
    end: isize,
) ?isize {
    const Counts = @Vector(lanes, u8);
    const Mask = @Int(.unsigned, lanes);
    const one: Counts = @splat(1);
    const zero: Counts = @splat(0);
    var offset = first;
    while (offset < end) : (offset += lanes) {
        const text = window.at(offset);
        const need = rule.need(lanes, offset);
        var matches = zero;
        var column: usize = 0;
        const reachable = while (column < splats.len) {
            const stop = @min(column + PRUNE_COLUMNS, splats.len);
            for (splats[column..stop], column..) |splat, index| {
                const bases: Counts = text[index..][0..lanes].*;
                matches += @select(u8, bases == splat, one, zero);
            }
            column = stop;
            const remaining: Counts = @splat(@intCast(splats.len - column));
            if (!@reduce(.Or, matches +| remaining >= need)) break false;
        } else true;
        if (!reachable) continue;
        const valid = std.math.shr(Mask, std.math.maxInt(Mask), lanes - @as(usize, @intCast(@min(end - offset, lanes))));
        const accepted: Mask = @bitCast(matches >= need);
        if (accepted & valid != 0) return offset + @ctz(accepted & valid);
    }
    return null;
}

fn longestGapMatch(inserted: []const u8, normal: []const u8, longest: usize) usize {
    std.debug.assert(longest <= MAX_ADAPTER_BASES and inserted.len > longest and normal.len >= longest);
    var left: i16 = 0;
    var shifted: i16 = 0;
    var lowest: i16 = std.math.maxInt(i16);
    var found: usize = 0;
    const most_allowed = @as(i16, @intCast(longest / BASES_PER_MISMATCH)) - 1;
    for (1..longest + 1) |length| {
        left += @intFromBool(inserted[length - 1] != normal[length - 1]);
        shifted += @intFromBool(inserted[length] != normal[length - 1]);
        const allowed: i16 = @intCast(length / BASES_PER_MISMATCH);
        const score = lowest +| shifted;
        if (allowed >= 1 and score <= allowed - 1) found = length;
        if (left > most_allowed and score > most_allowed) break;
        lowest = @min(lowest, left - shifted);
    }
    return found;
}

fn firstAcceptedScalar(rule: anytype, text: []const u8, pattern: []const u8, first: isize, end: isize) ?isize {
    var offset = first;
    while (offset < end) : (offset += 1) {
        var matches: usize = 0;
        for (pattern, 0..) |base, index| {
            const at = offset + @as(isize, @intCast(index));
            if (at >= 0 and at < text.len and text[@intCast(at)] == base) matches += 1;
        }
        if (matches >= rule.needAt(offset)) return offset;
    }
    return null;
}

test "[property] - [preprocess]: block searches find the first offset the scalar reference finds" {
    var text: [2600]u8 = undefined;
    var state: u32 = 1;
    for (&text) |*base| {
        state = state *% 1103515245 +% 12345;
        base.* = "ACGT"[(state >> 16) & 3];
    }
    var pattern: [MAX_ADAPTER_BASES]u8 = undefined;
    const Case = struct { length: usize, from: usize, mutations: usize, size: usize, complement: bool };
    const cases = [_]Case{
        .{ .length = 50, .from = 0, .mutations = 0, .size = 200, .complement = false },
        .{ .length = 50, .from = 117, .mutations = 5, .size = 200, .complement = false },
        .{ .length = 50, .from = 117, .mutations = 6, .size = 200, .complement = false },
        .{ .length = 24, .from = 160, .mutations = 2, .size = 200, .complement = true },
        .{ .length = 50, .from = 2400, .mutations = 3, .size = 2600, .complement = false },
        .{ .length = 33, .from = 190, .mutations = 1, .size = 200, .complement = false },
        .{ .length = 128, .from = 1500, .mutations = 16, .size = 2600, .complement = true },
        .{ .length = 19, .from = 0, .mutations = 0, .size = 7, .complement = false },
    };
    for (cases, 0..) |case, case_index| {
        errdefer std.debug.print("case {d}\n", .{case_index});
        const source = text[0..case.size];
        var window_text: [2600]u8 = undefined;
        for (window_text[0..case.size], 0..) |*byte, index| {
            byte.* = if (case.complement) COMPLEMENT[source[case.size - 1 - index]] else source[index];
        }
        const searched = window_text[0..case.size];
        for (pattern[0..case.length], 0..) |*base, index| {
            const at = case.from + index;
            base.* = if (at < searched.len) searched[at] else 'A';
            if (index * 7 % case.length < case.mutations) base.* = if (base.* == 'A') 'C' else 'A';
        }
        const overlap: OverlapRule = .{ .text_len = case.size, .other_len = case.length };
        const adapter: AdapterRule = .{ .read_len = case.size, .adapter_len = case.length };
        const overlap_end: isize = @as(isize, @intCast(case.size)) - OVERLAP_REQUIRE;
        const adapter_end: isize = @as(isize, @intCast(case.size)) - ADAPTER_MATCH_REQUIRE;
        const expected_overlap = firstAcceptedScalar(overlap, searched, pattern[0..@min(case.length, OVERLAP_COMPARE_LIMIT)], 0, overlap_end);
        const expected_adapter = firstAcceptedScalar(adapter, searched, pattern[0..case.length], -4, adapter_end);
        inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |lanes| {
            var splats: [MAX_ADAPTER_BASES]@Vector(lanes, u8) = undefined;
            var window: Window = .{ .source = source, .complement = case.complement };
            splatPattern(lanes, pattern[0..case.length], &splats);
            try std.testing.expectEqual(expected_adapter, firstAccepted(lanes, adapter, &window, splats[0..case.length], -4, adapter_end));
            const overlap_length = @min(case.length, OVERLAP_COMPARE_LIMIT);
            try std.testing.expectEqual(expected_overlap, firstAccepted(lanes, overlap, &window, splats[0..overlap_length], 0, overlap_end));
        }
    }
}
