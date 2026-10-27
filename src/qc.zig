//! Per-mate FASTQ quality distributions for `process`, computed as FastQC 0.13.0 computes its core
//! modules: per-cycle quality and base content, per-read quality and GC, read lengths, terminal
//! homopolymer runs, and Adapter Content for FastQC's bundled adapter list, where a read exactly
//! one adapter long counts as FastQC-Rust counts it. The retained population of a trimmed or
//! filtered run counts only the removed 3' ends and dropped reads per cycle and takes its
//! per-cycle tables as raw minus removed, until raw drops its tables for a read past the cycle
//! limit; from then on it counts its reads in full. Adapter matches are searched once per raw
//! read and count for the retained read when they lie inside its kept bases.

const std = @import("std");
const fastq = @import("fastq.zig");
const simd = @import("simd.zig");
const statistics = @import("stats.zig");

pub const MAX_CYCLES = 1024;
pub const TERMINAL_RUN_MIN = 10;
pub const TOO_LONG = "a read is longer than 1024 bases";

const PHRED_OFFSET = fastq.PHRED_OFFSET;
const QUALITY_BINS = fastq.PHRED_MAX_SCORE + 1;
const GC_BINS = 101;
const PERCENTILE_MIN_READS = 100;
const GROUPING_MIN_LENGTH = 75;
const GROUP_COUNT_LIMIT = 75;
const SINGLE_BASE_GROUPS = 9;
const LENGTH_BIN_LIMIT = 50;

const BaseKind = enum(u3) { a, c, g, t, n, other };
const COUNTED_LETTERS = [_]u8{ 'A', 'C', 'G', 'T', 'N' };
const COUNTED_KINDS = COUNTED_LETTERS.len;
const QualityCounts = [QUALITY_BINS]u64;
const BaseCounts = [COUNTED_KINDS]u64;
const PendingBases = [COUNTED_KINDS][MAX_CYCLES]u8;
const PENDING_READS_MAX = std.math.maxInt(u8);
const CASE_FOLD = 0xdf;

const BASE_KIND: [256]BaseKind = table: {
    var kinds: [256]BaseKind = @splat(.other);
    for ([_]struct { u8, BaseKind }{ .{ 'A', .a }, .{ 'C', .c }, .{ 'G', .g }, .{ 'T', .t }, .{ 'N', .n } }) |entry| {
        kinds[entry[0]] = entry[1];
        kinds[std.ascii.toLower(entry[0])] = entry[1];
    }
    break :table kinds;
};

const Adapter = struct {
    name: []const u8,
    sequence: *const [ADAPTER_LENGTH]u8,
};

const ADAPTER_LENGTH = 12;

const ADAPTERS = [_]Adapter{
    .{ .name = "Illumina Universal Adapter", .sequence = "AGATCGGAAGAG" },
    .{ .name = "Illumina Small RNA 3' Adapter", .sequence = "TGGAATTCTCGG" },
    .{ .name = "Illumina Small RNA 5' Adapter", .sequence = "GATCGTCGGACT" },
    .{ .name = "Nextera Transposase Sequence", .sequence = "CTGTCTCTTATA" },
    .{ .name = "PolyA", .sequence = "AAAAAAAAAAAA" },
    .{ .name = "PolyG", .sequence = "GGGGGGGGGGGG" },
};

const ADAPTER_COUNT = ADAPTERS.len;
const NO_MATCH = std.math.maxInt(u16);
const AdapterStarts = [ADAPTER_COUNT]u16;

const ReadFacts = struct {
    quality_sum: u64 = 0,
    gc_count: usize = 0,
    adapter_starts: AdapterStarts = @splat(NO_MATCH),
};

pub const MatePopulations = struct {
    raw: MateQc = .{},
    retained: MateQc = .{},
    derived: bool = true,
    retained_cycles: usize = 0,
    last: ReadFacts = .{},

    pub fn deinit(self: *MatePopulations, allocator: std.mem.Allocator) void {
        self.raw.deinit(allocator);
        self.retained.deinit(allocator);
    }

    pub fn addRead(
        self: *MatePopulations,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
    ) (error{OutOfMemory} || statistics.StatsError)!void {
        if (self.derived and !self.raw.too_long and sequence.len > MAX_CYCLES) self.materialize();
        self.last = try self.raw.collect(allocator, sequence, quality);
    }

    pub fn addRetained(
        self: *MatePopulations,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
        kept: usize,
    ) (error{OutOfMemory} || statistics.StatsError)!void {
        std.debug.assert(kept <= sequence.len and sequence.len == quality.len);
        if (!self.derived) return self.retained.addRead(allocator, sequence[0..kept], quality[0..kept]);
        const retained = &self.retained;
        var facts = self.last;
        if (kept < sequence.len) {
            const removed_sum = try self.removeCycles(allocator, sequence, quality, kept);
            facts.quality_sum -= removed_sum;
            facts.gc_count = gcCountScalar(sequence[0..gcKeptLength(kept)]);
        } else if (sequence.len > retained.quality.items.len) {
            try retained.growCycles(allocator, sequence.len);
        }
        retained.lengths[kept] += 1;
        self.retained_cycles = @max(self.retained_cycles, kept);
        for (facts.adapter_starts, &retained.adapter_starts) |start, *starts| {
            if (start != NO_MATCH and start + ADAPTER_LENGTH <= kept) starts[start] += 1;
        }
        if (kept != 0) retained.addSummary(sequence[0..kept], facts);
    }

    pub fn dropRetained(
        self: *MatePopulations,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
    ) error{OutOfMemory}!void {
        std.debug.assert(sequence.len == quality.len);
        if (self.derived) _ = try self.removeCycles(allocator, sequence, quality, 0);
    }

    fn removeCycles(
        self: *MatePopulations,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
        from: usize,
    ) error{OutOfMemory}!u64 {
        const retained = &self.retained;
        if (sequence.len > retained.quality.items.len) try retained.growCycles(allocator, sequence.len);
        var removed_sum: u64 = 0;
        for (
            sequence[from..],
            quality[from..],
            retained.quality.items[from..sequence.len],
            retained.bases.items[from..sequence.len],
        ) |base, byte, *qualities, *bases| {
            qualities[byte - PHRED_OFFSET] += 1;
            removed_sum += byte;
            const kind = BASE_KIND[base];
            if (kind != .other) bases[@intFromEnum(kind)] += 1;
        }
        return removed_sum;
    }

    pub fn finish(self: *MatePopulations) void {
        if (self.derived) self.materialize();
        self.raw.finish();
        self.retained.finish();
    }

    fn materialize(self: *MatePopulations) void {
        std.debug.assert(self.derived and !self.raw.too_long);
        self.raw.finish();
        const retained = &self.retained;
        const cycles = self.retained_cycles;
        std.debug.assert(cycles <= retained.quality.items.len);
        retained.quality.shrinkRetainingCapacity(cycles);
        retained.bases.shrinkRetainingCapacity(cycles);
        for (
            retained.quality.items,
            retained.bases.items,
            self.raw.quality.items[0..cycles],
            self.raw.bases.items[0..cycles],
        ) |*qualities, *bases, raw_qualities, raw_bases| {
            for (qualities, raw_qualities) |*count, total| count.* = total - count.*;
            for (bases, raw_bases) |*count, total| count.* = total - count.*;
        }
        self.derived = false;
    }
};

pub const MateQc = struct {
    quality: std.ArrayList(QualityCounts) = .empty,
    bases: std.ArrayList(BaseCounts) = .empty,
    pending_bases: PendingBases = @splat(@splat(0)),
    pending_reads: u8 = 0,
    lengths: [MAX_CYCLES + 1]u64 = @splat(0),
    read_quality: QualityCounts = @splat(0),
    gc: [GC_BINS]f64 = @splat(0),
    gc_model: GcModel = .{},
    terminal_g: u64 = 0,
    terminal_a: u64 = 0,
    adapter_starts: [ADAPTER_COUNT][MAX_CYCLES]u64 = @splat(@splat(0)),
    too_long: bool = false,
    long_stats: statistics.Stats = .{},

    pub fn deinit(self: *MateQc, allocator: std.mem.Allocator) void {
        self.quality.deinit(allocator);
        self.bases.deinit(allocator);
    }

    pub fn finish(self: *MateQc) void {
        for (self.bases.items, 0..) |*totals, cycle| {
            for (totals, &self.pending_bases) |*total, *pending| {
                total.* += pending[cycle];
                pending[cycle] = 0;
            }
        }
        self.pending_reads = 0;
    }

    pub fn addRead(
        self: *MateQc,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
    ) (error{OutOfMemory} || statistics.StatsError)!void {
        _ = try self.collect(allocator, sequence, quality);
    }

    fn collect(
        self: *MateQc,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
    ) (error{OutOfMemory} || statistics.StatsError)!ReadFacts {
        std.debug.assert(sequence.len == quality.len);
        if (!self.too_long and sequence.len > MAX_CYCLES) {
            self.finish();
            self.long_stats = try self.tableStats();
            self.too_long = true;
            self.quality.clearAndFree(allocator);
            self.bases.clearAndFree(allocator);
        }
        if (self.too_long) try self.long_stats.addRecord(.{ .header = "", .id = "", .sequence = sequence, .plus = "", .quality = quality });
        const gc_length = gcKeptLength(sequence.len);
        var quality_sum: u64 = 0;
        const gc_count = if (self.too_long) too_long: {
            for (quality) |byte| quality_sum += byte;
            break :too_long gcCountScalar(sequence[0..gc_length]);
        } else try self.addCycles(allocator, sequence, quality, gc_length, &quality_sum);
        var facts: ReadFacts = .{ .quality_sum = quality_sum, .gc_count = gc_count };
        if (!self.too_long) {
            facts.adapter_starts = firstAdapterMatches(sequence);
            for (facts.adapter_starts, &self.adapter_starts) |start, *starts| {
                if (start != NO_MATCH) starts[start] += 1;
            }
        }
        if (sequence.len != 0) self.addSummary(sequence, facts);
        return facts;
    }

    fn addSummary(self: *MateQc, sequence: []const u8, facts: ReadFacts) void {
        std.debug.assert(sequence.len != 0);
        const gc_length = gcKeptLength(sequence.len);
        self.read_quality[facts.quality_sum / sequence.len - PHRED_OFFSET] += 1;
        const weights = self.gc_model.weights(gc_length);
        const bins = gcBins(gc_length, facts.gc_count);
        for (bins[0]..bins[1] + 1) |bin| self.gc[bin] += weights[bin];
        self.terminal_g += @intFromBool(endsInRun(sequence, 'G'));
        self.terminal_a += @intFromBool(endsInRun(sequence, 'A'));
    }

    fn addCycles(
        self: *MateQc,
        allocator: std.mem.Allocator,
        sequence: []const u8,
        quality: []const u8,
        gc_length: usize,
        quality_sum: *u64,
    ) error{OutOfMemory}!usize {
        if (sequence.len > self.quality.items.len) try self.growCycles(allocator, sequence.len);
        self.lengths[sequence.len] += 1;
        var sum: u64 = 0;
        for (quality, self.quality.items[0..sequence.len]) |byte, *qualities| {
            qualities[byte - PHRED_OFFSET] += 1;
            sum += byte;
        }
        quality_sum.* = sum;
        const gc_count = countCycleBases(simd.INLINE_LANES, &self.pending_bases, sequence, gc_length);
        self.pending_reads += 1;
        if (self.pending_reads == PENDING_READS_MAX) self.finish();
        return gc_count;
    }

    fn growCycles(self: *MateQc, allocator: std.mem.Allocator, cycles: usize) error{OutOfMemory}!void {
        const capacity = @min(MAX_CYCLES, @max(cycles, self.quality.items.len * 2));
        try self.quality.ensureTotalCapacityPrecise(allocator, capacity);
        try self.bases.ensureTotalCapacityPrecise(allocator, capacity);
        self.quality.appendNTimesAssumeCapacity(@splat(0), cycles - self.quality.items.len);
        self.bases.appendNTimesAssumeCapacity(@splat(0), cycles - self.bases.items.len);
    }

    pub fn stats(self: *const MateQc) statistics.StatsError!statistics.Stats {
        std.debug.assert(self.pending_reads == 0);
        if (self.too_long) return self.long_stats;
        return self.tableStats();
    }

    fn tableStats(self: *const MateQc) statistics.StatsError!statistics.Stats {
        const add = std.math.add;
        var result: statistics.Stats = .{};
        for (self.lengths, 0..) |reads, length| {
            if (reads == 0) continue;
            if (result.reads == 0) result.min_length = length;
            result.max_length = length;
            result.reads = try add(u64, result.reads, reads);
            result.bases = try add(u64, result.bases, try std.math.mul(u64, reads, length));
        }
        var counted: u64 = 0;
        for (self.bases.items) |counts| {
            for ([_]*u64{ &result.a, &result.c, &result.g, &result.t, &result.n }, counts) |total, count| {
                total.* = try add(u64, total.*, count);
                counted = try add(u64, counted, count);
            }
        }
        result.other_bases = result.bases - counted;
        for (self.quality.items) |counts| {
            for (counts, 0..) |count, score| {
                result.quality_sum = try add(u64, result.quality_sum, try std.math.mul(u64, count, score));
                if (score >= statistics.Q20_SCORE) result.q20_bases = try add(u64, result.q20_bases, count);
                if (score >= statistics.Q30_SCORE) result.q30_bases = try add(u64, result.q30_bases, count);
            }
        }
        return result;
    }

    pub fn writeJson(self: *const MateQc, json: *std.json.Stringify) !void {
        std.debug.assert(self.pending_reads == 0);
        try json.beginObject();
        const cycles_known = !self.too_long;
        const groups = Groups.init(self.quality.items.len);

        try json.objectField("per_base_quality");
        if (cycles_known) try self.writeBaseQuality(json, groups) else try json.write(null);
        try json.objectField("per_base_content");
        if (cycles_known) try self.writeBaseContent(json, groups) else try json.write(null);
        try json.objectField("per_base_n");
        if (cycles_known) try self.writeBaseN(json, groups) else try json.write(null);
        try json.objectField("length_distribution");
        if (cycles_known) try self.writeLengths(json) else try json.write(null);
        try json.objectField("per_sequence_quality");
        try self.writeReadQuality(json);
        try json.objectField("per_sequence_gc");
        try json.beginArray();
        for (self.gc) |value| try json.write(value);
        try json.endArray();
        try json.objectField("terminal_runs");
        try json.beginObject();
        try json.objectField("min_length");
        try json.write(TERMINAL_RUN_MIN);
        try json.objectField("g");
        try json.write(self.terminal_g);
        try json.objectField("a");
        try json.write(self.terminal_a);
        try json.endObject();
        try json.objectField("adapter_content");
        if (cycles_known) try self.writeAdapterContent(json) else try json.write(null);
        try json.objectField("unavailable");
        if (cycles_known) try json.write(null) else try json.write(TOO_LONG);
        try json.endObject();
    }

    fn writeAdapterContent(self: *const MateQc, json: *std.json.Stringify) !void {
        try json.beginObject();
        try json.objectField("adapters");
        try json.beginArray();
        for (ADAPTERS) |adapter| try json.write(adapter.name);
        try json.endArray();
        try json.objectField("rows");
        try json.beginArray();
        var reads: u64 = 0;
        for (self.lengths) |count| reads += count;
        const longest = std.mem.lastIndexOfNone(u64, &self.lengths, &.{0}) orelse 0;
        if (reads != 0 and longest > ADAPTER_LENGTH) {
            const total: f64 = @floatFromInt(reads);
            var running: [ADAPTER_COUNT]u64 = @splat(0);
            var position: usize = 0;
            var groups = Groups.init(longest - ADAPTER_LENGTH + 1);
            while (groups.next()) |group| {
                var sums: [ADAPTER_COUNT]f64 = @splat(0);
                while (position < group.last) : (position += 1) {
                    for (&running, &sums, &self.adapter_starts) |*count, *sum, *starts| {
                        count.* += starts[position];
                        sum.* += @as(f64, @floatFromInt(count.*)) * 100.0 / total;
                    }
                }
                const width: f64 = @floatFromInt(group.last - group.first);
                try json.beginObject();
                try json.objectField("position");
                try writeRangeLabel(json, group.first + 1, group.last);
                try json.objectField("values");
                try json.beginArray();
                for (sums) |sum| try json.write(sum / width);
                try json.endArray();
                try json.endObject();
            }
        }
        try json.endArray();
        try json.endObject();
    }

    fn writeBaseQuality(self: *const MateQc, json: *std.json.Stringify, groups: Groups) !void {
        try json.beginArray();
        var iterator = groups;
        while (iterator.next()) |group| {
            try json.beginObject();
            try writeGroupLabel(json, group);
            try json.objectField("mean");
            try json.write(self.groupMean(group));
            for ([_]struct { []const u8, u7 }{
                .{ "median", 50 },
                .{ "lower_quartile", 25 },
                .{ "upper_quartile", 75 },
                .{ "percentile_10", 10 },
                .{ "percentile_90", 90 },
            }) |entry| {
                try json.objectField(entry[0]);
                try json.write(self.groupPercentile(group, entry[1]));
            }
            try json.endObject();
        }
        try json.endArray();
    }

    fn groupMean(self: *const MateQc, group: Group) f64 {
        var sum: f64 = 0;
        var cycles: u64 = 0;
        for (self.quality.items[group.first..group.last]) |counts| {
            var reads: u64 = 0;
            var scores: u64 = 0;
            for (counts, 0..) |count, score| {
                reads += count;
                scores += count * score;
            }
            if (reads == 0) continue;
            sum += @as(f64, @floatFromInt(scores)) / @as(f64, @floatFromInt(reads));
            cycles += 1;
        }
        if (cycles == 0) return 0;
        return sum / @as(f64, @floatFromInt(cycles));
    }

    fn groupPercentile(self: *const MateQc, group: Group, percent: u7) ?f64 {
        var sum: f64 = 0;
        var cycles: u64 = 0;
        for (self.quality.items[group.first..group.last]) |counts| {
            var reads: u64 = 0;
            for (counts) |count| reads += count;
            if (reads <= PERCENTILE_MIN_READS) continue;
            sum += @floatFromInt(cyclePercentile(counts, reads, percent));
            cycles += 1;
        }
        if (cycles == 0) return null;
        return sum / @as(f64, @floatFromInt(cycles));
    }

    fn pooledBases(self: *const MateQc, group: Group) BaseCounts {
        var pooled: BaseCounts = @splat(0);
        for (self.bases.items[group.first..group.last]) |counts| {
            for (&pooled, counts) |*total, count| total.* += count;
        }
        return pooled;
    }

    fn pooledReads(self: *const MateQc, group: Group) u64 {
        var reads: u64 = 0;
        for (self.quality.items[group.first..group.last]) |counts| {
            for (counts) |count| reads += count;
        }
        return reads;
    }

    fn writeBaseContent(self: *const MateQc, json: *std.json.Stringify, groups: Groups) !void {
        try json.beginArray();
        var iterator = groups;
        while (iterator.next()) |group| {
            const pooled = self.pooledBases(group);
            const total = pooled[@intFromEnum(BaseKind.a)] + pooled[@intFromEnum(BaseKind.c)] +
                pooled[@intFromEnum(BaseKind.g)] + pooled[@intFromEnum(BaseKind.t)];
            try json.beginObject();
            try writeGroupLabel(json, group);
            for ([_]struct { []const u8, BaseKind }{ .{ "g", .g }, .{ "a", .a }, .{ "t", .t }, .{ "c", .c } }) |entry| {
                try json.objectField(entry[0]);
                try json.write(percentOf(pooled[@intFromEnum(entry[1])], total));
            }
            try json.endObject();
        }
        try json.endArray();
    }

    fn writeBaseN(self: *const MateQc, json: *std.json.Stringify, groups: Groups) !void {
        try json.beginArray();
        var iterator = groups;
        while (iterator.next()) |group| {
            const pooled = self.pooledBases(group);
            const total = self.pooledReads(group);
            try json.beginObject();
            try writeGroupLabel(json, group);
            try json.objectField("n");
            try json.write(percentOf(pooled[@intFromEnum(BaseKind.n)], total));
            try json.endObject();
        }
        try json.endArray();
    }

    fn writeLengths(self: *const MateQc, json: *std.json.Stringify) !void {
        try json.beginArray();
        const first = std.mem.indexOfNone(u64, &self.lengths, &.{0}) orelse 0;
        const last = std.mem.lastIndexOfNone(u64, &self.lengths, &.{0}) orelse 0;
        const low = first -| 1;
        const high = last + 1;
        const size = lengthBinSize(high - low);
        var start = low / size * size;
        while (start <= high) : (start += size) {
            const end = @min(start + size - 1, high);
            var count: u64 = 0;
            for (start..@min(end, MAX_CYCLES) + 1) |length| count += self.lengths[length];
            const edge = start == low / size * size or start + size > high;
            if (edge and count == 0) continue;
            try json.beginObject();
            try json.objectField("length");
            try writeRangeLabel(json, start, if (size == 1) start else end);
            try json.objectField("count");
            try json.write(count);
            try json.endObject();
        }
        try json.endArray();
    }

    fn writeReadQuality(self: *const MateQc, json: *std.json.Stringify) !void {
        try json.beginArray();
        if (std.mem.indexOfNone(u64, &self.read_quality, &.{0})) |lowest| {
            const highest = std.mem.lastIndexOfNone(u64, &self.read_quality, &.{0}).?;
            for (lowest..highest + 1) |score| {
                try json.beginObject();
                try json.objectField("quality");
                try json.write(score);
                try json.objectField("count");
                try json.write(self.read_quality[score]);
                try json.endObject();
            }
        }
        try json.endArray();
    }
};

fn countCycleBases(comptime lanes: comptime_int, pending: *PendingBases, sequence: []const u8, gc_length: usize) usize {
    var gc_count: usize = 0;
    var cycle: usize = 0;
    while (cycle + lanes <= sequence.len) : (cycle += lanes) {
        gc_count += countBlock(lanes, pending, sequence, cycle, gc_length);
    }
    for (sequence[cycle..], cycle..) |base, at| {
        const kind = BASE_KIND[base];
        if (kind == .other) continue;
        pending[@intFromEnum(kind)][at] += 1;
        gc_count += @intFromBool(at < gc_length and (kind == .g or kind == .c));
    }
    return gc_count;
}

inline fn countBlock(
    comptime lanes: comptime_int,
    pending: *PendingBases,
    sequence: []const u8,
    start: usize,
    gc_length: usize,
) usize {
    const Bytes = @Vector(lanes, u8);
    const ones: Bytes = @splat(1);
    const zeros: Bytes = @splat(0);
    const upper = @as(Bytes, sequence[start..][0..lanes].*) & @as(Bytes, @splat(CASE_FOLD));
    var gc_count: usize = 0;
    inline for (COUNTED_LETTERS, pending) |letter, *row| {
        const hits = upper == @as(Bytes, @splat(letter));
        const counts = row[start..][0..lanes];
        counts.* = @as(Bytes, counts.*) + @select(u8, hits, ones, zeros);
        if (letter == 'G' or letter == 'C') gc_count += keptHits(lanes, hits, start, gc_length);
    }
    return gc_count;
}

fn keptHits(comptime lanes: comptime_int, hits: @Vector(lanes, bool), cycle: usize, gc_length: usize) usize {
    const Mask = @Int(.unsigned, lanes);
    if (cycle >= gc_length) return 0;
    const mask: Mask = @bitCast(hits);
    const kept = gc_length - cycle;
    if (kept >= lanes) return @popCount(mask);
    return @popCount(mask & ((@as(Mask, 1) << @intCast(kept)) - 1));
}

fn gcCountScalar(sequence: []const u8) usize {
    var gc_count: usize = 0;
    for (sequence) |base| {
        const kind = BASE_KIND[base];
        gc_count += @intFromBool(kind == .g or kind == .c);
    }
    return gc_count;
}

fn cyclePercentile(counts: QualityCounts, reads: u64, percent: u7) usize {
    const target = reads * percent / 100;
    var seen: u64 = 0;
    for (counts, 0..) |count, score| {
        seen += count;
        if (seen >= target) return score;
    }
    return QUALITY_BINS - 1;
}

fn percentOf(count: u64, total: u64) ?f64 {
    if (total == 0) return null;
    return @as(f64, @floatFromInt(count)) / @as(f64, @floatFromInt(total)) * 100.0;
}

const WINDOW_BYTES = 8;
const WINDOW_STRIDE = ADAPTER_LENGTH - WINDOW_BYTES;
const WINDOW_FOLD: u64 = 0x0101_0101_0101_0101 * CASE_FOLD;
const WINDOW_SLOT_BITS = 7;

const WindowSlot = struct {
    key: u64 = 0,
    adapters: u8 = 0,
    offsets: [ADAPTER_COUNT]u8 = @splat(0),
};

fn windowKey(bytes: *const [WINDOW_BYTES]u8) u64 {
    return std.mem.readInt(u64, bytes, .little) & WINDOW_FOLD;
}

fn windowSlot(key: u64, multiplier: u64) usize {
    return @intCast((key *% multiplier) >> (64 - WINDOW_SLOT_BITS));
}

const WindowTable = struct {
    multiplier: u64,
    slots: [1 << WINDOW_SLOT_BITS]WindowSlot,
};

const WINDOW_TABLE: WindowTable = table: {
    @setEvalBranchQuota(1_000_000);
    var multiplier: u64 = 0x9e37_79b9_7f4a_7c15;
    while (true) : (multiplier +%= 2) {
        var slots: [1 << WINDOW_SLOT_BITS]WindowSlot = @splat(.{});
        var clean = true;
        for (ADAPTERS, 0..) |adapter, index| {
            for (0..WINDOW_STRIDE) |offset| {
                const key = windowKey(adapter.sequence[offset..][0..WINDOW_BYTES]);
                const slot = &slots[windowSlot(key, multiplier)];
                if (slot.key != 0 and slot.key != key) clean = false;
                slot.key = key;
                slot.adapters |= 1 << index;
                slot.offsets[index] |= 1 << offset;
            }
        }
        if (clean) {
            const result: WindowTable = .{ .multiplier = multiplier, .slots = slots };
            break :table result;
        }
    }
};

fn firstAdapterMatches(sequence: []const u8) AdapterStarts {
    var starts: AdapterStarts = @splat(NO_MATCH);
    if (sequence.len < ADAPTER_LENGTH) return starts;
    var found: u8 = 0;
    var window: usize = 0;
    while (window + WINDOW_BYTES <= sequence.len) : (window += WINDOW_STRIDE) {
        const key = windowKey(sequence[window..][0..WINDOW_BYTES]);
        const slot = &WINDOW_TABLE.slots[windowSlot(key, WINDOW_TABLE.multiplier)];
        if (slot.key != key or slot.adapters & ~found == 0) continue;
        for (slot.offsets, 0..) |offsets, adapter| {
            const bit = @as(u8, 1) << @intCast(adapter);
            if (offsets == 0 or found & bit != 0) continue;
            var offset: usize = WINDOW_STRIDE;
            while (offset > 0) {
                offset -= 1;
                if (offsets & (@as(u8, 1) << @intCast(offset)) == 0 or offset > window) continue;
                const at = window - offset;
                if (at + ADAPTER_LENGTH <= sequence.len and adapterAt(sequence, at, ADAPTERS[adapter].sequence)) {
                    starts[adapter] = @intCast(at);
                    found |= bit;
                    break;
                }
            }
        }
        if (found == (1 << ADAPTER_COUNT) - 1) break;
    }
    return starts;
}

fn adapterAt(sequence: []const u8, at: usize, adapter: *const [ADAPTER_LENGTH]u8) bool {
    const head = std.mem.readInt(u64, sequence[at..][0..8], .little) & WINDOW_FOLD;
    const tail = std.mem.readInt(u32, sequence[at + 8 ..][0..4], .little) & @as(u32, @truncate(WINDOW_FOLD));
    return head == std.mem.readInt(u64, adapter[0..8], .little) and tail == std.mem.readInt(u32, adapter[8..12], .little);
}

fn endsInRun(sequence: []const u8, base: u8) bool {
    if (sequence.len < TERMINAL_RUN_MIN) return false;
    for (sequence[sequence.len - TERMINAL_RUN_MIN ..]) |byte| {
        if (std.ascii.toUpper(byte) != base) return false;
    }
    return true;
}

fn lengthBinSize(span: usize) usize {
    var scale: usize = 1;
    while (true) : (scale *= 10) {
        for ([_]usize{ 1, 2, 5 }) |step| {
            if (span / (step * scale) <= LENGTH_BIN_LIMIT) return step * scale;
        }
    }
}

const Group = struct {
    first: usize,
    last: usize,
};

const Groups = struct {
    cycles: usize,
    interval: usize,
    next_first: usize = 0,

    fn init(cycles: usize) Groups {
        return .{ .cycles = cycles, .interval = groupInterval(cycles) };
    }

    fn next(self: *Groups) ?Group {
        if (self.next_first >= self.cycles) return null;
        const first = self.next_first;
        const width = if (self.interval == 1 or first < SINGLE_BASE_GROUPS)
            1
        else if (first == SINGLE_BASE_GROUPS and self.interval > SINGLE_BASE_GROUPS + 1)
            self.interval - SINGLE_BASE_GROUPS - 1
        else
            self.interval;
        self.next_first = @min(first + width, self.cycles);
        return .{ .first = first, .last = self.next_first };
    }
};

fn groupInterval(cycles: usize) usize {
    if (cycles <= GROUPING_MIN_LENGTH) return 1;
    var scale: usize = 1;
    while (true) : (scale *= 10) {
        for ([_]usize{ 2, 5, 10 }) |step| {
            const interval = step * scale;
            const groups = SINGLE_BASE_GROUPS + (cycles - SINGLE_BASE_GROUPS + interval - 1) / interval;
            if (groups < GROUP_COUNT_LIMIT) return interval;
        }
    }
}

fn writeGroupLabel(json: *std.json.Stringify, group: Group) !void {
    try json.objectField("base");
    try writeRangeLabel(json, group.first + 1, group.last);
}

fn writeRangeLabel(json: *std.json.Stringify, first: usize, last: usize) !void {
    var buffer: [48]u8 = undefined;
    const label = if (first == last)
        std.fmt.bufPrint(&buffer, "{d}", .{first}) catch unreachable
    else
        std.fmt.bufPrint(&buffer, "{d}-{d}", .{ first, last }) catch unreachable;
    try json.write(label);
}

fn gcKeptLength(length: usize) usize {
    if (length > 1000) return length / 1000 * 1000;
    if (length > 100) return length / 100 * 100;
    return length;
}

fn gcBins(length: usize, gc_count: usize) [2]usize {
    const n: f64 = @floatFromInt(length);
    const k: f64 = @floatFromInt(gc_count);
    const low = @max(k - 0.5, 0.0);
    const high = @min(k + 0.5, n);
    return .{ javaRound(low * 100.0 / n), javaRound(high * 100.0 / n) };
}

fn javaRound(value: f64) usize {
    return @intFromFloat(value + 0.5);
}

const GC_HUNDREDS_MAX = 1000;
const GC_CACHED_LENGTHS = 100 + GC_HUNDREDS_MAX / 100;
const GcWeights = [GC_BINS]f64;

const GcModel = struct {
    filled: std.StaticBitSet(GC_CACHED_LENGTHS) = .initEmpty(),
    cached: [GC_CACHED_LENGTHS]GcWeights = undefined,
    long_length: usize = 0,
    long: GcWeights = undefined,

    fn weights(self: *GcModel, length: usize) *const GcWeights {
        if (cacheSlot(length)) |slot| {
            if (!self.filled.isSet(slot)) {
                fillWeights(&self.cached[slot], length);
                self.filled.set(slot);
            }
            return &self.cached[slot];
        }
        if (self.long_length != length) {
            fillWeights(&self.long, length);
            self.long_length = length;
        }
        return &self.long;
    }
};

fn cacheSlot(length: usize) ?usize {
    if (length <= 100) return length;
    if (length <= GC_HUNDREDS_MAX) return 100 + length / 100 - 1;
    return null;
}

fn fillWeights(weights: *GcWeights, length: usize) void {
    var coverage: [GC_BINS]u64 = @splat(0);
    for (0..length + 1) |gc_count| {
        const bins = gcBins(length, gc_count);
        for (bins[0]..bins[1] + 1) |bin| coverage[bin] += 1;
    }
    for (weights, coverage) |*weight, count| {
        weight.* = if (count == 0) 0 else 1.0 / @as(f64, @floatFromInt(count));
    }
}

const TestRead = struct { []const u8, []const u8 };

fn collect(allocator: std.mem.Allocator, reads: []const TestRead) !MateQc {
    var mate: MateQc = .{};
    errdefer mate.deinit(allocator);
    for (reads) |read| try mate.addRead(allocator, read[0], read[1]);
    return mate;
}

fn reportOf(arena: std.mem.Allocator, mate: *MateQc) !std.json.ObjectMap {
    mate.finish();
    var output: std.Io.Writer.Allocating = .init(arena);
    var json: std.json.Stringify = .{ .writer = &output.writer };
    try mate.writeJson(&json);
    const value = try std.json.parseFromSliceLeaky(std.json.Value, arena, output.written(), .{});
    return value.object;
}

fn numberOf(value: std.json.Value) ?f64 {
    return switch (value) {
        .null => null,
        .integer => |integer| @floatFromInt(integer),
        .float => |float| float,
        else => unreachable,
    };
}

fn expectColumn(rows: std.json.Value, key: []const u8, expected: []const ?f64) !void {
    try std.testing.expectEqual(expected.len, rows.array.items.len);
    for (rows.array.items, expected) |row, want| {
        try std.testing.expectEqual(want, numberOf(row.object.get(key).?));
    }
}

const PROBE_P1 = [_]TestRead{
    .{ "acgtNNa", "IIII##!" },
    .{ "GGGCCCAT", "JJJJJJJJ" },
    .{ "", "" },
    .{ "NNNNNNNNNNNN", "!!!!!!!!!!!!" },
};

test "[property] - [qc]: FastQC's probe P1 gives FastQC's values in every module" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var mate = try collect(std.testing.allocator, &PROBE_P1);
    defer mate.deinit(std.testing.allocator);
    const report = try reportOf(arena.allocator(), &mate);

    const third = 14.333333333333334;
    try expectColumn(report.get("per_base_quality").?, "mean", &.{ 27, 27, 27, 27, third, third, 13.666666666666666, 20.5, 0, 0, 0, 0 });
    try expectColumn(report.get("per_base_quality").?, "median", &(.{null} ** 12));
    try expectColumn(report.get("per_base_content").?, "g", &(.{ 50, 50, 100, 0, 0, 0, 0, 0 } ++ .{null} ** 4));
    try expectColumn(report.get("per_base_content").?, "c", &(.{ 0, 50, 0, 50, 100, 100, 0, 0 } ++ .{null} ** 4));
    const n_third = 33.33333333333333;
    try expectColumn(report.get("per_base_n").?, "n", &.{ n_third, n_third, n_third, n_third, 66.66666666666666, 66.66666666666666, n_third, 50, 100, 100, 100, 100 });
    try expectColumn(report.get("length_distribution").?, "count", &.{ 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1 });

    const reads_by_quality = report.get("per_sequence_quality").?.array.items;
    try std.testing.expectEqual(@as(usize, 42), reads_by_quality.len);
    for (reads_by_quality, 0..) |row, score| {
        const expected: f64 = if (score == 0 or score == 23 or score == 41) 1 else 0;
        try std.testing.expectEqual(expected, numberOf(row.object.get("count").?).?);
    }

    var gc: [GC_BINS]f64 = @splat(0);
    for (0..4) |bin| gc[bin] = 1;
    for (22..36) |bin| gc[bin] = 1;
    for (70..81) |bin| gc[bin] = 1;
    for ([_]usize{ 4, 21, 36, 69, 81 }) |bin| gc[bin] = 0.5;
    for (report.get("per_sequence_gc").?.array.items, gc) |value, want| {
        try std.testing.expectEqual(want, numberOf(value).?);
    }
}

test "[edge] - [qc]: percentiles need more than 100 reads and take the first score at the rank" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const Case = struct { q10: usize, q30: usize, mean: f64, median: ?f64, upper_quartile: ?f64 };
    const cases = [_]Case{
        .{ .q10 = 50, .q30 = 50, .mean = 20, .median = null, .upper_quartile = null },
        .{ .q10 = 50, .q30 = 51, .mean = 20.099009900990097, .median = 10, .upper_quartile = 30 },
        .{ .q10 = 50, .q30 = 52, .mean = 20.19607843137255, .median = 30, .upper_quartile = 30 },
    };
    for (cases) |case| {
        var mate: MateQc = .{};
        defer mate.deinit(std.testing.allocator);
        for (0..case.q10) |_| try mate.addRead(std.testing.allocator, "A", "+");
        for (0..case.q30) |_| try mate.addRead(std.testing.allocator, "A", "?");
        const quality = (try reportOf(arena.allocator(), &mate)).get("per_base_quality").?;
        try expectColumn(quality, "mean", &.{case.mean});
        try expectColumn(quality, "median", &.{case.median});
        try expectColumn(quality, "upper_quartile", &.{case.upper_quartile});
    }
}

test "[unit] - [qc]: cycles group as FastQC groups them at its interval boundaries" {
    const Case = struct { cycles: usize, groups: usize, tenth: Group, last: Group };
    const cases = [_]Case{
        .{ .cycles = 75, .groups = 75, .tenth = .{ .first = 9, .last = 10 }, .last = .{ .first = 74, .last = 75 } },
        .{ .cycles = 76, .groups = 43, .tenth = .{ .first = 9, .last = 11 }, .last = .{ .first = 75, .last = 76 } },
        .{ .cycles = 151, .groups = 38, .tenth = .{ .first = 9, .last = 14 }, .last = .{ .first = 149, .last = 151 } },
        .{ .cycles = 335, .groups = 42, .tenth = .{ .first = 9, .last = 19 }, .last = .{ .first = 329, .last = 335 } },
        .{ .cycles = 1310, .groups = 36, .tenth = .{ .first = 9, .last = 49 }, .last = .{ .first = 1299, .last = 1310 } },
    };
    for (cases) |case| {
        var groups = Groups.init(case.cycles);
        var count: usize = 0;
        var last: Group = undefined;
        while (groups.next()) |group| : (count += 1) {
            if (count == 9) try std.testing.expectEqual(case.tenth, group);
            last = group;
        }
        try std.testing.expectEqual(case.groups, count);
        try std.testing.expectEqual(case.last, last);
    }
}

test "[unit] - [qc]: GC keeps whole hundreds and thousands of bases and splits a count over shared bins" {
    for ([_][2]usize{ .{ 100, 100 }, .{ 101, 100 }, .{ 150, 100 }, .{ 250, 200 }, .{ 1000, 1000 }, .{ 1001, 1000 }, .{ 1500, 1000 } }) |case| {
        try std.testing.expectEqual(case[1], gcKeptLength(case[0]));
    }
    var mate: MateQc = .{};
    defer mate.deinit(std.testing.allocator);
    try mate.addRead(std.testing.allocator, "GCGAAAAAAA", "IIIIIIIIII");
    for (mate.gc, 0..) |value, bin| {
        const expected: f64 = if (bin == 25 or bin == 35) 0.5 else if (bin > 25 and bin < 35) 1 else 0;
        try std.testing.expectEqual(expected, value);
    }
}

test "[edge] - [qc]: a terminal run needs the last ten bases to match" {
    const Case = struct { sequence: []const u8, g: u64, a: u64 };
    const cases = [_]Case{
        .{ .sequence = "ACGGGGGGGGG", .g = 0, .a = 0 },
        .{ .sequence = "AGGGGGGGGGG", .g = 1, .a = 0 },
        .{ .sequence = "gggggggggggg", .g = 1, .a = 0 },
        .{ .sequence = "GGGGGNGGGGG", .g = 0, .a = 0 },
        .{ .sequence = "CAAAAAAAAAA", .g = 0, .a = 1 },
    };
    for (cases) |case| {
        var mate: MateQc = .{};
        defer mate.deinit(std.testing.allocator);
        var quality: [12]u8 = @splat('I');
        try mate.addRead(std.testing.allocator, case.sequence, quality[0..case.sequence.len]);
        try std.testing.expectEqual(case.g, mate.terminal_g);
        try std.testing.expectEqual(case.a, mate.terminal_a);
    }
}

test "[edge] - [qc]: a read past the cycle limit makes only the per-cycle modules unavailable, one at it keeps them" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const long_sequence = try arena.allocator().alloc(u8, MAX_CYCLES + 1);
    @memset(long_sequence, 'A');
    const long_quality = try arena.allocator().alloc(u8, MAX_CYCLES + 1);
    @memset(long_quality, 'I');
    var mate = try collect(std.testing.allocator, &.{ .{ "ACGT", "IIII" }, .{ long_sequence, long_quality } });
    defer mate.deinit(std.testing.allocator);

    const report = try reportOf(arena.allocator(), &mate);
    for ([_][]const u8{ "per_base_quality", "per_base_content", "per_base_n", "length_distribution", "adapter_content" }) |key| {
        try std.testing.expectEqual(std.json.Value.null, report.get(key).?);
    }
    try std.testing.expectEqualStrings(TOO_LONG, report.get("unavailable").?.string);
    try std.testing.expectEqual(@as(usize, 1), report.get("per_sequence_quality").?.array.items.len);
    try std.testing.expectEqual(@as(usize, 0), mate.quality.capacity);

    const limit_sequence = try arena.allocator().dupe(u8, long_sequence[0..MAX_CYCLES]);
    @memcpy(limit_sequence[MAX_CYCLES - ADAPTER_LENGTH ..], ADAPTERS[0].sequence);
    var limit = try collect(std.testing.allocator, &.{.{ limit_sequence, long_quality[0..MAX_CYCLES] }});
    defer limit.deinit(std.testing.allocator);
    const limit_report = try reportOf(arena.allocator(), &limit);
    try std.testing.expectEqual(std.json.Value.null, limit_report.get("unavailable").?);
    const rows = limit_report.get("adapter_content").?.object.get("rows").?.array.items;
    try std.testing.expectEqualStrings("1000-1013", rows[rows.len - 1].object.get("position").?.string);
    try std.testing.expectEqual(@as(?f64, 100.0 / 14.0), numberOf(rows[rows.len - 1].object.get("values").?.array.items[0]));
}

test "[failure] - [qc]: a failed cycle allocation leaves the counts unchanged" {
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    var mate: MateQc = .{};
    defer mate.deinit(failing.allocator());
    try std.testing.expectError(error.OutOfMemory, mate.addRead(failing.allocator(), "ACGT", "IIII"));
    try std.testing.expectEqual(@as(u64, 0), mate.lengths[4]);
    try std.testing.expectEqual(@as(usize, 0), mate.quality.items.len);
}

test "[failure] - [qc]: a failed retained allocation leaves the retained counts unchanged" {
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
    var populations: MatePopulations = .{};
    defer populations.deinit(failing.allocator());
    try populations.addRead(failing.allocator(), "ACGT", "IIII");
    try std.testing.expectError(error.OutOfMemory, populations.addRetained(failing.allocator(), "ACGT", "IIII", 2));
    try std.testing.expectEqual(@as(u64, 0), populations.retained.lengths[2]);
    try std.testing.expectEqual(@as(usize, 0), populations.retained_cycles);
    try std.testing.expectEqual(@as(u64, 1), populations.raw.lengths[4]);
}

test "[property] - [qc]: vector cycle counts match a scalar count at every tier and lane boundary" {
    const pattern = "ACGTNacgtnRYgc";
    var sequence: [151]u8 = undefined;
    for (&sequence, 0..) |*base, index| base.* = pattern[index % pattern.len];
    inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |lanes| {
        for ([_]usize{ 0, 1, lanes - 1, lanes, lanes + 1, 2 * lanes + 1, 100, 151 }) |length| {
            for ([_]usize{ 0, length / 2, length }) |gc_length| {
                var vector: PendingBases = @splat(@splat(0));
                var scalar: PendingBases = @splat(@splat(0));
                var scalar_gc: usize = 0;
                for (sequence[0..length], 0..) |base, cycle| {
                    const upper = std.ascii.toUpper(base);
                    for (COUNTED_LETTERS, &scalar) |letter, *row| row[cycle] += @intFromBool(upper == letter);
                    scalar_gc += @intFromBool(cycle < gc_length and (upper == 'G' or upper == 'C'));
                }
                const vector_gc = countCycleBases(lanes, &vector, sequence[0..length], gc_length);
                try std.testing.expectEqual(scalar_gc, vector_gc);
                try std.testing.expectEqualDeep(scalar, vector);
            }
        }
    }
}

test "[property] - [qc]: statistics from the tables equal Stats.addRecord, also past the cycle limit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const long_sequence = try arena.allocator().alloc(u8, MAX_CYCLES + 1);
    for (long_sequence, 0..) |*base, index| base.* = "ACGTNacgtR"[index % 10];
    const long_quality = try arena.allocator().alloc(u8, MAX_CYCLES + 1);
    for (long_quality, 0..) |*byte, index| byte.* = @intCast(PHRED_OFFSET + index % (QUALITY_BINS - 1));
    const short_reads = PROBE_P1;
    for ([_][]const TestRead{ &short_reads, &(short_reads ++ .{.{ long_sequence, long_quality }} ++ short_reads) }) |reads| {
        var expected: statistics.Stats = .{};
        for (reads) |read| try expected.addRecord(.{ .header = "", .id = "", .sequence = read[0], .plus = "", .quality = read[1] });
        var mate = try collect(std.testing.allocator, reads);
        defer mate.deinit(std.testing.allocator);
        mate.finish();
        try std.testing.expectEqualDeep(expected.result(), (try mate.stats()).result());
    }
}

test "[property] - [qc]: a retained population derived from raw equals one collected from the kept reads" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = std.testing.allocator;
    var long_sequence: [MAX_CYCLES + 5]u8 = undefined;
    for (&long_sequence, 0..) |*base, index| base.* = "ACGTGGCAN"[index % 9];
    const long_quality: [MAX_CYCLES + 5]u8 = @splat('5');
    const Read = struct { sequence: []const u8, quality: []const u8, kept: usize, dropped: bool = false };
    const reads = [_]Read{
        .{ .sequence = "ACGTNNacgtGGGG", .quality = "II#5!IIII?5+#2", .kept = 9 },
        .{ .sequence = "CCAGATCGGAAGAGTTCTGTCTCTTATACC", .quality = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII", .kept = 25 },
        .{ .sequence = "AGATCGGAAGAGTTAAAAAAAAAAAAAAAG", .quality = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII", .kept = 30 },
        .{ .sequence = "GGGCCCAT", .quality = "JJJJJJJJ", .kept = 8 },
        .{ .sequence = "AACCGGTTNNAC", .quality = "IIII#####III", .kept = 12, .dropped = true },
        .{ .sequence = "TTTTAAAAGGGGGGGGGGGG", .quality = "IIIIIIIIIIIIIIIIIII#", .kept = 0 },
        .{ .sequence = "", .quality = "", .kept = 0 },
        .{ .sequence = "NNNNACGTACGT", .quality = "!!!!IIIIIIII", .kept = 3 },
        .{ .sequence = &long_sequence, .quality = &long_quality, .kept = 700 },
        .{ .sequence = "CCCCAAAAAAAAAAAA", .quality = "IIIIIIIIIIIIIIII", .kept = 16 },
        .{ .sequence = "GGGGGGGGGGGGGGGGGG", .quality = "!!!!!!!!!!!!!!!!!!", .kept = 4, .dropped = true },
        .{ .sequence = "ACGT", .quality = "IIII", .kept = 1 },
    };
    for ([_]usize{ 6, reads.len }) |count| {
        errdefer std.debug.print("{d} reads\n", .{count});
        var populations: MatePopulations = .{};
        defer populations.deinit(allocator);
        var expected: MateQc = .{};
        defer expected.deinit(allocator);
        for (reads[0..count]) |read| {
            try populations.addRead(allocator, read.sequence, read.quality);
            if (read.dropped) {
                try populations.dropRetained(allocator, read.sequence, read.quality);
                continue;
            }
            try populations.addRetained(allocator, read.sequence, read.quality, read.kept);
            try expected.addRead(allocator, read.sequence[0..read.kept], read.quality[0..read.kept]);
        }
        populations.finish();
        const want = try reportOf(arena.allocator(), &expected);
        const got = try reportOf(arena.allocator(), &populations.retained);
        inline for (.{ "per_base_quality", "per_base_content", "per_base_n", "length_distribution", "per_sequence_quality", "per_sequence_gc", "terminal_runs", "adapter_content", "unavailable" }) |module| {
            errdefer std.debug.print("{s}\n", .{module});
            const want_text = try std.json.Stringify.valueAlloc(arena.allocator(), want.get(module).?, .{});
            const got_text = try std.json.Stringify.valueAlloc(arena.allocator(), got.get(module).?, .{});
            try std.testing.expectEqualStrings(want_text, got_text);
        }
        try std.testing.expectEqual(try expected.stats(), try populations.retained.stats());
    }
}

fn adapterRows(arena: std.mem.Allocator, reads: []const []const u8) ![]const std.json.Value {
    var mate: MateQc = .{};
    defer mate.deinit(std.testing.allocator);
    for (reads) |read| {
        const quality = try arena.alloc(u8, read.len);
        @memset(quality, 'I');
        try mate.addRead(std.testing.allocator, read, quality);
    }
    const report = try reportOf(arena, &mate);
    return report.get("adapter_content").?.object.get("rows").?.array.items;
}

fn expectAdapterRow(row: std.json.Value, label: []const u8, values: [ADAPTER_COUNT]f64) !void {
    try std.testing.expectEqualStrings(label, row.object.get("position").?.string);
    for (row.object.get("values").?.array.items, values) |got, want| try std.testing.expectEqual(@as(?f64, want), numberOf(got));
}

test "[property] - [qc]: adapter content counts first exact matches cumulatively over all reads, as FastQC-Rust does" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const a = arena.allocator();
    const u = "AGATCGGAAGAG";
    const third = 33.333333333333336;
    const two_thirds = 66.66666666666667;

    const two_starts = try adapterRows(a, &.{ u ++ "CCCCCCCC", "TTTTT" ++ u ++ "GGGG", "ACACACACACACAC" });
    try std.testing.expectEqual(@as(usize, 10), two_starts.len);
    try expectAdapterRow(two_starts[0], "1", .{ third, 0, 0, 0, 0, 0 });
    try expectAdapterRow(two_starts[4], "5", .{ third, 0, 0, 0, 0, 0 });
    try expectAdapterRow(two_starts[5], "6", .{ two_thirds, 0, 0, 0, 0, 0 });
    try expectAdapterRow(two_starts[9], "10", .{ two_thirds, 0, 0, 0, 0, 0 });

    const short_first = try adapterRows(a, &.{ "ACGTA", "ACGTACGTACG", u, "C" ++ u });
    try std.testing.expectEqual(@as(usize, 2), short_first.len);
    try expectAdapterRow(short_first[0], "1", .{ 25, 0, 0, 0, 0, 0 });
    try expectAdapterRow(short_first[1], "2", .{ 50, 0, 0, 0, 0, 0 });

    const empty_reads = try adapterRows(a, &.{ "", u ++ "CCCCCCCCC", "", "ACGT" ** 5 });
    try expectAdapterRow(empty_reads[0], "1", .{ 25, 0, 0, 0, 0, 0 });

    const several_adapters = try adapterRows(a, &.{ u ++ "CCCCC" ++ "TGGAATTCTCGG" ++ "A", "GGGG" ++ "CTGTCTCTTATA" });
    try expectAdapterRow(several_adapters[several_adapters.len - 1], "19", .{ 50, 50, 0, 50, 0, 0 });

    const repeated = try adapterRows(a, &.{ u ++ "AAAA" ++ u, "C" ** 29 });
    try expectAdapterRow(repeated[repeated.len - 1], "18", .{ 50, 0, 0, 0, 0, 0 });

    const unknown_bases = try adapterRows(a, &.{ "AGATNGGAAGAG" ++ "CCCCCCCC", "AGATCGGAAGNG" ++ "CCCC", u ++ "CCCC" });
    try expectAdapterRow(unknown_bases[0], "1", .{ third, 0, 0, 0, 0, 0 });

    const grouped = try adapterRows(a, &.{u ++ "C" ** 75});
    try std.testing.expectEqual(@as(usize, 43), grouped.len);
    try expectAdapterRow(grouped[9], "10-11", .{ 100, 0, 0, 0, 0, 0 });
    try std.testing.expectEqualStrings("76", grouped[42].object.get("position").?.string);

    const short = try adapterRows(a, &.{ u, "ACGTACG" });
    try std.testing.expectEqual(@as(usize, 0), short.len);
    const none = try adapterRows(a, &.{});
    try std.testing.expectEqual(@as(usize, 0), none.len);
}

fn firstAdapterMatchesScalar(sequence: []const u8) AdapterStarts {
    var starts: AdapterStarts = @splat(NO_MATCH);
    for (ADAPTERS, &starts) |adapter, *start| {
        for (0..sequence.len -| (ADAPTER_LENGTH - 1)) |at| {
            if (std.ascii.eqlIgnoreCase(sequence[at..][0..ADAPTER_LENGTH], adapter.sequence)) {
                start.* = @intCast(at);
                break;
            }
        }
    }
    return starts;
}

test "[property] - [qc]: the strided adapter search finds the first match a direct search finds" {
    const reads = [_][]const u8{
        "",                         "AGATCGGAAGA",                          "AGATCGGAAGAG",      "agatcggaagagNAGATCGGAAGAG", "GGGGGGGGGGGGGG", "AAAAAAAAAAAAAAAAAGATCGGAAGAG",
        "CTGTCTCTTATACTGTCTCTTATA", "TGGAATTCTCGGATCGTCGGACTRAAAAAAAAAAAA", "ACGTNACGTACGTACGT", "GATCGTCGGACTGATCGTCGGACT",
    };
    for (reads, 0..) |text, index| {
        errdefer std.debug.print("read {d}\n", .{index});
        try std.testing.expectEqual(firstAdapterMatchesScalar(text), firstAdapterMatches(text));
    }
    var read: [80]u8 = undefined;
    var state: u32 = 3;
    for (&read) |*base| {
        state = state *% 1103515245 +% 12345;
        base.* = "ACGTacgt"[(state >> 16) & 7];
    }
    for (ADAPTERS, 0..) |adapter, adapter_index| {
        for ([_]usize{ 0, 11, 12, 13, 15, 16, 19, 20, 23, 31, 40, 79, 80 }) |length| {
            for (0..length -| (ADAPTER_LENGTH - 1)) |at| {
                errdefer std.debug.print("adapter {d} length {d} at {d}\n", .{ adapter_index, length, at });
                var copy = read;
                @memcpy(copy[at..][0..ADAPTER_LENGTH], adapter.sequence);
                if (at % 3 == 1) copy[at + 4] = std.ascii.toLower(copy[at + 4]);
                if (at % 7 == 2 and at + 20 < length) @memcpy(copy[at + 8 ..][0..ADAPTER_LENGTH], adapter.sequence);
                if (at % 5 == 3) copy[(at + 6) % length] = 'N';
                try std.testing.expectEqual(firstAdapterMatchesScalar(copy[0..length]), firstAdapterMatches(copy[0..length]));
            }
        }
    }
}
