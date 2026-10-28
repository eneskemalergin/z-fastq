//! Per-mate FASTQ quality distributions for `process`, computed as FastQC 0.13.0 computes its core
//! modules: per-cycle quality and base content, per-read quality and GC, read lengths, terminal
//! homopolymer runs, and Adapter Content for FastQC's bundled adapter list, where a read exactly
//! one adapter long counts as FastQC-Rust counts it. The retained population of a trimmed or
//! filtered run counts only the removed 3' ends and dropped reads per cycle and takes its
//! per-cycle tables as raw minus removed, until raw drops its tables for a read past the cycle
//! limit; from then on it counts its reads in full. Adapter matches are searched once per raw
//! read and count for the retained read when they lie inside its kept bases. Sequence Duplication
//! Levels and Overrepresented sequences share FastQC's sequence table, one per population: the
//! uppercase read cut to 50 bases is the key, at most 100,000 keys are stored, and later reads
//! count only toward stored keys and the total. A key is kept as 50 five-bit codes, the low five
//! bits of each IUPAC letter, which are distinct and nonzero. Overrepresented ties keep
//! first-appearance order, and Possible Source comes from a FastQC contaminant list matched by
//! FastQC's rules. A population also writes as a FastQC 0.13.0 fastqc_data.txt: the same tables
//! in FastQC's text layout and Java number format, each module graded pass, warn, or fail by
//! FastQC's rules at its default limits, without Kmer Content. Basic Statistics names the
//! Sanger / Illumina 1.9 encoding, the only one z-fastq reads. As a read section of fastp 1.1.0's
//! JSON it writes fastp's totals and per-cycle mean quality and base content curves.

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
    raw_sequences: SequenceCounts = .{},
    retained_sequences: SequenceCounts = .{},
    retained_enabled: bool = false,
    last_key: SequenceKey = .{},
    raw_tiles: TileQuality = .{},
    retained_tiles: TileQuality = .{},

    pub fn deinit(self: *MatePopulations, allocator: std.mem.Allocator) void {
        self.raw.deinit(allocator);
        self.retained.deinit(allocator);
        self.raw_sequences.deinit(allocator);
        self.retained_sequences.deinit(allocator);
        self.raw_tiles.deinit(allocator);
        self.retained_tiles.deinit(allocator);
    }

    pub fn enableRetained(self: *MatePopulations) void {
        std.debug.assert(self.raw_sequences.reads == 0);
        self.retained_enabled = true;
    }

    pub fn addRead(
        self: *MatePopulations,
        allocator: std.mem.Allocator,
        header: []const u8,
        sequence: []const u8,
        quality: []const u8,
    ) (error{OutOfMemory} || statistics.StatsError)!void {
        if (self.derived and !self.raw.too_long and sequence.len > MAX_CYCLES) self.materialize();
        const key: SequenceKey = .of(sequence);
        try self.raw_sequences.reserve(allocator);
        self.raw_sequences.prefetch(key);
        self.last = try self.raw.collect(allocator, sequence, quality);
        self.raw_sequences.add(key);
        try self.raw_tiles.add(allocator, header, quality);
        if (self.retained_enabled) {
            self.last_key = key;
            self.retained_sequences.prefetch(key);
        }
    }

    pub fn addRetained(
        self: *MatePopulations,
        allocator: std.mem.Allocator,
        header: []const u8,
        sequence: []const u8,
        quality: []const u8,
        kept: usize,
    ) (error{OutOfMemory} || statistics.StatsError)!void {
        std.debug.assert(kept <= sequence.len and sequence.len == quality.len);
        std.debug.assert(self.retained_enabled);
        const key = if (kept >= @min(sequence.len, KEY_LENGTH)) self.last_key else SequenceKey.of(sequence[0..kept]);
        try self.retained_sequences.reserve(allocator);
        try self.retained_tiles.add(allocator, header, quality[0..kept]);
        if (!self.derived) {
            try self.retained.addRead(allocator, sequence[0..kept], quality[0..kept]);
            self.retained_sequences.add(key);
            return;
        }
        const retained = &self.retained;
        if (sequence.len > retained.quality.items.len) try retained.growCycles(allocator, sequence.len);
        self.retained_sequences.add(key);
        var facts = self.last;
        if (kept < sequence.len) {
            const removed_sum = try self.removeCycles(allocator, sequence, quality, kept);
            facts.quality_sum -= removed_sum;
            facts.gc_count = gcCountScalar(sequence[0..gcKeptLength(kept)]);
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

    pub fn finish(self: *MatePopulations, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        if (self.derived) self.materialize();
        self.raw.finish();
        self.retained.finish();
        try self.raw_sequences.finish(allocator);
        if (self.retained_enabled) try self.retained_sequences.finish(allocator);
    }

    pub fn writeJson(self: *const MatePopulations, json: *std.json.Stringify, retained: bool, contaminants: ?*const ContaminantList) !void {
        try self.population(retained).writeJson(json, self.view(retained, contaminants));
    }

    pub fn tableUnavailable(self: *const MatePopulations, retained: bool) ?[]const u8 {
        return if (self.population(retained).too_long) TOO_LONG else null;
    }

    pub fn writeQcTable(
        self: *const MatePopulations,
        out: *std.Io.Writer,
        retained: bool,
        filename: []const u8,
        version: []const u8,
        contaminants: ?*const ContaminantList,
    ) !void {
        try self.population(retained).writeQcTable(out, filename, version, self.view(retained, contaminants));
    }

    pub fn writeSummaryRead(self: *const MatePopulations, json: *std.json.Stringify, retained: bool) !void {
        try self.population(retained).writeSummaryRead(json);
    }

    fn population(self: *const MatePopulations, retained: bool) *const MateQc {
        return if (retained) &self.retained else &self.raw;
    }

    fn view(self: *const MatePopulations, retained: bool, contaminants: ?*const ContaminantList) PopulationView {
        return .{
            .counts = if (retained) &self.retained_sequences else &self.raw_sequences,
            .tiles = if (retained) &self.retained_tiles else &self.raw_tiles,
            .contaminants = contaminants,
        };
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

    fn writeJson(self: *const MateQc, json: *std.json.Stringify, sequences: PopulationView) !void {
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
        try json.objectField("duplication");
        try sequences.writeDuplication(json);
        try json.objectField("overrepresented");
        try sequences.writeOverrepresented(json);
        try json.objectField("per_tile_quality");
        try sequences.tiles.writeJson(json);
        try json.objectField("per_tile_unavailable");
        try json.write(sequences.tiles.unavailable());
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
        var rows = self.adapterRows();
        while (rows.next()) |row| {
            try json.beginObject();
            try json.objectField("position");
            try writeRangeLabel(json, row.group.first + 1, row.group.last);
            try json.objectField("values");
            try json.beginArray();
            for (row.values) |value| try json.write(value);
            try json.endArray();
            try json.endObject();
        }
        try json.endArray();
        try json.endObject();
    }

    const AdapterRow = struct {
        group: Group,
        values: [ADAPTER_COUNT]f64,
    };

    const AdapterRows = struct {
        mate: *const MateQc,
        groups: ?Groups,
        total: f64,
        running: [ADAPTER_COUNT]u64 = @splat(0),
        position: usize = 0,

        fn next(self: *AdapterRows) ?AdapterRow {
            const groups = if (self.groups) |*groups| groups else return null;
            const group = groups.next() orelse return null;
            var sums: [ADAPTER_COUNT]f64 = @splat(0);
            while (self.position < group.last) : (self.position += 1) {
                for (&self.running, &sums, &self.mate.adapter_starts) |*count, *sum, *starts| {
                    count.* += starts[self.position];
                    sum.* += @as(f64, @floatFromInt(count.*)) * 100.0 / self.total;
                }
            }
            const width: f64 = @floatFromInt(group.last - group.first);
            var row: AdapterRow = .{ .group = group, .values = undefined };
            for (&row.values, sums) |*value, sum| value.* = sum / width;
            return row;
        }
    };

    fn adapterRows(self: *const MateQc) AdapterRows {
        var reads: u64 = 0;
        for (self.lengths) |count| reads += count;
        const longest = self.longestRead();
        const groups: ?Groups = if (reads != 0 and longest > ADAPTER_LENGTH) .init(longest - ADAPTER_LENGTH + 1) else null;
        return .{ .mate = self, .groups = groups, .total = @floatFromInt(reads) };
    }

    fn longestRead(self: *const MateQc) usize {
        return std.mem.lastIndexOfNone(u64, &self.lengths, &.{0}) orelse 0;
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
            try json.beginObject();
            try writeGroupLabel(json, group);
            for (CONTENT_NAMES, self.groupContent(group)) |name, percent| {
                try json.objectField(name);
                try json.write(percent);
            }
            try json.endObject();
        }
        try json.endArray();
    }

    const CONTENT_ORDER = [_]BaseKind{ .g, .a, .t, .c };
    const CONTENT_NAMES = [_][]const u8{ "g", "a", "t", "c" };

    fn groupContent(self: *const MateQc, group: Group) [CONTENT_ORDER.len]?f64 {
        const pooled = self.pooledBases(group);
        const total = pooled[@intFromEnum(BaseKind.a)] + pooled[@intFromEnum(BaseKind.c)] +
            pooled[@intFromEnum(BaseKind.g)] + pooled[@intFromEnum(BaseKind.t)];
        var percents: [CONTENT_ORDER.len]?f64 = undefined;
        for (&percents, CONTENT_ORDER) |*percent, kind| percent.* = percentOf(pooled[@intFromEnum(kind)], total);
        return percents;
    }

    fn writeBaseN(self: *const MateQc, json: *std.json.Stringify, groups: Groups) !void {
        try json.beginArray();
        var iterator = groups;
        while (iterator.next()) |group| {
            try json.beginObject();
            try writeGroupLabel(json, group);
            try json.objectField("n");
            try json.write(self.groupN(group));
            try json.endObject();
        }
        try json.endArray();
    }

    fn groupN(self: *const MateQc, group: Group) ?f64 {
        return percentOf(self.pooledBases(group)[@intFromEnum(BaseKind.n)], self.pooledReads(group));
    }

    fn writeLengths(self: *const MateQc, json: *std.json.Stringify) !void {
        try json.beginArray();
        var bins = self.lengthBins();
        while (bins.next()) |bin| {
            try json.beginObject();
            try json.objectField("length");
            try writeRangeLabel(json, bin.first, bin.last);
            try json.objectField("count");
            try json.write(bin.count);
            try json.endObject();
        }
        try json.endArray();
    }

    const LengthBin = struct {
        first: usize,
        last: usize,
        count: u64,
    };

    const LengthBins = struct {
        lengths: *const [MAX_CYCLES + 1]u64,
        low_start: usize,
        high: usize,
        size: usize,
        start: usize,

        fn next(self: *LengthBins) ?LengthBin {
            while (self.start <= self.high) {
                const start = self.start;
                self.start += self.size;
                const end = @min(start + self.size - 1, self.high);
                var count: u64 = 0;
                for (start..@min(end, MAX_CYCLES) + 1) |length| count += self.lengths[length];
                const edge = start == self.low_start or start + self.size > self.high;
                if (edge and count == 0) continue;
                return .{ .first = start, .last = if (self.size == 1) start else end, .count = count };
            }
            return null;
        }
    };

    fn lengthBins(self: *const MateQc) LengthBins {
        const first = std.mem.indexOfNone(u64, &self.lengths, &.{0}) orelse 0;
        const low = first -| 1;
        const high = self.longestRead() + 1;
        const size = lengthBinSize(high - low);
        return .{ .lengths = &self.lengths, .low_start = low / size * size, .high = high, .size = size, .start = low / size * size };
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

    fn writeSummaryRead(self: *const MateQc, json: *std.json.Stringify) !void {
        std.debug.assert(self.pending_reads == 0);
        const totals = try self.stats();
        try json.beginObject();
        for ([_][]const u8{ "total_reads", "total_bases", "q20_bases", "q30_bases" }, [_]u64{ totals.reads, totals.bases, totals.q20_bases, totals.q30_bases }) |name, value| {
            try json.objectField(name);
            try json.write(value);
        }
        if (!self.too_long) {
            const cycles = self.quality.items[0..self.longestRead()];
            var q40_bases: u64 = 0;
            for (cycles) |counts| {
                for (counts[Q40_SCORE..]) |count| q40_bases += count;
            }
            try json.objectField("q40_bases");
            try json.write(q40_bases);
            try json.objectField("total_cycles");
            try json.write(cycles.len);
            try json.objectField("quality_curves");
            try json.beginObject();
            try json.objectField("mean");
            try json.beginArray();
            for (cycles) |counts| {
                var reads: u64 = 0;
                var scores: u64 = 0;
                for (counts, 0..) |count, score| {
                    reads += count;
                    scores += count * score;
                }
                try json.write(@as(f64, @floatFromInt(scores)) / @as(f64, @floatFromInt(reads)));
            }
            try json.endArray();
            try json.endObject();
            try json.objectField("content_curves");
            try json.beginObject();
            for ([_][]const u8{ "A", "T", "C", "G", "N", "GC" }, [_][]const BaseKind{ &.{.a}, &.{.t}, &.{.c}, &.{.g}, &.{.n}, &.{ .g, .c } }) |name, kinds| {
                try json.objectField(name);
                try json.beginArray();
                for (cycles, self.bases.items[0..cycles.len]) |qualities, bases| {
                    var reads: u64 = 0;
                    for (qualities) |count| reads += count;
                    var matched: u64 = 0;
                    for (kinds) |kind| matched += bases[@intFromEnum(kind)];
                    try json.write(@as(f64, @floatFromInt(matched)) / @as(f64, @floatFromInt(reads)));
                }
                try json.endArray();
            }
            try json.endObject();
        }
        try json.endObject();
    }

    fn writeQcTable(self: *const MateQc, out: *std.Io.Writer, filename: []const u8, version: []const u8, sequences: PopulationView) !void {
        std.debug.assert(self.pending_reads == 0 and !self.too_long);
        const groups = Groups.init(self.quality.items.len);
        try out.print("##FastQC\t{s}\n", .{version});
        try self.writeBasicStatistics(out, filename);
        if (self.quality.items.len != 0) try self.writeTableBaseQuality(out, groups);
        if (sequences.tiles.unavailable() == null) try sequences.tiles.writeTable(out);
        if (std.mem.indexOfNone(u64, &self.read_quality, &.{0}) != null) try self.writeTableReadQuality(out);
        try self.writeTableBaseContent(out, groups);
        try self.writeTableGc(out);
        try self.writeTableBaseN(out, groups);
        try self.writeTableLengths(out);
        try sequences.counts.writeTableDuplication(out);
        try sequences.counts.writeTableOverrepresented(out, sequences.contaminants);
        try self.writeTableAdapters(out);
    }

    fn writeBasicStatistics(self: *const MateQc, out: *std.Io.Writer, filename: []const u8) !void {
        const totals = try self.tableStats();
        try beginModule(out, "Basic Statistics", .pass, "Measure\tValue");
        const file_type = if (totals.reads == 0) "null" else "Conventional base calls";
        try out.print("Filename\t{s}\nFile type\t{s}\nEncoding\tSanger / Illumina 1.9\n", .{ filename, file_type });
        var buffer: [LENGTH_TEXT_BYTES]u8 = undefined;
        try out.print("Total Sequences\t{d}\nTotal Bases\t{s}\n", .{ totals.reads, formatLength(&buffer, totals.bases) });
        if (totals.min_length == totals.max_length)
            try out.print("Sequence length\t{d}\n", .{totals.min_length})
        else
            try out.print("Sequence length\t{d}-{d}\n", .{ totals.min_length, totals.max_length });
        const counted: u128 = totals.a + totals.c + totals.g + totals.t;
        try out.print("Mean Length\t{d}\nMedian Length\t{d}\n%GC\t{d}\n", .{
            if (totals.reads == 0) 0 else totals.bases / totals.reads,
            self.medianLength(),
            if (counted == 0) 0 else (@as(u128, totals.g) + totals.c) * 100 / counted,
        });
        try endModule(out);
    }

    fn medianLength(self: *const MateQc) usize {
        var reads: u64 = 0;
        for (self.lengths) |count| reads += count;
        var seen: u64 = 0;
        for (self.lengths, 0..) |count, length| {
            if (seen + count > reads / 2) return length;
            seen += count;
        }
        return 0;
    }

    fn writeTableBaseQuality(self: *const MateQc, out: *std.Io.Writer, groups: Groups) !void {
        var grade: Grade = .pass;
        var scan = groups;
        while (scan.next()) |group| {
            const lower = self.groupPercentile(group, 25) orelse continue;
            const median = self.groupPercentile(group, 50).?;
            grade = grade.worst(.below(lower, QUALITY_LOWER_WARN, QUALITY_LOWER_FAIL));
            grade = grade.worst(.below(median, QUALITY_MEDIAN_WARN, QUALITY_MEDIAN_FAIL));
        }
        try beginModule(out, "Per base sequence quality", grade, "Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile");
        var rows = groups;
        while (rows.next()) |group| {
            try writeLabel(out, group.first + 1, group.last);
            try writeCell(out, self.groupMean(group));
            for ([_]u7{ 50, 25, 75, 10, 90 }) |percent| try writeCell(out, self.groupPercentile(group, percent));
            try out.writeByte('\n');
        }
        try endModule(out);
    }

    fn writeTableReadQuality(self: *const MateQc, out: *std.Io.Writer) !void {
        var mode: usize = 0;
        var most: u64 = 0;
        for (self.read_quality, 0..) |count, score| {
            if (count <= most) continue;
            most = count;
            mode = score;
        }
        try beginModule(out, "Per sequence quality scores", .atMost(@floatFromInt(mode), READ_QUALITY_WARN, READ_QUALITY_FAIL), "Quality\tCount");
        if (std.mem.indexOfNone(u64, &self.read_quality, &.{0})) |lowest| {
            const highest = std.mem.lastIndexOfNone(u64, &self.read_quality, &.{0}).?;
            for (lowest..highest + 1) |score| try writeCountRow(out, score, self.read_quality[score]);
        }
        try endModule(out);
    }

    fn writeTableBaseContent(self: *const MateQc, out: *std.Io.Writer, groups: Groups) !void {
        const g, const a, const t, const c = .{ 0, 1, 2, 3 };
        var grade: Grade = .pass;
        var scan = groups;
        while (scan.next()) |group| {
            const percents = self.groupContent(group);
            const at = if (percents[a] != null and percents[t] != null) @abs(percents[a].? - percents[t].?) else 0;
            const gc = if (percents[g] != null and percents[c] != null) @abs(percents[g].? - percents[c].?) else 0;
            grade = grade.worst(.above(@max(at, gc), CONTENT_WARN, CONTENT_FAIL));
        }
        try beginModule(out, "Per base sequence content", grade, "Base\tG\tA\tT\tC");
        var rows = groups;
        while (rows.next()) |group| {
            try writeLabel(out, group.first + 1, group.last);
            for (self.groupContent(group)) |percent| try writeCell(out, percent);
            try out.writeByte('\n');
        }
        try endModule(out);
    }

    fn writeTableGc(self: *const MateQc, out: *std.Io.Writer) !void {
        try beginModule(out, "Per sequence GC content", .above(gcDeviation(&self.gc), GC_WARN, GC_FAIL), "GC Content\tCount");
        for (self.gc, 0..) |count, bin| {
            try out.print("{d}", .{bin});
            try writeCell(out, count);
            try out.writeByte('\n');
        }
        try endModule(out);
    }

    fn writeTableBaseN(self: *const MateQc, out: *std.Io.Writer, groups: Groups) !void {
        var grade: Grade = .pass;
        var scan = groups;
        while (scan.next()) |group| grade = grade.worst(.above(self.groupN(group) orelse 0, N_WARN, N_FAIL));
        try beginModule(out, "Per base N content", grade, "Base\tN-Count");
        var rows = groups;
        while (rows.next()) |group| {
            try writeLabel(out, group.first + 1, group.last);
            try writeCell(out, self.groupN(group));
            try out.writeByte('\n');
        }
        try endModule(out);
    }

    fn writeTableLengths(self: *const MateQc, out: *std.Io.Writer) !void {
        var lengths: usize = 0;
        for (self.lengths) |count| lengths += @intFromBool(count != 0);
        const grade: Grade = if (self.lengths[0] != 0) .fail else if (lengths > 1) .warn else .pass;
        try beginModule(out, "Sequence Length Distribution", grade, "Length\tCount");
        var bins = self.lengthBins();
        while (bins.next()) |bin| {
            try writeLabel(out, bin.first, bin.last);
            try writeCell(out, @floatFromInt(bin.count));
            try out.writeByte('\n');
        }
        try endModule(out);
    }

    fn writeTableAdapters(self: *const MateQc, out: *std.Io.Writer) !void {
        var rows = self.adapterRows();
        if (rows.groups == null) {
            try beginModule(out, "Adapter Content", .warn, null);
            return endModule(out);
        }
        var grade: Grade = .pass;
        var scan = rows;
        while (scan.next()) |row| {
            for (row.values) |value| grade = grade.worst(.above(value, ADAPTER_WARN, ADAPTER_FAIL));
        }
        try beginModule(out, "Adapter Content", grade, ADAPTER_HEADER);
        while (rows.next()) |row| {
            try writeLabel(out, row.group.first + 1, row.group.last);
            for (row.values) |value| try writeCell(out, value);
            try out.writeByte('\n');
        }
        try endModule(out);
    }
};

const Q40_SCORE = 40;
const QUALITY_LOWER_WARN = 10;
const QUALITY_LOWER_FAIL = 5;
const QUALITY_MEDIAN_WARN = 25;
const QUALITY_MEDIAN_FAIL = 20;
const TILE_WARN = 5;
const TILE_FAIL = 10;
const READ_QUALITY_WARN = 27;
const READ_QUALITY_FAIL = 20;
const CONTENT_WARN = 10;
const CONTENT_FAIL = 20;
const GC_WARN = 15;
const GC_FAIL = 30;
const N_WARN = 5;
const N_FAIL = 20;
const DUPLICATION_WARN = 70;
const DUPLICATION_FAIL = 50;
const OVERREPRESENTED_FAIL = 1;
const ADAPTER_WARN = 5;
const ADAPTER_FAIL = 10;

const ADAPTER_HEADER = header: {
    var columns: []const u8 = "Position";
    for (ADAPTERS) |adapter| columns = columns ++ "\t" ++ adapter.name;
    break :header columns;
};

const Grade = enum {
    pass,
    warn,
    fail,

    fn above(value: f64, warn: f64, fail: f64) Grade {
        return if (value > fail) .fail else if (value > warn) .warn else .pass;
    }

    fn below(value: f64, warn: f64, fail: f64) Grade {
        return if (value < fail) .fail else if (value < warn) .warn else .pass;
    }

    fn atMost(value: f64, warn: f64, fail: f64) Grade {
        return if (value <= fail) .fail else if (value <= warn) .warn else .pass;
    }

    fn worst(self: Grade, other: Grade) Grade {
        return @enumFromInt(@max(@intFromEnum(self), @intFromEnum(other)));
    }
};

fn beginModule(out: *std.Io.Writer, name: []const u8, grade: Grade, header: ?[]const u8) !void {
    try out.print(">>{s}\t{t}\n", .{ name, grade });
    if (header) |columns| try out.print("#{s}\n", .{columns});
}

fn endModule(out: *std.Io.Writer) !void {
    try out.writeAll(">>END_MODULE\n");
}

fn writeLabel(out: *std.Io.Writer, first: usize, last: usize) !void {
    var buffer: [48]u8 = undefined;
    try out.writeAll(rangeLabel(&buffer, first, last));
}

fn writeCell(out: *std.Io.Writer, value: ?f64) !void {
    var buffer: [JAVA_DOUBLE_BYTES]u8 = undefined;
    try out.writeByte('\t');
    try out.writeAll(javaDouble(&buffer, value orelse std.math.nan(f64)));
}

fn writeCountRow(out: *std.Io.Writer, label: usize, count: u64) !void {
    try out.print("{d}", .{label});
    try writeCell(out, @floatFromInt(count));
    try out.writeByte('\n');
}

const JAVA_DOUBLE_BYTES = std.fmt.float.min_buffer_size;
const LENGTH_TEXT_BYTES = JAVA_DOUBLE_BYTES + 4;

fn javaDouble(buffer: *[JAVA_DOUBLE_BYTES]u8, value: f64) []const u8 {
    if (std.math.isNan(value)) return "NaN";
    if (std.math.isInf(value)) return if (value > 0) "Infinity" else "-Infinity";
    const magnitude = @abs(value);
    if (magnitude == 0 or (magnitude >= 1e-3 and magnitude < 1e7)) {
        const text = std.fmt.float.render(buffer, value, .{ .mode = .decimal }) catch unreachable;
        if (std.mem.indexOfScalar(u8, text, '.') != null) return text;
        buffer[text.len..][0..2].* = ".0".*;
        return buffer[0 .. text.len + 2];
    }
    var scientific: [JAVA_DOUBLE_BYTES]u8 = undefined;
    const text = std.fmt.float.render(&scientific, value, .{}) catch unreachable;
    const exponent = std.mem.indexOfScalar(u8, text, 'e').?;
    const fraction = if (std.mem.indexOfScalar(u8, text[0..exponent], '.') == null) ".0" else "";
    return std.fmt.bufPrint(buffer, "{s}{s}E{s}", .{ text[0..exponent], fraction, text[exponent + 1 ..] }) catch unreachable;
}

fn formatLength(buffer: *[LENGTH_TEXT_BYTES]u8, bases: u64) []const u8 {
    var length: f64 = @floatFromInt(bases);
    var unit: []const u8 = " bp";
    for ([_]struct { f64, []const u8 }{ .{ 1e9, " Gbp" }, .{ 1e6, " Mbp" }, .{ 1e3, " kbp" } }) |scale| {
        if (length < scale[0]) continue;
        length /= scale[0];
        unit = scale[1];
        break;
    }
    var digits: [JAVA_DOUBLE_BYTES]u8 = undefined;
    const text = javaDouble(&digits, length);
    var last: usize = 0;
    for (text, 0..) |char, index| {
        last = index;
        if (char == '.') break;
    }
    if (last + 1 < text.len and text[last + 1] != '0') {
        last += 1;
    } else if (last > 0 and text[last] == '.') {
        last -= 1;
    }
    return std.fmt.bufPrint(buffer, "{s}{s}", .{ text[0 .. last + 1], unit }) catch unreachable;
}

fn gcDeviation(gc: *const [GC_BINS]f64) f64 {
    var first_mode: usize = 0;
    var highest: f64 = 0;
    var total: f64 = 0;
    for (gc, 0..) |count, bin| {
        total += count;
        if (count <= highest) continue;
        highest = count;
        first_mode = bin;
    }
    const floor = gc[first_mode] - gc[first_mode] / 10.0;
    var mode_sum: f64 = 0;
    var mode_bins: f64 = 0;
    var bin = first_mode;
    const fell_off_top = while (bin < GC_BINS) : (bin += 1) {
        if (!(gc[bin] > floor)) break false;
        mode_sum += @floatFromInt(bin);
        mode_bins += 1;
    } else true;
    bin = first_mode;
    const fell_off_bottom = while (bin > 0) {
        bin -= 1;
        if (!(gc[bin] > floor)) break false;
        mode_sum += @floatFromInt(bin);
        mode_bins += 1;
    } else true;
    const mode: f64 = if (fell_off_top or fell_off_bottom) @floatFromInt(first_mode) else mode_sum / mode_bins;
    var variance: f64 = 0;
    for (gc, 0..) |count, index| variance += std.math.pow(f64, @as(f64, @floatFromInt(index)) - mode, 2) * count;
    const deviation = @sqrt(variance / (total - 1));
    const scale = 1.0 / @sqrt(2 * std.math.pi * deviation * deviation);
    var difference: f64 = 0;
    for (gc, 0..) |count, index| {
        const offset = @as(f64, @floatFromInt(index)) - mode;
        const density = scale * std.math.pow(f64, std.math.e, 0 - std.math.pow(f64, offset, 2) / (2 * deviation * deviation));
        difference += @abs(density * total - count);
    }
    return difference / total * 100;
}

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
    try json.write(rangeLabel(&buffer, first, last));
}

fn rangeLabel(buffer: *[48]u8, first: usize, last: usize) []const u8 {
    return if (first == last)
        std.fmt.bufPrint(buffer, "{d}", .{first}) catch unreachable
    else
        std.fmt.bufPrint(buffer, "{d}-{d}", .{ first, last }) catch unreachable;
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

const KEY_LIMIT = 100_000;
const KEY_LENGTH = 50;
const KEY_WORDS = 7;
const SLOT_LOAD_PERCENT = 80;
const SLOT_MIN = 1024;
const CODE_BITS = 5;
const CODE_MASK = (1 << CODE_BITS) - 1;
const PACKED_WORDS = 4;
const OVERREPRESENTED_PERCENT = 0.1;
const OVERREPRESENTED_MAX = 1000;
const LEVEL_LABELS = [_][]const u8{ "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+" };
pub const CONTAMINANT_FILE_LIMIT = 1 << 20;

const LEVEL_EDGES = [_]u64{ 9, 49, 99, 499, 999, 4999, 9999 };
const EXACT_MATCH_MIN = 8;
const EXACT_MATCH_END = 20;
const WINDOW_RUN_MIN = 21;
const KEY_SECRETS = [_]u64{
    0xa076_1d64_78bd_642f, 0xe703_7ed1_a0b4_28db, 0x8ebc_6af0_9c88_c6e3, 0x5899_65cc_7537_4cc3,
    0x1d8e_4e27_c47d_124f, 0x9e37_79b9_7f4a_7c15, 0xc2b2_ae3d_27d4_eb4f, 0x1656_67b1_9e37_79f9,
};

const SequenceKey = struct {
    words: [KEY_WORDS]u64 = @splat(0),
    hash: u64 = 0,

    const WORD_MASKS: [KEY_WORDS]u64 = masks: {
        var bytes: [KEY_WORDS * 8]u8 = @splat(CODE_MASK);
        for (bytes[KEY_LENGTH..]) |*byte| byte.* = 0;
        break :masks @bitCast(bytes);
    };

    fn of(sequence: []const u8) SequenceKey {
        var words: [KEY_WORDS]u64 = undefined;
        if (sequence.len >= KEY_WORDS * 8) {
            for (&words, WORD_MASKS, 0..) |*word, mask, index| {
                word.* = std.mem.readInt(u64, sequence[index * 8 ..][0..8], .little) & mask;
            }
        } else {
            var bytes: [KEY_WORDS * 8]u8 = @splat(0);
            const length = @min(sequence.len, KEY_LENGTH);
            @memcpy(bytes[0..length], sequence[0..length]);
            for (&words, WORD_MASKS, 0..) |*word, mask, index| {
                word.* = std.mem.readInt(u64, bytes[index * 8 ..][0..8], .little) & mask;
            }
        }
        return .{ .words = words, .hash = mix(
            mix(words[0] ^ KEY_SECRETS[0], words[1] ^ KEY_SECRETS[1]) ^
                mix(words[2] ^ KEY_SECRETS[2], words[3] ^ KEY_SECRETS[3]) ^
                mix(words[4] ^ KEY_SECRETS[4], words[5] ^ KEY_SECRETS[5]) ^ KEY_SECRETS[6],
            words[6] ^ KEY_SECRETS[7],
        ) };
    }

    fn mix(a: u64, b: u64) u64 {
        const product = @as(u128, a) * b;
        return @as(u64, @truncate(product)) ^ @as(u64, @truncate(product >> 64));
    }

    fn tag(self: SequenceKey) u32 {
        return @truncate(self.hash >> 32);
    }

    fn pack(self: SequenceKey) [PACKED_WORDS]u64 {
        var packed_words: [PACKED_WORDS + 1]u64 = @splat(0);
        for (self.words, 0..) |word, chunk| {
            var codes = word;
            codes = (codes & 0x001f_001f_001f_001f) | ((codes >> 3) & 0x03e0_03e0_03e0_03e0);
            codes = (codes & 0x0000_03ff_0000_03ff) | ((codes >> 6) & 0x000f_fc00_000f_fc00);
            codes = (codes & 0x0000_0000_000f_ffff) | ((codes >> 12) & 0x0000_00ff_fff0_0000);
            const bit = chunk * 8 * CODE_BITS;
            packed_words[bit / 64] |= codes << @intCast(bit % 64);
            if (bit % 64 > 64 - 8 * CODE_BITS) packed_words[bit / 64 + 1] |= codes >> @intCast(64 - bit % 64);
        }
        std.debug.assert(packed_words[PACKED_WORDS] == 0);
        return packed_words[0..PACKED_WORDS].*;
    }
};

const PackedKey = [PACKED_WORDS]u64;

fn decodeKey(words: PackedKey, buffer: *[KEY_LENGTH]u8) []const u8 {
    for (buffer, 0..) |*byte, index| {
        const bit = index * CODE_BITS;
        var code = words[bit / 64] >> @intCast(bit % 64);
        if (bit % 64 > 64 - CODE_BITS) code |= words[bit / 64 + 1] << @intCast(64 - bit % 64);
        code &= CODE_MASK;
        if (code == 0) return buffer[0..index];
        byte.* = 0x40 | @as(u8, @intCast(code));
    }
    return buffer;
}

const Levels = struct {
    bins: [LEVEL_LABELS.len]f64,
    deduplicated: f64,
    total: f64,
};

const SequenceCounts = struct {
    records: std.ArrayList(Record) = .empty,
    slots: []u64 = &.{},
    limit: u32 = KEY_LIMIT,
    reads: u64 = 0,
    reads_at_limit: u64 = 0,
    levels: ?Levels = null,

    const Record = struct {
        words: PackedKey,
        reads: u64,
    };

    fn deinit(self: *SequenceCounts, allocator: std.mem.Allocator) void {
        self.records.deinit(allocator);
        allocator.free(self.slots);
    }

    fn limited(self: *const SequenceCounts) bool {
        return self.records.items.len == self.limit;
    }

    inline fn reserve(self: *SequenceCounts, allocator: std.mem.Allocator) error{ OutOfMemory, Overflow }!void {
        std.debug.assert(self.levels == null);
        _ = try std.math.add(u64, self.reads, 1);
        if (self.limited()) return;
        const needed = self.records.items.len + 1;
        if (needed <= self.records.capacity and needed * 100 <= self.slots.len * SLOT_LOAD_PERCENT) return;
        return self.grow(allocator);
    }

    fn grow(self: *SequenceCounts, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        const needed = self.records.items.len + 1;
        if (needed > self.records.capacity) {
            try self.records.ensureTotalCapacityPrecise(allocator, @min(self.limit, @max(needed, self.records.capacity * 2, 256)));
        }
        if (needed * 100 > self.slots.len * SLOT_LOAD_PERCENT) try self.growSlots(allocator);
    }

    fn growSlots(self: *SequenceCounts, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        const slots = try allocator.alloc(u64, @max(SLOT_MIN, self.slots.len * 2));
        @memset(slots, 0);
        const mask = slots.len - 1;
        for (self.slots) |slot| {
            if (slot == 0) continue;
            var index: usize = @as(u32, @truncate(slot >> 32)) & mask;
            while (slots[index] != 0) index = (index + 1) & mask;
            slots[index] = slot;
        }
        allocator.free(self.slots);
        self.slots = slots;
    }

    fn prefetch(self: *const SequenceCounts, key: SequenceKey) void {
        if (self.slots.len == 0) return;
        @prefetch(&self.slots[key.tag() & (self.slots.len - 1)], .{ .rw = .write });
    }

    fn add(self: *SequenceCounts, key: SequenceKey) void {
        std.debug.assert(self.levels == null and self.slots.len != 0);
        const reads = self.reads + 1;
        const full = self.limited();
        const mask = self.slots.len - 1;
        const tag = key.tag();
        var packed_key: ?PackedKey = null;
        var index: usize = tag & mask;
        while (true) : (index = (index + 1) & mask) {
            const slot = self.slots[index];
            if (slot == 0) break;
            if (@as(u32, @truncate(slot >> 32)) != tag) continue;
            const record = &self.records.items[@as(u32, @truncate(slot)) - 1];
            const wanted = packed_key orelse key.pack();
            packed_key = wanted;
            if (!std.mem.eql(u64, &record.words, &wanted)) continue;
            record.reads += 1;
            if (!full) self.reads_at_limit = reads;
            self.reads = reads;
            return;
        }
        self.reads = reads;
        if (full) return;
        self.records.appendAssumeCapacity(.{ .words = packed_key orelse key.pack(), .reads = 1 });
        self.slots[index] = @as(u64, tag) << 32 | self.records.items.len;
        self.reads_at_limit = reads;
    }

    fn finish(self: *SequenceCounts, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        if (self.levels != null) return;
        const counts = try allocator.alloc(u64, self.records.items.len);
        defer allocator.free(counts);
        for (self.records.items, counts) |record, *count| count.* = record.reads;
        std.sort.pdq(u64, counts, {}, std.sort.asc(u64));
        var levels: Levels = .{ .bins = @splat(0), .deduplicated = 0, .total = 0 };
        var start: usize = 0;
        while (start < counts.len) {
            const level = counts[start];
            var end = start + 1;
            while (end < counts.len and counts[end] == level) end += 1;
            const corrected = correctedCount(self.reads_at_limit, self.reads, level, end - start);
            const weighted = corrected * @as(f64, @floatFromInt(level));
            levels.deduplicated += corrected;
            levels.total += weighted;
            levels.bins[levelBin(level)] += weighted;
            start = end;
        }
        self.levels = levels;
    }

    fn writeDuplication(self: *const SequenceCounts, json: *std.json.Stringify) !void {
        const levels = self.levels.?;
        try json.beginObject();
        try json.objectField("reads");
        try json.write(self.reads);
        try json.objectField("stored_sequences");
        try json.write(self.records.items.len);
        try json.objectField("reads_at_limit");
        try json.write(self.reads_at_limit);
        try json.objectField("total_deduplicated_percentage");
        try json.write(deduplicatedPercentage(levels));
        try json.objectField("levels");
        try json.beginArray();
        for (LEVEL_LABELS, levels.bins) |label, bin| {
            try json.beginObject();
            try json.objectField("level");
            try json.write(label);
            try json.objectField("percentage");
            try json.write(levelPercentage(levels, bin));
            try json.endObject();
        }
        try json.endArray();
        try json.endObject();
    }

    fn writeTableDuplication(self: *const SequenceCounts, out: *std.Io.Writer) !void {
        const levels = self.levels.?;
        const deduplicated = deduplicatedPercentage(levels);
        try beginModule(out, "Sequence Duplication Levels", .below(deduplicated, DUPLICATION_WARN, DUPLICATION_FAIL), null);
        try out.writeAll("#Total Deduplicated Percentage");
        try writeCell(out, deduplicated);
        try out.writeAll("\n#Duplication Level\tPercentage of total\n");
        for (LEVEL_LABELS, levels.bins) |label, bin| {
            try out.writeAll(label);
            try writeCell(out, levelPercentage(levels, bin));
            try out.writeByte('\n');
        }
        try endModule(out);
    }

    fn deduplicatedPercentage(levels: Levels) f64 {
        return if (levels.total == 0) 100.0 else levels.deduplicated / levels.total * 100.0;
    }

    fn levelPercentage(levels: Levels, bin: f64) ?f64 {
        return if (levels.total == 0) null else bin / levels.total * 100.0;
    }

    fn percentage(self: *const SequenceCounts, reads: u64) f64 {
        return @as(f64, @floatFromInt(reads)) / @as(f64, @floatFromInt(self.reads)) * 100.0;
    }

    fn moreReads(self: *const SequenceCounts, a: u32, b: u32) bool {
        const first = self.records.items[a].reads;
        const second = self.records.items[b].reads;
        return first > second or (first == second and a < b);
    }

    fn overrepresentedRows(self: *const SequenceCounts, rows: *[OVERREPRESENTED_MAX]u32) []u32 {
        var row_count: usize = 0;
        for (self.records.items, 0..) |record, index| {
            if (!(self.percentage(record.reads) > OVERREPRESENTED_PERCENT)) continue;
            rows[row_count] = @intCast(index);
            row_count += 1;
        }
        std.sort.pdq(u32, rows[0..row_count], self, moreReads);
        return rows[0..row_count];
    }

    fn writeOverrepresented(self: *const SequenceCounts, json: *std.json.Stringify, contaminants: ?*const ContaminantList) !void {
        var rows: [OVERREPRESENTED_MAX]u32 = undefined;
        try json.beginArray();
        for (self.overrepresentedRows(&rows)) |index| {
            const record = self.records.items[index];
            var buffer: [KEY_LENGTH]u8 = undefined;
            const key = decodeKey(record.words, &buffer);
            try json.beginObject();
            try json.objectField("sequence");
            try json.write(key);
            try json.objectField("count");
            try json.write(record.reads);
            try json.objectField("percentage");
            try json.write(self.percentage(record.reads));
            try json.objectField("possible_source");
            const hit = if (contaminants) |list| list.bestHit(key) else null;
            if (hit) |found| {
                try json.beginObject();
                try json.objectField("name");
                try json.write(found.name);
                try json.objectField("identity");
                try json.write(found.identity);
                try json.objectField("length");
                try json.write(found.length);
                try json.endObject();
            } else try json.write(null);
            try json.endObject();
        }
        try json.endArray();
    }

    fn writeTableOverrepresented(self: *const SequenceCounts, out: *std.Io.Writer, contaminants: ?*const ContaminantList) !void {
        var rows: [OVERREPRESENTED_MAX]u32 = undefined;
        const found = self.overrepresentedRows(&rows);
        const grade: Grade = if (found.len == 0)
            .pass
        else if (self.percentage(self.records.items[found[0]].reads) > OVERREPRESENTED_FAIL)
            .fail
        else
            .warn;
        try beginModule(out, "Overrepresented sequences", grade, if (found.len == 0) null else "Sequence\tCount\tPercentage\tPossible Source");
        for (found) |index| {
            const record = self.records.items[index];
            var buffer: [KEY_LENGTH]u8 = undefined;
            const key = decodeKey(record.words, &buffer);
            try out.print("{s}\t{d}", .{ key, record.reads });
            try writeCell(out, @round(self.percentage(record.reads) * 100) / 100);
            if (if (contaminants) |list| list.bestHit(key) else null) |hit| {
                try out.print("\t{s} ({d}% over {d}bp)\n", .{ hit.name, hit.identity, hit.length });
            } else try out.writeAll("\tNo Hit\n");
        }
        try endModule(out);
    }
};

const PopulationView = struct {
    counts: *const SequenceCounts,
    tiles: *const TileQuality,
    contaminants: ?*const ContaminantList,

    fn writeDuplication(self: PopulationView, json: *std.json.Stringify) !void {
        try self.counts.writeDuplication(json);
    }

    fn writeOverrepresented(self: PopulationView, json: *std.json.Stringify) !void {
        try self.counts.writeOverrepresented(json, self.contaminants);
    }
};

fn levelBin(level: u64) usize {
    const offset = level - 1;
    if (offset <= LEVEL_EDGES[0]) return @intCast(offset);
    for (LEVEL_EDGES[1..], LEVEL_EDGES[0]..) |edge, bin| {
        if (offset <= edge) return bin;
    }
    return LEVEL_LABELS.len - 1;
}

fn correctedCount(limit: u64, total: u64, level: u64, observations: usize) f64 {
    const count: f64 = @floatFromInt(observations);
    if (limit == total or total - observations < limit) return count;
    std.debug.assert(level <= total and total < 1 << 53);
    const floor = 1.0 - count / (count + 0.01);
    var unseen: f64 = 1.0;
    var numerator: f64 = @floatFromInt(total - level);
    var denominator: f64 = @floatFromInt(total);
    for (0..limit) |_| {
        unseen *= numerator / denominator;
        if (unseen < floor) {
            unseen = 0;
            break;
        }
        numerator -= 1;
        denominator -= 1;
    }
    return count / (1.0 - unseen);
}

const Hit = struct {
    name: []const u8,
    identity: u32,
    length: u32,
};

const Contaminant = struct {
    name: []const u8,
    forward: []const u8,
    reverse: []const u8,

    fn match(self: Contaminant, key: []const u8) ?Hit {
        if (key.len >= EXACT_MATCH_MIN and key.len < EXACT_MATCH_END) {
            if (std.mem.indexOf(u8, self.forward, key) != null or std.mem.indexOf(u8, self.reverse, key) != null) {
                return .{ .name = self.name, .identity = 100, .length = @intCast(key.len) };
            }
        }
        var best: ?Hit = null;
        for ([_][]const u8{ self.forward, self.reverse }) |sequence| {
            var offset: isize = @as(isize, EXACT_MATCH_END) - @as(isize, @intCast(sequence.len));
            while (offset < @as(isize, @intCast(key.len)) - EXACT_MATCH_END) : (offset += 1) {
                const hit = self.alignAt(sequence, key, offset) orelse continue;
                if (best == null or hit.length > best.?.length) best = hit;
            }
        }
        return best;
    }

    fn alignAt(self: Contaminant, sequence: []const u8, key: []const u8, offset: isize) ?Hit {
        var best: ?Hit = null;
        var mismatches: u32 = 0;
        var start: isize = 0;
        var end: isize = 0;
        for (sequence, 0..) |base, index| {
            const position = @as(isize, @intCast(index)) + offset;
            const signed_index: isize = @intCast(index);
            if (position < 0) {
                start = signed_index + 1;
                continue;
            }
            if (position >= key.len) break;
            if (base == key[@intCast(position)]) {
                end = signed_index;
                continue;
            }
            mismatches += 1;
            if (mismatches > 1) {
                best = better(best, self.run(start, end, mismatches - 1));
                start = signed_index + 1;
                end = signed_index + 1;
                mismatches = 0;
            }
        }
        return better(best, self.run(start, end, mismatches));
    }

    fn run(self: Contaminant, start: isize, end: isize, mismatches: u32) ?Hit {
        const length = 1 + end - start;
        if (length < WINDOW_RUN_MIN) return null;
        const bases: u32 = @intCast(length);
        return .{ .name = self.name, .identity = (bases - mismatches) * 100 / bases, .length = bases };
    }

    fn better(best: ?Hit, hit: ?Hit) ?Hit {
        const found = hit orelse return best;
        const current = best orelse return found;
        if (current.length < found.length) return found;
        if (current.length == found.length and current.identity < found.identity) return found;
        return current;
    }
};

pub const ContaminantList = struct {
    arena: std.heap.ArenaAllocator,
    items: []const Contaminant,

    pub const ParseError = error{ OutOfMemory, InvalidBase };

    pub fn parse(allocator: std.mem.Allocator, text: []const u8, failed_line: *usize) ParseError!ContaminantList {
        var arena = std.heap.ArenaAllocator.init(allocator);
        errdefer arena.deinit();
        const owned = arena.allocator();
        var items: std.ArrayList(Contaminant) = .empty;
        var lines = std.mem.splitScalar(u8, text, '\n');
        var line_number: usize = 0;
        while (lines.next()) |raw_line| {
            line_number += 1;
            const line = std.mem.trimEnd(u8, raw_line, "\r");
            if (std.mem.startsWith(u8, line, "#")) continue;
            if (isBlank(line)) continue;
            var fields: [2][]const u8 = undefined;
            const count = splitTabs(line, &fields) orelse continue;
            if (count != 2) continue;
            const forward = try owned.alloc(u8, fields[1].len);
            const reverse = try owned.alloc(u8, fields[1].len);
            for (fields[1], forward, 0..) |byte, *upper, index| {
                upper.* = std.ascii.toUpper(byte);
                reverse[reverse.len - 1 - index] = switch (upper.*) {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    else => {
                        failed_line.* = line_number;
                        return error.InvalidBase;
                    },
                };
            }
            try items.append(owned, .{ .name = try owned.dupe(u8, fields[0]), .forward = forward, .reverse = reverse });
        }
        return .{ .arena = arena, .items = items.items };
    }

    fn isBlank(line: []const u8) bool {
        for (line) |byte| {
            if (byte > ' ') return false;
        }
        return true;
    }

    fn splitTabs(line: []const u8, fields: *[2][]const u8) ?usize {
        var count: usize = 0;
        var parts = std.mem.splitScalar(u8, std.mem.trimEnd(u8, line, "\t"), '\t');
        var first = true;
        while (parts.next()) |part| {
            if (part.len == 0 and !first) continue;
            first = false;
            if (count == fields.len) return null;
            fields[count] = part;
            count += 1;
        }
        return count;
    }

    pub fn deinit(self: *ContaminantList) void {
        self.arena.deinit();
    }

    pub fn bestHit(self: *const ContaminantList, key: []const u8) ?Hit {
        var best: ?Hit = null;
        for (self.items) |contaminant| {
            const hit = contaminant.match(key) orelse continue;
            if (best == null or hit.length > best.?.length or
                (hit.length == best.?.length and hit.identity > best.?.identity)) best = hit;
        }
        return best;
    }
};

const TILE_LIMIT = 2500;
const TILE_SAMPLE_ALL = 10_000;
const TILE_SAMPLE_EVERY = 10;
const TILE_FLUSH_READS = 1 << 24;
const NO_TILE_FIELD = "the read headers have no tile field";
const BAD_TILE = "a read header has no valid tile number";
const TOO_MANY_TILES = "more than 2500 tiles";
const NO_TILE_READS = "no read has a quality string";

const TileQuality = struct {
    tiles: std.ArrayList(Tile) = .empty,
    index: std.AutoHashMapUnmanaged(i32, u32) = .empty,
    reads: u64 = 0,
    countdown: u8 = TILE_SAMPLE_EVERY,
    last: ?u32 = null,
    longest: usize = 0,
    field: ?u3 = null,
    stopped: ?[]const u8 = null,

    const Tile = struct {
        number: i32,
        pending: u32 = 0,
        sums: std.ArrayList(u32) = .empty,
        lengths: std.ArrayList(u32) = .empty,
        carried_sums: std.ArrayList(u64) = .empty,
        carried_lengths: std.ArrayList(u64) = .empty,

        fn deinit(self: *Tile, allocator: std.mem.Allocator) void {
            self.sums.deinit(allocator);
            self.lengths.deinit(allocator);
            self.carried_sums.deinit(allocator);
            self.carried_lengths.deinit(allocator);
        }

        fn flush(self: *Tile, allocator: std.mem.Allocator) error{OutOfMemory}!void {
            try self.carried_sums.ensureTotalCapacity(allocator, self.sums.items.len);
            try self.carried_lengths.ensureTotalCapacity(allocator, self.lengths.items.len);
            self.carried_sums.appendNTimesAssumeCapacity(0, self.sums.items.len - self.carried_sums.items.len);
            self.carried_lengths.appendNTimesAssumeCapacity(0, self.lengths.items.len - self.carried_lengths.items.len);
            for (self.carried_sums.items, self.sums.items) |*total, *partial| {
                total.* += partial.*;
                partial.* = 0;
            }
            for (self.carried_lengths.items, self.lengths.items) |*total, *count| {
                total.* += count.*;
                count.* = 0;
            }
            self.pending = 0;
        }

        fn qualitySum(self: *const Tile, position: usize) u64 {
            const carried = if (position < self.carried_sums.items.len) self.carried_sums.items[position] else 0;
            return carried + self.sums.items[position];
        }

        fn readsOfLength(self: *const Tile, length: usize) u64 {
            const carried = if (length < self.carried_lengths.items.len) self.carried_lengths.items[length] else 0;
            return carried + self.lengths.items[length];
        }
    };

    fn deinit(self: *TileQuality, allocator: std.mem.Allocator) void {
        for (self.tiles.items) |*tile| tile.deinit(allocator);
        self.tiles.deinit(allocator);
        self.index.deinit(allocator);
    }

    fn stop(self: *TileQuality, allocator: std.mem.Allocator, reason: []const u8) void {
        self.deinit(allocator);
        self.* = .{ .reads = self.reads, .stopped = reason };
    }

    inline fn add(self: *TileQuality, allocator: std.mem.Allocator, header: []const u8, quality: []const u8) error{OutOfMemory}!void {
        if (self.stopped != null or quality.len == 0) return;
        self.reads += 1;
        if (self.reads > TILE_SAMPLE_ALL) {
            self.countdown -= 1;
            if (self.countdown != 0) return;
            self.countdown = TILE_SAMPLE_EVERY;
        }
        return self.addSampled(allocator, header, quality);
    }

    fn addSampled(self: *TileQuality, allocator: std.mem.Allocator, header: []const u8, quality: []const u8) error{OutOfMemory}!void {
        if (quality.len > MAX_CYCLES) return self.stop(allocator, TOO_LONG);
        const number = self.tileNumber(header) catch |err| return self.stop(allocator, switch (err) {
            error.NoTileField => NO_TILE_FIELD,
            error.BadTile => BAD_TILE,
        });
        const cached = if (self.last) |position| self.tiles.items[position].number == number else false;
        const found = if (cached) self.last else self.index.get(number);
        const tile = if (found) |position| tile: {
            self.last = position;
            break :tile &self.tiles.items[position];
        } else created: {
            if (self.tiles.items.len > TILE_LIMIT) return self.stop(allocator, TOO_MANY_TILES);
            try self.tiles.ensureUnusedCapacity(allocator, 1);
            try self.index.ensureUnusedCapacity(allocator, 1);
            self.index.putAssumeCapacityNoClobber(number, @intCast(self.tiles.items.len));
            self.tiles.appendAssumeCapacity(.{ .number = number });
            self.last = @intCast(self.tiles.items.len - 1);
            break :created &self.tiles.items[self.tiles.items.len - 1];
        };
        if (tile.sums.items.len < quality.len) {
            try tile.sums.ensureTotalCapacity(allocator, quality.len);
            try tile.lengths.ensureTotalCapacity(allocator, quality.len + 1);
            tile.sums.appendNTimesAssumeCapacity(0, quality.len - tile.sums.items.len);
            tile.lengths.appendNTimesAssumeCapacity(0, quality.len + 1 - tile.lengths.items.len);
        }
        if (tile.pending == TILE_FLUSH_READS) try tile.flush(allocator);
        for (quality, tile.sums.items[0..quality.len]) |byte, *sum| sum.* += byte - PHRED_OFFSET;
        tile.lengths.items[quality.len] += 1;
        tile.pending += 1;
        self.longest = @max(self.longest, quality.len);
    }

    fn tileNumber(self: *TileQuality, header: []const u8) error{ NoTileField, BadTile }!i32 {
        if (self.field) |field| {
            const text = headerField(header, field) orelse return error.BadTile;
            return javaInt(text) orelse error.BadTile;
        }
        var fields: [7][]const u8 = undefined;
        var count: usize = 0;
        var rest = header;
        var trailing_empty: usize = 0;
        while (true) {
            const end = std.mem.indexOfScalar(u8, rest, ':') orelse rest.len;
            const field = rest[0..end];
            if (count < fields.len) fields[count] = field;
            trailing_empty = if (count > 0 and field.len == 0) trailing_empty + 1 else 0;
            count += 1;
            if (end == rest.len) break;
            rest = rest[end + 1 ..];
        }
        count -= trailing_empty;
        const field: u3 = if (count >= 7) 4 else if (count >= 5) 2 else return error.NoTileField;
        self.field = field;
        return javaInt(fields[field]) orelse error.BadTile;
    }

    fn headerField(header: []const u8, field: usize) ?[]const u8 {
        const Chunk = @Vector(16, u8);
        var seen: usize = 0;
        var start: usize = 0;
        var offset: usize = 0;
        while (offset + 16 <= header.len) : (offset += 16) {
            const chunk: Chunk = header[offset..][0..16].*;
            var colons: u16 = @bitCast(chunk == @as(Chunk, @splat(':')));
            while (colons != 0) : (colons &= colons - 1) {
                const at = offset + @ctz(colons);
                if (seen == field) return header[start..at];
                seen += 1;
                start = at + 1;
            }
        }
        for (header[offset..], offset..) |byte, at| {
            if (byte != ':') continue;
            if (seen == field) return header[start..at];
            seen += 1;
            start = at + 1;
        }
        return if (seen == field) header[start..] else null;
    }

    fn javaInt(text: []const u8) ?i32 {
        const digits = if (text.len > 0 and (text[0] == '+' or text[0] == '-')) text[1..] else text;
        if (digits.len == 0) return null;
        var value: i64 = 0;
        for (digits) |byte| {
            if (byte < '0' or byte > '9') return null;
            value = value * 10 + (byte - '0');
            if (value > 1 << 31) return null;
        }
        if (text[0] == '-') value = -value;
        return std.math.cast(i32, value);
    }

    fn unavailable(self: *const TileQuality) ?[]const u8 {
        if (self.stopped) |reason| return reason;
        if (self.longest == 0) return NO_TILE_READS;
        return null;
    }

    fn mean(tile: *const Tile, group: Group, counts: []const u64) f64 {
        var total: f64 = 0;
        var positions: usize = 0;
        for (group.first..group.last) |position| {
            if (position >= tile.sums.items.len or counts[position] == 0) continue;
            total += @as(f64, @floatFromInt(tile.qualitySum(position))) / @as(f64, @floatFromInt(counts[position]));
            positions += 1;
        }
        return if (positions == 0) 0 else total / @as(f64, @floatFromInt(positions));
    }

    fn positionCounts(tile: *const Tile, counts: *[MAX_CYCLES]u64) void {
        var reads: u64 = 0;
        var position = tile.lengths.items.len;
        while (position > 0) {
            position -= 1;
            if (position < MAX_CYCLES) counts[position] = reads;
            reads += tile.readsOfLength(position);
        }
    }

    fn lessNumber(self: *const TileQuality, a: u16, b: u16) bool {
        return self.tiles.items[a].number < self.tiles.items[b].number;
    }

    fn writeJson(self: *const TileQuality, json: *std.json.Stringify) !void {
        if (self.unavailable() != null) return json.write(null);
        const table: TileTable = .init(self);
        try json.beginObject();
        try json.objectField("positions");
        try json.beginArray();
        var labels = table.groups;
        while (labels.next()) |group| try writeRangeLabel(json, group.first + 1, group.last);
        try json.endArray();
        try json.objectField("tiles");
        try json.beginArray();
        var deviations: [MAX_CYCLES]f64 = undefined;
        for (0..table.count) |rank| {
            try json.beginObject();
            try json.objectField("tile");
            try json.write(table.deviations(rank, &deviations));
            try json.objectField("deviations");
            try json.beginArray();
            for (deviations[0..table.group_count]) |deviation| try json.write(deviation);
            try json.endArray();
            try json.endObject();
        }
        try json.endArray();
        try json.endObject();
    }

    fn writeTable(self: *const TileQuality, out: *std.Io.Writer) !void {
        const table: TileTable = .init(self);
        try beginModule(out, "Per tile sequence quality", .above(table.maxDeviation(), TILE_WARN, TILE_FAIL), "Tile\tBase\tMean");
        var deviations: [MAX_CYCLES]f64 = undefined;
        for (0..table.count) |rank| {
            const tile = table.deviations(rank, &deviations);
            var labels = table.groups;
            var index: usize = 0;
            while (labels.next()) |group| : (index += 1) {
                try out.print("{d}\t", .{tile});
                try writeLabel(out, group.first + 1, group.last);
                try writeCell(out, deviations[index]);
                try out.writeByte('\n');
            }
        }
        try endModule(out);
    }
};

const TileTable = struct {
    quality: *const TileQuality,
    order: [TILE_LIMIT + 1]u16,
    count: usize,
    groups: Groups,
    group_count: usize,
    averages: [MAX_CYCLES]f64,

    fn init(quality: *const TileQuality) TileTable {
        std.debug.assert(quality.unavailable() == null);
        var table: TileTable = .{
            .quality = quality,
            .order = undefined,
            .count = quality.tiles.items.len,
            .groups = .init(quality.longest),
            .group_count = 0,
            .averages = @splat(0),
        };
        const sorted = table.order[0..table.count];
        for (sorted, 0..) |*entry, index| entry.* = @intCast(index);
        std.sort.pdq(u16, sorted, quality, TileQuality.lessNumber);
        var counts: [MAX_CYCLES]u64 = undefined;
        for (sorted) |entry| {
            const tile = &quality.tiles.items[entry];
            TileQuality.positionCounts(tile, &counts);
            var iterator = table.groups;
            var index: usize = 0;
            while (iterator.next()) |group| : (index += 1) table.averages[index] += TileQuality.mean(tile, group, &counts);
            table.group_count = index;
        }
        const tile_count: f64 = @floatFromInt(table.count);
        for (table.averages[0..table.group_count]) |*average| average.* /= tile_count;
        return table;
    }

    fn deviations(self: *const TileTable, rank: usize, out: *[MAX_CYCLES]f64) i32 {
        const tile = &self.quality.tiles.items[self.order[rank]];
        var counts: [MAX_CYCLES]u64 = undefined;
        TileQuality.positionCounts(tile, &counts);
        var iterator = self.groups;
        var index: usize = 0;
        while (iterator.next()) |group| : (index += 1) out[index] = TileQuality.mean(tile, group, &counts) - self.averages[index];
        return tile.number;
    }

    fn maxDeviation(self: *const TileTable) f64 {
        var largest: f64 = 0;
        var deviations_buffer: [MAX_CYCLES]f64 = undefined;
        for (0..self.count) |rank| {
            _ = self.deviations(rank, &deviations_buffer);
            for (deviations_buffer[0..self.group_count]) |deviation| {
                if (@abs(deviation) > largest) largest = @abs(deviation);
            }
        }
        return largest;
    }
};

const TestRead = struct { []const u8, []const u8 };

fn collect(allocator: std.mem.Allocator, reads: []const TestRead) !MateQc {
    var mate: MateQc = .{};
    errdefer mate.deinit(allocator);
    for (reads) |read| try mate.addRead(allocator, read[0], read[1]);
    return mate;
}

fn reportOf(arena: std.mem.Allocator, mate: *MateQc) !std.json.ObjectMap {
    mate.finish();
    var sequences: SequenceCounts = .{};
    try sequences.finish(arena);
    var output: std.Io.Writer.Allocating = .init(arena);
    var json: std.json.Stringify = .{ .writer = &output.writer };
    const tiles: TileQuality = .{};
    try mate.writeJson(&json, .{ .counts = &sequences, .tiles = &tiles, .contaminants = null });
    const value = try std.json.parseFromSliceLeaky(std.json.Value, arena, output.written(), .{});
    return value.object;
}

fn populationReportOf(arena: std.mem.Allocator, populations: *const MatePopulations, retained: bool) !std.json.ObjectMap {
    var output: std.Io.Writer.Allocating = .init(arena);
    var json: std.json.Stringify = .{ .writer = &output.writer };
    try populations.writeJson(&json, retained, null);
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

test "[property] - [qc]: FastQC data for probe P1 grades and summarizes it as FastQC-Rust does" {
    const allocator = std.testing.allocator;
    var populations: MatePopulations = .{};
    defer populations.deinit(allocator);
    for (PROBE_P1) |read| try populations.addRead(allocator, "@r", read[0], read[1]);
    try populations.finish(allocator);
    var output: std.Io.Writer.Allocating = .init(allocator);
    defer output.deinit();
    try populations.writeQcTable(&output.writer, false, "p1.fastq", "z-fastq test", null);

    var kept: std.Io.Writer.Allocating = .init(allocator);
    defer kept.deinit();
    var basic = false;
    var lines = std.mem.splitScalar(u8, output.written(), '\n');
    while (lines.next()) |line| {
        if (std.mem.eql(u8, line, ">>END_MODULE")) {
            basic = false;
            continue;
        }
        if (std.mem.startsWith(u8, line, ">>Basic")) basic = true;
        if (basic or std.mem.startsWith(u8, line, ">>") or std.mem.startsWith(u8, line, "##")) try kept.writer.print("{s}\n", .{line});
    }
    try std.testing.expectEqualStrings("##FastQC\tz-fastq test\n" ++
        ">>Basic Statistics\tpass\n#Measure\tValue\nFilename\tp1.fastq\nFile type\tConventional base calls\n" ++
        "Encoding\tSanger / Illumina 1.9\nTotal Sequences\t4\nTotal Bases\t27 bp\nSequence length\t0-12\n" ++
        "Mean Length\t6\nMedian Length\t8\n%GC\t61\n" ++
        ">>Per base sequence quality\tpass\n>>Per sequence quality scores\tfail\n>>Per base sequence content\tfail\n" ++
        ">>Per sequence GC content\tfail\n>>Per base N content\tfail\n>>Sequence Length Distribution\tfail\n" ++
        ">>Sequence Duplication Levels\tpass\n>>Overrepresented sequences\tfail\n>>Adapter Content\twarn\n", kept.written());
}

test "[unit] - [qc]: doubles print as Java prints them and base totals as FastQC shortens them" {
    var buffer: [JAVA_DOUBLE_BYTES]u8 = undefined;
    for ([_]struct { f64, []const u8 }{
        .{ 0, "0.0" },                                 .{ -0.0, "-0.0" },                    .{ 1, "1.0" },
        .{ 31.398000000000003, "31.398000000000003" }, .{ 0.001, "0.001" },                  .{ 9.99e-4, "9.99E-4" },
        .{ 9999999.5, "9999999.5" },                   .{ 1e7, "1.0E7" },                    .{ 12345678, "1.2345678E7" },
        .{ std.math.nan(f64), "NaN" },                 .{ -std.math.inf(f64), "-Infinity" },
    }) |case| try std.testing.expectEqualStrings(case[1], javaDouble(&buffer, case[0]));
    var text: [LENGTH_TEXT_BYTES]u8 = undefined;
    for ([_]struct { u64, []const u8 }{
        .{ 0, "0 bp" },          .{ 999, "999 bp" },          .{ 14_300, "14.3 kbp" },         .{ 14_350, "14.3 kbp" },
        .{ 2_000_000, "2 Mbp" }, .{ 1_050_000_000, "1 Gbp" }, .{ 15_600_000_000, "15.6 Gbp" },
    }) |case| try std.testing.expectEqualStrings(case[1], formatLength(&text, case[0]));
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
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    var populations: MatePopulations = .{};
    defer populations.deinit(failing.allocator());
    populations.enableRetained();
    try populations.addRead(failing.allocator(), "r", "ACGT", "IIII");
    failing.fail_index = failing.alloc_index;
    try std.testing.expectError(error.OutOfMemory, populations.addRetained(failing.allocator(), "r", "ACGT", "IIII", 2));
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
        populations.enableRetained();
        var expected: MatePopulations = .{};
        defer expected.deinit(allocator);
        for (reads[0..count]) |read| {
            try populations.addRead(allocator, "A00123:8:HFLKW:1:1101:1:2 1:N:0:ACGT", read.sequence, read.quality);
            if (read.dropped) {
                try populations.dropRetained(allocator, read.sequence, read.quality);
                continue;
            }
            try populations.addRetained(allocator, "A00123:8:HFLKW:1:1101:1:2 1:N:0:ACGT", read.sequence, read.quality, read.kept);
            try expected.addRead(allocator, "A00123:8:HFLKW:1:1101:1:2 1:N:0:ACGT", read.sequence[0..read.kept], read.quality[0..read.kept]);
        }
        try populations.finish(allocator);
        try expected.finish(allocator);
        const want = try populationReportOf(arena.allocator(), &expected, false);
        const got = try populationReportOf(arena.allocator(), &populations, true);
        inline for (.{ "per_base_quality", "per_base_content", "per_base_n", "length_distribution", "per_sequence_quality", "per_sequence_gc", "terminal_runs", "adapter_content", "duplication", "overrepresented", "per_tile_quality", "per_tile_unavailable", "unavailable" }) |module| {
            errdefer std.debug.print("{s}\n", .{module});
            const want_text = try std.json.Stringify.valueAlloc(arena.allocator(), want.get(module).?, .{});
            const got_text = try std.json.Stringify.valueAlloc(arena.allocator(), got.get(module).?, .{});
            try std.testing.expectEqualStrings(want_text, got_text);
        }
        try std.testing.expectEqual(try expected.raw.stats(), try populations.retained.stats());
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

fn countsOf(allocator: std.mem.Allocator, reads: []const []const u8) !SequenceCounts {
    var counts: SequenceCounts = .{};
    errdefer counts.deinit(allocator);
    for (reads) |read| {
        try counts.reserve(allocator);
        counts.add(.of(read));
    }
    try counts.finish(allocator);
    return counts;
}

fn sequenceReportOf(arena: std.mem.Allocator, counts: *const SequenceCounts, contaminants: ?*const ContaminantList) !std.json.ObjectMap {
    var output: std.Io.Writer.Allocating = .init(arena);
    var json: std.json.Stringify = .{ .writer = &output.writer };
    try json.beginObject();
    try json.objectField("duplication");
    try counts.writeDuplication(&json);
    try json.objectField("overrepresented");
    try counts.writeOverrepresented(&json, contaminants);
    try json.endObject();
    const value = try std.json.parseFromSliceLeaky(std.json.Value, arena, output.written(), .{});
    return value.object;
}

fn jsonNumber(value: std.json.Value) f64 {
    return switch (value) {
        .integer => |integer| @floatFromInt(integer),
        .float => |float| float,
        else => std.math.nan(f64),
    };
}

fn sharePercent(part: f64, whole: f64) f64 {
    return part / whole * 100.0;
}

fn levelValue(report: std.json.ObjectMap, bin: usize) f64 {
    return jsonNumber(report.get("duplication").?.object.get("levels").?.array.items[bin].object.get("percentage").?);
}

fn deduplicatedPercent(report: std.json.ObjectMap) f64 {
    return jsonNumber(report.get("duplication").?.object.get("total_deduplicated_percentage").?);
}

fn indexedKey(buffer: *[KEY_LENGTH]u8, index: usize) []const u8 {
    var value = index;
    for (buffer) |*base| {
        base.* = "ACGT"[value & 3];
        value >>= 2;
    }
    return buffer;
}

test "[unit] - [qc]: duplication levels count uppercase keys cut to 50 bases" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const a = arena.allocator();
    const long = "ACGTACGTAC" ** 5;

    var mixed_reads: std.ArrayList([]const u8) = .empty;
    for ([_]usize{ 7, 5, 3, 2, 2, 1, 1, 1 }, 0..) |copies, index| {
        const read = try std.fmt.allocPrint(a, "{s}{d}", .{ "GATTACA", index });
        for (0..copies) |_| try mixed_reads.append(a, read);
    }
    var mixed = try countsOf(std.testing.allocator, mixed_reads.items);
    defer mixed.deinit(std.testing.allocator);
    const mixed_report = try sequenceReportOf(a, &mixed, null);
    try std.testing.expectEqual(@as(f64, 36.36363636363637), deduplicatedPercent(mixed_report));
    try std.testing.expectEqual(sharePercent(3, 22), levelValue(mixed_report, 0));
    try std.testing.expectEqual(sharePercent(4, 22), levelValue(mixed_report, 1));
    try std.testing.expectEqual(sharePercent(7, 22), levelValue(mixed_report, 6));

    var folded = try countsOf(std.testing.allocator, &.{ long ++ "A", "acgtacgtac" ** 5 ++ "TTT", long, "", "" });
    defer folded.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(usize, 2), folded.records.items.len);
    for ([_]usize{ 50, 55, 56, 60 }) |length| {
        const lower = ("acgtacgtac" ** 6)[0..length];
        try std.testing.expectEqual(SequenceKey.of(long), SequenceKey.of(lower));
    }
    try std.testing.expect(SequenceKey.of(long[0..49]).hash != SequenceKey.of(long).hash);
    try std.testing.expectEqual(@as(u64, 3), folded.records.items[0].reads);
    for ([_][]const u8{ "acgturyswkmbdhvn" ++ "N" ** 40, "ACGTURYSWKMBDHVN", "", "T", "GATTACA" ** 7 }) |sequence| {
        const length = @min(sequence.len, KEY_LENGTH);
        var buffer: [KEY_LENGTH]u8 = undefined;
        var upper: [KEY_LENGTH]u8 = undefined;
        try std.testing.expectEqualStrings(std.ascii.upperString(upper[0..length], sequence[0..length]), decodeKey(SequenceKey.of(sequence).pack(), &buffer));
    }
    try std.testing.expect(!std.mem.eql(u64, &SequenceKey.of("ACGT").pack(), &SequenceKey.of("ACGTA").pack()));
    try std.testing.expectEqual(@as(f64, 40), deduplicatedPercent(try sequenceReportOf(a, &folded, null)));

    var empty = try countsOf(std.testing.allocator, &.{});
    defer empty.deinit(std.testing.allocator);
    const empty_report = try sequenceReportOf(a, &empty, null);
    try std.testing.expectEqual(@as(f64, 100), deduplicatedPercent(empty_report));
    try std.testing.expectEqual(std.json.Value.null, empty_report.get("duplication").?.object.get("levels").?.array.items[0].object.get("percentage").?);
    try std.testing.expectEqual(@as(usize, 0), empty_report.get("overrepresented").?.array.items.len);

    try std.testing.expectEqual(@as(usize, 8), levelBin(9));
    try std.testing.expectEqual(@as(usize, 9), levelBin(10));
    try std.testing.expectEqual(@as(usize, 9), levelBin(50));
    try std.testing.expectEqual(@as(usize, 10), levelBin(51));
    try std.testing.expectEqual(@as(usize, 14), levelBin(10000));
    try std.testing.expectEqual(@as(usize, 15), levelBin(10001));
}

test "[edge] - [qc]: past 100,000 keys new reads count only in the total and the levels are corrected" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var counts: SequenceCounts = .{};
    defer counts.deinit(std.testing.allocator);
    var buffer: [KEY_LENGTH]u8 = undefined;
    for (0..KEY_LIMIT) |index| {
        try counts.reserve(std.testing.allocator);
        counts.add(.of(indexedKey(&buffer, index)));
    }
    try std.testing.expect(counts.limited());
    try std.testing.expectEqual(@as(u64, KEY_LIMIT), counts.reads_at_limit);
    try counts.reserve(std.testing.allocator);
    counts.add(.of(indexedKey(&buffer, KEY_LIMIT)));
    try std.testing.expectEqual(@as(usize, KEY_LIMIT), counts.records.items.len);
    try std.testing.expectEqual(@as(u64, KEY_LIMIT + 1), counts.reads);
    for (0..KEY_LIMIT / 2) |index| {
        try counts.reserve(std.testing.allocator);
        counts.add(.of(indexedKey(&buffer, index)));
    }
    try std.testing.expectEqual(@as(u64, KEY_LIMIT), counts.reads_at_limit);
    try std.testing.expectEqual(@as(usize, KEY_LIMIT), counts.records.capacity);
    try std.testing.expectEqual(@as(usize, 131072), counts.slots.len);

    var probe: SequenceCounts = .{};
    defer probe.deinit(std.testing.allocator);
    for (0..KEY_LIMIT) |index| {
        try probe.reserve(std.testing.allocator);
        probe.add(.of(indexedKey(&buffer, index)));
    }
    for (0..KEY_LIMIT / 2) |index| {
        try probe.reserve(std.testing.allocator);
        probe.add(.of(indexedKey(&buffer, index)));
    }
    try probe.finish(std.testing.allocator);
    const report = try sequenceReportOf(arena.allocator(), &probe, null);
    try std.testing.expectEqual(@as(f64, 70.00002000011865), deduplicatedPercent(report));
    try std.testing.expectEqual(@as(f64, 40.0000400002373), levelValue(report, 0));
}

test "[unit] - [qc]: overrepresented keys are above 0.1% of all reads, by count then first appearance" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const a = arena.allocator();
    var reads: std.ArrayList([]const u8) = .empty;
    try reads.appendSlice(a, &.{ "TTTT", "CCCC", "CCCC", "GGGG", "TTTT", "GGGG", "GGGG" });
    for (0..993) |index| try reads.append(a, try std.fmt.allocPrint(a, "A{d}", .{index}));
    var counts = try countsOf(std.testing.allocator, reads.items);
    defer counts.deinit(std.testing.allocator);
    const rows = (try sequenceReportOf(a, &counts, null)).get("overrepresented").?.array.items;
    try std.testing.expectEqual(@as(usize, 3), rows.len);
    for (rows, [_][]const u8{ "GGGG", "TTTT", "CCCC" }) |row, sequence| {
        try std.testing.expectEqualStrings(sequence, row.object.get("sequence").?.string);
        try std.testing.expectEqual(std.json.Value.null, row.object.get("possible_source").?);
    }
    try std.testing.expectEqual(sharePercent(3, 1000), jsonNumber(rows[0].object.get("percentage").?));
}

const ADAPTER_1 = "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG";
const ADAPTER_2 = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT";

test "[unit] - [qc]: contaminant hits follow FastQC's exact and one-mismatch window rules" {
    var failed_line: usize = 0;
    const text = "# comment\n\n  \t \nIllumina Single End Adapter 1\t\t\t" ++ ADAPTER_1 ++ "\r\nthree\tfields\there\n" ++
        "Illumina Single End Adapter 2\t" ++ "caagcagaagacggcatacgagctcttccgatct" ++ "\t\n";
    var list = try ContaminantList.parse(std.testing.allocator, text, &failed_line);
    defer list.deinit();
    try std.testing.expectEqual(@as(usize, 2), list.items.len);
    try std.testing.expectEqualStrings(ADAPTER_2, list.items[1].forward);

    const Case = struct { key: []const u8, name: ?[]const u8, identity: u32 = 0, length: u32 = 0 };
    const one = ADAPTER_1[0..15] ++ "A" ++ ADAPTER_1[16..30];
    const two = ADAPTER_1[0..10] ++ "C" ++ ADAPTER_1[11..20] ++ "A" ++ ADAPTER_1[21..30];
    const reverse = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATC";
    const cases = [_]Case{
        .{ .key = ADAPTER_1[0..7], .name = null },
        .{ .key = ADAPTER_1[0..8], .name = "Illumina Single End Adapter 1", .identity = 100, .length = 8 },
        .{ .key = ADAPTER_1[0..19], .name = "Illumina Single End Adapter 1", .identity = 100, .length = 19 },
        .{ .key = ADAPTER_1[0..20], .name = null },
        .{ .key = ADAPTER_1[0..21], .name = "Illumina Single End Adapter 1", .identity = 100, .length = 21 },
        .{ .key = one, .name = "Illumina Single End Adapter 1", .identity = 96, .length = 30 },
        .{ .key = two, .name = null },
        .{ .key = reverse, .name = "Illumina Single End Adapter 1", .identity = 100, .length = 33 },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        const hit = list.bestHit(case.key);
        if (case.name) |name| {
            try std.testing.expectEqualStrings(name, hit.?.name);
            try std.testing.expectEqual(case.identity, hit.?.identity);
            try std.testing.expectEqual(case.length, hit.?.length);
        } else try std.testing.expectEqual(@as(?Hit, null), hit);
    }

    try std.testing.expectError(error.InvalidBase, ContaminantList.parse(std.testing.allocator, "ok\tACGT\nbad\tACNT\n", &failed_line));
    try std.testing.expectEqual(@as(usize, 2), failed_line);
}

fn casavaHeader(buffer: []u8, tile: i64) []const u8 {
    return std.fmt.bufPrint(buffer, "A00123:8:HFLKW:1:{d}:1000:2000 1:N:0:ACGT", .{tile}) catch unreachable;
}

fn tileRows(arena: std.mem.Allocator, tiles: *const TileQuality) !std.json.Value {
    var output: std.Io.Writer.Allocating = .init(arena);
    var json: std.json.Stringify = .{ .writer = &output.writer };
    try tiles.writeJson(&json);
    return std.json.parseFromSliceLeaky(std.json.Value, arena, output.written(), .{});
}

test "[unit] - [qc]: the tile is header field 4 of 7 or more, or field 2 of 5 or 6, as a Java int" {
    const Case = struct { headers: []const []const u8, tile: ?i32, reason: ?[]const u8 = null };
    const cases = [_]Case{
        .{ .headers = &.{"A00123:8:HFLKW:1:1101:1000:2000 1:N:0:ACGT"}, .tile = 1101 },
        .{ .headers = &.{"HWUSI-EAS100R:6:73:941:1973#0/1"}, .tile = 73 },
        .{ .headers = &.{"A:B:C:D:+0012:E:F"}, .tile = 12 },
        .{ .headers = &.{"A:B:C:D:-2147483648:E:F"}, .tile = -2147483648 },
        .{ .headers = &.{"A:B:C:D:2147483648:E:F"}, .tile = null, .reason = BAD_TILE },
        .{ .headers = &.{"A:B:C:D:+:E:F"}, .tile = null, .reason = BAD_TILE },
        .{ .headers = &.{"A:B:12:D:"}, .tile = null, .reason = NO_TILE_FIELD },
        .{ .headers = &.{"A:B:12:D:E::"}, .tile = 12 },
        .{ .headers = &.{"read1/1"}, .tile = null, .reason = NO_TILE_FIELD },
        .{ .headers = &.{ "HWUSI:6:73:941:1973", "A00123:8:HFLKW:1:1101:1000:2000" }, .tile = null, .reason = BAD_TILE },
        .{ .headers = &.{ "A00123:8:HFLKW:1:1101:1000:2000", "A:B:C" }, .tile = null, .reason = BAD_TILE },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        var tiles: TileQuality = .{};
        defer tiles.deinit(std.testing.allocator);
        for (case.headers) |header| try tiles.add(std.testing.allocator, header, "IIII");
        try std.testing.expectEqual(case.reason, tiles.stopped);
        if (case.tile) |tile| try std.testing.expectEqual(tile, tiles.tiles.items[0].number);
    }
}

test "[unit] - [qc]: tile deviations are tile means minus the average of tile means, 0 for an uncovered group" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var tiles: TileQuality = .{};
    defer tiles.deinit(std.testing.allocator);
    var buffer: [64]u8 = undefined;
    for (0..3) |_| try tiles.add(std.testing.allocator, casavaHeader(&buffer, 2101), "??");
    for (0..2) |_| try tiles.add(std.testing.allocator, casavaHeader(&buffer, 1101), "5555");
    try tiles.add(std.testing.allocator, casavaHeader(&buffer, 1101), "+");
    const value = try tileRows(arena.allocator(), &tiles);
    const rows = value.object.get("tiles").?.array.items;
    try std.testing.expectEqual(@as(usize, 4), value.object.get("positions").?.array.items.len);
    try std.testing.expectEqual(@as(i64, 1101), rows[0].object.get("tile").?.integer);
    const low = rows[0].object.get("deviations").?.array.items;
    const high = rows[1].object.get("deviations").?.array.items;
    const first_low: f64 = (20.0 + 20.0 + 10.0) / 3.0;
    try std.testing.expectEqual(first_low - (first_low + 30.0) / 2.0, numberOf(low[0]).?);
    try std.testing.expectEqual(20.0 - (20.0 + 30.0) / 2.0, numberOf(low[1]).?);
    try std.testing.expectEqual(20.0 - (20.0 + 0.0) / 2.0, numberOf(low[2]).?);
    try std.testing.expectEqual(0.0 - (20.0 + 0.0) / 2.0, numberOf(high[3]).?);
}

test "[edge] - [qc]: after 10,000 reads every tenth read counts, and a 2,502nd tile stops the module" {
    var buffer: [64]u8 = undefined;
    var tiles: TileQuality = .{};
    defer tiles.deinit(std.testing.allocator);
    for (0..10_000) |_| try tiles.add(std.testing.allocator, casavaHeader(&buffer, 1101), "I");
    try tiles.add(std.testing.allocator, casavaHeader(&buffer, 1101), "");
    for (0..9) |_| try tiles.add(std.testing.allocator, casavaHeader(&buffer, 7), "I");
    try std.testing.expectEqual(@as(usize, 1), tiles.tiles.items.len);
    try tiles.add(std.testing.allocator, casavaHeader(&buffer, 7), "I");
    try std.testing.expectEqual(@as(usize, 2), tiles.tiles.items.len);

    var many: TileQuality = .{};
    defer many.deinit(std.testing.allocator);
    for (0..TILE_LIMIT + 1) |tile| try many.add(std.testing.allocator, casavaHeader(&buffer, @intCast(tile)), "I");
    try std.testing.expectEqual(@as(?[]const u8, null), many.stopped);
    try many.add(std.testing.allocator, casavaHeader(&buffer, TILE_LIMIT + 1), "I");
    try std.testing.expectEqual(@as(?[]const u8, TOO_MANY_TILES), many.stopped);
    try std.testing.expectEqual(@as(usize, 0), many.tiles.items.len);
}

test "[edge] - [qc]: per-tile counters flushed into wide carries give the same values" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    var buffer: [64]u8 = undefined;
    var flushed: TileQuality = .{};
    defer flushed.deinit(std.testing.allocator);
    var direct: TileQuality = .{};
    defer direct.deinit(std.testing.allocator);
    const reads = [_][]const u8{ "IIII", "#", "??55", "+++++++", "I" };
    for (reads, 0..) |quality, index| {
        for ([_]*TileQuality{ &flushed, &direct }) |tiles| {
            try tiles.add(std.testing.allocator, casavaHeader(&buffer, 1101 + @as(i64, @intCast(index % 2))), quality);
        }
        for (flushed.tiles.items) |*tile| try tile.flush(std.testing.allocator);
    }
    try std.testing.expect(flushed.tiles.items[0].carried_sums.items.len != 0);
    const want = try std.json.Stringify.valueAlloc(arena.allocator(), try tileRows(arena.allocator(), &direct), .{});
    const got = try std.json.Stringify.valueAlloc(arena.allocator(), try tileRows(arena.allocator(), &flushed), .{});
    try std.testing.expectEqualStrings(want, got);
    try std.testing.expectEqualStrings("1101", TileQuality.headerField("A00123:8:HFLKW:1:1101:1000:2000", 4).?);
    try std.testing.expectEqualStrings("", TileQuality.headerField("A:B:C:D:", 4).?);
    try std.testing.expectEqual(@as(?[]const u8, null), TileQuality.headerField("A:B:C", 4));
    try std.testing.expectEqualStrings("12", TileQuality.headerField("LONG-INSTRUMENT-NAME-0123:run:flowcell:lane:12", 4).?);
}
