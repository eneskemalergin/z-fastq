//! Checked aggregate length, base-composition, and Phred+33 statistics.

const std = @import("std");
const fastq = @import("fastq.zig");
const simd = @import("simd.zig");

pub const StatsError = error{
    S005LengthMismatch,
    S006InvalidQuality,
    Overflow,
};

pub const QualityError = struct {
    byte_index: usize,
    byte: u8,
};

pub const Ratio = struct {
    numerator: u128,
    denominator: u128,

    pub fn value(self: Ratio) ?f64 {
        if (self.denominator == 0) return null;
        return @as(f64, @floatFromInt(self.numerator)) / @as(f64, @floatFromInt(self.denominator));
    }
};

pub const Ratios = struct {
    mean_length: Ratio,
    gc_fraction: Ratio,
    mean_quality: Ratio,
    q20_fraction: Ratio,
    q30_fraction: Ratio,
};

pub const StatsResult = struct {
    reads: u64,
    bases: u64,
    min_length: ?u64,
    max_length: ?u64,
    mean_length: ?f64,
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
    other_bases: u64,
    gc_fraction: ?f64,
    quality_sum: u64,
    mean_quality: ?f64,
    q20_bases: u64,
    q20_fraction: ?f64,
    q30_bases: u64,
    q30_fraction: ?f64,
};

pub const Stats = struct {
    reads: u64 = 0,
    bases: u64 = 0,
    min_length: u64 = 0,
    max_length: u64 = 0,
    a: u64 = 0,
    c: u64 = 0,
    g: u64 = 0,
    t: u64 = 0,
    n: u64 = 0,
    other_bases: u64 = 0,
    quality_sum: u64 = 0,
    q20_bases: u64 = 0,
    q30_bases: u64 = 0,
    last_quality_error: ?QualityError = null,

    pub fn addRecord(self: *Stats, record: fastq.Record) StatsError!void {
        self.last_quality_error = null;
        var quality_error: ?QualityError = null;
        const delta = recordTotals(record, &quality_error) catch |err| {
            self.last_quality_error = quality_error;
            return err;
        };
        try self.applyTotals(delta);
    }

    pub fn addStreamed(self: *Stats, streamed: *const StreamedTotals) StatsError!void {
        self.last_quality_error = null;
        if (streamed.quality_error) |quality_error| {
            self.last_quality_error = quality_error;
            return error.S006InvalidQuality;
        }
        if (streamed.overflow) return error.Overflow;
        try self.applyTotals(try totalsFromPayload(streamed.totals, streamed.sequence_len));
    }

    inline fn applyTotals(self: *Stats, delta: RecordTotals) StatsError!void {
        var updated = self.*;
        updated.reads = try add(self.reads, 1);
        updated.bases = try add(self.bases, delta.length);
        updated.a = try add(self.a, delta.a);
        updated.c = try add(self.c, delta.c);
        updated.g = try add(self.g, delta.g);
        updated.t = try add(self.t, delta.t);
        updated.n = try add(self.n, delta.n);
        updated.other_bases = try add(self.other_bases, delta.other_bases);
        updated.quality_sum = try add(self.quality_sum, delta.quality_sum);
        updated.q20_bases = try add(self.q20_bases, delta.q20_bases);
        updated.q30_bases = try add(self.q30_bases, delta.q30_bases);
        if (self.reads == 0) {
            updated.min_length = delta.length;
            updated.max_length = delta.length;
        } else {
            updated.min_length = @min(self.min_length, delta.length);
            updated.max_length = @max(self.max_length, delta.length);
        }
        self.* = updated;
    }

    pub fn takeLastQualityError(self: *Stats) ?QualityError {
        const quality_error = self.last_quality_error;
        self.last_quality_error = null;
        return quality_error;
    }

    pub fn result(self: *const Stats) StatsResult {
        const exact = self.ratios();
        return .{
            .reads = self.reads,
            .bases = self.bases,
            .min_length = if (self.reads == 0) null else self.min_length,
            .max_length = if (self.reads == 0) null else self.max_length,
            .mean_length = exact.mean_length.value(),
            .a = self.a,
            .c = self.c,
            .g = self.g,
            .t = self.t,
            .n = self.n,
            .other_bases = self.other_bases,
            .gc_fraction = exact.gc_fraction.value(),
            .quality_sum = self.quality_sum,
            .mean_quality = exact.mean_quality.value(),
            .q20_bases = self.q20_bases,
            .q20_fraction = exact.q20_fraction.value(),
            .q30_bases = self.q30_bases,
            .q30_fraction = exact.q30_fraction.value(),
        };
    }

    pub fn ratios(self: *const Stats) Ratios {
        return .{
            .mean_length = .{ .numerator = self.bases, .denominator = self.reads },
            .gc_fraction = .{
                .numerator = @as(u128, self.g) + self.c,
                .denominator = @as(u128, self.a) + self.c + self.g + self.t,
            },
            .mean_quality = .{ .numerator = self.quality_sum, .denominator = self.bases },
            .q20_fraction = .{ .numerator = self.q20_bases, .denominator = self.bases },
            .q30_fraction = .{ .numerator = self.q30_bases, .denominator = self.bases },
        };
    }
};

const RecordTotals = struct {
    length: u64,
    a: u64,
    c: u64,
    g: u64,
    t: u64,
    n: u64,
    other_bases: u64,
    quality_sum: u64,
    q20_bases: u64,
    q30_bases: u64,
};

const Q20_SCORE = 20;
const Q30_SCORE = 30;

const PayloadTotals = extern struct {
    base_counts: [6]usize = .{0} ** 6,
    quality_sum: u64 = 0,
    q20_bases: usize = 0,
    q30_bases: usize = 0,
};

fn recordTotals(
    record: fastq.Record,
    quality_error: *?QualityError,
) StatsError!RecordTotals {
    if (record.sequence.len != record.quality.len) return error.S005LengthMismatch;

    const payload_totals = if (maxQualitySum(record.quality.len) != null)
        try payloadTotals(u64, record.sequence, record.quality, quality_error)
    else
        try payloadTotals(u128, record.sequence, record.quality, quality_error);
    return totalsFromPayload(payload_totals, record.sequence.len);
}

fn totalsFromPayload(payload_totals: PayloadTotals, length: usize) StatsError!RecordTotals {
    return .{
        .length = std.math.cast(u64, length) orelse return error.Overflow,
        .a = std.math.cast(u64, payload_totals.base_counts[0]) orelse return error.Overflow,
        .c = std.math.cast(u64, payload_totals.base_counts[1]) orelse return error.Overflow,
        .g = std.math.cast(u64, payload_totals.base_counts[2]) orelse return error.Overflow,
        .t = std.math.cast(u64, payload_totals.base_counts[3]) orelse return error.Overflow,
        .n = std.math.cast(u64, payload_totals.base_counts[4]) orelse return error.Overflow,
        .other_bases = std.math.cast(u64, payload_totals.base_counts[5]) orelse
            return error.Overflow,
        .quality_sum = payload_totals.quality_sum,
        .q20_bases = std.math.cast(u64, payload_totals.q20_bases) orelse return error.Overflow,
        .q30_bases = std.math.cast(u64, payload_totals.q30_bases) orelse return error.Overflow,
    };
}

pub const StreamedTotals = struct {
    totals: PayloadTotals = .{},
    sequence_len: usize = 0,
    quality_error: ?QualityError = null,
    overflow: bool = false,

    pub fn consume(
        self: *StreamedTotals,
        line_kind: fastq.ExpectedLine,
        bytes: []const u8,
        start_index: usize,
    ) fastq.ReaderError!void {
        switch (line_kind) {
            .sequence => self.consumeField(.sequence, bytes, start_index),
            .quality => self.consumeField(.quality, bytes, start_index),
            .header, .plus => {},
        }
    }

    fn consumeField(self: *StreamedTotals, comptime field: fastq.SemanticField, bytes: []const u8, start_index: usize) void {
        self.add(field == .quality, bytes, start_index);
    }

    fn add(self: *StreamedTotals, is_quality: bool, bytes: []const u8, start_index: usize) void {
        if (bytes.len == 0) return;
        if (!is_quality) {
            self.sequence_len += bytes.len;
            var no_quality_error: ?QualityError = null;
            const chunk = chunkTotals(.sequence, bytes, &no_quality_error) catch unreachable;
            for (&self.totals.base_counts, chunk.base_counts) |*total, count| total.* += count;
            return;
        }
        if (self.quality_error != null or self.overflow) return;
        var chunk_error: ?QualityError = null;
        const chunk = chunkTotals(.quality, bytes, &chunk_error) catch |err| switch (err) {
            error.S006InvalidQuality => {
                self.quality_error = .{ .byte_index = start_index + chunk_error.?.byte_index, .byte = chunk_error.?.byte };
                return;
            },
            else => {
                self.overflow = true;
                return;
            },
        };
        self.totals.quality_sum = std.math.add(u64, self.totals.quality_sum, chunk.quality_sum) catch {
            self.overflow = true;
            return;
        };
        self.totals.q20_bases += chunk.q20_bases;
        self.totals.q30_bases += chunk.q30_bases;
    }
};

fn maxQualitySum(length: usize) ?u64 {
    const length_u64 = std.math.cast(u64, length) orelse return null;
    return std.math.mul(u64, length_u64, fastq.PHRED_MAX_SCORE) catch null;
}

fn payloadTotals(
    comptime sum_type: type,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    if (sum_type == u64 and simd.useAvx2(quality.len)) return avx2Totals(.both, sequence, quality, quality_error);
    return payloadTotalsInline(sum_type, .both, sequence, quality, quality_error);
}

fn chunkTotals(comptime parts: Parts, bytes: []const u8, quality_error: *?QualityError) StatsError!PayloadTotals {
    const sequence = if (parts == .sequence) bytes else "";
    const quality = if (parts == .quality) bytes else "";
    if (simd.useAvx2(bytes.len)) return avx2Totals(parts, sequence, quality, quality_error);
    return payloadTotalsInline(u64, parts, sequence, quality, quality_error);
}

fn avx2Totals(parts: Parts, sequence: []const u8, quality: []const u8, quality_error: *?QualityError) StatsError!PayloadTotals {
    var totals: PayloadTotals = undefined;
    var error_index: usize = undefined;
    var error_byte: u8 = undefined;
    switch (simd.avx2Kernel("payload_totals", &payloadTotalsExport)(@intFromEnum(parts), sequence.ptr, sequence.len, quality.ptr, quality.len, &totals, &error_index, &error_byte)) {
        0 => return totals,
        1 => {
            quality_error.* = .{ .byte_index = error_index, .byte = error_byte };
            return error.S006InvalidQuality;
        },
        else => return error.Overflow,
    }
}

comptime {
    simd.exportKernel("payload_totals", &payloadTotalsExport);
}

fn payloadTotalsExport(
    parts: u8,
    sequence: [*]const u8,
    sequence_len: usize,
    quality: [*]const u8,
    quality_len: usize,
    totals: *PayloadTotals,
    error_index: *usize,
    error_byte: *u8,
) callconv(.c) u8 {
    var quality_error: ?QualityError = null;
    totals.* = switch (@as(Parts, @enumFromInt(parts))) {
        inline else => |known| payloadTotalsInline(u64, known, sequence[0..sequence_len], quality[0..quality_len], &quality_error),
    } catch |err| switch (err) {
        error.S006InvalidQuality => {
            error_index.* = quality_error.?.byte_index;
            error_byte.* = quality_error.?.byte;
            return 1;
        },
        else => return 2,
    };
    return 0;
}

const Parts = enum { both, sequence, quality };

fn payloadTotalsInline(
    comptime sum_type: type,
    comptime parts: Parts,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    const len = if (parts == .sequence) sequence.len else quality.len;
    if (sum_type == u64 and len >= simd.INLINE_LANES) {
        return payloadTotalsVector(simd.INLINE_LANES, parts, sequence, quality, quality_error);
    }
    return payloadTotalsScalar(sum_type, parts, sequence, quality, quality_error);
}

fn payloadTotalsScalar(
    comptime sum_type: type,
    comptime parts: Parts,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    var base_counts = [_]usize{0} ** 6;
    var quality_sum: sum_type = 0;
    var q20_bases: usize = 0;
    var q30_bases: usize = 0;
    if (parts != .sequence) {
        for (quality, 0..) |quality_byte, byte_index| {
            const score = fastq.decodePhred33(quality_byte) catch {
                quality_error.* = .{ .byte_index = byte_index, .byte = quality_byte };
                return error.S006InvalidQuality;
            };
            quality_sum += score;
            if (score >= Q20_SCORE) q20_bases += 1;
            if (score >= Q30_SCORE) q30_bases += 1;
        }
    }
    if (parts != .quality) {
        for (sequence) |base| {
            const index: usize = switch (base) {
                'A', 'a' => 0,
                'C', 'c' => 1,
                'G', 'g' => 2,
                'T', 't' => 3,
                'N', 'n' => 4,
                else => 5,
            };
            base_counts[index] += 1;
        }
    }
    return .{
        .base_counts = base_counts,
        .quality_sum = std.math.cast(u64, quality_sum) orelse return error.Overflow,
        .q20_bases = q20_bases,
        .q30_bases = q30_bases,
    };
}

fn payloadTotalsVector(
    comptime vector_len: comptime_int,
    comptime parts: Parts,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    comptime std.debug.assert(vector_len >= 2 and vector_len % 2 == 0 and
        vector_len <= std.math.maxInt(u16) / std.math.maxInt(u8));
    const half_len = vector_len / 2;
    const vectors_per_quality_block: usize = std.math.maxInt(u8) - 1;
    const vectors_per_base_block: usize = std.math.maxInt(u8) - 1;
    const Bytes = @Vector(vector_len, u8);
    const Mask = @Vector(vector_len, bool);
    const Lanes = @Vector(half_len, u16);
    const ReductionLanes = @Vector(half_len, u32);
    const phred_offset: Bytes = @splat(fastq.PHRED_OFFSET);
    const maximum_score: Bytes = @splat(fastq.PHRED_MAX_SCORE);
    const q20_score: Bytes = @splat(Q20_SCORE);
    const q30_score: Bytes = @splat(Q30_SCORE);
    const case_mask: Bytes = @splat(0xdf);
    const a_base: Bytes = @splat('A');
    const c_base: Bytes = @splat('C');
    const g_base: Bytes = @splat('G');
    const t_base: Bytes = @splat('T');
    const n_base: Bytes = @splat('N');
    const ones: Bytes = @splat(1);
    const zeros: Bytes = @splat(0);
    const lane_indices = std.simd.iota(u8, vector_len);
    const tail_len = (if (parts == .sequence) sequence.len else quality.len) % vector_len;

    var quality_sum: u64 = 0;
    var q20_bases: usize = 0;
    var q30_bases: usize = 0;
    var byte_index: usize = 0;
    while (parts != .sequence and quality.len - byte_index >= vector_len) {
        const remaining_vectors = (quality.len - byte_index) / vector_len;
        const block_vectors = @min(remaining_vectors, vectors_per_quality_block);
        const block_start = byte_index;
        var sum_lanes: Lanes = @splat(0);
        var q20_lanes: Bytes = @splat(0);
        var q30_lanes: Bytes = @splat(0);
        var invalid_lanes: Mask = @splat(false);

        var vectors_left = block_vectors;
        while (vectors_left >= 2) : (vectors_left -= 2) {
            const first_encoded: Bytes = quality[byte_index..][0..vector_len].*;
            const second_encoded: Bytes = quality[byte_index + vector_len ..][0..vector_len].*;
            const first_scores = first_encoded -% phred_offset;
            const second_scores = second_encoded -% phred_offset;
            invalid_lanes |= (first_scores > maximum_score) | (second_scores > maximum_score);

            const score_pairs = first_scores +% second_scores;
            const scores = std.simd.deinterlace(2, score_pairs);
            sum_lanes += @as(Lanes, @intCast(scores[0])) + @as(Lanes, @intCast(scores[1]));

            q20_lanes += @select(u8, first_scores >= q20_score, ones, zeros) +
                @select(u8, second_scores >= q20_score, ones, zeros);
            q30_lanes += @select(u8, first_scores >= q30_score, ones, zeros) +
                @select(u8, second_scores >= q30_score, ones, zeros);
            byte_index += 2 * vector_len;
        }
        if (vectors_left == 1) {
            const encoded: Bytes = quality[byte_index..][0..vector_len].*;
            const score_bytes = encoded -% phred_offset;
            invalid_lanes |= score_bytes > maximum_score;

            const scores = std.simd.deinterlace(2, score_bytes);
            sum_lanes += @as(Lanes, @intCast(scores[0])) + @as(Lanes, @intCast(scores[1]));
            q20_lanes += @select(u8, score_bytes >= q20_score, ones, zeros);
            q30_lanes += @select(u8, score_bytes >= q30_score, ones, zeros);
            byte_index += vector_len;
        }
        if (block_vectors == remaining_vectors and tail_len != 0) {
            const tail_start = quality.len - vector_len;
            const overlap: u8 = @intCast(byte_index - tail_start);
            const active_lanes = lane_indices >= @as(Bytes, @splat(overlap));
            const encoded: Bytes = quality[tail_start..][0..vector_len].*;
            const score_bytes = encoded -% phred_offset;
            invalid_lanes |= active_lanes & (score_bytes > maximum_score);

            const active_scores = @select(u8, active_lanes, score_bytes, zeros);
            const scores = std.simd.deinterlace(2, active_scores);
            sum_lanes += @as(Lanes, @intCast(scores[0])) + @as(Lanes, @intCast(scores[1]));
            q20_lanes += @select(
                u8,
                active_lanes & (score_bytes >= q20_score),
                ones,
                zeros,
            );
            q30_lanes += @select(
                u8,
                active_lanes & (score_bytes >= q30_score),
                ones,
                zeros,
            );
            byte_index = quality.len;
        }

        if (@reduce(.Or, invalid_lanes)) {
            return rejectInvalidQualityBlock(
                quality[block_start..byte_index],
                block_start,
                quality_error,
            );
        }

        quality_sum += @reduce(.Add, @as(ReductionLanes, @intCast(sum_lanes)));
        q20_bases += sumBaseLanes(vector_len, q20_lanes);
        q30_bases += sumBaseLanes(vector_len, q30_lanes);
    }

    var base_counts = [_]usize{0} ** 6;
    var base_index: usize = 0;
    while (parts != .quality and sequence.len - base_index >= vector_len) {
        const remaining_vectors = (sequence.len - base_index) / vector_len;
        const block_vectors = @min(remaining_vectors, vectors_per_base_block);
        var a_lanes: Bytes = @splat(0);
        var c_lanes: Bytes = @splat(0);
        var g_lanes: Bytes = @splat(0);
        var t_lanes: Bytes = @splat(0);
        var n_lanes: Bytes = @splat(0);

        var vectors_left = block_vectors;
        while (vectors_left >= 2) : (vectors_left -= 2) {
            const first_bases: Bytes = sequence[base_index..][0..vector_len].*;
            const second_bases: Bytes = sequence[base_index + vector_len ..][0..vector_len].*;
            const first_normalized = first_bases & case_mask;
            const second_normalized = second_bases & case_mask;
            a_lanes += @select(u8, first_normalized == a_base, ones, zeros) +
                @select(u8, second_normalized == a_base, ones, zeros);
            c_lanes += @select(u8, first_normalized == c_base, ones, zeros) +
                @select(u8, second_normalized == c_base, ones, zeros);
            g_lanes += @select(u8, first_normalized == g_base, ones, zeros) +
                @select(u8, second_normalized == g_base, ones, zeros);
            t_lanes += @select(u8, first_normalized == t_base, ones, zeros) +
                @select(u8, second_normalized == t_base, ones, zeros);
            n_lanes += @select(u8, first_normalized == n_base, ones, zeros) +
                @select(u8, second_normalized == n_base, ones, zeros);
            base_index += 2 * vector_len;
        }
        if (vectors_left == 1) {
            const bases: Bytes = sequence[base_index..][0..vector_len].*;
            const normalized = bases & case_mask;
            a_lanes += @select(u8, normalized == a_base, ones, zeros);
            c_lanes += @select(u8, normalized == c_base, ones, zeros);
            g_lanes += @select(u8, normalized == g_base, ones, zeros);
            t_lanes += @select(u8, normalized == t_base, ones, zeros);
            n_lanes += @select(u8, normalized == n_base, ones, zeros);
            base_index += vector_len;
        }
        if (block_vectors == remaining_vectors and tail_len != 0) {
            const tail_start = sequence.len - vector_len;
            const overlap: u8 = @intCast(base_index - tail_start);
            const active_lanes = lane_indices >= @as(Bytes, @splat(overlap));
            const bases: Bytes = sequence[tail_start..][0..vector_len].*;
            const normalized = bases & case_mask;
            a_lanes += @select(u8, active_lanes & (normalized == a_base), ones, zeros);
            c_lanes += @select(u8, active_lanes & (normalized == c_base), ones, zeros);
            g_lanes += @select(u8, active_lanes & (normalized == g_base), ones, zeros);
            t_lanes += @select(u8, active_lanes & (normalized == t_base), ones, zeros);
            n_lanes += @select(u8, active_lanes & (normalized == n_base), ones, zeros);
            base_index = sequence.len;
        }

        base_counts[0] += sumBaseLanes(vector_len, a_lanes);
        base_counts[1] += sumBaseLanes(vector_len, c_lanes);
        base_counts[2] += sumBaseLanes(vector_len, g_lanes);
        base_counts[3] += sumBaseLanes(vector_len, t_lanes);
        base_counts[4] += sumBaseLanes(vector_len, n_lanes);
    }
    base_counts[5] = sequence.len -
        (base_counts[0] + base_counts[1] + base_counts[2] + base_counts[3] + base_counts[4]);
    return .{
        .base_counts = base_counts,
        .quality_sum = quality_sum,
        .q20_bases = q20_bases,
        .q30_bases = q30_bases,
    };
}

fn rejectInvalidQualityBlock(
    quality: []const u8,
    base_index: usize,
    quality_error: *?QualityError,
) error{S006InvalidQuality} {
    for (quality, base_index..) |quality_byte, byte_index| {
        _ = fastq.decodePhred33(quality_byte) catch {
            quality_error.* = .{ .byte_index = byte_index, .byte = quality_byte };
            return error.S006InvalidQuality;
        };
    }
    unreachable;
}

fn sumBaseLanes(
    comptime vector_len: comptime_int,
    lanes: @Vector(vector_len, u8),
) usize {
    const WideLanes = @Vector(vector_len, u16);
    return @reduce(.Add, @as(WideLanes, @intCast(lanes)));
}

fn add(left: u64, right: u64) error{Overflow}!u64 {
    return std.math.add(u64, left, right);
}

test "[unit] - [statistics]: quality accumulator bound is the largest score sum that fits 64 bits" {
    const common_limit = std.math.maxInt(u64) / fastq.PHRED_MAX_SCORE;
    const boundary: usize = @intCast(common_limit);
    try std.testing.expectEqual(common_limit * fastq.PHRED_MAX_SCORE, maxQualitySum(boundary).?);
    try std.testing.expect(maxQualitySum(boundary + 1) == null);
}

fn blockBoundaryLengths(comptime vector_len: comptime_int) [14]usize {
    const quality_block_len = vector_len * (std.math.maxInt(u8) - 1);
    const base_block_len = vector_len * std.math.maxInt(u8);
    return .{
        0,
        1,
        vector_len - 1,
        vector_len,
        vector_len + 1,
        2 * vector_len + 1,
        base_block_len - 1,
        base_block_len,
        base_block_len + 1,
        quality_block_len - vector_len,
        quality_block_len - 1,
        quality_block_len,
        quality_block_len + 1,
        quality_block_len + vector_len,
    };
}

fn fillPayload(sequence: []u8, quality: []u8) void {
    const qualities = [_]u8{ 33, 52, 53, 62, 63, 126 };
    for (sequence, quality, 0..) |*base, *quality_byte, index| {
        base.* = @intCast(index % 256);
        quality_byte.* = qualities[index % qualities.len];
    }
}

test "[property] - [statistics]: vector payload totals match the scalar path at block boundaries" {
    inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |vector_len| {
        const lengths = blockBoundaryLengths(vector_len);
        const max_len = std.mem.max(usize, &lengths);
        const sequence = try std.testing.allocator.alloc(u8, max_len);
        defer std.testing.allocator.free(sequence);
        const quality = try std.testing.allocator.alloc(u8, max_len);
        defer std.testing.allocator.free(quality);
        fillPayload(sequence, quality);

        for (lengths) |length| {
            const sequence_part = sequence[0..length];
            const quality_part = quality[0..length];
            var quality_error: ?QualityError = null;
            try std.testing.expectEqualDeep(
                try payloadTotalsScalar(u64, .both, sequence_part, quality_part, &quality_error),
                try payloadTotals(u64, sequence_part, quality_part, &quality_error),
            );
            if (length < vector_len) continue;
            inline for (.{ Parts.both, Parts.sequence, Parts.quality }) |parts| {
                const sequence_arg = if (parts == .quality) sequence_part[0..0] else sequence_part;
                const quality_arg = if (parts == .sequence) quality_part[0..0] else quality_part;
                try std.testing.expectEqualDeep(
                    try payloadTotalsScalar(u64, parts, sequence_arg, quality_arg, &quality_error),
                    try payloadTotalsVector(vector_len, parts, sequence_arg, quality_arg, &quality_error),
                );
            }
        }
    }
}

fn expectFirstInvalidQuality(
    comptime vector_len: comptime_int,
    sequence: []const u8,
    quality: []const u8,
    index: usize,
) !void {
    const expected: ?QualityError = .{ .byte_index = index, .byte = quality[index] };
    var scalar_error: ?QualityError = null;
    var vector_error: ?QualityError = null;
    try std.testing.expectError(
        error.S006InvalidQuality,
        payloadTotalsScalar(u64, .both, sequence, quality, &scalar_error),
    );
    try std.testing.expectError(
        error.S006InvalidQuality,
        payloadTotalsVector(vector_len, .both, sequence, quality, &vector_error),
    );
    try std.testing.expectEqualDeep(expected, scalar_error);
    try std.testing.expectEqualDeep(expected, vector_error);
}

test "[failure] - [statistics]: vector payload totals report the first invalid quality byte" {
    inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |vector_len| {
        const quality_block_len = vector_len * (std.math.maxInt(u8) - 1);
        const block_len = quality_block_len + vector_len;
        const sequence = try std.testing.allocator.alloc(u8, block_len);
        defer std.testing.allocator.free(sequence);
        const quality = try std.testing.allocator.alloc(u8, block_len);
        defer std.testing.allocator.free(quality);

        const short_len = vector_len + 1;
        for ([_]usize{ vector_len - 1, vector_len }) |invalid_index| {
            fillPayload(sequence, quality);
            quality[invalid_index] = 127;
            try expectFirstInvalidQuality(vector_len, sequence[0..short_len], quality[0..short_len], invalid_index);
        }

        const invalid_block = [_]struct { index: usize, byte: u8 }{
            .{ .index = 0, .byte = 32 },
            .{ .index = quality_block_len - 1, .byte = 127 },
            .{ .index = quality_block_len, .byte = 32 },
            .{ .index = block_len - 1, .byte = 127 },
        };
        for (invalid_block) |invalid| {
            fillPayload(sequence, quality);
            quality[invalid.index] = invalid.byte;
            try expectFirstInvalidQuality(vector_len, sequence, quality, invalid.index);
        }

        fillPayload(sequence, quality);
        quality[block_len - 1] = 127;
        quality[3 * vector_len + 1] = 32;
        try expectFirstInvalidQuality(vector_len, sequence, quality, 3 * vector_len + 1);
    }
}

test "[edge] - [statistics]: vector payload totals at the maximum quality block" {
    inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |vector_len| {
        const block_len = vector_len * (std.math.maxInt(u8) - 1);
        const sequence = try std.testing.allocator.alloc(u8, block_len);
        defer std.testing.allocator.free(sequence);
        const quality = try std.testing.allocator.alloc(u8, block_len);
        defer std.testing.allocator.free(quality);
        @memset(sequence, 'A');
        @memset(quality, fastq.PHRED_MAX_BYTE);

        var quality_error: ?QualityError = null;
        const totals = try payloadTotalsVector(vector_len, .both, sequence, quality, &quality_error);
        try std.testing.expectEqual(@as(usize, block_len), totals.base_counts[0]);
        try std.testing.expectEqual(@as(u64, block_len) * fastq.PHRED_MAX_SCORE, totals.quality_sum);
        try std.testing.expectEqual(@as(usize, block_len), totals.q30_bases);
        try std.testing.expectEqual(@as(?QualityError, null), quality_error);
    }
}
