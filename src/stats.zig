//! Checked aggregate length, base-composition, and Phred+33 statistics.

const std = @import("std");
const fastq = @import("fastq.zig");

pub const StatsError = error{
    S005LengthMismatch,
    S006InvalidQuality,
    Overflow,
};

pub const QualityError = struct {
    byte_index: usize,
    byte: u8,
};

/// Materialized counters and ratios; optional values are null for an empty population.
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

/// Allocation-free accumulator whose totals change only after a complete valid update.
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

    /// Validates and applies one record without partially changing any counter.
    pub fn addRecord(self: *Stats, record: fastq.Record) StatsError!void {
        self.last_quality_error = null;
        var quality_error: ?QualityError = null;
        const delta = recordTotals(record, &quality_error) catch |err| {
            self.last_quality_error = quality_error;
            return err;
        };

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

    /// Returns and clears details retained after S006.
    pub fn takeLastQualityError(self: *Stats) ?QualityError {
        const quality_error = self.last_quality_error;
        self.last_quality_error = null;
        return quality_error;
    }

    /// Materializes ratios without changing the accumulator.
    pub fn result(self: *const Stats) StatsResult {
        const gc_bases = @as(u128, self.a) + self.c + self.g + self.t;
        return .{
            .reads = self.reads,
            .bases = self.bases,
            .min_length = if (self.reads == 0) null else self.min_length,
            .max_length = if (self.reads == 0) null else self.max_length,
            .mean_length = ratio(self.bases, self.reads),
            .a = self.a,
            .c = self.c,
            .g = self.g,
            .t = self.t,
            .n = self.n,
            .other_bases = self.other_bases,
            .gc_fraction = ratioWide(@as(u128, self.g) + self.c, gc_bases),
            .quality_sum = self.quality_sum,
            .mean_quality = ratio(self.quality_sum, self.bases),
            .q20_bases = self.q20_bases,
            .q20_fraction = ratio(self.q20_bases, self.bases),
            .q30_bases = self.q30_bases,
            .q30_fraction = ratio(self.q30_bases, self.bases),
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

const PayloadTotals = struct {
    base_counts: [6]usize,
    quality_sum: u64,
    q20_bases: usize,
    q30_bases: usize,
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

    return .{
        .length = std.math.cast(u64, record.sequence.len) orelse return error.Overflow,
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

fn maxQualitySum(length: usize) ?u64 {
    const length_u64 = std.math.cast(u64, length) orelse return null;
    return std.math.mul(u64, length_u64, 93) catch null;
}

fn payloadTotals(
    comptime sum_type: type,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    if (sum_type == u64) {
        if (std.simd.suggestVectorLength(u8)) |vector_len| {
            if (vector_len >= 2 and
                vector_len % 2 == 0 and
                vector_len <= std.math.maxInt(u16) / std.math.maxInt(u8) and
                quality.len >= vector_len)
            {
                return payloadTotalsVector(vector_len, sequence, quality, quality_error);
            }
        }
    }
    return payloadTotalsScalar(sum_type, sequence, quality, quality_error);
}

fn payloadTotalsScalar(
    comptime sum_type: type,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    var base_counts = [_]usize{0} ** 6;
    var quality_sum: sum_type = 0;
    var q20_bases: usize = 0;
    var q30_bases: usize = 0;
    for (sequence, quality, 0..) |base, quality_byte, byte_index| {
        const score = fastq.decodePhred33(quality_byte) catch {
            quality_error.* = .{ .byte_index = byte_index, .byte = quality_byte };
            return error.S006InvalidQuality;
        };
        quality_sum += score;
        if (score >= 20) q20_bases += 1;
        if (score >= 30) q30_bases += 1;

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
    return .{
        .base_counts = base_counts,
        .quality_sum = std.math.cast(u64, quality_sum) orelse return error.Overflow,
        .q20_bases = q20_bases,
        .q30_bases = q30_bases,
    };
}

fn payloadTotalsVector(
    comptime vector_len: comptime_int,
    sequence: []const u8,
    quality: []const u8,
    quality_error: *?QualityError,
) StatsError!PayloadTotals {
    const half_len = vector_len / 2;
    const vectors_per_block = std.math.maxInt(u16) / (2 * 93);
    const Bytes = @Vector(vector_len, u8);
    const Lanes = @Vector(half_len, u16);
    const ReductionLanes = @Vector(half_len, u32);
    const minimum: Bytes = @splat(33);
    const maximum_score: Bytes = @splat(93);
    const q20_minimum: Bytes = @splat(53);
    const q30_minimum: Bytes = @splat(63);
    const case_mask: Bytes = @splat(0xdf);
    const a_base: Bytes = @splat('A');
    const c_base: Bytes = @splat('C');
    const g_base: Bytes = @splat('G');
    const t_base: Bytes = @splat('T');
    const n_base: Bytes = @splat('N');
    const ones: Bytes = @splat(1);
    const zeros: Bytes = @splat(0);

    var quality_sum: u64 = 0;
    var q20_bases: usize = 0;
    var q30_bases: usize = 0;
    var byte_index: usize = 0;
    while (quality.len - byte_index >= vector_len) {
        const remaining_vectors = (quality.len - byte_index) / vector_len;
        const block_vectors = @min(remaining_vectors, vectors_per_block);
        var sum_lanes: Lanes = @splat(0);
        var q20_lanes: Lanes = @splat(0);
        var q30_lanes: Lanes = @splat(0);

        var vector_index: usize = 0;
        while (vector_index < block_vectors) : (vector_index += 1) {
            const encoded: Bytes = quality[byte_index..][0..vector_len].*;
            const score_bytes = encoded -% minimum;
            const invalid = score_bytes > maximum_score;
            if (@reduce(.Or, invalid)) {
                const lane = std.simd.firstTrue(invalid).?;
                quality_error.* = .{
                    .byte_index = byte_index + lane,
                    .byte = quality[byte_index + lane],
                };
                return error.S006InvalidQuality;
            }

            const scores = std.simd.deinterlace(2, score_bytes);
            sum_lanes += @as(Lanes, @intCast(scores[0])) + @as(Lanes, @intCast(scores[1]));

            const q20 = std.simd.deinterlace(2, @select(u8, encoded >= q20_minimum, ones, zeros));
            q20_lanes += @as(Lanes, @intCast(q20[0])) + @as(Lanes, @intCast(q20[1]));

            const q30 = std.simd.deinterlace(2, @select(u8, encoded >= q30_minimum, ones, zeros));
            q30_lanes += @as(Lanes, @intCast(q30[0])) + @as(Lanes, @intCast(q30[1]));
            byte_index += vector_len;
        }

        quality_sum += @reduce(.Add, @as(ReductionLanes, @intCast(sum_lanes)));
        q20_bases += @reduce(.Add, @as(ReductionLanes, @intCast(q20_lanes)));
        q30_bases += @reduce(.Add, @as(ReductionLanes, @intCast(q30_lanes)));
    }

    var base_counts = [_]usize{0} ** 6;
    var base_index: usize = 0;
    while (byte_index - base_index >= vector_len) {
        const remaining_vectors = (byte_index - base_index) / vector_len;
        const block_vectors = @min(remaining_vectors, std.math.maxInt(u8));
        var a_lanes: Bytes = @splat(0);
        var c_lanes: Bytes = @splat(0);
        var g_lanes: Bytes = @splat(0);
        var t_lanes: Bytes = @splat(0);
        var n_lanes: Bytes = @splat(0);

        var vector_index: usize = 0;
        while (vector_index < block_vectors) : (vector_index += 1) {
            const bases: Bytes = sequence[base_index..][0..vector_len].*;
            const normalized = bases & case_mask;
            a_lanes += @select(u8, normalized == a_base, ones, zeros);
            c_lanes += @select(u8, normalized == c_base, ones, zeros);
            g_lanes += @select(u8, normalized == g_base, ones, zeros);
            t_lanes += @select(u8, normalized == t_base, ones, zeros);
            n_lanes += @select(u8, normalized == n_base, ones, zeros);
            base_index += vector_len;
        }

        base_counts[0] += sumBaseLanes(vector_len, a_lanes);
        base_counts[1] += sumBaseLanes(vector_len, c_lanes);
        base_counts[2] += sumBaseLanes(vector_len, g_lanes);
        base_counts[3] += sumBaseLanes(vector_len, t_lanes);
        base_counts[4] += sumBaseLanes(vector_len, n_lanes);
    }
    base_counts[5] = byte_index -
        (base_counts[0] + base_counts[1] + base_counts[2] + base_counts[3] + base_counts[4]);

    for (
        sequence[byte_index..],
        quality[byte_index..],
        byte_index..,
    ) |base, quality_byte, tail_index| {
        const score = fastq.decodePhred33(quality_byte) catch {
            quality_error.* = .{ .byte_index = tail_index, .byte = quality_byte };
            return error.S006InvalidQuality;
        };
        quality_sum += score;
        if (score >= 20) q20_bases += 1;
        if (score >= 30) q30_bases += 1;

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
    return .{
        .base_counts = base_counts,
        .quality_sum = quality_sum,
        .q20_bases = q20_bases,
        .q30_bases = q30_bases,
    };
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

fn ratio(numerator: u64, denominator: u64) ?f64 {
    if (denominator == 0) return null;
    return @as(f64, @floatFromInt(numerator)) / @as(f64, @floatFromInt(denominator));
}

fn ratioWide(numerator: u128, denominator: u128) ?f64 {
    if (denominator == 0) return null;
    return @as(f64, @floatFromInt(numerator)) / @as(f64, @floatFromInt(denominator));
}

test "[unit] - [statistics]: quality accumulator width uses the checked bound" {
    if (comptime @bitSizeOf(usize) < 64) {
        try std.testing.expect(maxQualitySum(std.math.maxInt(usize)) != null);
    } else {
        const common_limit = std.math.maxInt(u64) / 93;
        const boundary: usize = @intCast(common_limit);
        try std.testing.expectEqual(common_limit * 93, maxQualitySum(boundary).?);
        try std.testing.expect(maxQualitySum(boundary + 1) == null);
    }
}

test "[property] - [statistics]: wrapped quality scores preserve the valid range" {
    for (0..std.math.maxInt(u8) + 1) |value| {
        const encoded: u8 = @intCast(value);
        const wrapped_score = encoded -% 33;
        if (fastq.decodePhred33(encoded)) |score| {
            try std.testing.expect(wrapped_score <= 93);
            try std.testing.expectEqual(score, wrapped_score);
        } else |err| {
            try std.testing.expectEqual(error.InvalidQuality, err);
            try std.testing.expect(wrapped_score > 93);
        }
    }
}

test "[property] - [statistics]: vector payload kernels match the scalar path" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse return;
    if (vector_len < 2 or
        vector_len % 2 != 0 or
        vector_len > std.math.maxInt(u16) / std.math.maxInt(u8)) return;

    const quality_block_len = vector_len * (std.math.maxInt(u16) / (2 * 93));
    const base_block_len = vector_len * std.math.maxInt(u8);
    const max_len = @max(quality_block_len, base_block_len) + vector_len + 1;
    const lengths = [_]usize{
        0,
        1,
        vector_len - 1,
        vector_len,
        vector_len + 1,
        2 * vector_len - 1,
        2 * vector_len,
        2 * vector_len + 1,
        base_block_len - 1,
        base_block_len,
        base_block_len + 1,
        quality_block_len - 1,
        quality_block_len,
        quality_block_len + 1,
    };
    const qualities = [_]u8{ 33, 52, 53, 62, 63, 126 };
    const sequence = try std.testing.allocator.alloc(u8, max_len);
    defer std.testing.allocator.free(sequence);
    const quality = try std.testing.allocator.alloc(u8, max_len);
    defer std.testing.allocator.free(quality);
    for (sequence, quality, 0..) |*base, *quality_byte, index| {
        base.* = @intCast(index % 256);
        quality_byte.* = qualities[index % qualities.len];
    }

    for (lengths) |length| {
        var scalar_error: ?QualityError = null;
        var vector_error: ?QualityError = null;
        const scalar = try payloadTotalsScalar(
            u64,
            sequence[0..length],
            quality[0..length],
            &scalar_error,
        );
        const vector = try payloadTotals(
            u64,
            sequence[0..length],
            quality[0..length],
            &vector_error,
        );
        try std.testing.expectEqualDeep(scalar, vector);
        try std.testing.expectEqualDeep(scalar_error, vector_error);
    }

    const invalid_len = 2 * vector_len + 3;
    const invalid_indices = [_]usize{ 0, vector_len - 1, vector_len, invalid_len - 1 };
    for (invalid_indices, 0..) |invalid_index, case_index| {
        const original = quality[invalid_index];
        quality[invalid_index] = if (case_index % 2 == 0) 32 else 127;
        var scalar_error: ?QualityError = null;
        var vector_error: ?QualityError = null;
        try std.testing.expectError(
            error.S006InvalidQuality,
            payloadTotalsScalar(
                u64,
                sequence[0..invalid_len],
                quality[0..invalid_len],
                &scalar_error,
            ),
        );
        try std.testing.expectError(
            error.S006InvalidQuality,
            payloadTotals(u64, sequence[0..invalid_len], quality[0..invalid_len], &vector_error),
        );
        try std.testing.expectEqualDeep(scalar_error, vector_error);
        quality[invalid_index] = original;
    }

    quality[1] = 127;
    quality[vector_len - 1] = 32;
    var scalar_error: ?QualityError = null;
    var vector_error: ?QualityError = null;
    try std.testing.expectError(
        error.S006InvalidQuality,
        payloadTotalsScalar(
            u64,
            sequence[0..invalid_len],
            quality[0..invalid_len],
            &scalar_error,
        ),
    );
    try std.testing.expectError(
        error.S006InvalidQuality,
        payloadTotals(
            u64,
            sequence[0..invalid_len],
            quality[0..invalid_len],
            &vector_error,
        ),
    );
    try std.testing.expectEqualDeep(scalar_error, vector_error);
}
