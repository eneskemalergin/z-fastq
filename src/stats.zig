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

/// Decodes one Phred+33 byte and rejects values outside ASCII 33 through 126.
pub fn decodePhred33(quality_byte: u8) error{InvalidQuality}!u8 {
    if (quality_byte < 33 or quality_byte > 126) return error.InvalidQuality;
    return quality_byte - 33;
}

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

fn recordTotals(
    record: fastq.Record,
    quality_error: *?QualityError,
) StatsError!RecordTotals {
    if (record.sequence.len != record.quality.len) return error.S005LengthMismatch;

    var base_counts = [_]usize{0} ** 6;
    for (record.sequence) |base| {
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

    var quality_sum: u128 = 0;
    var q20_bases: usize = 0;
    var q30_bases: usize = 0;
    for (record.quality, 0..) |quality_byte, byte_index| {
        const score = decodePhred33(quality_byte) catch {
            quality_error.* = .{ .byte_index = byte_index, .byte = quality_byte };
            return error.S006InvalidQuality;
        };
        quality_sum += score;
        if (score >= 20) q20_bases += 1;
        if (score >= 30) q30_bases += 1;
    }

    return .{
        .length = std.math.cast(u64, record.sequence.len) orelse return error.Overflow,
        .a = std.math.cast(u64, base_counts[0]) orelse return error.Overflow,
        .c = std.math.cast(u64, base_counts[1]) orelse return error.Overflow,
        .g = std.math.cast(u64, base_counts[2]) orelse return error.Overflow,
        .t = std.math.cast(u64, base_counts[3]) orelse return error.Overflow,
        .n = std.math.cast(u64, base_counts[4]) orelse return error.Overflow,
        .other_bases = std.math.cast(u64, base_counts[5]) orelse return error.Overflow,
        .quality_sum = std.math.cast(u64, quality_sum) orelse return error.Overflow,
        .q20_bases = std.math.cast(u64, q20_bases) orelse return error.Overflow,
        .q30_bases = std.math.cast(u64, q30_bases) orelse return error.Overflow,
    };
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
