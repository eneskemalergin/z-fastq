//! Deterministic numeric parsing and record selection for the private sample command.

const std = @import("std");

pub const FractionError = error{InvalidFraction};
pub const CountError = error{ InvalidCount, Overflow };
pub const SeedError = error{ InvalidSeed, Overflow };
pub const ReservoirError = std.mem.Allocator.Error || error{Overflow};

pub const Fraction = union(enum) {
    none,
    all,
    probability: f64,

    pub fn parse(text: []const u8) FractionError!Fraction {
        if (std.mem.eql(u8, text, "0")) return .none;
        if (std.mem.eql(u8, text, "1")) return .all;
        if (text.len < 3 or text[1] != '.') return error.InvalidFraction;

        if (text[0] == '1') {
            for (text[2..]) |byte| {
                if (byte != '0') return error.InvalidFraction;
            }
            return .all;
        }
        if (text[0] != '0') return error.InvalidFraction;
        var nonzero = false;
        for (text[2..]) |byte| {
            if (byte < '0' or byte > '9') return error.InvalidFraction;
            nonzero = nonzero or byte != '0';
        }
        if (!nonzero) return .none;
        const probability = std.fmt.parseFloat(f64, text) catch
            return error.InvalidFraction;
        if (probability == 0.0) return .none;
        if (probability == 1.0) return .all;
        return .{ .probability = probability };
    }
};

pub fn parseSeed(text: []const u8) SeedError!u64 {
    if (text.len == 0) return error.InvalidSeed;
    for (text) |byte| {
        if (byte < '0' or byte > '9') return error.InvalidSeed;
    }
    return std.fmt.parseInt(u64, text, 10) catch |err| switch (err) {
        error.Overflow => error.Overflow,
        error.InvalidCharacter => error.InvalidSeed,
    };
}

pub fn parseCount(text: []const u8) CountError!u64 {
    return parseSeed(text) catch |err| switch (err) {
        error.InvalidSeed => error.InvalidCount,
        error.Overflow => error.Overflow,
    };
}

pub const Selector = union(enum) {
    none,
    all,
    probability: Probability,

    const Probability = struct {
        fraction: f64,
        generator: Mt19937_64,
    };

    pub fn init(fraction: Fraction, seed: u64) Selector {
        return switch (fraction) {
            .none => .none,
            .all => .all,
            .probability => |probability| .{ .probability = .{
                .fraction = probability,
                .generator = Mt19937_64.init(seed),
            } },
        };
    }

    pub fn selectRecord(self: *Selector) bool {
        return switch (self.*) {
            .none => false,
            .all => true,
            .probability => |*state| state.generator.nextUnitFloat() < state.fraction,
        };
    }
};

pub const ExactSelector = struct {
    target: u64,
    generator: Mt19937_64,
    indexes: std.ArrayList(u32) = .empty,
    index_high_words: std.ArrayList(u32) = .empty,
    upper_bit_count: u8 = 0,
    unordered: bool = false,
    record_count: u64 = 0,

    pub fn init(target: u64, seed: u64) ExactSelector {
        return .{
            .target = target,
            .generator = Mt19937_64.init(seed),
        };
    }

    pub fn deinit(self: *ExactSelector, allocator: std.mem.Allocator) void {
        self.indexes.deinit(allocator);
        self.index_high_words.deinit(allocator);
        self.* = undefined;
    }

    pub fn considerRecord(
        self: *ExactSelector,
        allocator: std.mem.Allocator,
        record_number: u64,
    ) ReservoirError!void {
        std.debug.assert(record_number != 0);
        std.debug.assert(record_number - 1 == self.record_count);
        if (self.indexes.capacity != 0 and self.upper_bit_count < 8 and
            record_number == @as(u64, 1) << @intCast(24 + self.upper_bit_count))
        {
            self.activateUpperBits(self.upper_bit_count + 1);
        }
        try self.considerPreparedRecord(allocator, record_number);
    }

    pub fn considerRecordsThrough(
        self: *ExactSelector,
        allocator: std.mem.Allocator,
        completed_records: u64,
    ) ReservoirError!void {
        std.debug.assert(completed_records >= self.record_count);
        if (self.target == 0) {
            self.record_count = completed_records;
            return;
        }
        while (self.record_count < completed_records) {
            if (self.indexes.capacity == 0) {
                const allocation_record = if (self.target == std.math.maxInt(u64))
                    completed_records
                else
                    @min(completed_records, self.target + 1);
                while (self.record_count < allocation_record) {
                    try self.considerPreparedRecord(allocator, self.record_count + 1);
                }
                continue;
            }

            if (self.upper_bit_count < 8) {
                const boundary = @as(u64, 1) << @intCast(24 + self.upper_bit_count);
                if (self.record_count + 1 == boundary) {
                    self.activateUpperBits(self.upper_bit_count + 1);
                    continue;
                }
                const range_end = @min(completed_records, boundary - 1);
                while (self.record_count < range_end) {
                    try self.considerPreparedRecord(allocator, self.record_count + 1);
                }
                continue;
            }

            while (self.record_count < completed_records) {
                try self.considerPreparedRecord(allocator, self.record_count + 1);
            }
        }
    }

    fn considerPreparedRecord(
        self: *ExactSelector,
        allocator: std.mem.Allocator,
        record_number: u64,
    ) ReservoirError!void {
        self.record_count = record_number;
        if (self.target == 0) return;

        const draw = self.generator.nextUnitFloat();
        if (record_number <= self.target) return;

        if (self.indexes.capacity == 0) {
            const target_len = std.math.cast(usize, self.target) orelse
                return error.Overflow;
            const storage_words = std.math.add(usize, target_len, 2) catch
                return error.Overflow;
            _ = std.math.mul(usize, storage_words, @sizeOf(u32)) catch
                return error.Overflow;
            try self.indexes.ensureTotalCapacityPrecise(allocator, storage_words);
            self.indexes.items.len = target_len;
            if (self.target > std.math.maxInt(u32)) {
                try self.index_high_words.ensureTotalCapacityPrecise(
                    allocator,
                    target_len,
                );
                self.index_high_words.expandToCapacity();
            }
            self.activateUpperBits(requiredUpperBits(record_number));
            for (0..self.indexes.items.len) |slot| {
                const initial_record = slot + 1;
                self.writeLowIndex(slot, @truncate(initial_record));
                if (self.index_high_words.items.len != 0) {
                    self.index_high_words.items[slot] =
                        @truncate(initial_record >> 32);
                }
            }
        }

        const candidate_float = draw * @as(f64, @floatFromInt(record_number));
        const candidate = if (candidate_float >= @as(f64, @floatFromInt(record_number)))
            record_number - 1
        else
            @as(u64, @intFromFloat(candidate_float));
        if (candidate >= self.target) return;

        const slot = std.math.cast(usize, candidate) orelse return error.Overflow;
        try self.storeIndex(allocator, slot, record_number);
        self.unordered = self.unordered or slot + 1 != self.indexes.items.len;
    }

    fn storeIndex(
        self: *ExactSelector,
        allocator: std.mem.Allocator,
        slot: usize,
        record_number: u64,
    ) std.mem.Allocator.Error!void {
        const high_word: u32 = @truncate(record_number >> 32);
        if (high_word != 0 and self.index_high_words.items.len == 0) {
            try self.index_high_words.ensureTotalCapacityPrecise(
                allocator,
                self.indexes.items.len,
            );
            self.index_high_words.expandToCapacity();
            @memset(self.index_high_words.items, 0);
        }
        self.writeLowIndex(slot, @truncate(record_number));
        if (self.index_high_words.items.len != 0) {
            self.index_high_words.items[slot] = high_word;
        }
    }

    fn writeLowIndex(self: *ExactSelector, slot: usize, value: u32) void {
        self.lowWords()[slot] = @truncate(value);
        self.middleBytes()[slot] = @truncate(value >> 16);
        writeUpperIndex(
            self.upperBytes(),
            self.indexes.items.len,
            self.upper_bit_count,
            slot,
            @truncate(value >> 24),
        );
    }

    fn activateUpperBits(self: *ExactSelector, needed: u8) void {
        std.debug.assert(needed <= 8);
        if (needed <= self.upper_bit_count) return;

        const plane_bytes = std.math.divCeil(usize, self.indexes.items.len, 8) catch
            unreachable;
        const upper_bytes = self.upperBytes();
        while (self.upper_bit_count < needed) {
            const start = @as(usize, self.upper_bit_count) * plane_bytes;
            @memset(upper_bytes[start .. start + plane_bytes], 0);
            self.upper_bit_count += 1;
        }
    }

    fn lowWords(self: *ExactSelector) []u16 {
        const words: [*]u16 = @ptrCast(self.indexes.items.ptr);
        return words[0..self.indexes.items.len];
    }

    fn middleBytes(self: *ExactSelector) []u8 {
        const bytes: [*]u8 = @ptrCast(self.indexes.items.ptr);
        const start = self.indexes.items.len * @sizeOf(u16);
        return bytes[start .. start + self.indexes.items.len];
    }

    fn upperBytes(self: *ExactSelector) []u8 {
        const bytes: [*]u8 = @ptrCast(self.indexes.items.ptr);
        const start = self.indexes.items.len * 3;
        const plane_bytes = std.math.divCeil(usize, self.indexes.items.len, 8) catch
            unreachable;
        std.debug.assert(
            start + plane_bytes * 8 <= self.indexes.capacity * @sizeOf(u32),
        );
        return bytes[start .. start + plane_bytes * 8];
    }

    pub fn finish(self: *ExactSelector) ExactSelection {
        if (self.target == 0 or self.record_count == 0) return .none;
        if (self.record_count <= self.target) return .all;
        if (self.unordered) {
            if (self.index_high_words.items.len == 0) {
                sortSplitIndexes(
                    self.lowWords(),
                    self.middleBytes(),
                    self.upperBytes(),
                    self.upper_bit_count,
                );
            } else {
                std.sort.pdqContext(0, self.indexes.items.len, IndexSortContext{
                    .low_words = self.lowWords(),
                    .middle_bytes = self.middleBytes(),
                    .upper_bytes = self.upperBytes(),
                    .upper_bit_count = self.upper_bit_count,
                    .high_words = self.index_high_words.items,
                });
            }
        }
        return .{ .indexes = .{
            .low_words = self.lowWords(),
            .middle_bytes = self.middleBytes(),
            .upper_bytes = self.upperBytes(),
            .upper_bit_count = self.upper_bit_count,
            .high_words = self.index_high_words.items,
        } };
    }

    const IndexSortContext = struct {
        low_words: []u16,
        middle_bytes: []u8,
        upper_bytes: []u8,
        upper_bit_count: u8,
        high_words: []u32,

        pub fn lessThan(self: IndexSortContext, left: usize, right: usize) bool {
            if (self.high_words[left] != self.high_words[right]) {
                return self.high_words[left] < self.high_words[right];
            }
            return readSplitIndex(
                self.low_words,
                self.middle_bytes,
                self.upper_bytes,
                self.upper_bit_count,
                left,
            ) < readSplitIndex(
                self.low_words,
                self.middle_bytes,
                self.upper_bytes,
                self.upper_bit_count,
                right,
            );
        }

        pub fn swap(self: IndexSortContext, left: usize, right: usize) void {
            std.mem.swap(u16, &self.low_words[left], &self.low_words[right]);
            std.mem.swap(u8, &self.middle_bytes[left], &self.middle_bytes[right]);
            swapUpperIndexes(
                self.upper_bytes,
                self.low_words.len,
                self.upper_bit_count,
                left,
                right,
            );
            std.mem.swap(u32, &self.high_words[left], &self.high_words[right]);
        }
    };
};

pub const ExactIndexes = struct {
    low_words: []const u16 = &.{},
    middle_bytes: []const u8 = &.{},
    upper_bytes: []const u8 = &.{},
    upper_bit_count: u8 = 0,
    high_words: []const u32 = &.{},

    pub const empty: ExactIndexes = .{};

    pub fn len(self: ExactIndexes) usize {
        return self.low_words.len;
    }

    pub fn at(self: ExactIndexes, index: usize) u64 {
        const low = readSplitIndex(
            self.low_words,
            self.middle_bytes,
            self.upper_bytes,
            self.upper_bit_count,
            index,
        );
        const high_word = if (self.high_words.len == 0)
            0
        else
            @as(u64, self.high_words[index]);
        return high_word << 32 | low;
    }
};

pub const ExactSelection = union(enum) {
    none,
    all,
    indexes: ExactIndexes,
};

fn requiredUpperBits(record_number: u64) u8 {
    if (record_number <= 0xff_ffff) return 0;
    if (record_number > std.math.maxInt(u32)) return 8;
    return @intCast(std.math.log2_int(u32, @truncate(record_number >> 24)) + 1);
}

fn readSplitIndex(
    low_words: []const u16,
    middle_bytes: []const u8,
    upper_bytes: []const u8,
    upper_bit_count: u8,
    index: usize,
) u32 {
    return @as(u32, readUpperIndex(
        upper_bytes,
        low_words.len,
        upper_bit_count,
        index,
    )) << 24 |
        @as(u32, middle_bytes[index]) << 16 |
        @as(u32, low_words[index]);
}

fn readUpperIndex(
    storage: []const u8,
    len: usize,
    bit_count: u8,
    index: usize,
) u8 {
    if (bit_count == 0) return 0;
    const plane_bytes = std.math.divCeil(usize, len, 8) catch unreachable;

    const mask = @as(u8, 1) << @intCast(index % 8);
    const byte_index = index / 8;
    if (bit_count == 1) {
        return @intFromBool(storage[byte_index] & mask != 0);
    }
    var value: u8 = 0;
    for (0..bit_count) |bit| {
        if (storage[bit * plane_bytes + byte_index] & mask != 0) {
            value |= @as(u8, 1) << @intCast(bit);
        }
    }
    return value;
}

fn writeUpperIndex(
    storage: []u8,
    len: usize,
    bit_count: u8,
    index: usize,
    value: u8,
) void {
    if (bit_count == 0) return;
    const plane_bytes = std.math.divCeil(usize, len, 8) catch unreachable;

    const mask = @as(u8, 1) << @intCast(index % 8);
    const byte_index = index / 8;
    if (bit_count == 1) {
        if (value == 0) {
            storage[byte_index] &= ~mask;
        } else {
            storage[byte_index] |= mask;
        }
        return;
    }
    for (0..bit_count) |bit| {
        const byte = &storage[bit * plane_bytes + byte_index];
        if (value & (@as(u8, 1) << @intCast(bit)) == 0) {
            byte.* &= ~mask;
        } else {
            byte.* |= mask;
        }
    }
}

fn swapUpperIndexes(
    storage: []u8,
    len: usize,
    bit_count: u8,
    left: usize,
    right: usize,
) void {
    if (bit_count == 0 or left == right) return;
    const left_value = readUpperIndex(storage, len, bit_count, left);
    const right_value = readUpperIndex(storage, len, bit_count, right);
    writeUpperIndex(storage, len, bit_count, left, right_value);
    writeUpperIndex(storage, len, bit_count, right, left_value);
}

fn upperIndexBit(
    storage: []const u8,
    len: usize,
    index: usize,
    bit: u8,
) bool {
    const plane_bytes = std.math.divCeil(usize, len, 8) catch unreachable;
    return storage[@as(usize, bit) * plane_bytes + index / 8] &
        (@as(u8, 1) << @intCast(index % 8)) != 0;
}

fn partitionUpperBit(
    low_words: []u16,
    middle_bytes: []u8,
    upper_bytes: []u8,
    upper_bit_count: u8,
    start: usize,
    end: usize,
    bit: u8,
) usize {
    var left = start;
    var right = end;
    while (left < right) {
        while (left < right and !upperIndexBit(upper_bytes, low_words.len, left, bit)) {
            left += 1;
        }
        while (left < right and upperIndexBit(upper_bytes, low_words.len, right - 1, bit)) {
            right -= 1;
        }
        if (left == right) break;

        const swap_right = right - 1;
        std.mem.swap(u16, &low_words[left], &low_words[swap_right]);
        std.mem.swap(u8, &middle_bytes[left], &middle_bytes[swap_right]);
        swapUpperIndexes(
            upper_bytes,
            low_words.len,
            upper_bit_count,
            left,
            swap_right,
        );
        left += 1;
        right = swap_right;
    }
    return left;
}

fn appendUpperBoundaries(
    low_words: []u16,
    middle_bytes: []u8,
    upper_bytes: []u8,
    upper_bit_count: u8,
    start: usize,
    end: usize,
    remaining_bits: u8,
    boundaries: *[257]usize,
    boundary_count: *usize,
) void {
    if (remaining_bits == 0) {
        boundaries[boundary_count.*] = end;
        boundary_count.* += 1;
        return;
    }
    const bit = remaining_bits - 1;
    const split = partitionUpperBit(
        low_words,
        middle_bytes,
        upper_bytes,
        upper_bit_count,
        start,
        end,
        bit,
    );
    appendUpperBoundaries(
        low_words,
        middle_bytes,
        upper_bytes,
        upper_bit_count,
        start,
        split,
        bit,
        boundaries,
        boundary_count,
    );
    appendUpperBoundaries(
        low_words,
        middle_bytes,
        upper_bytes,
        upper_bit_count,
        split,
        end,
        bit,
        boundaries,
        boundary_count,
    );
}

fn sortSplitIndexes(
    low_words: []u16,
    middle_bytes: []u8,
    upper_bytes: []u8,
    upper_bit_count: u8,
) void {
    var upper_boundaries: [257]usize = undefined;
    upper_boundaries[0] = 0;
    var upper_boundary_count: usize = 1;
    appendUpperBoundaries(
        low_words,
        middle_bytes,
        upper_bytes,
        upper_bit_count,
        0,
        low_words.len,
        upper_bit_count,
        &upper_boundaries,
        &upper_boundary_count,
    );

    var counts = [_]usize{0} ** 256;
    var next: [256]usize = undefined;
    var middle_boundaries: [257]usize = undefined;
    var present = [_]u64{0} ** 1024;
    for (0..upper_boundary_count - 1) |upper| {
        const upper_start = upper_boundaries[upper];
        const upper_end = upper_boundaries[upper + 1];
        if (upper_start == upper_end) continue;

        @memset(&counts, 0);
        for (middle_bytes[upper_start..upper_end]) |middle| counts[middle] += 1;
        middle_boundaries[0] = upper_start;
        for (counts, 0..) |count, bucket| {
            middle_boundaries[bucket + 1] = middle_boundaries[bucket] + count;
        }
        @memcpy(next[0..], middle_boundaries[0..256]);

        for (0..256) |bucket| {
            while (next[bucket] < middle_boundaries[bucket + 1]) {
                const slot = next[bucket];
                const owner = middle_bytes[slot];
                if (owner == bucket) {
                    next[bucket] += 1;
                } else {
                    std.mem.swap(u16, &low_words[slot], &low_words[next[owner]]);
                    std.mem.swap(u8, &middle_bytes[slot], &middle_bytes[next[owner]]);
                    next[owner] += 1;
                }
            }
        }

        for (0..256) |middle| {
            const start = middle_boundaries[middle];
            const end = middle_boundaries[middle + 1];
            if (start == end) continue;

            @memset(&present, 0);
            for (low_words[start..end]) |low| {
                present[low / 64] |= @as(u64, 1) << @intCast(low % 64);
            }
            var output = start;
            for (present, 0..) |word, word_index| {
                var remaining = word;
                while (remaining != 0) {
                    low_words[output] = @intCast(word_index * 64 + @ctz(remaining));
                    output += 1;
                    remaining &= remaining - 1;
                }
            }
            std.debug.assert(output == end);
        }
    }
}

fn exerciseWideExactIndexes(allocator: std.mem.Allocator) !void {
    var selector = ExactSelector.init(2, 11);
    defer selector.deinit(allocator);
    for (1..4) |record_number| {
        try selector.considerRecord(allocator, record_number);
    }
    selector.activateUpperBits(8);
    const wide_index = @as(u64, 1) << 32 | 0xab_00_0007;
    try selector.storeIndex(allocator, 0, wide_index);
    try selector.storeIndex(allocator, 1, 4);
    selector.unordered = true;

    const selection = selector.finish();
    try std.testing.expect(selection == .indexes);
    try std.testing.expectEqual(@as(usize, 2), selection.indexes.len());
    try std.testing.expectEqual(@as(u64, 4), selection.indexes.at(0));
    try std.testing.expectEqual(wide_index, selection.indexes.at(1));
}

fn exerciseBatchedUpperBoundary(allocator: std.mem.Allocator) !void {
    var direct = ExactSelector.init(2, 11);
    defer direct.deinit(allocator);
    var batched = ExactSelector.init(2, 11);
    defer batched.deinit(allocator);
    for (1..4) |record_number| {
        try direct.considerRecord(allocator, record_number);
        try batched.considerRecord(allocator, record_number);
    }

    const before_boundary = (@as(u64, 1) << 24) - 2;
    direct.record_count = before_boundary;
    batched.record_count = before_boundary;
    for (before_boundary + 1..before_boundary + 4) |record_number| {
        try direct.considerRecord(allocator, record_number);
    }
    try batched.considerRecordsThrough(allocator, before_boundary + 3);

    const direct_indexes = direct.finish().indexes;
    const batched_indexes = batched.finish().indexes;
    try std.testing.expectEqual(direct_indexes.len(), batched_indexes.len());
    for (0..direct_indexes.len()) |index| {
        try std.testing.expectEqual(
            direct_indexes.at(index),
            batched_indexes.at(index),
        );
    }
}

fn exerciseUpperBitIndexSort(allocator: std.mem.Allocator) !void {
    var selector = ExactSelector.init(10, 11);
    defer selector.deinit(allocator);
    for (1..12) |record_number| {
        try selector.considerRecord(allocator, record_number);
    }

    selector.activateUpperBits(8);
    const unsorted = [_]u32{
        0x80_00_0001,
        0x00_ff_0002,
        0x01_00_ffff,
        0x00_00_0003,
        0xff_00_0000,
        0x7f_ff_0000,
        0x00_00_0001,
        0x03_00_0000,
        0x00_ff_0001,
        0x01_00_0000,
    };
    for (unsorted, 0..) |value, slot| selector.writeLowIndex(slot, value);
    selector.unordered = true;

    const expected = [_]u64{
        0x00_00_0001,
        0x00_00_0003,
        0x00_ff_0001,
        0x00_ff_0002,
        0x01_00_0000,
        0x01_00_ffff,
        0x03_00_0000,
        0x7f_ff_0000,
        0x80_00_0001,
        0xff_00_0000,
    };
    const selection = selector.finish();
    try std.testing.expect(selection == .indexes);
    for (expected, 0..) |value, index| {
        try std.testing.expectEqual(value, selection.indexes.at(index));
    }
}

test "[unit] - [exact selector]: compact indexes retain wide record numbers" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseWideExactIndexes,
        .{},
    );
}

test "[unit] - [exact selector]: split indexes sort upper planes" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseUpperBitIndexSort,
        .{},
    );
}

test "[unit] - [exact selector]: batched records preserve an upper-bit boundary" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseBatchedUpperBoundary,
        .{},
    );
}

pub const Mt19937_64 = struct {
    state: [STATE_WORDS]u64,
    index: usize,

    const STATE_WORDS = 312;
    const MIDDLE_WORD = 156;
    const MATRIX_A: u64 = 0xb5026f5aa96619e9;
    const UPPER_MASK: u64 = 0xffffffff80000000;
    const LOWER_MASK: u64 = 0x000000007fffffff;

    pub fn init(seed: u64) Mt19937_64 {
        var generator: Mt19937_64 = undefined;
        generator.state[0] = seed;
        for (1..STATE_WORDS) |index| {
            const previous = generator.state[index - 1];
            generator.state[index] = 6364136223846793005 *%
                (previous ^ (previous >> 62)) +% index;
        }
        generator.index = STATE_WORDS;
        return generator;
    }

    pub fn nextU64(self: *Mt19937_64) u64 {
        if (self.index == STATE_WORDS) self.twist();

        var value = self.state[self.index];
        self.index += 1;
        value ^= (value >> 29) & 0x5555555555555555;
        value ^= (value << 17) & 0x71d67fffeda60000;
        value ^= (value << 37) & 0xfff7eee000000000;
        value ^= value >> 43;
        return value;
    }

    pub fn nextUnitFloat(self: *Mt19937_64) f64 {
        const top_53 = self.nextU64() >> 11;
        return @as(f64, @floatFromInt(top_53)) * (1.0 / 9007199254740992.0);
    }

    fn twist(self: *Mt19937_64) void {
        for (0..STATE_WORDS - MIDDLE_WORD) |index| {
            const joined = (self.state[index] & UPPER_MASK) |
                (self.state[index + 1] & LOWER_MASK);
            self.state[index] = self.state[index + MIDDLE_WORD] ^
                (joined >> 1) ^ matrixTerm(joined);
        }
        for (STATE_WORDS - MIDDLE_WORD..STATE_WORDS - 1) |index| {
            const joined = (self.state[index] & UPPER_MASK) |
                (self.state[index + 1] & LOWER_MASK);
            self.state[index] = self.state[index + MIDDLE_WORD - STATE_WORDS] ^
                (joined >> 1) ^ matrixTerm(joined);
        }
        const joined = (self.state[STATE_WORDS - 1] & UPPER_MASK) |
            (self.state[0] & LOWER_MASK);
        self.state[STATE_WORDS - 1] = self.state[MIDDLE_WORD - 1] ^
            (joined >> 1) ^ matrixTerm(joined);
        self.index = 0;
    }

    fn matrixTerm(value: u64) u64 {
        return if (value & 1 == 0) 0 else MATRIX_A;
    }
};
