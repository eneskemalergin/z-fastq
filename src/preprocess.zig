//! Paired read trimming and pair retention with fastp 1.1.0's rules: the 3' quality window,
//! poly-G tails, adapters found by mate overlap or by adapter sequence, then the empty-read,
//! adapter-dimer, quality, N, and length checks. Every trim keeps a prefix of the read, so a
//! result is one kept length per mate; sequence and quality are cut together by the caller.
//! The overlap and adapter searches score a block of consecutive offsets per pass over the
//! pattern, reading the read through a fixed stack window padded with zero bytes, which no read
//! or adapter base equals. The block kernel, the reverse complement, and the low-quality and N
//! count use the inline vector width only, with no AVX2 object for builds whose target lacks
//! AVX2: each costs a small share of a trimmed or filtered run at every read length. Insert sizes
//! follow fastp's overlap rule. The pair duplication estimate is fastp's: the bases of mate 1, then
//! mate 2, hash to one position in each of two 2^32-bit arrays, and a pair whose two positions are
//! both set already counts as a duplicate, hash collisions included. Only the set positions are
//! kept, per 65,536-bit block as a sorted list of offsets or, past 4,096 of them, as a bitmap, so
//! memory follows the pairs read instead of fastp's fixed 1 GiB. The pair sum has an AVX2 kernel.

const std = @import("std");

const fastq = @import("fastq.zig");
const simd = @import("simd.zig");

pub const MIN_ADAPTER_BASES = 4;
pub const MAX_ADAPTER_BASES = 128;
pub const MAX_TAIL_WINDOW = 1000;
pub const MAX_TAIL_QUALITY = 30;
pub const MAX_N_LIMIT = 50;
pub const MAX_LOW_QUALITY = fastq.PHRED_MAX_SCORE;
pub const MAX_LOW_QUALITY_PERCENT = 100;

const PHRED_OFFSET = fastq.PHRED_OFFSET;
const OVERLAP_REQUIRE = 30;
const OVERLAP_DIFF_LIMIT = 5;
const OVERLAP_DIFF_DIVISOR = 5;
const OVERLAP_COMPARE_LIMIT = 50;
const ADAPTER_MATCH_REQUIRE = 4;
const ADAPTER_LEAD = 4;
const BASES_PER_MISMATCH = 8;
const POLY_G_MAX_MISMATCHES = 5;
const DIMER_MAX_BASES = 2;
const LANES = simd.INLINE_LANES;
const MAX_LANES = simd.tierWidth(.avx512);
const WINDOW_OFFSETS = 1024;
const WINDOW_SPAN = MAX_ADAPTER_BASES + MAX_LANES;
const WINDOW_BYTES = ADAPTER_LEAD + WINDOW_OFFSETS + WINDOW_SPAN;
const THRESHOLD_CLAMP = 1024;
const PRUNE_COLUMNS = 16;
const CASE_FOLD = 0xdf;

const COMPLEMENT: [256]u8 = table: {
    var table = [_]u8{'N'} ** 256;
    for ("ACGTacgt", "TGCATGCA") |base, complement| table[base] = complement;
    break :table table;
};

pub const Tail = struct {
    window: u16,
    quality: u8,
};

pub const LowQuality = struct {
    quality: u8,
    percent: u8,
};

pub const Options = struct {
    adapters: ?[2][]const u8 = null,
    poly_g: ?u32 = null,
    tail: ?Tail = null,
    min_length: ?u32 = null,
    max_n: ?u8 = null,
    low_quality: ?LowQuality = null,

    pub fn trims(self: Options) bool {
        return self.adapters != null or self.poly_g != null or self.tail != null;
    }

    pub fn filters(self: Options) bool {
        return self.min_length != null or self.max_n != null or self.low_quality != null;
    }

    pub fn changesReads(self: Options) bool {
        return self.trims() or self.filters();
    }
};

pub const INSERT_BINS = 512;

pub const InsertCounts = struct {
    pairs: u64 = 0,
    overlapping: u64 = 0,
    histogram: [INSERT_BINS + 1]u64 = @splat(0),

    pub fn add(self: *InsertCounts, result: *const Result) error{Overflow}!void {
        const lengths = result.overlap_lengths orelse return;
        self.pairs = try std.math.add(u64, self.pairs, 1);
        var bin: usize = INSERT_BINS;
        if (result.overlap) |overlap| {
            self.overlapping += 1;
            const size = if (overlap.offset > 0) lengths[0] + lengths[1] - overlap.length else overlap.length;
            bin = @min(size, INSERT_BINS);
        }
        self.histogram[bin] += 1;
    }

    pub fn unknown(self: *const InsertCounts) u64 {
        return self.histogram[INSERT_BINS];
    }

    pub fn peak(self: *const InsertCounts) usize {
        var best: usize = 0;
        for (self.histogram[0..INSERT_BINS], 0..) |count, size| {
            if (count > self.histogram[best]) best = size;
        }
        return best;
    }

    pub fn median(self: *const InsertCounts) ?usize {
        const known = self.pairs - self.unknown();
        if (known == 0) return null;
        const half = known - known / 2;
        var seen: u64 = 0;
        for (self.histogram[0..INSERT_BINS], 0..) |count, size| {
            seen += count;
            if (seen >= half) return size;
        }
        unreachable;
    }
};

pub const PairDuplication = struct {
    sets: [ARRAY_COUNT]PositionSet = .{ .{}, .{} },
    primes: Primes = undefined,
    pending: ?Hash = null,
    pairs: u64 = 0,
    duplicates: u64 = 0,

    const ARRAY_COUNT = 2;
    const PRIME_ROW = 512;
    const PRIME_START = 10_000;
    const PRIME_GAP = 10_000;
    const OTHER_VALUE = 13;

    const Hash = [ARRAY_COUNT]u32;
    const Sums = [ARRAY_COUNT]u64;
    const SUM_LANES = 8;
    const Primes = [ARRAY_COUNT][PRIME_ROW + SUM_LANES]u32;
    const VECTOR_LENGTH_LIMIT = 1 << 31;

    const BASE_VALUE: [256]u64 = table: {
        var values: [256]u64 = @splat(OTHER_VALUE);
        values['A'] = 7;
        values['C'] = 74;
        values['G'] = 31;
        values['T'] = 222;
        break :table values;
    };

    const DIVISORS = divisors: {
        @setEvalBranchQuota(1_000_000);
        var found: [600]u32 = undefined;
        var count: usize = 0;
        var candidate: u32 = 2;
        while (candidate < 4096) : (candidate += 1) {
            var prime = true;
            for (found[0..count]) |divisor| {
                if (divisor * divisor > candidate) break;
                if (candidate % divisor == 0) {
                    prime = false;
                    break;
                }
            }
            if (prime) {
                found[count] = candidate;
                count += 1;
            }
        }
        const result = found[0..count].*;
        break :divisors result;
    };

    pub fn init(self: *PairDuplication, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        self.* = .{};
        errdefer self.deinit(allocator);
        for (&self.sets) |*set| try set.init(allocator);
        var number: u32 = PRIME_START;
        for (0..PRIME_ROW * ARRAY_COUNT) |index| {
            number += 1;
            while (!isPrime(number)) number += 1;
            self.primes[index % ARRAY_COUNT][index / ARRAY_COUNT] = number;
            number += PRIME_GAP;
        }
        for (&self.primes) |*row| @memcpy(row[PRIME_ROW..], row[0..SUM_LANES]);
    }

    fn isPrime(number: u32) bool {
        std.debug.assert(number < 4096 * 4096);
        for (DIVISORS) |divisor| {
            if (divisor * divisor > number) return true;
            if (number % divisor == 0) return false;
        }
        return true;
    }

    pub fn deinit(self: *PairDuplication, allocator: std.mem.Allocator) void {
        for (&self.sets) |*set| set.deinit(allocator);
    }

    pub fn hash(self: *const PairDuplication, read1: []const u8, read2: []const u8) Hash {
        var sums: Sums = @splat(0);
        if (read1.len + read2.len >= VECTOR_LENGTH_LIMIT) {
            sumsScalar(&self.primes, read1, 0, &sums);
            sumsScalar(&self.primes, read2, read1.len, &sums);
        } else if (simd.INLINE_TIER == .baseline and simd.avx2Selected()) {
            simd.avx2Kernel("pair_sums", &pairSumsExport)(&self.primes, read1.ptr, read1.len, read2.ptr, read2.len, &sums);
        } else if (simd.INLINE_TIER == .avx2) {
            pairSumsVector(&self.primes, read1, read2, &sums);
        } else {
            sumsScalar(&self.primes, read1, 0, &sums);
            sumsScalar(&self.primes, read2, read1.len, &sums);
        }
        var bits: Hash = undefined;
        for (sums, &bits, &self.sets) |sum, *bit, *set| {
            bit.* = @truncate(sum);
            set.prefetch(bit.*);
        }
        return bits;
    }

    fn sumsScalar(primes: *const Primes, read: []const u8, start: usize, sums: *Sums) void {
        for (read, start..) |base, position| {
            const value = BASE_VALUE[base] +% position;
            for (sums, primes) |*sum, *row| sum.* +%= @as(u64, row[position % PRIME_ROW]) *% value;
        }
    }

    fn pairSumsVector(primes: *const Primes, read1: []const u8, read2: []const u8, sums: *Sums) void {
        sumsVector(primes, read1, 0, sums);
        sumsVector(primes, read2, read1.len, sums);
    }

    fn pairSumsExport(primes: *const Primes, read1: [*]const u8, len1: usize, read2: [*]const u8, len2: usize, sums: *Sums) callconv(.c) void {
        pairSumsVector(primes, read1[0..len1], read2[0..len2], sums);
    }

    fn sumsVector(primes: *const Primes, read: []const u8, start: usize, sums: *Sums) void {
        const Lanes32 = @Vector(SUM_LANES, u32);
        const Lanes64 = @Vector(SUM_LANES, u64);
        var totals: [ARRAY_COUNT]Lanes64 = @splat(@splat(0));
        const steps: Lanes32 = std.simd.iota(u32, SUM_LANES);
        var index: usize = 0;
        while (index + SUM_LANES <= read.len) : (index += SUM_LANES) {
            const bytes: @Vector(SUM_LANES, u8) = read[index..][0..SUM_LANES].*;
            var values: Lanes32 = @splat(OTHER_VALUE);
            inline for (.{ 'A', 'C', 'G', 'T' }) |base| {
                const matches = bytes == @as(@Vector(SUM_LANES, u8), @splat(base));
                values = @select(u32, matches, @as(Lanes32, @splat(@intCast(BASE_VALUE[base]))), values);
            }
            const position = start + index;
            const terms: Lanes64 = values + @as(Lanes32, @splat(@intCast(position))) + steps;
            const row = position % PRIME_ROW;
            for (&totals, primes) |*total, *factors| {
                const wide: Lanes64 = @as(Lanes32, factors[row..][0..SUM_LANES].*);
                total.* +%= wide * terms;
            }
        }
        for (&totals, sums) |total, *sum| {
            for (@as([SUM_LANES]u64, total)) |lane| sum.* +%= lane;
        }
        sumsScalar(primes, read[index..], start + index, sums);
    }

    pub fn push(self: *PairDuplication, allocator: std.mem.Allocator, bits: Hash) error{ OutOfMemory, Overflow }!void {
        if (self.pending) |previous| try self.record(allocator, previous);
        for (&self.sets, bits) |*set, bit| set.prefetchBlock(bit);
        self.pending = bits;
    }

    pub fn flush(self: *PairDuplication, allocator: std.mem.Allocator) error{ OutOfMemory, Overflow }!void {
        const previous = self.pending orelse return;
        try self.record(allocator, previous);
        self.pending = null;
    }

    fn record(self: *PairDuplication, allocator: std.mem.Allocator, bits: Hash) error{ OutOfMemory, Overflow }!void {
        const pairs = try std.math.add(u64, self.pairs, 1);
        var probes: [ARRAY_COUNT]PositionSet.Probe = undefined;
        for (&probes, &self.sets, bits) |*probe, *set, bit| probe.* = set.find(bit);
        for (&probes, &self.sets) |*probe, *set| try set.reserve(allocator, probe);
        var seen = true;
        for (probes) |probe| {
            if (probe.present) continue;
            seen = false;
            PositionSet.insert(probe);
        }
        self.pairs = pairs;
        if (seen) self.duplicates += 1;
    }

    pub fn rate(self: *const PairDuplication) f64 {
        if (self.pairs == 0) return 0;
        return @as(f64, @floatFromInt(self.duplicates)) / @as(f64, @floatFromInt(self.pairs));
    }

    pub fn writeJson(self: *const PairDuplication, json: *std.json.Stringify) !void {
        try json.beginObject();
        try json.objectField("pairs");
        try json.write(self.pairs);
        try json.objectField("duplicates");
        try json.write(self.duplicates);
        try json.objectField("rate");
        try json.write(self.rate());
        try json.endObject();
    }
};

comptime {
    simd.exportKernel("pair_sums", &PairDuplication.pairSumsExport);
}

const PositionSet = struct {
    blocks: []Block = &.{},

    const BLOCK_COUNT = 1 << 16;
    const LIST_MAX = 4096;
    const BITMAP_WORDS = BLOCK_COUNT / 64;
    const BITMAP = std.math.maxInt(u32);

    const Block = struct {
        data: ?*anyopaque = null,
        count: u32 = 0,
        capacity: u32 = 0,

        fn list(self: Block) []u16 {
            const items: [*]u16 = @ptrCast(@alignCast(self.data.?));
            return items[0..self.count];
        }

        fn bitmap(self: Block) *[BITMAP_WORDS]u64 {
            return @ptrCast(@alignCast(self.data.?));
        }
    };

    const Probe = struct {
        block: *Block,
        low: u16,
        index: u32,
        present: bool,
    };

    fn init(self: *PositionSet, allocator: std.mem.Allocator) error{OutOfMemory}!void {
        self.blocks = try allocator.alloc(Block, BLOCK_COUNT);
        @memset(self.blocks, .{});
    }

    fn deinit(self: *PositionSet, allocator: std.mem.Allocator) void {
        for (self.blocks) |block| {
            const data = block.data orelse continue;
            if (block.count == BITMAP) {
                allocator.destroy(block.bitmap());
            } else {
                const items: [*]u16 = @ptrCast(@alignCast(data));
                allocator.free(items[0..block.capacity]);
            }
        }
        allocator.free(self.blocks);
        self.blocks = &.{};
    }

    fn prefetch(self: *const PositionSet, position: u32) void {
        @prefetch(&self.blocks[position >> 16], .{ .rw = .write });
    }

    fn prefetchBlock(self: *const PositionSet, position: u32) void {
        const block = self.blocks[position >> 16];
        if (block.data == null) return;
        const low: u16 = @truncate(position);
        if (block.count == BITMAP) {
            @prefetch(&block.bitmap()[low >> 6], .{ .rw = .write });
        } else {
            const items = block.list();
            @prefetch(items.ptr, .{ .rw = .write });
            if (items.len > 32) @prefetch(&items[items.len / 2], .{});
        }
    }

    fn find(self: *PositionSet, position: u32) Probe {
        const block = &self.blocks[position >> 16];
        const low: u16 = @truncate(position);
        if (block.count == BITMAP) {
            const present = block.bitmap()[low >> 6] & (@as(u64, 1) << @intCast(low & 63)) != 0;
            return .{ .block = block, .low = low, .index = 0, .present = present };
        }
        const index = if (block.data == null) 0 else std.sort.lowerBound(u16, block.list(), low, orderU16);
        const present = index < block.count and block.list()[index] == low;
        return .{ .block = block, .low = low, .index = @intCast(index), .present = present };
    }

    fn orderU16(wanted: u16, item: u16) std.math.Order {
        return std.math.order(wanted, item);
    }

    fn reserve(self: *PositionSet, allocator: std.mem.Allocator, probe: *Probe) error{OutOfMemory}!void {
        _ = self;
        const block = probe.block;
        if (probe.present or block.count == BITMAP or block.count < block.capacity) return;
        if (block.count == LIST_MAX) {
            const words = try allocator.create([BITMAP_WORDS]u64);
            @memset(words, 0);
            const items = block.list();
            for (items) |item| words[item >> 6] |= @as(u64, 1) << @intCast(item & 63);
            allocator.free(@as([*]u16, @ptrCast(@alignCast(block.data.?)))[0..block.capacity]);
            block.* = .{ .data = words, .count = BITMAP, .capacity = 0 };
            return;
        }
        const capacity: u32 = @min(LIST_MAX, @max(4, block.capacity * 2));
        const old: []u16 = if (block.data) |data| @as([*]u16, @ptrCast(@alignCast(data)))[0..block.capacity] else &.{};
        const items = try allocator.realloc(old, capacity);
        block.data = items.ptr;
        block.capacity = capacity;
    }

    fn insert(probe: Probe) void {
        std.debug.assert(!probe.present);
        const block = probe.block;
        if (block.count == BITMAP) {
            block.bitmap()[probe.low >> 6] |= @as(u64, 1) << @intCast(probe.low & 63);
            return;
        }
        std.debug.assert(block.count < block.capacity);
        const items: [*]u16 = @ptrCast(@alignCast(block.data.?));
        var index = block.count;
        while (index > probe.index) : (index -= 1) items[index] = items[index - 1];
        items[probe.index] = probe.low;
        block.count += 1;
    }
};

pub const Reason = enum { too_many_n, too_short, low_quality, adapter_dimer };

pub const REASON_COUNT = @typeInfo(Reason).@"enum".fields.len;

pub const Decision = struct {
    reads: [2]?Reason,
    pair: ?Reason,
};

pub const FilterCounts = struct {
    dropped_pairs: [REASON_COUNT]u64 = @splat(0),
    failed_reads: [2][REASON_COUNT]u64 = @splat(@splat(0)),

    pub fn add(self: *FilterCounts, decision: Decision) error{Overflow}!void {
        const pair = decision.pair orelse return;
        const dropped = &self.dropped_pairs[@intFromEnum(pair)];
        dropped.* = try std.math.add(u64, dropped.*, 1);
        for (&self.failed_reads, decision.reads) |*mate, reason| {
            const failed = &mate[@intFromEnum(reason orelse continue)];
            failed.* = try std.math.add(u64, failed.*, 1);
        }
    }
};

pub fn decide(options: Options, result: *const Result, sequences: [2][]const u8, qualities: [2][]const u8) Decision {
    std.debug.assert(sequences[0].len == result.kept[0] and sequences[1].len == result.kept[1]);
    var reads = [2]?Reason{ readReason(options, sequences[0], qualities[0]), readReason(options, sequences[1], qualities[1]) };
    const cut = result.adapter_cut[0] or result.adapter_cut[1];
    if (cut and sequences[0].len <= DIMER_MAX_BASES and sequences[1].len <= DIMER_MAX_BASES) {
        reads = .{ .adapter_dimer, .adapter_dimer };
    }
    return .{ .reads = reads, .pair = higherReason(reads[0], reads[1]) };
}

fn higherReason(first: ?Reason, second: ?Reason) ?Reason {
    const a = first orelse return second;
    const b = second orelse return first;
    return if (@intFromEnum(a) >= @intFromEnum(b)) a else b;
}

const LowAndN = struct {
    low: usize,
    n: usize,
};

fn countLowAndN(comptime lanes: comptime_int, sequence: []const u8, quality: []const u8, threshold: u8) LowAndN {
    std.debug.assert(sequence.len == quality.len);
    const Bytes = @Vector(lanes, u8);
    const Wide = @Vector(lanes, u16);
    const one: Bytes = @splat(1);
    const zero: Bytes = @splat(0);
    var counts: LowAndN = .{ .low = 0, .n = 0 };
    const whole = sequence.len - sequence.len % lanes;
    var offset: usize = 0;
    while (offset < whole) {
        const stop = @min(whole, offset + std.math.maxInt(u8) * lanes);
        var low = zero;
        var n_bases = zero;
        while (offset < stop) : (offset += lanes) {
            const bytes: Bytes = quality[offset..][0..lanes].*;
            const bases: Bytes = sequence[offset..][0..lanes].*;
            low += @select(u8, bytes < @as(Bytes, @splat(threshold)), one, zero);
            n_bases += @select(u8, bases == @as(Bytes, @splat('N')), one, zero);
        }
        counts.low += @reduce(.Add, @as(Wide, low));
        counts.n += @reduce(.Add, @as(Wide, n_bases));
    }
    for (sequence[whole..], quality[whole..]) |base, byte| {
        counts.low += @intFromBool(byte < threshold);
        counts.n += @intFromBool(base == 'N');
    }
    return counts;
}

fn readReason(options: Options, sequence: []const u8, quality: []const u8) ?Reason {
    std.debug.assert(sequence.len == quality.len);
    if (sequence.len == 0) return .too_short;
    if (options.low_quality != null or options.max_n != null) {
        const threshold: u8 = if (options.low_quality) |low_quality| PHRED_OFFSET + low_quality.quality else 0;
        const counts = countLowAndN(LANES, sequence, quality, threshold);
        if (options.low_quality) |low_quality| {
            if (counts.low * 100 > @as(usize, low_quality.percent) * sequence.len) return .low_quality;
        }
        if (options.max_n) |max_n| {
            if (counts.n > max_n) return .too_many_n;
        }
    }
    if (options.min_length) |min_length| {
        if (sequence.len < min_length) return .too_short;
    }
    return null;
}

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
    overlap_lengths: ?[2]usize = null,
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
        result.overlap_lengths = .{ reads[0].len, reads[1].len };
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
            const end = self.source.len - begin;
            reverseComplement(LANES, target, self.source[end - count .. end]);
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
    reverseComplement(LANES, head2[0..head2_len], read2[read2.len - head2_len ..]);
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

fn reverseComplement(comptime lanes: comptime_int, target: []u8, source: []const u8) void {
    std.debug.assert(target.len == source.len);
    const Bytes = @Vector(lanes, u8);
    const reverse = comptime mask: {
        var indexes: [lanes]i32 = undefined;
        for (&indexes, 0..) |*index, lane| index.* = lanes - 1 - lane;
        break :mask indexes;
    };
    var done: usize = 0;
    while (done + lanes <= target.len) : (done += lanes) {
        const block: Bytes = source[source.len - done - lanes ..][0..lanes].*;
        const folded = @shuffle(u8, block, undefined, reverse) & @as(Bytes, @splat(CASE_FOLD));
        var complement: Bytes = @splat('N');
        inline for ("ACGT", "TGCA") |base, paired| {
            complement = @select(u8, folded == @as(Bytes, @splat(base)), @as(Bytes, @splat(paired)), complement);
        }
        target[done..][0..lanes].* = complement;
    }
    for (target[done..], done..) |*byte, index| byte.* = COMPLEMENT[source[source.len - 1 - index]];
}

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

fn decideReads(options: Options, reads: [2][]const u8, qualities: [2][]const u8, adapter_cut: [2]bool) Decision {
    const result: Result = .{ .kept = .{ reads[0].len, reads[1].len }, .adapter_cut = adapter_cut };
    return decide(options, &result, reads, qualities);
}

test "[unit] - [preprocess]: each read check fails on one side of its threshold and passes on the other" {
    const quality: Options = .{ .low_quality = .{ .quality = 15, .percent = 40 } };
    const Case = struct { options: Options, sequence: []const u8, quality: []const u8, expected: ?Reason };
    const cases = [_]Case{
        .{ .options = .{}, .sequence = "", .quality = "", .expected = .too_short },
        .{ .options = quality, .sequence = "ACGTA", .quality = "//000", .expected = null },
        .{ .options = quality, .sequence = "ACGTA", .quality = "///00", .expected = .low_quality },
        .{ .options = quality, .sequence = "ACG", .quality = "/00", .expected = null },
        .{ .options = quality, .sequence = "ACG", .quality = "//0", .expected = .low_quality },
        .{ .options = .{ .low_quality = .{ .quality = 15, .percent = 0 } }, .sequence = "ACG", .quality = "000", .expected = null },
        .{ .options = .{ .low_quality = .{ .quality = 15, .percent = 0 } }, .sequence = "ACG", .quality = "00/", .expected = .low_quality },
        .{ .options = .{ .low_quality = .{ .quality = 0, .percent = 0 } }, .sequence = "A", .quality = "!", .expected = null },
        .{ .options = .{ .max_n = 5 }, .sequence = "NNNNNA", .quality = "IIIIII", .expected = null },
        .{ .options = .{ .max_n = 5 }, .sequence = "NNNNNN", .quality = "IIIIII", .expected = .too_many_n },
        .{ .options = .{ .max_n = 0 }, .sequence = "AnnnA", .quality = "IIIII", .expected = null },
        .{ .options = .{ .max_n = 0 }, .sequence = "ANA", .quality = "III", .expected = .too_many_n },
        .{ .options = .{ .min_length = 15 }, .sequence = "A" ** 14, .quality = "I" ** 14, .expected = .too_short },
        .{ .options = .{ .min_length = 15 }, .sequence = "A" ** 15, .quality = "I" ** 15, .expected = null },
        .{ .options = .{ .min_length = 30 }, .sequence = "A" ** 29, .quality = "I" ** 29, .expected = .too_short },
        .{ .options = .{ .min_length = 30 }, .sequence = "A" ** 30, .quality = "I" ** 30, .expected = null },
        .{ .options = .{ .min_length = 30 }, .sequence = "A" ** 31, .quality = "I" ** 31, .expected = null },
        .{ .options = .{ .min_length = 0 }, .sequence = "", .quality = "", .expected = .too_short },
        .{ .options = .{ .low_quality = .{ .quality = 15, .percent = 40 }, .max_n = 0, .min_length = 9 }, .sequence = "NNNN", .quality = "!!!!", .expected = .low_quality },
        .{ .options = .{ .max_n = 0, .min_length = 9 }, .sequence = "NNNN", .quality = "!!!!", .expected = .too_many_n },
        .{ .options = .{ .min_length = 9 }, .sequence = "NNNN", .quality = "!!!!", .expected = .too_short },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        try std.testing.expectEqual(case.expected, readReason(case.options, case.sequence, case.quality));
    }
}

test "[unit] - [preprocess]: a pair keeps both mates or drops both, under the higher reason" {
    const options: Options = .{ .adapters = .{ "ACGT", "ACGT" }, .low_quality = .{ .quality = 15, .percent = 40 }, .max_n = 1, .min_length = 4 };
    const good = "ACGTA";
    const good_quality = "IIIII";
    const Case = struct { reads: [2][]const u8, qualities: [2][]const u8, cut: [2]bool = .{ false, false }, expected: Decision };
    const cases = [_]Case{
        .{ .reads = .{ good, good }, .qualities = .{ good_quality, good_quality }, .expected = .{ .reads = .{ null, null }, .pair = null } },
        .{ .reads = .{ "NNACG", good }, .qualities = .{ good_quality, good_quality }, .expected = .{ .reads = .{ .too_many_n, null }, .pair = .too_many_n } },
        .{ .reads = .{ good, "ACG" }, .qualities = .{ good_quality, "III" }, .expected = .{ .reads = .{ null, .too_short }, .pair = .too_short } },
        .{ .reads = .{ "NNACG", "ACG" }, .qualities = .{ good_quality, "III" }, .expected = .{ .reads = .{ .too_many_n, .too_short }, .pair = .too_short } },
        .{ .reads = .{ "ACG", good }, .qualities = .{ "III", "!!!II" }, .expected = .{ .reads = .{ .too_short, .low_quality }, .pair = .low_quality } },
        .{ .reads = .{ "", good }, .qualities = .{ "", good_quality }, .expected = .{ .reads = .{ .too_short, null }, .pair = .too_short } },
        .{ .reads = .{ "AC", "G" }, .qualities = .{ "II", "I" }, .cut = .{ true, false }, .expected = .{ .reads = .{ .adapter_dimer, .adapter_dimer }, .pair = .adapter_dimer } },
        .{ .reads = .{ "", "" }, .qualities = .{ "", "" }, .cut = .{ false, true }, .expected = .{ .reads = .{ .adapter_dimer, .adapter_dimer }, .pair = .adapter_dimer } },
        .{ .reads = .{ "AC", "G" }, .qualities = .{ "II", "I" }, .expected = .{ .reads = .{ .too_short, .too_short }, .pair = .too_short } },
        .{ .reads = .{ "ACG", "G" }, .qualities = .{ "III", "I" }, .cut = .{ true, true }, .expected = .{ .reads = .{ .too_short, .too_short }, .pair = .too_short } },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        try std.testing.expectEqual(case.expected, decideReads(options, case.reads, case.qualities, case.cut));
    }
    var counts: FilterCounts = .{};
    for (cases) |case| try counts.add(decideReads(options, case.reads, case.qualities, case.cut));
    try std.testing.expectEqual([REASON_COUNT]u64{ 1, 5, 1, 2 }, counts.dropped_pairs);
    try std.testing.expectEqual([REASON_COUNT]u64{ 2, 4, 0, 2 }, counts.failed_reads[0]);
    try std.testing.expectEqual([REASON_COUNT]u64{ 0, 4, 1, 2 }, counts.failed_reads[1]);
}

test "[property] - [preprocess]: vector low-quality and N counts match a scalar count at every tier width" {
    var sequence: [255 * 64 + 65]u8 = undefined;
    var quality: [255 * 64 + 65]u8 = undefined;
    for (&sequence, &quality, 0..) |*base, *byte, index| {
        base.* = "ACGTNnRN"[index * 7 % 8];
        byte.* = @intCast(33 + (index * 13) % 94);
    }
    for ([_]usize{ 0, 1, 15, 16, 17, 63, 64, 65, 255 * 16 + 1, 255 * 32 + 1, 255 * 64 + 1, 255 * 64 + 65 }) |length| {
        for ([_]u8{ 0, 48, 126, 127 }) |threshold| {
            var expected: LowAndN = .{ .low = 0, .n = 0 };
            for (sequence[0..length], quality[0..length]) |base, byte| {
                expected.low += @intFromBool(byte < threshold);
                expected.n += @intFromBool(base == 'N');
            }
            inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |lanes| {
                try std.testing.expectEqual(expected, countLowAndN(lanes, sequence[0..length], quality[0..length], threshold));
            }
        }
    }
}

test "[property] - [preprocess]: vector reverse complement matches the scalar table at every tier width" {
    var source: [150]u8 = undefined;
    for (&source, 0..) |*byte, index| byte.* = "ACGTacgtNnRYSWKMBDHV.-*XxIi"[index % 27];
    for ([_]usize{ 0, 15, 16, 17, 63, 64, 65, 150 }) |length| {
        var expected: [150]u8 = undefined;
        for (expected[0..length], 0..) |*byte, index| byte.* = COMPLEMENT[source[length - 1 - index]];
        inline for (.{ simd.tierWidth(.baseline), simd.tierWidth(.avx2), simd.tierWidth(.avx512) }) |lanes| {
            var target: [150]u8 = undefined;
            reverseComplement(lanes, target[0..length], source[0..length]);
            try std.testing.expectEqualSlices(u8, expected[0..length], target[0..length]);
        }
    }
}

fn insertOf(options: Options, read1: []const u8, read2: []const u8) InsertCounts {
    var counts: InsertCounts = .{};
    const quality1 = ("I" ** 200)[0..read1.len];
    const quality2 = ("I" ** 200)[0..read2.len];
    const trimmer: Trimmer = .init(options);
    const result = trimmer.trimPair(.{ read1, read2 }, .{ quality1, quality2 });
    counts.add(&result) catch unreachable;
    return counts;
}

fn reverseComplementOf(comptime text: []const u8) [text.len]u8 {
    var out: [text.len]u8 = undefined;
    for (&out, 0..) |*byte, index| byte.* = COMPLEMENT[text[text.len - 1 - index]];
    return out;
}

test "[unit] - [preprocess]: insert sizes follow fastp's overlap formula and its unknown bin" {
    const fragment = "ACCTTAATGGGCGTGGACGTCTCTTGCTTGGAGGCTGGTGTGCACCCCTCATTATGACATCTCAGGAGAGGTATAGAGGAG";
    const tail = "AAAGCGGCACTTGTGAAGTGTTCCCCACGCCGCTTGGGTCTTCTGTGTTGTTC";
    const options: Options = .{ .adapters = .{ "CTGTCTCTTATACACATCT", "CTGTCTCTTATACACATCT" } };
    const Case = struct { insert: usize, length1: usize, length2: usize, expected: ?usize };
    const cases = [_]Case{
        .{ .insert = 40, .length1 = 50, .length2 = 50, .expected = 40 },
        .{ .insert = 31, .length1 = 50, .length2 = 50, .expected = 31 },
        .{ .insert = 30, .length1 = 50, .length2 = 50, .expected = null },
        .{ .insert = 50, .length1 = 50, .length2 = 50, .expected = 50 },
        .{ .insert = 60, .length1 = 50, .length2 = 50, .expected = 60 },
        .{ .insert = 69, .length1 = 50, .length2 = 50, .expected = 69 },
        .{ .insert = 70, .length1 = 50, .length2 = 50, .expected = null },
        .{ .insert = 45, .length1 = 50, .length2 = 40, .expected = 50 },
        .{ .insert = 45, .length1 = 40, .length2 = 50, .expected = 40 },
    };
    inline for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        const insert = fragment[0..case.insert];
        var read1: [100]u8 = undefined;
        var read2: [100]u8 = undefined;
        const source1 = insert ++ tail;
        @memcpy(read1[0..case.length1], source1[0..case.length1]);
        const complement = reverseComplementOf(insert);
        var source2: [complement.len + tail.len]u8 = undefined;
        @memcpy(source2[0..complement.len], &complement);
        @memcpy(source2[complement.len..], tail);
        @memcpy(read2[0..case.length2], source2[0..case.length2]);
        const counts = insertOf(options, read1[0..case.length1], read2[0..case.length2]);
        try std.testing.expectEqual(@as(u64, 1), counts.pairs);
        try std.testing.expectEqual(@as(u64, @intFromBool(case.expected != null)), counts.overlapping);
        try std.testing.expectEqual(@as(u64, @intFromBool(case.expected == null)), counts.unknown());
        if (case.expected) |size| try std.testing.expectEqual(@as(u64, 1), counts.histogram[size]);
    }
    const no_adapters = insertOf(.{ .poly_g = 10 }, fragment[0..40], fragment[0..40]);
    try std.testing.expectEqual(@as(u64, 0), no_adapters.pairs);
    const emptied = insertOf(.{ .adapters = options.adapters, .tail = .{ .window = 4, .quality = 20 } }, fragment[0..40], "ACGT");
    try std.testing.expectEqual(@as(u64, 0), emptied.pairs);
}

test "[unit] - [preprocess]: the insert peak takes the lowest tied size and the median the lower middle" {
    var counts: InsertCounts = .{};
    try std.testing.expectEqual(@as(usize, 0), counts.peak());
    try std.testing.expectEqual(@as(?usize, null), counts.median());
    for ([_]usize{ 40, 40, 60, 60, 100 }) |size| {
        const result: Result = .{ .kept = .{ size, size }, .overlap = .{ .offset = -1, .length = size }, .overlap_lengths = .{ size, size } };
        try counts.add(&result);
    }
    const unknown: Result = .{ .kept = .{ 50, 50 }, .overlap_lengths = .{ 50, 50 } };
    try counts.add(&unknown);
    const large: Result = .{ .kept = .{ 600, 600 }, .overlap = .{ .offset = -1, .length = 512 }, .overlap_lengths = .{ 600, 600 } };
    try counts.add(&large);
    try std.testing.expectEqual(@as(u64, 7), counts.pairs);
    try std.testing.expectEqual(@as(u64, 6), counts.overlapping);
    try std.testing.expectEqual(@as(u64, 2), counts.unknown());
    try std.testing.expectEqual(@as(usize, 40), counts.peak());
    try std.testing.expectEqual(@as(?usize, 60), counts.median());
}

test "[unit] - [preprocess]: pair duplicates follow fastp's hash of R1 then R2" {
    var pairs: PairDuplication = undefined;
    try pairs.init(std.testing.allocator);
    defer pairs.deinit(std.testing.allocator);
    try std.testing.expectEqual(@as(u32, 10007), pairs.primes[0][0]);
    const Case = struct { []const u8, []const u8, bool };
    const cases = [_]Case{
        .{ "ACGTACGT", "TTGGCCAA", false },
        .{ "ACGTACGT", "TTGGCCAA", true },
        .{ "TTGGCCAA", "ACGTACGT", false },
        .{ "ACGTACGT", "TTGGCCAT", false },
        .{ "ACGTAC", "GTTTGGCCAA", true },
        .{ "acgtNNNN", "TTGGCCAA", false },
        .{ "NNNNnnnn", "TTGGCCAA", true },
        .{ "", "", false },
        .{ "", "", true },
    };
    var duplicates: u64 = 0;
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        try pairs.record(std.testing.allocator, pairs.hash(case[0], case[1]));
        duplicates += @intFromBool(case[2]);
        try std.testing.expectEqual(duplicates, pairs.duplicates);
    }
    try std.testing.expectEqual(@as(u64, cases.len), pairs.pairs);
    try std.testing.expectEqual(4.0 / 9.0, pairs.rate());

    var streamed: PairDuplication = undefined;
    try streamed.init(std.testing.allocator);
    defer streamed.deinit(std.testing.allocator);
    for (cases) |case| try streamed.push(std.testing.allocator, streamed.hash(case[0], case[1]));
    try std.testing.expectEqual(@as(u64, cases.len - 1), streamed.pairs);
    try streamed.flush(std.testing.allocator);
    try std.testing.expectEqual(pairs.pairs, streamed.pairs);
    try std.testing.expectEqual(pairs.duplicates, streamed.duplicates);
}

test "[property] - [preprocess]: the position set answers membership as a full bitmap does" {
    var set: PositionSet = .{};
    try set.init(std.testing.allocator);
    defer set.deinit(std.testing.allocator);
    var reference = std.AutoHashMap(u32, void).init(std.testing.allocator);
    defer reference.deinit();
    var state: u64 = 5;
    for (0..30_000) |step| {
        state = state *% 6364136223846793005 +% 1442695040888963407;
        const random: u32 = @truncate(state >> 32);
        const position = switch (step % 3) {
            0 => random,
            1 => 0x0007_0000 | (random & 0x3fff),
            else => 0xffff_0000 | (random & 0x00ff),
        };
        var probe = set.find(position);
        try std.testing.expectEqual(reference.contains(position), probe.present);
        try set.reserve(std.testing.allocator, &probe);
        if (!probe.present) PositionSet.insert(probe);
        try reference.put(position, {});
    }
    try std.testing.expectEqual(PositionSet.BITMAP, set.blocks[7].count);
    try std.testing.expect(set.blocks[0xffff].count <= 256);
}

test "[failure] - [preprocess]: a failed position allocation leaves the pair counts unchanged" {
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    var pairs: PairDuplication = undefined;
    try pairs.init(failing.allocator());
    defer pairs.deinit(failing.allocator());
    failing.fail_index = failing.alloc_index;
    try std.testing.expectError(error.OutOfMemory, pairs.record(failing.allocator(), pairs.hash("ACGT", "TTGG")));
    try std.testing.expectEqual(@as(u64, 0), pairs.pairs);
    failing.fail_index = std.math.maxInt(usize);
    try pairs.record(failing.allocator(), pairs.hash("ACGT", "TTGG"));
    try pairs.record(failing.allocator(), pairs.hash("ACGT", "TTGG"));
    try std.testing.expectEqual(@as(u64, 2), pairs.pairs);
    try std.testing.expectEqual(@as(u64, 1), pairs.duplicates);
}

test "[property] - [preprocess]: the vector pair sums equal fastp's per-base sums" {
    var pairs: PairDuplication = undefined;
    try pairs.init(std.testing.allocator);
    defer pairs.deinit(std.testing.allocator);
    var read: [1500]u8 = undefined;
    var state: u32 = 9;
    for (&read) |*base| {
        state = state *% 1103515245 +% 12345;
        base.* = "ACGTNacgtR"[(state >> 16) % 10];
    }
    for ([_]usize{ 0, 1, 7, 8, 9, 151, 300, 505, 511, 512, 513, 1500 }) |length1| {
        for ([_]usize{ 0, 1, 8, 151, 512, 1300 }) |length2| {
            errdefer std.debug.print("lengths {d} {d}\n", .{ length1, length2 });
            const read2 = read[(length1 * 7) % 200 ..][0..length2];
            var scalar: PairDuplication.Sums = @splat(0);
            PairDuplication.sumsScalar(&pairs.primes, read[0..length1], 0, &scalar);
            PairDuplication.sumsScalar(&pairs.primes, read2, length1, &scalar);
            var vector: PairDuplication.Sums = @splat(0);
            PairDuplication.pairSumsVector(&pairs.primes, read[0..length1], read2, &vector);
            try std.testing.expectEqual(scalar, vector);
        }
    }
}
