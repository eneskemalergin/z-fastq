//! Decodes checked RFC 1951 streams with a bounded dynamic-block fast path.

const std = @import("std");
const builtin = @import("builtin");
const assert = std.debug.assert;
const flate = std.compress.flate;
const testing = std.testing;
const Writer = std.Io.Writer;
const Reader = std.Io.Reader;
const Container = flate.Container;

// --- Token codes ---

const MIN_LENGTH = 3;
const MAX_LENGTH = 258;

const MIN_DISTANCE = 1;
const MAX_DISTANCE = std.compress.flate.history_len;

const CODEGEN_ORDER: [19]u8 = .{
    16, 17, 18,
    0,  8,  7,
    9,  6,  10,
    5,  11, 4,
    12, 3,  13,
    2,  14, 1,
    15,
};

const LenCode = if (builtin.mode != .ReleaseSmall) LookupLenCode else ShortLenCode;
const DistCode = if (builtin.mode != .ReleaseSmall) LookupDistCode else ShortDistCode;
const ShortLenCode = ShortCode(u8, u2, u3, true);
const ShortDistCode = ShortCode(u15, u1, u4, false);

// This encoding keeps each length and distance base plus its scale in five bits.
fn ShortCode(Value: type, HighBits: type, HighLog2: type, len_special: bool) type {
    return packed struct(u5) {
        high_bits: HighBits,
        high_log2: HighLog2,

        pub fn fromVal(v: Value) @This() {
            if (len_special and v == 255) return .fromInt(28);
            const high_bits = @bitSizeOf(HighBits) + 1;
            const bits = @bitSizeOf(Value) - @clz(v);
            if (bits <= high_bits) return @bitCast(@as(u5, @intCast(v)));
            const high = v >> @intCast(bits - high_bits);
            return .{ .high_bits = @truncate(high), .high_log2 = @intCast(bits - high_bits + 1) };
        }

        pub fn base(c: @This()) Value {
            if (len_special and c.toInt() == 28) return 255;
            if (c.high_log2 <= 1) return @as(u5, @bitCast(c));
            const high_value = (@as(Value, @intFromBool(c.high_log2 != 0)) << @bitSizeOf(HighBits)) | c.high_bits;
            const high_start = @as(std.math.Log2Int(Value), c.high_log2 - 1);
            return @shlExact(high_value, high_start);
        }

        const max_extra = @bitSizeOf(Value) - (1 + @bitSizeOf(HighLog2));
        pub fn extraBits(c: @This()) std.math.IntFittingRange(0, max_extra) {
            if (len_special and c.toInt() == 28) return 0;
            return @intCast(c.high_log2 -| 1);
        }

        pub fn toInt(c: @This()) u5 {
            return @bitCast(c);
        }

        pub fn fromInt(x: u5) @This() {
            return @bitCast(x);
        }
    };
}

const LookupLenCode = packed struct(u5) {
    code: ShortLenCode,

    const code_table = table: {
        var codes: [256]ShortLenCode = undefined;
        for (0.., &codes) |v, *c| {
            c.* = .fromVal(v);
        }
        break :table codes;
    };

    const base_table = table: {
        var bases: [29]u8 = undefined;
        for (0.., &bases) |c, *b| {
            b.* = ShortLenCode.fromInt(c).base();
        }
        break :table bases;
    };

    pub fn fromVal(v: u8) LookupLenCode {
        return .{ .code = code_table[v] };
    }

    pub fn base(c: LookupLenCode) u8 {
        return base_table[c.toInt()];
    }

    pub fn extraBits(c: LookupLenCode) u3 {
        return c.code.extraBits();
    }

    pub fn toInt(c: LookupLenCode) u5 {
        return @bitCast(c);
    }

    pub fn fromInt(x: u5) LookupLenCode {
        return @bitCast(x);
    }
};

const LookupDistCode = packed struct(u5) {
    code: ShortDistCode,

    const base_table = table: {
        var bases: [30]u15 = undefined;
        for (0.., &bases) |c, *b| {
            b.* = ShortDistCode.fromInt(c).base();
        }
        break :table bases;
    };

    pub fn fromVal(v: u15) LookupDistCode {
        return .{ .code = .fromVal(v) };
    }

    pub fn base(c: LookupDistCode) u15 {
        return base_table[c.toInt()];
    }

    pub fn extraBits(c: LookupDistCode) u4 {
        return c.code.extraBits();
    }

    pub fn toInt(c: LookupDistCode) u5 {
        return @bitCast(c);
    }

    pub fn fromInt(x: u5) LookupDistCode {
        return @bitCast(x);
    }
};

// --- Stream decoder ---

const Decompress = @This();

const DYNAMIC_FAST_INPUT_BYTES: usize = 16;

input: *Reader,
consumed_bits: u3,

reader: Reader,

container_metadata: Container.Metadata,

lit_dec: LiteralDecoder,
dst_dec: DistanceDecoder,

final_block: bool,
state: State,

err: ?Error,
fast_iterations: if (builtin.is_test) usize else void,

const BlockType = enum(u2) {
    stored = 0,
    fixed = 1,
    dynamic = 2,
    invalid = 3,
};

const State = union(enum) {
    protocol_header,
    block_header,
    stored_block: u16,
    fixed_block,
    fixed_block_literal: u8,
    fixed_block_match: u16,
    dynamic_block,
    dynamic_block_literal: u8,
    dynamic_block_match: u16,
    protocol_footer,
    end,
};

const DynamicFastResult = enum {
    input_boundary,
    output_boundary,
    end_block,
};

pub const Error = Container.Error || error{
    InvalidCode,
    InvalidMatch,
    WrongStoredBlockNlen,
    InvalidBlockType,
    InvalidDynamicBlockHeader,
    ReadFailed,
    OversubscribedHuffmanTree,
    IncompleteHuffmanTree,
    MissingEndOfBlockCode,
    EndOfStream,
};

const direct_vtable: Reader.VTable = .{
    .stream = streamDirect,
    .rebase = rebaseFallible,
    .discard = discardDirect,
    .readVec = readVec,
};

const indirect_vtable: Reader.VTable = .{
    .stream = streamIndirect,
    .rebase = rebaseFallible,
    .discard = discardIndirect,
    .readVec = readVec,
};

/// Initializes a decoder borrowing `input`; a nonempty history buffer must hold one DEFLATE window.
pub fn init(input: *Reader, container: Container, buffer: []u8) Decompress {
    if (buffer.len != 0) assert(buffer.len >= flate.max_window_len);
    return .{
        .reader = .{
            .vtable = if (buffer.len == 0) &direct_vtable else &indirect_vtable,
            .buffer = buffer,
            .seek = 0,
            .end = 0,
        },
        .input = input,
        .consumed_bits = 0,
        .container_metadata = .init(container),
        .lit_dec = .{},
        .dst_dec = .{},
        .final_block = false,
        .state = .protocol_header,
        .err = null,
        .fast_iterations = if (builtin.is_test) 0 else {},
    };
}

fn rebaseFallible(r: *Reader, capacity: usize) Reader.RebaseError!void {
    rebase(r, capacity);
}

fn rebase(r: *Reader, capacity: usize) void {
    assert(capacity <= r.buffer.len - flate.history_len);
    assert(r.end + capacity > r.buffer.len);
    const discard_n = @min(r.seek, r.end - flate.history_len);
    const keep = r.buffer[discard_n..r.end];
    @memmove(r.buffer[0..keep.len], keep);
    r.end = keep.len;
    r.seek -= discard_n;
}

fn discardDirect(r: *Reader, limit: std.Io.Limit) Reader.Error!usize {
    if (r.end + flate.history_len > r.buffer.len) rebase(r, flate.history_len);
    var writer: Writer = .{
        .vtable = &.{
            .drain = std.Io.Writer.Discarding.drain,
            .sendFile = std.Io.Writer.Discarding.sendFile,
        },
        .buffer = r.buffer,
        .end = r.end,
    };
    defer {
        assert(writer.end != 0);
        r.end = writer.end;
        r.seek = r.end;
    }
    const n = r.stream(&writer, limit) catch |err| switch (err) {
        error.WriteFailed => unreachable,
        error.ReadFailed => return error.ReadFailed,
        error.EndOfStream => return error.EndOfStream,
    };
    assert(n <= @intFromEnum(limit));
    return n;
}

fn discardIndirect(r: *Reader, limit: std.Io.Limit) Reader.Error!usize {
    const d: *Decompress = @alignCast(@fieldParentPtr("reader", r));
    if (r.end + flate.history_len > r.buffer.len) rebase(r, flate.history_len);
    var writer: Writer = .{
        .buffer = r.buffer,
        .end = r.end,
        .vtable = &.{ .drain = Writer.unreachableDrain },
    };
    {
        defer r.end = writer.end;
        _ = streamFallible(d, &writer, .limited(writer.buffer.len - writer.end)) catch |err| switch (err) {
            error.WriteFailed => unreachable,
            else => |e| return e,
        };
    }
    const n = limit.minInt(r.end - r.seek);
    r.seek += n;
    return n;
}

fn readVec(r: *Reader, data: [][]u8) Reader.Error!usize {
    _ = data;
    const d: *Decompress = @alignCast(@fieldParentPtr("reader", r));
    return streamIndirectInner(d);
}

fn streamIndirectInner(d: *Decompress) Reader.Error!usize {
    const r = &d.reader;
    if (r.buffer.len - r.end < flate.history_len) rebase(r, flate.history_len);
    var writer: Writer = .{
        .buffer = r.buffer,
        .end = r.end,
        .vtable = &.{
            .drain = Writer.unreachableDrain,
            .rebase = Writer.unreachableRebase,
        },
    };
    defer r.end = writer.end;
    _ = streamFallible(d, &writer, .limited(writer.buffer.len - writer.end)) catch |err| switch (err) {
        error.WriteFailed => unreachable,
        else => |e| return e,
    };
    return 0;
}

fn decodeLength(self: *Decompress, code_int: u5) !u16 {
    if (code_int > 28) return error.InvalidCode;
    const l: LenCode = .fromInt(code_int);
    const base = l.base();
    const extra = l.extraBits();
    return MIN_LENGTH + (base | try self.takeBits(extra));
}

fn decodeDistance(self: *Decompress, code_int: u5) !u16 {
    if (code_int > 29) return error.InvalidCode;
    const d: DistCode = .fromInt(code_int);
    const base = d.base();
    const extra = d.extraBits();
    return MIN_DISTANCE + (base | try self.takeBits(extra));
}

fn dynamicCodeLength(self: *Decompress, code: u16, lens: []u4, pos: usize) !usize {
    if (pos >= lens.len)
        return error.InvalidDynamicBlockHeader;

    switch (code) {
        0...15 => {
            lens[pos] = @intCast(code);
            return 1;
        },
        16 => {
            const n: u8 = @as(u8, try self.takeIntBits(u2)) + 3;
            if (pos == 0 or pos + n > lens.len)
                return error.InvalidDynamicBlockHeader;
            for (0..n) |i| {
                lens[pos + i] = lens[pos + i - 1];
            }
            return n;
        },
        17 => return @as(u8, try self.takeIntBits(u3)) + 3,
        18 => return @as(u8, try self.takeIntBits(u7)) + 11,
        else => return error.InvalidDynamicBlockHeader,
    }
}

fn streamDirect(r: *Reader, w: *Writer, limit: std.Io.Limit) Reader.StreamError!usize {
    const d: *Decompress = @alignCast(@fieldParentPtr("reader", r));
    return streamFallible(d, w, limit);
}

fn streamIndirect(r: *Reader, w: *Writer, limit: std.Io.Limit) Reader.StreamError!usize {
    const d: *Decompress = @alignCast(@fieldParentPtr("reader", r));
    _ = limit;
    _ = w;
    return streamIndirectInner(d);
}

fn streamFallible(d: *Decompress, w: *Writer, limit: std.Io.Limit) Reader.StreamError!usize {
    return streamInner(d, w, limit) catch |err| switch (err) {
        error.EndOfStream => {
            if (d.state == .end) {
                return error.EndOfStream;
            } else {
                d.err = error.EndOfStream;
                return error.ReadFailed;
            }
        },
        error.WriteFailed => return error.WriteFailed,
        else => |e| {
            // In the event of an error, state is unmodified so that it can be
            // better used to diagnose the failure.
            d.err = e;
            return error.ReadFailed;
        },
    };
}

fn streamInner(d: *Decompress, w: *Writer, limit: std.Io.Limit) (Error || Reader.StreamError)!usize {
    var remaining = @intFromEnum(limit);
    const in = d.input;
    sw: switch (d.state) {
        .protocol_header => switch (d.container_metadata.container()) {
            .gzip => {
                const Header = extern struct {
                    magic: u16 align(1),
                    method: u8,
                    flags: packed struct(u8) {
                        text: bool,
                        hcrc: bool,
                        extra: bool,
                        name: bool,
                        comment: bool,
                        reserved: u3,
                    },
                    mtime: u32 align(1),
                    xfl: u8,
                    os: u8,
                };
                const header = try in.takeStruct(Header, .little);
                if (header.magic != 0x8b1f or header.method != 0x08)
                    return error.BadGzipHeader;
                if (header.flags.extra) {
                    const extra_len = try in.takeInt(u16, .little);
                    try in.discardAll(extra_len);
                }
                if (header.flags.name) {
                    _ = try in.discardDelimiterInclusive(0);
                }
                if (header.flags.comment) {
                    _ = try in.discardDelimiterInclusive(0);
                }
                if (header.flags.hcrc) {
                    try in.discardAll(2);
                }
                continue :sw .block_header;
            },
            .zlib => {
                const header = try in.takeArray(2);
                const cmf: packed struct(u8) { cm: u4, cinfo: u4 } = @bitCast(header[0]);
                if (cmf.cm != 8 or cmf.cinfo > 7) return error.BadZlibHeader;
                continue :sw .block_header;
            },
            .raw => continue :sw .block_header,
        },
        .block_header => {
            d.final_block = (try d.takeIntBits(u1)) != 0;
            const block_type: BlockType = @enumFromInt(try d.takeIntBits(u2));
            switch (block_type) {
                .stored => {
                    d.alignBitsForward();
                    const len = try in.takeInt(u16, .little);
                    const nlen = try in.takeInt(u16, .little);
                    if (len != ~nlen) return error.WrongStoredBlockNlen;
                    continue :sw .{ .stored_block = len };
                },
                .fixed => continue :sw .fixed_block,
                .dynamic => {
                    const literal_count: u16 = @as(u16, try d.takeIntBits(u5)) + 257;
                    const distance_count: u16 = @as(u16, try d.takeIntBits(u5)) + 1;
                    const code_length_count: u8 = @as(u8, try d.takeIntBits(u4)) + 4;

                    if (literal_count > 286 or distance_count > 30)
                        return error.InvalidDynamicBlockHeader;

                    var code_length_alphabet: [19]u4 = @splat(0);
                    for (CODEGEN_ORDER[0..code_length_count]) |i| {
                        code_length_alphabet[i] = try d.takeIntBits(u3);
                    }
                    var code_length_decoder: CodegenDecoder = .{};
                    try code_length_decoder.generate(&code_length_alphabet);

                    var symbol_lengths: [286 + 30]u4 = @splat(0);
                    var length_index: usize = 0;
                    while (length_index < literal_count + distance_count) {
                        const peeked = try d.peekIntBitsShort(u7);
                        const sym = try code_length_decoder.find(peeked);
                        try d.tossBitsShort(sym.code_bits);
                        length_index += try d.dynamicCodeLength(
                            sym.value,
                            &symbol_lengths,
                            length_index,
                        );
                    }
                    if (length_index > literal_count + distance_count) {
                        return error.InvalidDynamicBlockHeader;
                    }

                    try d.lit_dec.generate(symbol_lengths[0..literal_count]);
                    try d.dst_dec.generate(symbol_lengths[literal_count..][0..distance_count]);

                    continue :sw .dynamic_block;
                },
                .invalid => return error.InvalidBlockType,
            }
        },
        .stored_block => |remaining_len| {
            const out: []u8 = if (remaining != 0)
                try w.writableSliceGreedyPreserve(flate.history_len, 1)
            else
                &.{};
            var limited_out: [1][]u8 = .{limit.min(.limited(remaining_len)).slice(out)};
            const n = try in.readVec(&limited_out);
            if (remaining_len - n == 0) {
                d.state = if (d.final_block) .protocol_footer else .block_header;
            } else {
                d.state = .{ .stored_block = @intCast(remaining_len - n) };
            }
            w.advance(n);
            return @intFromEnum(limit) - remaining + n;
        },
        .fixed_block => while (true) {
            const sym = try d.readFixedCode();

            if (sym >= 256) {
                @branchHint(.unlikely);

                if (sym == 256) {
                    @branchHint(.unlikely);
                    d.state = if (d.final_block) .protocol_footer else .block_header;
                    continue :sw d.state;
                }

                const length = try d.decodeLength(@intCast(sym - 257));
                continue :sw .{ .fixed_block_match = length };
            }

            const byte: u8 = @intCast(sym);
            if (remaining != 0) {
                @branchHint(.likely);
                remaining -= 1;
                try w.writeBytePreserve(flate.history_len, byte);
            } else {
                d.state = .{ .fixed_block_literal = byte };
                return @intFromEnum(limit) - remaining;
            }
        },
        .fixed_block_literal => |symbol| {
            assert(remaining != 0);
            remaining -= 1;
            try w.writeBytePreserve(flate.history_len, symbol);
            continue :sw .fixed_block;
        },
        .fixed_block_match => |length| {
            if (remaining >= length) {
                @branchHint(.likely);
                const distance = try d.decodeDistance(@bitReverse(try d.takeIntBits(u5)));
                try writeMatch(w, length, distance);
                remaining -= length;
                continue :sw .fixed_block;
            } else {
                d.state = .{ .fixed_block_match = length };
                return @intFromEnum(limit) - remaining;
            }
        },
        // Corpus profiling found only dynamic blocks, so this path owns the bounded fast loop.
        .dynamic_block => {
            const fast_result = try decodeDynamicFast(d, w, &remaining);
            if (fast_result == .end_block) {
                d.state = if (d.final_block) .protocol_footer else .block_header;
                continue :sw d.state;
            }

            while (true) {
                const token = try d.lit_dec.find(try d.peekIntBitsShort(u15));
                try d.tossBitsShort(token.code_bits);

                if (token.operation != .literal) {
                    @branchHint(.unlikely);

                    if (token.operation == .end) {
                        @branchHint(.unlikely);
                        d.state = if (d.final_block) .protocol_footer else .block_header;
                        continue :sw d.state;
                    }

                    const length = token.value + try d.takeBits(token.extraBits());
                    continue :sw .{ .dynamic_block_match = length };
                }

                const byte: u8 = @intCast(token.value);
                if (remaining != 0) {
                    @branchHint(.likely);
                    remaining -= 1;
                    try w.writeBytePreserve(flate.history_len, byte);
                    if (fast_result == .input_boundary) continue :sw .dynamic_block;
                } else {
                    d.state = .{ .dynamic_block_literal = byte };
                    return @intFromEnum(limit) - remaining;
                }
            }
        },
        .dynamic_block_literal => |symbol| {
            assert(remaining != 0);
            remaining -= 1;
            try w.writeBytePreserve(flate.history_len, symbol);
            continue :sw .dynamic_block;
        },
        .dynamic_block_match => |length| {
            if (remaining >= length) {
                @branchHint(.likely);
                remaining -= length;
                const token = try d.dst_dec.find(try d.peekIntBitsShort(u15));
                try d.tossBitsShort(token.code_bits);
                const distance = token.value + try d.takeBits(token.extra_bits);
                try writeMatch(w, length, distance);
                continue :sw .dynamic_block;
            } else {
                d.state = .{ .dynamic_block_match = length };
                return @intFromEnum(limit) - remaining;
            }
        },
        .protocol_footer => {
            d.alignBitsForward();
            switch (d.container_metadata) {
                .gzip => |*gzip| {
                    gzip.crc = try in.takeInt(u32, .little);
                    gzip.count = try in.takeInt(u32, .little);
                },
                .zlib => |*zlib| {
                    zlib.adler = try in.takeInt(u32, .big);
                },
                .raw => {},
            }
            d.state = .end;
            return @intFromEnum(limit) - remaining;
        },
        .end => return error.EndOfStream,
    }
}

// --- Dynamic block fast path ---

fn decodeDynamicFast(
    d: *Decompress,
    w: *Writer,
    remaining: *usize,
) Error!DynamicFastResult {
    if (remaining.* < MAX_LENGTH or w.buffer.len - w.end < MAX_LENGTH) {
        return .output_boundary;
    }

    const input = d.input.buffered();
    if (input.len < DYNAMIC_FAST_INPUT_BYTES) return .input_boundary;

    var input_pos: usize = @sizeOf(u64);
    var bit_buffer = std.mem.readInt(u64, input[0..@sizeOf(u64)], .little) >> d.consumed_bits;
    var bit_count: u7 = @as(u7, 64) - d.consumed_bits;
    var output_pos = w.end;
    var output_remaining = remaining.*;
    const lit_dec = &d.lit_dec;
    const dst_dec = &d.dst_dec;

    defer {
        const consumed_bit_count = input_pos * 8 - @as(usize, bit_count);
        d.input.toss(consumed_bit_count / 8);
        d.consumed_bits = @truncate(consumed_bit_count);
        w.end = output_pos;
        remaining.* = output_remaining;
    }

    while (output_remaining >= MAX_LENGTH and
        w.buffer.len - output_pos >= MAX_LENGTH and
        input.len - input_pos >= 8)
    {
        if (builtin.is_test) d.fast_iterations += 1;
        refillDynamicBits(input, &input_pos, &bit_buffer, &bit_count, 15);
        const literal_token = try lit_dec.find(@truncate(bit_buffer));
        bit_buffer >>= @intCast(literal_token.code_bits);
        bit_count -= literal_token.code_bits;

        if (literal_token.operation == .literal) {
            w.buffer[output_pos] = @intCast(literal_token.value);
            output_pos += 1;
            output_remaining -= 1;
            continue;
        }
        if (literal_token.operation == .end) return .end_block;

        const length = literal_token.value + takeDynamicBits(
            input,
            &input_pos,
            &bit_buffer,
            &bit_count,
            literal_token.extraBits(),
        );

        refillDynamicBits(input, &input_pos, &bit_buffer, &bit_count, 15);
        const distance_token = try dst_dec.find(@truncate(bit_buffer));
        bit_buffer >>= @intCast(distance_token.code_bits);
        bit_count -= distance_token.code_bits;

        const distance = distance_token.value + takeDynamicBits(
            input,
            &input_pos,
            &bit_buffer,
            &bit_count,
            distance_token.extra_bits,
        );
        if (output_pos < distance) return error.InvalidMatch;

        copyMatch(w.buffer, output_pos, length, distance);
        output_pos += length;
        output_remaining -= length;
    }

    return if (output_remaining < MAX_LENGTH or
        w.buffer.len - output_pos < MAX_LENGTH)
        .output_boundary
    else
        .input_boundary;
}

inline fn refillDynamicBits(
    input: []const u8,
    input_pos: *usize,
    bit_buffer: *u64,
    bit_count: *u7,
    needed: u4,
) void {
    if (bit_count.* >= needed) return;
    assert(input.len - input_pos.* >= @sizeOf(u32));
    const word = std.mem.readInt(u32, input[input_pos.*..][0..@sizeOf(u32)], .little);
    bit_buffer.* |= @as(u64, word) << @intCast(bit_count.*);
    bit_count.* += 32;
    input_pos.* += @sizeOf(u32);
}

inline fn takeDynamicBits(
    input: []const u8,
    input_pos: *usize,
    bit_buffer: *u64,
    bit_count: *u7,
    count: u4,
) u16 {
    refillDynamicBits(input, input_pos, bit_buffer, bit_count, count);
    const mask = @shlExact(@as(u16, 1), count) - 1;
    const value: u16 = @intCast(bit_buffer.* & mask);
    bit_buffer.* >>= @intCast(count);
    bit_count.* -= count;
    return value;
}

fn writeMatch(w: *Writer, length: u16, distance: u16) !void {
    if (w.end < distance) return error.InvalidMatch;
    assert(length >= MIN_LENGTH);
    assert(length <= MAX_LENGTH);
    assert(distance >= MIN_DISTANCE);
    assert(distance <= MAX_DISTANCE);

    // This is not a @memmove; it intentionally repeats patterns caused by
    // iterating one byte at a time.
    const dest = try w.writableSlicePreserve(flate.history_len, length);
    const end = dest.ptr - w.buffer.ptr;
    copyMatch(w.buffer, end, length, distance);
}

fn copyMatch(buffer: []u8, end: usize, length: u16, distance: u16) void {
    const dest = buffer[end..][0..length];
    const src = buffer[end - distance ..][0..length];
    if (distance >= length) {
        @memcpy(dest, src);
    } else if (distance == 1) {
        @memset(dest, src[0]);
    } else {
        for (dest, src) |*d, s| d.* = s;
    }
}

fn peekBits(d: *Decompress, n: u4) !u16 {
    const bits = d.input.peekInt(u32, .little) catch |e| return switch (e) {
        error.ReadFailed => error.ReadFailed,
        error.EndOfStream => d.peekBitsEnding(n),
    };
    const mask = @shlExact(@as(u16, 1), n) - 1;
    return @intCast((bits >> d.consumed_bits) & mask);
}

fn peekBitsEnding(d: *Decompress, n: u4) !u16 {
    @branchHint(.unlikely);

    const left = d.input.buffered();
    if (left.len * 8 < @as(usize, n) + d.consumed_bits) return error.EndOfStream;
    const bits = std.mem.readVarInt(u32, left, .little);
    const mask = @shlExact(@as(u16, 1), n) - 1;
    return @intCast((bits >> d.consumed_bits) & mask);
}

// Callers already proved these bits are buffered, so advancing cannot refill or fail.
fn tossBits(d: *Decompress, n: u4) void {
    d.input.toss((@as(u8, n) + d.consumed_bits) / 8);
    d.consumed_bits +%= @truncate(n);
}

fn takeBits(d: *Decompress, n: u4) !u16 {
    const bits = try d.peekBits(n);
    d.tossBits(n);
    return bits;
}

fn alignBitsForward(d: *Decompress) void {
    d.input.toss(@intFromBool(d.consumed_bits != 0));
    d.consumed_bits = 0;
}

fn peekBitsShort(d: *Decompress, n: u4) !u16 {
    const bits = d.input.peekInt(u32, .little) catch |e| return switch (e) {
        error.ReadFailed => error.ReadFailed,
        error.EndOfStream => d.peekBitsShortEnding(n),
    };
    const mask = @shlExact(@as(u16, 1), n) - 1;
    return @intCast((bits >> d.consumed_bits) & mask);
}

fn peekBitsShortEnding(d: *Decompress, n: u4) !u16 {
    @branchHint(.unlikely);

    const left = d.input.buffered();
    const bits = std.mem.readVarInt(u32, left, .little);
    const mask = @shlExact(@as(u16, 1), n) - 1;
    return @intCast((bits >> d.consumed_bits) & mask);
}

fn tossBitsShort(d: *Decompress, n: u4) !void {
    if (d.input.bufferedLen() * 8 < @as(usize, n) + d.consumed_bits) {
        return error.EndOfStream;
    }
    d.tossBits(n);
}

fn takeIntBits(d: *Decompress, T: type) !T {
    return @intCast(try d.takeBits(@bitSizeOf(T)));
}

fn peekIntBitsShort(d: *Decompress, T: type) !T {
    return @intCast(try d.peekBitsShort(@bitSizeOf(T)));
}

// Fixed literals share seven prefix bits and need at most two more bits to select a symbol.
fn readFixedCode(d: *Decompress) !u16 {
    const code7 = @bitReverse(try d.takeIntBits(u7));
    return switch (code7) {
        0...0b0010_111 => @as(u16, code7) + 256,
        0b0010_111 + 1...0b1011_111 => (@as(u16, code7) << 1) + @as(u16, try d.takeIntBits(u1)) - 0b0011_0000,
        0b1011_111 + 1...0b1100_011 => (@as(u16, code7 - 0b1100000) << 1) + try d.takeIntBits(u1) + 280,
        else => (@as(u16, code7 - 0b1100_100) << 2) + @as(u16, @bitReverse(try d.takeIntBits(u2))) + 144,
    };
}

// --- Huffman tables ---

const Symbol = packed struct(u16) {
    value: u12 = 0,
    code_bits: u4 = 0,

    fn fromSymbol(value: u16, code_bits: u4) @This() {
        return .{ .value = @intCast(value), .code_bits = code_bits };
    }

    fn invalid() @This() {
        return .{ .value = 0xfff };
    }

    fn isInvalid(self: @This()) bool {
        return self.value == 0xfff;
    }

    fn secondary(group: u16) @This() {
        return .{ .value = @intCast(group) };
    }

    fn secondaryGroup(self: @This()) u16 {
        return self.value;
    }
};

const LiteralOperation = enum(u3) {
    literal,
    end,
    length_0,
    length_1,
    length_2,
    length_3,
    length_4,
    length_5,
};

const LiteralToken = packed struct(u16) {
    value: u9 = 0,
    code_bits: u4 = 0,
    operation: LiteralOperation = .literal,

    fn fromSymbol(value: u16, code_bits: u4) @This() {
        if (value < 256) return .{ .value = @intCast(value), .code_bits = code_bits };
        if (value == 256) return .{ .value = 256, .code_bits = code_bits, .operation = .end };
        const length_code: LenCode = .fromInt(@intCast(value - 257));
        return .{
            .value = @intCast(@as(u16, MIN_LENGTH) + length_code.base()),
            .code_bits = code_bits,
            .operation = @enumFromInt(2 + length_code.extraBits()),
        };
    }

    fn invalid() @This() {
        return .{ .value = std.math.maxInt(u9) };
    }

    fn isInvalid(self: @This()) bool {
        return self.value == std.math.maxInt(u9);
    }

    fn secondary(group: u16) @This() {
        return .{ .value = @intCast(group) };
    }

    fn secondaryGroup(self: @This()) u16 {
        return self.value;
    }

    fn extraBits(self: @This()) u3 {
        return @intCast(@intFromEnum(self.operation) - 2);
    }
};

const DistanceToken = packed struct(u32) {
    value: u15 = 0,
    code_bits: u4 = 0,
    extra_bits: u4 = 0,
    padding: u9 = 0,

    fn fromSymbol(value: u16, code_bits: u4) @This() {
        const distance_code: DistCode = .fromInt(@intCast(value));
        return .{
            .value = MIN_DISTANCE + distance_code.base(),
            .code_bits = code_bits,
            .extra_bits = distance_code.extraBits(),
        };
    }

    fn invalid() @This() {
        return .{ .value = std.math.maxInt(u15) };
    }

    fn isInvalid(self: @This()) bool {
        return self.value == std.math.maxInt(u15);
    }

    fn secondary(group: u16) @This() {
        return .{ .value = @intCast(group) };
    }

    fn secondaryGroup(self: @This()) u16 {
        return self.value;
    }
};

const LiteralDecoder = HuffmanDecoder(286, 15, 11, LiteralToken);
const DistanceDecoder = HuffmanDecoder(30, 15, 9, DistanceToken);
const CodegenDecoder = HuffmanDecoder(19, 7, 7, Symbol);

// Two levels avoid a complete 15-bit table while keeping misses to one bounded lookup.
fn HuffmanDecoder(
    comptime alphabet_size: u16,
    comptime max_code_bits: u4,
    comptime lookup_bits: u4,
    comptime Entry: type,
) type {
    const lookup_shift = max_code_bits - lookup_bits;
    const lookup_mask = (1 << lookup_bits) - 1;
    const secondary_bits = max_code_bits - lookup_bits;
    const secondary_group_size = @as(usize, 1) << secondary_bits;

    return struct {
        lookup: [1 << lookup_bits]Entry = undefined,
        secondary: if (lookup_bits == max_code_bits)
            void
        else
            [@as(usize, alphabet_size) * secondary_group_size]Entry = undefined,

        const Self = @This();

        fn reverseIdx(idx: usize) u16 {
            return @bitReverse(@as(@Int(.unsigned, lookup_bits), @intCast(idx)));
        }

        pub fn generate(self: *Self, lens: []const u4) !void {
            try checkCompleteness(lens);
            self.lookup = @splat(Entry.invalid());

            var buckets: [1 + @as(usize, max_code_bits)][alphabet_size]Entry = undefined;
            var bucket_len: [buckets.len]u16 = @splat(0);
            for (0.., lens) |symbol, bits| {
                buckets[bits][bucket_len[bits]] = Entry.fromSymbol(@intCast(symbol), bits);
                bucket_len[bits] += 1;
            }

            var code: u16 = 0;
            var idx: u16 = 0;
            var secondary_group_count: u16 = 0;
            for (1..lookup_bits + 1) |bits| {
                const inc = @as(u16, 1) << @intCast(max_code_bits - bits);
                for (buckets[bits][0..bucket_len[bits]]) |lookup_sym| {
                    const next_code = code + inc;
                    const next_idx = next_code >> lookup_shift;
                    for (idx..next_idx) |i| {
                        self.lookup[reverseIdx(i)] = lookup_sym;
                    }
                    code = next_code;
                    idx = next_idx;
                }
            }
            for (lookup_bits + 1..buckets.len) |bits| {
                const inc = @as(u16, 1) << @intCast(max_code_bits - bits);
                for (buckets[bits][0..bucket_len[bits]]) |secondary_sym| {
                    const next_code = code + inc;
                    const reversed = @bitReverse(
                        @as(@Int(.unsigned, max_code_bits), @intCast(code)),
                    );
                    const primary = &self.lookup[reversed & lookup_mask];
                    const group = if (primary.isInvalid()) group: {
                        const group_index = secondary_group_count;
                        secondary_group_count += 1;
                        primary.* = Entry.secondary(group_index);
                        const start = @as(usize, group_index) * secondary_group_size;
                        @memset(
                            self.secondary[start..][0..secondary_group_size],
                            Entry.invalid(),
                        );
                        break :group group_index;
                    } else primary.secondaryGroup();

                    const suffix_bits: u4 = @intCast(bits - lookup_bits);
                    const suffix_mask = (@as(u16, 1) << suffix_bits) - 1;
                    const suffix = (reversed >> lookup_bits) & suffix_mask;
                    const start = @as(usize, group) * secondary_group_size;
                    for (0..@as(usize, 1) << (secondary_bits - suffix_bits)) |repeat| {
                        const secondary_index = suffix + (repeat << suffix_bits);
                        self.secondary[start + secondary_index] = secondary_sym;
                    }

                    code = next_code;
                }
            }
        }

        // This is puff's canonical-code left-space completeness check.
        fn checkCompleteness(lens: []const u4) !void {
            if (alphabet_size == 286)
                if (lens[256] == 0) return error.MissingEndOfBlockCode;

            var count = [_]u16{0} ** (@as(usize, max_code_bits) + 1);
            var max: usize = 0;
            for (lens) |n| {
                if (n == 0) continue;
                if (n > max) max = n;
                count[n] += 1;
            }
            if (max == 0) return;

            var left: usize = 1;
            for (1..count.len) |len| {
                left <<= 1;
                if (count[len] > left)
                    return error.OversubscribedHuffmanTree;
                left -= count[len];
            }
            if (left > 0) {
                // RFC 1951 permits an incomplete tree only for one single-bit code.
                if (max_code_bits > 7 and max == count[0] + count[1]) return;
                return error.IncompleteHuffmanTree;
            }
        }

        pub fn find(self: *Self, code: u16) !Entry {
            const idx = code & lookup_mask;
            const sym = self.lookup[idx];
            if (sym.code_bits != 0) return sym;
            if (sym.isInvalid() or lookup_bits == max_code_bits) return error.InvalidCode;
            const suffix = (code >> lookup_bits) & (secondary_group_size - 1);
            const secondary_idx = @as(usize, sym.secondaryGroup()) * secondary_group_size + suffix;
            const secondary_sym = self.secondary[secondary_idx];
            if (secondary_sym.code_bits == 0) return error.InvalidCode;
            return secondary_sym;
        }
    };
}

fn testFailure(container: Container, in: []const u8, expected_err: anyerror) !void {
    var reader: Reader = .fixed(in);
    var aw: Writer.Allocating = .init(testing.allocator);
    defer aw.deinit();

    var decompress: Decompress = .init(&reader, container, &.{});
    try testing.expectError(error.ReadFailed, decompress.reader.streamRemaining(&aw.writer));
    try testing.expectEqual(expected_err, decompress.err orelse return error.TestFailed);
}

fn testDecompress(container: Container, compressed: []const u8, expected_plain: []const u8) !void {
    var in: std.Io.Reader = .fixed(compressed);
    var aw: std.Io.Writer.Allocating = .init(testing.allocator);
    defer aw.deinit();

    var decompress: Decompress = .init(&in, container, &.{});
    const decompressed_len = try decompress.reader.streamRemaining(&aw.writer);
    try testing.expectEqual(expected_plain.len, decompressed_len);
    try testing.expectEqualSlices(u8, expected_plain, aw.written());
}

const DYNAMIC_FAST_FRAGMENTED = [_]u8{
    0xed, 0xce, 0x67, 0x02, 0x81, 0x00, 0x18, 0x00,
    0xd0, 0xab, 0x64, 0x24, 0x23, 0x32, 0x33, 0x22,
    0x91, 0x91, 0x8c, 0xec, 0x44, 0x46, 0x22, 0x91,
    0x99, 0x9d, 0x71, 0x76, 0xa7, 0xf0, 0xef, 0x7b,
    0x27, 0x78, 0x48, 0x71, 0x13, 0xec, 0x59, 0x8c,
    0x86, 0xb7, 0x6e, 0x69, 0x05, 0xab, 0x9b, 0x89,
    0xb1, 0x83, 0xdb, 0x45, 0xa5, 0x6f, 0x61, 0x4d,
    0x74, 0x2d, 0x66, 0xe9, 0x6b, 0x5e, 0x29, 0x05,
    0xab, 0x9d, 0xe2, 0xb2, 0x9d, 0xdb, 0x45, 0x06,
    0x1f, 0x56, 0x27, 0xba, 0xcf, 0xec, 0xc2, 0x2b,
    0x5c, 0xa9, 0x99, 0xab, 0x7a, 0x24, 0x65, 0x7b,
    0x79, 0x1b, 0x16, 0xdf, 0xac, 0x1e, 0xe8, 0x3c,
    0x68, 0xd5, 0x2b, 0x5c, 0x52, 0x53, 0x94, 0x3f,
    0x92, 0x23, 0x5b, 0xc9, 0x08, 0x89, 0xef, 0xfc,
    0xca, 0xdf, 0xbe, 0xd3, 0xaa, 0xa7, 0x71, 0x4e,
    0x4e, 0x50, 0xfe, 0x10, 0x1b, 0x22, 0x45, 0x23,
    0xd4, 0x7f, 0xe5, 0x34, 0xbc, 0x7d, 0xcf, 0xcc,
    0xdd, 0x75, 0x33, 0x39, 0x71, 0x56, 0xf6, 0x51,
    0x09, 0x81, 0x04, 0x24, 0x20, 0x01, 0x09, 0x48,
    0x40, 0x02, 0x12, 0x90, 0x80, 0x04, 0x24, 0x20,
    0x01, 0x09, 0x48, 0x40, 0x02, 0x12, 0xff, 0x4c,
    0xfc, 0x00,
};

// --- Huffman verification ---

test "[property] - [inflate length codes]: match RFC 1951 tables" {
    inline for ([_]type{ ShortLenCode, LookupLenCode }) |Code| {
        for (0.., [_]struct {
            base: u8,
            extra_bits: u4,
        }{
            // zig fmt: off
            .{ .base = 3   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 4   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 5   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 6   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 7   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 8   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 9   - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 10  - MIN_LENGTH, .extra_bits = 0 },
            .{ .base = 11  - MIN_LENGTH, .extra_bits = 1 },
            .{ .base = 13  - MIN_LENGTH, .extra_bits = 1 },
            .{ .base = 15  - MIN_LENGTH, .extra_bits = 1 },
            .{ .base = 17  - MIN_LENGTH, .extra_bits = 1 },
            .{ .base = 19  - MIN_LENGTH, .extra_bits = 2 },
            .{ .base = 23  - MIN_LENGTH, .extra_bits = 2 },
            .{ .base = 27  - MIN_LENGTH, .extra_bits = 2 },
            .{ .base = 31  - MIN_LENGTH, .extra_bits = 2 },
            .{ .base = 35  - MIN_LENGTH, .extra_bits = 3 },
            .{ .base = 43  - MIN_LENGTH, .extra_bits = 3 },
            .{ .base = 51  - MIN_LENGTH, .extra_bits = 3 },
            .{ .base = 59  - MIN_LENGTH, .extra_bits = 3 },
            .{ .base = 67  - MIN_LENGTH, .extra_bits = 4 },
            .{ .base = 83  - MIN_LENGTH, .extra_bits = 4 },
            .{ .base = 99  - MIN_LENGTH, .extra_bits = 4 },
            .{ .base = 115 - MIN_LENGTH, .extra_bits = 4 },
            .{ .base = 131 - MIN_LENGTH, .extra_bits = 5 },
            .{ .base = 163 - MIN_LENGTH, .extra_bits = 5 },
            .{ .base = 195 - MIN_LENGTH, .extra_bits = 5 },
            .{ .base = 227 - MIN_LENGTH, .extra_bits = 5 },
            .{ .base = 258 - MIN_LENGTH, .extra_bits = 0 },
        }) |code, params| {
            // zig fmt: on
            const c: u5 = @intCast(code);
            try std.testing.expectEqual(params.extra_bits, Code.extraBits(.fromInt(@intCast(c))));
            try std.testing.expectEqual(params.base, Code.base(.fromInt(@intCast(c))));
            for (params.base..params.base + @shlExact(@as(u16, 1), params.extra_bits) -
                @intFromBool(c == 27)) |v|
            {
                try std.testing.expectEqual(c, Code.fromVal(@intCast(v)).toInt());
            }
        }
    }
}

test "[property] - [inflate distance codes]: match RFC 1951 tables" {
    inline for ([_]type{ ShortDistCode, LookupDistCode }) |Code| {
        for (0.., [_]struct {
            base: u15,
            extra_bits: u4,
        }{
            // zig fmt: off
            .{ .base = 1     - MIN_DISTANCE, .extra_bits =  0 },
            .{ .base = 2     - MIN_DISTANCE, .extra_bits =  0 },
            .{ .base = 3     - MIN_DISTANCE, .extra_bits =  0 },
            .{ .base = 4     - MIN_DISTANCE, .extra_bits =  0 },
            .{ .base = 5     - MIN_DISTANCE, .extra_bits =  1 },
            .{ .base = 7     - MIN_DISTANCE, .extra_bits =  1 },
            .{ .base = 9     - MIN_DISTANCE, .extra_bits =  2 },
            .{ .base = 13    - MIN_DISTANCE, .extra_bits =  2 },
            .{ .base = 17    - MIN_DISTANCE, .extra_bits =  3 },
            .{ .base = 25    - MIN_DISTANCE, .extra_bits =  3 },
            .{ .base = 33    - MIN_DISTANCE, .extra_bits =  4 },
            .{ .base = 49    - MIN_DISTANCE, .extra_bits =  4 },
            .{ .base = 65    - MIN_DISTANCE, .extra_bits =  5 },
            .{ .base = 97    - MIN_DISTANCE, .extra_bits =  5 },
            .{ .base = 129   - MIN_DISTANCE, .extra_bits =  6 },
            .{ .base = 193   - MIN_DISTANCE, .extra_bits =  6 },
            .{ .base = 257   - MIN_DISTANCE, .extra_bits =  7 },
            .{ .base = 385   - MIN_DISTANCE, .extra_bits =  7 },
            .{ .base = 513   - MIN_DISTANCE, .extra_bits =  8 },
            .{ .base = 769   - MIN_DISTANCE, .extra_bits =  8 },
            .{ .base = 1025  - MIN_DISTANCE, .extra_bits =  9 },
            .{ .base = 1537  - MIN_DISTANCE, .extra_bits =  9 },
            .{ .base = 2049  - MIN_DISTANCE, .extra_bits = 10 },
            .{ .base = 3073  - MIN_DISTANCE, .extra_bits = 10 },
            .{ .base = 4097  - MIN_DISTANCE, .extra_bits = 11 },
            .{ .base = 6145  - MIN_DISTANCE, .extra_bits = 11 },
            .{ .base = 8193  - MIN_DISTANCE, .extra_bits = 12 },
            .{ .base = 12289 - MIN_DISTANCE, .extra_bits = 12 },
            .{ .base = 16385 - MIN_DISTANCE, .extra_bits = 13 },
            .{ .base = 24577 - MIN_DISTANCE, .extra_bits = 13 },
        }) |code, params| {
            // zig fmt: on
            const c: u5 = @intCast(code);
            try std.testing.expectEqual(params.extra_bits, Code.extraBits(.fromInt(@intCast(c))));
            try std.testing.expectEqual(params.base, Code.base(.fromInt(@intCast(c))));
            for (params.base..params.base + @shlExact(@as(u16, 1), params.extra_bits)) |v| {
                try std.testing.expectEqual(c, Code.fromVal(@intCast(v)).toInt());
            }
        }
    }
}

test "[property] - [inflate token metadata]: matches RFC 1951 symbols" {
    try testing.expectEqual(@as(usize, 2), @sizeOf(LiteralToken));
    try testing.expectEqual(@as(usize, 4), @sizeOf(DistanceToken));

    const literal = LiteralToken.fromSymbol('A', 8);
    try testing.expectEqual(LiteralOperation.literal, literal.operation);
    try testing.expectEqual(@as(u9, 'A'), literal.value);
    try testing.expectEqual(@as(u4, 8), literal.code_bits);

    const end = LiteralToken.fromSymbol(256, 7);
    try testing.expectEqual(LiteralOperation.end, end.operation);
    try testing.expectEqual(@as(u9, 256), end.value);
    try testing.expectEqual(@as(u4, 7), end.code_bits);

    for (0..29) |code| {
        const length_code: LenCode = .fromInt(@intCast(code));
        const token = LiteralToken.fromSymbol(@intCast(257 + code), 15);
        try testing.expectEqual(
            @as(u16, MIN_LENGTH) + length_code.base(),
            @as(u16, token.value),
        );
        try testing.expectEqual(length_code.extraBits(), token.extraBits());
    }

    for (0..30) |code| {
        const distance_code: DistCode = .fromInt(@intCast(code));
        const token = DistanceToken.fromSymbol(@intCast(code), 15);
        try testing.expectEqual(MIN_DISTANCE + distance_code.base(), token.value);
        try testing.expectEqual(distance_code.extraBits(), token.extra_bits);
    }
}

test "[property] - [inflate Huffman table]: resolves every canonical prefix" {
    const code_lens = [_]u4{ 4, 3, 0, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 2 };
    var h: CodegenDecoder = .{};
    try h.generate(&code_lens);

    for (0b0000_000..0b0100_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 3, .code_bits = 2 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );

    for (0b0100_000..0b1000_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 18, .code_bits = 2 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );

    for (0b1000_000..0b1010_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 1, .code_bits = 3 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );

    for (0b1010_000..0b1100_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 4, .code_bits = 3 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );

    for (0b1100_000..0b1110_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 17, .code_bits = 3 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );

    for (0b1110_000..0b1111_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 0, .code_bits = 4 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );

    for (0b1111_000..0b1_0000_000) |c|
        try testing.expectEqual(
            Symbol{ .value = 16, .code_bits = 4 },
            try h.find(@bitReverse(@as(u7, @intCast(c)))),
        );
}

test "[property] - [inflate Huffman table]: decodes RFC 1951 literal codes" {
    const max_bits = 5;
    var decoder: HuffmanDecoder(16, max_bits, 3, Symbol) = .{};
    try decoder.generate(&.{ 3, 3, 3, 3, 0, 0, 3, 2, 4, 4 });

    inline for (0.., .{
        @as(u3, 0b010),
        @as(u3, 0b011),
        @as(u3, 0b100),
        @as(u3, 0b101),
        @as(u0, 0),
        @as(u0, 0),
        @as(u3, 0b110),
        @as(u2, 0b00),
        @as(u4, 0b1110),
        @as(u4, 0b1111),
    }) |i, code| {
        const bits = @bitSizeOf(@TypeOf(code));
        if (bits == 0) continue;
        for (0..1 << (max_bits - bits)) |extra| {
            const full = (@as(u16, code) << (max_bits - bits)) | @as(u16, @intCast(extra));
            const symbol = try decoder.find(@bitReverse(@as(u5, @intCast(full))));
            try testing.expectEqual(i, symbol.value);
            try testing.expectEqual(bits, symbol.code_bits);
        }
    }
}

test "[failure] - [inflate Huffman table]: rejects an unused incomplete prefix" {
    var decoder: HuffmanDecoder(2, 15, 11, Symbol) = .{};
    try decoder.generate(&.{ 1, 0 });

    try testing.expectEqual(
        Symbol{ .value = 0, .code_bits = 1 },
        try decoder.find(0),
    );
    try testing.expectError(error.InvalidCode, decoder.find(1));
}

// --- Stream verification ---

test "[unit] - [inflate raw stream]: decodes a stored block" {
    try testDecompress(.raw, &[_]u8{
        0b0000_0001, 0b0000_1100, 0x00, 0b1111_0011, 0xff,
        'H',         'e',         'l',  'l',         'o',
        ' ',         'w',         'o',  'r',         'l',
        'd',         0x0a,
    }, "Hello world\n");
}

test "[unit] - [inflate raw stream]: decodes a fixed block" {
    try testDecompress(.raw, &[_]u8{
        0xf3, 0x48, 0xcd, 0xc9, 0xc9, 0x57, 0x28, 0xcf,
        0x2f, 0xca, 0x49, 0xe1, 0x02, 0x00,
    }, "Hello world\n");
}

test "[unit] - [inflate raw stream]: decodes a dynamic block" {
    try testDecompress(.raw, &[_]u8{
        0x3d, 0xc6, 0x39, 0x11, 0x00, 0x00, 0x0c, 0x02,
        0x30, 0x2b, 0xb5, 0x52, 0x1e, 0xff, 0x96, 0x38,
        0x16, 0x96, 0x5c, 0x1e, 0x94, 0xcb, 0x6d, 0x01,
    }, "ABCDEABCD ABCDEABCD");
}

test "[unit] - [dynamic fast loop]: decodes through the bounded path" {
    const compressed = [_]u8{
        0x3d, 0xc6, 0x39, 0x11, 0x00, 0x00, 0x0c, 0x02,
        0x30, 0x2b, 0xb5, 0x52, 0x1e, 0xff, 0x96, 0x38,
        0x16, 0x96, 0x5c, 0x1e, 0x94, 0xcb, 0x6d, 0x01,
    } ++ [_]u8{0} ** 64;
    var input: Reader = .fixed(&compressed);
    var decompressor: Decompress = .init(&input, .raw, &.{});
    var output: [512]u8 = undefined;
    var writer: Writer = .fixed(&output);

    const written = try decompressor.reader.streamRemaining(&writer);

    try testing.expect(decompressor.fast_iterations > 0);
    try testing.expectEqualStrings("ABCDEABCD ABCDEABCD", writer.buffered());
    try testing.expectEqual(@as(usize, 19), written);
    try testing.expectEqual(@as(usize, 64), input.bufferedLen());
}

test "[unit] - [dynamic fast loop]: resumes across fragmented input" {
    var input_buffer: [64]u8 = undefined;
    var input = testing.Reader.init(&input_buffer, &.{.{ .buffer = &DYNAMIC_FAST_FRAGMENTED }});
    input.artificial_limit = .limited(40);
    var decompressor: Decompress = .init(&input.interface, .raw, &.{});
    var output: [4096]u8 = undefined;
    var writer: Writer = .fixed(&output);

    const written = try decompressor.reader.streamRemaining(&writer);

    try testing.expect(decompressor.fast_iterations > 0);
    try testing.expectEqual(output.len, written);
    for (writer.buffered(), 0..) |byte, index| {
        try testing.expectEqual(@as(u8, @intCast(32 + (index * 37 + index / 7) % 95)), byte);
    }
}

test "[property] - [inflate dynamic stream]: every truncation fails without a panic" {
    for (0..DYNAMIC_FAST_FRAGMENTED.len) |end| {
        var input: Reader = .fixed(DYNAMIC_FAST_FRAGMENTED[0..end]);
        var decompressor: Decompress = .init(&input, .raw, &.{});
        var output: [4096]u8 = undefined;
        var writer: Writer = .fixed(&output);

        try testing.expectError(error.ReadFailed, decompressor.reader.streamRemaining(&writer));
    }
}

test "[unit] - [inflate gzip stream]: decodes a stored block" {
    try testDecompress(.gzip, &[_]u8{
        0x1f,        0x8b,        0x08, 0x00,        0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0b0000_0001, 0b0000_1100, 0x00, 0b1111_0011, 0xff, 'H',  'e',  'l',  'l',  'o',
        ' ',         'w',         'o',  'r',         'l',  'd',  0x0a, 0xd5, 0xe0, 0x39,
        0xb7,        0x0c,        0x00, 0x00,        0x00,
    }, "Hello world\n");
}

test "[unit] - [inflate gzip stream]: decodes a fixed block" {
    try testDecompress(.gzip, &[_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x03,
        0xf3, 0x48, 0xcd, 0xc9, 0xc9, 0x57, 0x28, 0xcf, 0x2f, 0xca,
        0x49, 0xe1, 0x02, 0x00, 0xd5, 0xe0, 0x39, 0xb7, 0x0c, 0x00,
        0x00, 0x00,
    }, "Hello world\n");
}

test "[unit] - [inflate gzip stream]: decodes a dynamic block" {
    try testDecompress(.gzip, &[_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x3d, 0xc6, 0x39, 0x11, 0x00, 0x00, 0x0c, 0x02, 0x30, 0x2b,
        0xb5, 0x52, 0x1e, 0xff, 0x96, 0x38, 0x16, 0x96, 0x5c, 0x1e,
        0x94, 0xcb, 0x6d, 0x01, 0x17, 0x1c, 0x39, 0xb4, 0x13, 0x00,
        0x00, 0x00,
    }, "ABCDEABCD ABCDEABCD");
}

test "[unit] - [inflate gzip stream]: accepts a named member" {
    try testDecompress(.gzip, &[_]u8{
        0x1f, 0x8b, 0x08, 0x08, 0xe5, 0x70, 0xb1, 0x65, 0x00, 0x03, 0x68, 0x65, 0x6c, 0x6c, 0x6f, 0x2e,
        0x74, 0x78, 0x74, 0x00, 0xf3, 0x48, 0xcd, 0xc9, 0xc9, 0x57, 0x28, 0xcf, 0x2f, 0xca, 0x49, 0xe1,
        0x02, 0x00, 0xd5, 0xe0, 0x39, 0xb7, 0x0c, 0x00, 0x00, 0x00,
    }, "Hello world\n");
}

test "[unit] - [inflate zlib stream]: decodes a stored block" {
    try testDecompress(.zlib, &[_]u8{
        0x78,        0b10_0_11100,
        0b0000_0001, 0b0000_1100,
        0x00,        0b1111_0011,
        0xff,        'H',
        'e',         'l',
        'l',         'o',
        ' ',         'w',
        'o',         'r',
        'l',         'd',
        0x0a,        0x1c,
        0xf2,        0x04,
        0x47,
    }, "Hello world\n");
}

test "[failure] - [inflate raw stream]: rejects a reserved block type" {
    try testFailure(.raw, &[_]u8{0b110}, error.InvalidBlockType);
}

test "[edge] - [inflate raw stream]: accepts an empty destination" {
    const input = &[_]u8{
        0b0000_0001, 0b0000_1100, 0x00, 0b1111_0011, 0xff,
        'H',         'e',         'l',  'l',         'o',
        ' ',         'w',         'o',  'r',         'l',
        'd',         0x0a,
    };
    var in: Reader = .fixed(input);
    var decomp: Decompress = .init(&in, .raw, &.{});
    const r = &decomp.reader;
    var bufs: [1][]u8 = .{&.{}};
    try testing.expectEqual(0, try r.readVec(&bufs));
}

test "[failure] - [inflate zlib stream]: rejects malformed framing" {
    try testFailure(.zlib, &[_]u8{0x78}, error.EndOfStream);

    try testFailure(.zlib, &[_]u8{ 0x79, 0x94 }, error.BadZlibHeader);

    try testFailure(.zlib, &[_]u8{ 0x88, 0x98 }, error.BadZlibHeader);

    try testFailure(.zlib, &[_]u8{ 0x78, 0xda, 0x03, 0x00, 0x00 }, error.EndOfStream);
}

test "[failure] - [inflate gzip stream]: rejects malformed framing" {
    try testFailure(.gzip, &[_]u8{ 0x1f, 0x8B }, error.EndOfStream);

    try testFailure(.gzip, &[_]u8{
        0x1f, 0x8b, 0x09, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x03,
    }, error.BadGzipHeader);

    try testFailure(.gzip, &[_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x03, 0x03, 0x00, 0x00, 0x00, 0x00,
    }, error.EndOfStream);

    try testFailure(.gzip, &[_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x03, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00,
    }, error.EndOfStream);

    try testDecompress(.gzip, &[_]u8{
        0x1f, 0x8b, 0x08, 0x12, 0x00, 0x09, 0x6e, 0x88, 0x00, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x00,
        0x99, 0xd6, 0x01, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    }, "");
}

test "[edge] - [inflate zlib stream]: leaves trailing bytes unread" {
    const data = [_]u8{
        0x78, 0x9c, 0x73, 0xce, 0x2f, 0xa8, 0x2c, 0xca, 0x4c, 0xcf, 0x28, 0x51, 0x08, 0xcf, 0xcc, 0xc9,
        0x49, 0xcd, 0x55, 0x28, 0x4b, 0xcc, 0x53, 0x08, 0x4e, 0xce, 0x48, 0xcc, 0xcc, 0xd6, 0x51, 0x08,
        0xce, 0xcc, 0x4b, 0x4f, 0x2c, 0xc8, 0x2f, 0x4a, 0x55, 0x30, 0xb4, 0xb4, 0x34, 0xd5, 0xb5, 0x34,
        0x03, 0x00, 0x8b, 0x61, 0x0f, 0xa4, 0x52, 0x5a, 0x94, 0x12,
    };

    var reader: std.Io.Reader = .fixed(&data);

    var decompress_buffer: [flate.max_window_len]u8 = undefined;
    var decompress: Decompress = .init(&reader, .zlib, &decompress_buffer);
    var out: [128]u8 = undefined;

    {
        const n = try decompress.reader.readSliceShort(&out);
        try std.testing.expectEqual(46, n);
        try std.testing.expectEqualStrings("Copyright Willem van Schaik, Singapore 1995-96", out[0..n]);
    }

    const n = try reader.readSliceShort(&out);
    try std.testing.expectEqual(n, 4);
    try std.testing.expectEqualSlices(u8, data[data.len - 4 .. data.len], out[0..n]);
}
