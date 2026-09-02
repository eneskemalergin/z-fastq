//! Single-pass plain FASTQ record counting.
//!
//! Zero heap allocations with incremental byte-slice input.

const std = @import("std");
const io_layer = @import("io.zig");
const fastq = @import("fastq.zig");

pub const Options = struct {
    /// Maximum content bytes in each logical FASTQ line, excluding CRLF or LF.
    max_line_bytes: usize = io_layer.DEFAULT_MAX_LINE_BYTES,
};

// Bare plus lines let the fast path derive every tail offset from the sequence length.
const DenseLayout = struct {
    sequence_start: usize,
    sequence_end: usize,
    plus_start: usize,
    quality_start: usize,
    quality_end: usize,
    record_stride: usize,

    fn fromRecordEnd(record_bytes: usize, seq_len: usize) ?DenseLayout {
        const doubled = std.math.mul(usize, seq_len, 2) catch return null;
        const tail = std.math.add(usize, doubled, 4) catch return null;
        if (record_bytes <= tail) return null;
        const sequence_start = record_bytes - tail;
        const sequence_end = std.math.add(usize, sequence_start, seq_len) catch return null;
        const plus_start = std.math.add(usize, sequence_end, 1) catch return null;
        const quality_start = std.math.add(usize, plus_start, 2) catch return null;
        const quality_end = std.math.add(usize, quality_start, seq_len) catch return null;
        const derived_record_bytes = std.math.add(usize, quality_end, 1) catch return null;
        if (derived_record_bytes != record_bytes) return null;
        return .{
            .sequence_start = sequence_start,
            .sequence_end = sequence_end,
            .plus_start = plus_start,
            .quality_start = quality_start,
            .quality_end = quality_end,
            .record_stride = record_bytes,
        };
    }

    fn validateRecordAt(
        data: []const u8,
        off: usize,
        layout: DenseLayout,
        max_line_bytes: usize,
    ) bool {
        const end = std.math.add(usize, off, layout.record_stride) catch return false;
        if (end > data.len) return false;
        const record = data[off..end];
        const header = record[0 .. layout.sequence_start - 1];
        const sequence = record[layout.sequence_start..layout.sequence_end];
        const quality = record[layout.quality_start..layout.quality_end];
        const sequence_len = lineContentLen(sequence);
        const quality_len = lineContentLen(quality);
        return fastq.headerPrefixIsValid(record[0], record[1]) and
            record[layout.sequence_start - 1] == '\n' and
            !@call(.always_inline, containsNewline, .{header}) and
            lineContentLen(header) <= max_line_bytes and
            record[layout.sequence_end] == '\n' and
            !@call(.always_inline, containsNewlinePair, .{ sequence, quality }) and
            sequence_len <= max_line_bytes and
            record[layout.plus_start] == '+' and
            record[layout.plus_start + 1] == '\n' and
            record[layout.quality_end] == '\n' and
            quality_len <= max_line_bytes and
            quality_len == sequence_len;
    }

    fn eql(a: DenseLayout, b: DenseLayout) bool {
        return a.sequence_start == b.sequence_start and
            a.sequence_end == b.sequence_end and
            a.plus_start == b.plus_start and
            a.quality_start == b.quality_start and
            a.quality_end == b.quality_end and
            a.record_stride == b.record_stride;
    }
};

const StrideRun = struct {
    count: usize = 0,
    stride: usize = 0,
};

const FastProgress = struct {
    bytes: usize = 0,
    records: usize = 0,
};

/// Allocation-free scanner; only `record_index` and `byte_offset` are caller-readable state.
pub const Scanner = struct {
    options: Options,
    machine: fastq.Machine = .{},
    fast_path_enabled: bool = false,
    layout: ?DenseLayout = null,
    layout_confirmed: bool = false,
    last_header_line_bytes: usize = 0,
    current_record_minimal_plus: bool = false,
    current_record_dense_eligible: bool = false,
    line_raw_len: usize = 0,
    line_first_byte: u8 = 0,
    line_second_byte: u8 = 0,
    line_last_byte: u8 = 0,
    line_start_offset: u64 = 0,
    record_index: u64 = 0,
    byte_offset: u64 = 0,
    last_error: ?fastq.ParseError = null,

    pub fn init(options: Options) Scanner {
        return .{ .options = options };
    }

    /// Returns and clears the structural diagnostic retained after a parse error.
    pub fn takeLastError(self: *Scanner) ?fastq.ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    // Keep this hot function cache-line aligned; half-line placement raised cycles.
    /// Consumes as many bytes as possible from one input chunk.
    ///
    /// On success the returned count equals `data.len`. Structural failures leave
    /// details for `takeLastError`.
    pub fn feed(self: *Scanner, data: []const u8) align(64) fastq.ReaderError!usize {
        var pos: usize = 0;
        while (pos < data.len) {
            if (self.fast_path_enabled and
                self.machine.expected == .header and
                self.line_raw_len == 0)
            {
                const progress = try self.feedFast(data[pos..]);
                if (progress.bytes > 0) {
                    const next_offset = try progressAfter(self.byte_offset, progress.bytes);
                    const next_record_index = try progressAfter(
                        self.record_index,
                        progress.records,
                    );
                    self.byte_offset = next_offset;
                    self.line_start_offset = next_offset;
                    self.record_index = next_record_index;
                    pos += progress.bytes;
                    continue;
                }
            }

            const advanced = try self.feedSlow(data[pos..]);
            if (advanced == 0) break;
            pos += advanced;
        }
        return pos;
    }

    /// Finalizes an unterminated last line and rejects an incomplete record.
    pub fn finishEof(self: *Scanner) fastq.ReaderError!void {
        if (self.line_raw_len > 0) try self.finishLine(false);
        const missing_line = self.machine.missingLine() orelse return;
        self.storeError(
            .s004_truncated_record,
            fastq.truncatedMessage(missing_line),
            missing_line,
            self.byte_offset,
        );
        return error.S004TruncatedRecord;
    }

    fn feedFast(self: *Scanner, data: []const u8) fastq.ReaderError!FastProgress {
        var cursor: usize = 0;
        var records: usize = 0;

        while (cursor < data.len) {
            const remaining = data[cursor..];
            if (self.layout_confirmed) {
                if (self.layout) |active| {
                    const run = self.tryStrideBlock(remaining, active);
                    if (run.count > 0) {
                        records = std.math.add(usize, records, run.count) catch
                            return error.ArithmeticLimit;
                        cursor += run.count * run.stride;
                        continue;
                    }
                    if (remaining.len >= active.record_stride) self.observeDenseLayout(null);
                }
            }

            const consumed = try self.tryFastRecord(remaining);
            if (consumed == 0) return .{ .bytes = cursor, .records = records };
            records = std.math.add(usize, records, 1) catch return error.ArithmeticLimit;
            cursor += consumed;
        }
        return .{ .bytes = cursor, .records = records };
    }

    fn tryStrideBlock(self: *Scanner, data: []const u8, layout: DenseLayout) StrideRun {
        if (data.len < layout.record_stride) return .{};
        const count = consumeStrideBlock(data, layout, self.options.max_line_bytes);
        if (count > 0) {
            return .{ .count = count, .stride = layout.record_stride };
        }
        return .{};
    }

    fn consumeStrideBlock(data: []const u8, layout: DenseLayout, max_line_bytes: usize) usize {
        if (data.len < layout.record_stride or data[0] != '@') return 0;

        const stride = layout.record_stride;
        const max_records = data.len / stride;
        var accepted: usize = 0;

        while (accepted < max_records) {
            const off = accepted * stride;
            if (!DenseLayout.validateRecordAt(data, off, layout, max_line_bytes)) break;
            accepted += 1;
        }

        return accepted;
    }

    fn tryFastRecord(self: *Scanner, data: []const u8) fastq.ReaderError!usize {
        if (data.len == 0) return 0;
        if (data[0] != '@') {
            self.disableFast();
            return 0;
        }

        var line_ends: [4]usize = undefined;
        if (!findRecordNewlines(data, &line_ends)) return 0;

        const header_end = line_ends[0];
        const sequence_start = header_end + 1;
        const sequence_end = line_ends[1];
        const sequence_len = sequence_end - sequence_start;
        const plus_start = sequence_end + 1;
        const plus_end = line_ends[2];
        const plus_len = plus_end - plus_start;
        const quality_start = plus_end + 1;
        const quality_end = line_ends[3];
        const quality_len = quality_end - quality_start;

        if ((header_end > 0 and data[header_end - 1] == '\r') or
            (sequence_len > 0 and data[sequence_end - 1] == '\r') or
            (plus_len > 0 and data[plus_end - 1] == '\r') or
            (quality_len > 0 and data[quality_end - 1] == '\r'))
        {
            return 0;
        }
        if (header_end < 2 or !fastq.headerPrefixIsValid(data[0], data[1])) {
            self.disableFast();
            return 0;
        }
        if (header_end > self.options.max_line_bytes) return error.LineTooLong;
        if (sequence_len > self.options.max_line_bytes) return error.LineTooLong;
        if (plus_len == 0 or data[plus_start] != '+') return 0;
        if (plus_len > self.options.max_line_bytes) return error.LineTooLong;
        if (quality_len > self.options.max_line_bytes) return error.LineTooLong;
        if (quality_len != sequence_len) return 0;

        const record_end = quality_end + 1;
        if (plus_len == 1) {
            self.observeDenseLayout(.{
                .sequence_start = sequence_start,
                .sequence_end = sequence_end,
                .plus_start = plus_start,
                .quality_start = quality_start,
                .quality_end = quality_end,
                .record_stride = record_end,
            });
        } else {
            self.observeDenseLayout(null);
        }
        self.fast_path_enabled = true;
        return record_end;
    }

    fn feedSlow(self: *Scanner, data: []const u8) fastq.ReaderError!usize {
        var pos: usize = 0;

        while (pos < data.len) {
            const rel = std.mem.indexOfScalar(u8, data[pos..], '\n');
            if (rel == null) {
                try self.consumeLineBytes(data[pos..]);
                self.byte_offset = try progressAfter(self.byte_offset, data.len - pos);
                return data.len;
            }

            const segment = data[pos .. pos + rel.?];
            try self.consumeLineBytes(segment);
            self.byte_offset = try progressAfter(self.byte_offset, segment.len + 1);
            pos += segment.len + 1;
            try self.finishLine(true);

            if (self.fast_path_enabled and self.machine.expected == .header) {
                return pos;
            }
        }
        return pos;
    }

    fn disableFast(self: *Scanner) void {
        self.fast_path_enabled = false;
        self.layout = null;
        self.layout_confirmed = false;
        self.last_header_line_bytes = 0;
        self.current_record_minimal_plus = false;
        self.current_record_dense_eligible = false;
    }

    fn learnDenseLayout(self: *Scanner, seq_len: usize, header_line_bytes: usize) void {
        const doubled = std.math.mul(usize, seq_len, 2) catch return;
        const with_header = std.math.add(usize, header_line_bytes, doubled) catch return;
        const record_bytes = std.math.add(usize, with_header, 4) catch return;
        if (DenseLayout.fromRecordEnd(record_bytes, seq_len)) |layout| {
            self.observeDenseLayout(layout);
            self.fast_path_enabled = true;
        }
    }

    fn observeDenseLayout(self: *Scanner, next: ?DenseLayout) void {
        const candidate = next orelse {
            self.layout = null;
            self.layout_confirmed = false;
            return;
        };
        self.layout_confirmed = if (self.layout) |previous| previous.eql(candidate) else false;
        self.layout = candidate;
    }

    fn consumeLineBytes(self: *Scanner, bytes: []const u8) fastq.ReaderError!void {
        if (bytes.len == 0) return;
        const previous_len = self.line_raw_len;
        if (previous_len == 0) self.line_first_byte = bytes[0];
        self.line_raw_len = std.math.add(usize, self.line_raw_len, bytes.len) catch
            return error.LineTooLong;
        if (previous_len < 2 and self.line_raw_len >= 2) {
            self.line_second_byte = bytes[1 - previous_len];
        }
        self.line_last_byte = bytes[bytes.len - 1];
        if (self.line_raw_len > self.options.max_line_bytes) {
            const excess = self.line_raw_len - self.options.max_line_bytes;
            if (excess > 1 or self.line_last_byte != '\r') return error.LineTooLong;
        }
    }

    fn finishLine(self: *Scanner, terminated_by_lf: bool) fastq.ReaderError!void {
        const had_cr = terminated_by_lf and
            self.line_raw_len > 0 and
            self.line_last_byte == '\r';
        const content_len = self.line_raw_len - @intFromBool(had_cr);
        if (content_len > self.options.max_line_bytes) return error.LineTooLong;

        const line_kind = self.machine.expected;
        switch (line_kind) {
            .header => {
                self.last_header_line_bytes = std.math.add(usize, self.line_raw_len, 1) catch
                    return error.LineTooLong;
                self.current_record_dense_eligible = !had_cr;
            },
            .sequence => {
                self.current_record_dense_eligible = self.current_record_dense_eligible and !had_cr;
            },
            .plus => {
                self.current_record_minimal_plus = content_len == 1;
                self.current_record_dense_eligible = self.current_record_dense_eligible and !had_cr;
            },
            .quality => {},
        }

        const first_byte = if (content_len == 0) null else self.line_first_byte;
        const second_byte = if (content_len < 2) null else self.line_second_byte;
        var machine = self.machine;
        const record_ready = machine.push(content_len, first_byte, second_byte) catch |err| {
            return self.structuralError(err, self.line_start_offset);
        };
        const next_record_index = if (record_ready)
            try progressAfter(self.record_index, 1)
        else
            self.record_index;
        self.machine = machine;
        if (record_ready) {
            self.record_index = next_record_index;
            if (self.current_record_dense_eligible and !had_cr) {
                self.fast_path_enabled = true;
                if (self.current_record_minimal_plus) {
                    self.learnDenseLayout(self.machine.sequence_len, self.last_header_line_bytes);
                } else {
                    self.observeDenseLayout(null);
                }
            } else {
                self.observeDenseLayout(null);
            }
        }
        self.line_raw_len = 0;
        self.line_first_byte = 0;
        self.line_second_byte = 0;
        self.line_last_byte = 0;
        self.line_start_offset = self.byte_offset;
    }

    fn structuralError(
        self: *Scanner,
        err: fastq.Error,
        offset: u64,
    ) fastq.ReaderError {
        const details = fastq.diagnostic(err);
        self.storeError(details.code, details.message, details.line, offset);
        return err;
    }

    fn storeError(
        self: *Scanner,
        code: fastq.LintCode,
        message: []const u8,
        line: u3,
        offset: u64,
    ) void {
        self.last_error = .{
            .code = code,
            .message = message,
            .record_index = self.record_index,
            .byte_offset = offset,
            .line_in_record = line,
        };
    }
};

fn progressAfter(current: u64, amount: usize) fastq.ReaderError!u64 {
    const amount_u64 = std.math.cast(u64, amount) orelse return error.ArithmeticLimit;
    return std.math.add(u64, current, amount_u64) catch error.ArithmeticLimit;
}

fn lineContentLen(line: []const u8) usize {
    if (line.len > 0 and line[line.len - 1] == '\r') return line.len - 1;
    return line.len;
}

fn containsNewline(bytes: []const u8) bool {
    const vector_len = std.simd.suggestVectorLength(u8) orelse @sizeOf(usize);
    const Vector = @Vector(vector_len, u8);
    const newline: Vector = @splat('\n');

    if (bytes.len >= vector_len) {
        const full_end = bytes.len - bytes.len % vector_len;
        var pos: usize = 0;
        while (pos < full_end) : (pos += vector_len) {
            const chunk: Vector = bytes[pos..][0..vector_len].*;
            if (@reduce(.Or, chunk == newline)) return true;
        }
        if (pos < bytes.len) {
            const chunk: Vector = bytes[bytes.len - vector_len ..][0..vector_len].*;
            if (@reduce(.Or, chunk == newline)) return true;
        }
        return false;
    }

    for (bytes) |byte| {
        if (byte == '\n') return true;
    }
    return false;
}

fn containsNewlinePair(first: []const u8, second: []const u8) bool {
    std.debug.assert(first.len == second.len);

    const vector_len = std.simd.suggestVectorLength(u8) orelse @sizeOf(usize);
    const Vector = @Vector(vector_len, u8);
    const newline: Vector = @splat('\n');

    if (first.len >= vector_len) {
        const full_end = first.len - first.len % vector_len;
        var pos: usize = 0;
        while (pos < full_end) : (pos += vector_len) {
            const first_chunk: Vector = first[pos..][0..vector_len].*;
            const second_chunk: Vector = second[pos..][0..vector_len].*;
            const first_has_newline = @reduce(.Or, first_chunk == newline);
            const second_has_newline = @reduce(.Or, second_chunk == newline);
            if (first_has_newline or second_has_newline) return true;
        }
        if (pos < first.len) {
            const first_chunk: Vector = first[first.len - vector_len ..][0..vector_len].*;
            const second_chunk: Vector = second[second.len - vector_len ..][0..vector_len].*;
            const first_has_newline = @reduce(.Or, first_chunk == newline);
            const second_has_newline = @reduce(.Or, second_chunk == newline);
            if (first_has_newline or second_has_newline) return true;
        }
        return false;
    }

    for (first, second) |first_byte, second_byte| {
        if (first_byte == '\n' or second_byte == '\n') return true;
    }
    return false;
}

fn findRecordNewlines(bytes: []const u8, line_ends: *[4]usize) bool {
    const Vector = @Vector(16, u8);
    const newline: Vector = @splat('\n');

    var found: usize = 0;
    var pos: usize = 0;
    while (bytes.len - pos >= 32) : (pos += 32) {
        const first: Vector = bytes[pos..][0..16].*;
        const second: Vector = bytes[pos + 16 ..][0..16].*;
        const first_mask: u16 = @bitCast(first == newline);
        const second_mask: u16 = @bitCast(second == newline);
        var mask = @as(u32, first_mask) | @as(u32, second_mask) << 16;
        while (mask != 0) {
            line_ends[found] = pos + @as(usize, @intCast(@ctz(mask)));
            found += 1;
            if (found == line_ends.len) return true;
            mask &= mask - 1;
        }
    }
    if (bytes.len - pos >= 16) {
        const block: Vector = bytes[pos..][0..16].*;
        var mask: u16 = @bitCast(block == newline);
        while (mask != 0) {
            line_ends[found] = pos + @as(usize, @intCast(@ctz(mask)));
            found += 1;
            if (found == line_ends.len) return true;
            mask &= mask - 1;
        }
        pos += 16;
    }
    for (bytes[pos..], pos..) |byte, index| {
        if (byte != '\n') continue;
        line_ends[found] = index;
        found += 1;
        if (found == line_ends.len) return true;
    }
    return false;
}

/// Counts all records in one slice and replaces `scanner` with the final state.
///
/// Structural failures leave details in `scanner` for `takeLastError`.
pub fn countSlice(
    data: []const u8,
    options: Options,
    scanner: *Scanner,
) fastq.ReaderError!u64 {
    scanner.* = Scanner.init(options);
    _ = try scanner.feed(data);
    try scanner.finishEof();
    return scanner.record_index;
}

test "[unit] - [dense layout]: construction derives complete record offsets" {
    const layout = DenseLayout.fromRecordEnd(11, 2).?;

    try std.testing.expect(layout.sequence_start == 3);
    try std.testing.expect(layout.sequence_end == 5);
    try std.testing.expect(layout.plus_start == 6);
    try std.testing.expect(layout.quality_start == 8);
    try std.testing.expect(layout.quality_end == 10);
    try std.testing.expect(layout.record_stride == 11);
}

test "[failure] - [dense layout]: construction rejects impossible geometry" {
    const max = std.math.maxInt(usize);

    try std.testing.expect(DenseLayout.fromRecordEnd(0, 0) == null);
    try std.testing.expect(DenseLayout.fromRecordEnd(8, 2) == null);
    try std.testing.expect(DenseLayout.fromRecordEnd(max, max) == null);
}

test "[edge] - [record newline scan]: vector boundaries and incomplete input are exact" {
    const cases = [_][4]usize{
        .{ 14, 15, 30, 31 },
        .{ 15, 16, 31, 32 },
        .{ 30, 31, 62, 63 },
        .{ 31, 32, 63, 64 },
        .{ 126, 127, 128, 129 },
    };

    for (cases) |expected| {
        var data: [131]u8 = @splat('x');
        for (expected) |index| data[index] = '\n';
        data[130] = '\n';
        var actual: [4]usize = undefined;

        try std.testing.expect(findRecordNewlines(&data, &actual));
        try std.testing.expectEqual(expected, actual);
        try std.testing.expect(!findRecordNewlines(data[0..expected[3]], &actual));
    }
}

test "[unit] - [derived record]: annotated plus lines stay on the complete-record path" {
    const data = "@a\nAC\n+first\n!!\n@longer\nG\n+second note\n#\n";
    var scanner = Scanner.init(.{ .max_line_bytes = 64 });
    scanner.fast_path_enabled = true;

    try std.testing.expectEqual(data.len, try scanner.feed(data));
    try std.testing.expectEqual(@as(u64, 2), scanner.record_index);
    try std.testing.expect(scanner.layout == null);
}

test "[edge] - [count scanner]: fast progress rejects maximum offset and record count" {
    const data = "@r\nA\n+\n!\n";

    var offset = Scanner.init(.{});
    offset.fast_path_enabled = true;
    offset.byte_offset = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, offset.feed(data));
    try std.testing.expectEqual(std.math.maxInt(u64), offset.byte_offset);
    try std.testing.expectEqual(@as(u64, 0), offset.record_index);

    var records = Scanner.init(.{});
    records.fast_path_enabled = true;
    records.record_index = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, records.feed(data));
    try std.testing.expectEqual(std.math.maxInt(u64), records.record_index);
    try std.testing.expectEqual(@as(u64, 0), records.byte_offset);
}

test "[edge] - [count scanner]: incremental progress rejects maximum offset and record count" {
    const data = "@r\nA\n+\n!\n";

    var offset = Scanner.init(.{});
    offset.byte_offset = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, offset.feed(data));
    try std.testing.expectEqual(std.math.maxInt(u64), offset.byte_offset);
    try std.testing.expectEqual(@as(u64, 0), offset.record_index);

    var records = Scanner.init(.{});
    records.record_index = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, records.feed(data));
    try std.testing.expectEqual(std.math.maxInt(u64), records.record_index);
    try std.testing.expectEqual(@as(u64, data.len), records.byte_offset);
}
