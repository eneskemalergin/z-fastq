//! Single-pass plain FASTQ record counting.
//!
//! Zero heap allocations with incremental byte-slice input.

const std = @import("std");
const limits = @import("../io/limits.zig");
const parse_error = @import("ParseError.zig");
const structure = @import("structure.zig");

pub const Options = struct {
    max_line_bytes: usize = limits.DEFAULT_MAX_LINE_BYTES,
};

const MAX_KNOWN_LAYOUTS = 8;
const MAX_HEADER_INDEX_BYTES = 128;
const NO_LAYOUT_INDEX = 255;

// Fixed tail geometry for `+\n` plus lines: `seq\n`, `+\n`, `qual\n`.
const DenseLayout = struct {
    seq_len: usize,
    record_stride: usize,
    header_line_bytes: usize,

    fn fromRecordEnd(record_bytes: usize, seq_len: usize) ?DenseLayout {
        const doubled = std.math.mul(usize, seq_len, 2) catch return null;
        const tail = std.math.add(usize, doubled, 4) catch return null;
        if (record_bytes <= tail) return null;
        return .{
            .seq_len = seq_len,
            .record_stride = record_bytes,
            .header_line_bytes = record_bytes - tail,
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
        const h = layout.header_line_bytes;
        const s = layout.seq_len;
        if (h == 0 or h > layout.record_stride) return false;
        const seq_start = std.math.add(usize, off, h) catch return false;
        const seq_end = std.math.add(usize, seq_start, s) catch return false;
        const plus_start = std.math.add(usize, seq_end, 1) catch return false;
        const qual_start = std.math.add(usize, plus_start, 2) catch return false;
        const qual_end = std.math.add(usize, qual_start, s) catch return false;
        if (qual_end >= end or qual_end + 1 != end) return false;
        const header = data[off .. seq_start - 1];
        const sequence = data[seq_start..seq_end];
        const quality = data[qual_start..qual_end];
        const sequence_len = lineContentLen(sequence);
        const quality_len = lineContentLen(quality);
        return data[off] == '@' and
            data[seq_start - 1] == '\n' and
            !containsNewline(header) and
            lineContentLen(header) <= max_line_bytes and
            data[seq_end] == '\n' and
            !containsNewline(sequence) and
            sequence_len <= max_line_bytes and
            data[plus_start] == '+' and
            data[plus_start + 1] == '\n' and
            data[qual_end] == '\n' and
            !containsNewline(quality) and
            quality_len <= max_line_bytes and
            quality_len == sequence_len;
    }
};

const StrideRun = struct {
    count: usize = 0,
    stride: usize = 0,
};

pub const Scanner = struct {
    options: Options,
    machine: structure.Machine = .{},
    fast_path_enabled: bool = false,
    layout: ?DenseLayout = null,
    known_layouts: [MAX_KNOWN_LAYOUTS]DenseLayout = undefined,
    known_layout_count: u8 = 0,
    layout_by_header: [MAX_HEADER_INDEX_BYTES]u8 = [_]u8{NO_LAYOUT_INDEX} ** MAX_HEADER_INDEX_BYTES,
    last_header_line_bytes: usize = 0,
    current_record_minimal_plus: bool = false,
    current_record_dense_eligible: bool = false,
    line_raw_len: usize = 0,
    line_first_byte: u8 = 0,
    line_last_byte: u8 = 0,
    line_start_offset: u64 = 0,
    record_index: u64 = 0,
    byte_offset: u64 = 0,
    last_error: ?parse_error.ParseError = null,

    pub fn init(options: Options) Scanner {
        return .{ .options = options };
    }

    pub fn takeLastError(self: *Scanner) ?parse_error.ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    pub fn feed(self: *Scanner, data: []const u8) parse_error.ReaderError!usize {
        var pos: usize = 0;
        while (pos < data.len) {
            if (self.fast_path_enabled and
                self.machine.expected == .header and
                self.line_raw_len == 0)
            {
                const advanced = try self.feedFast(data[pos..]);
                if (advanced > 0) {
                    self.byte_offset += @intCast(advanced);
                    self.line_start_offset = self.byte_offset;
                    pos += advanced;
                    continue;
                }
            }

            const advanced = try self.feedSlow(data[pos..]);
            if (advanced == 0) break;
            pos += advanced;
        }
        return pos;
    }

    pub fn finishEof(self: *Scanner) parse_error.ReaderError!void {
        if (self.line_raw_len > 0) try self.finishLine();
        const missing_line = self.machine.missingLine() orelse return;
        self.storeError(
            .s004_truncated_record,
            structure.truncatedMessage(missing_line),
            missing_line,
            self.byte_offset,
        );
        return error.S004TruncatedRecord;
    }

    fn feedFast(self: *Scanner, data: []const u8) parse_error.ReaderError!usize {
        var cursor: usize = 0;

        while (cursor < data.len) {
            const run = consumeWithKnownLayouts(self, data[cursor..]);
            if (run.count > 0) {
                self.record_index += run.count;
                cursor += run.count * run.stride;
                continue;
            }

            const consumed = try self.tryFastRecord(data[cursor..]);
            if (consumed == 0) return cursor;
            self.record_index += 1;
            cursor += consumed;
        }
        return cursor;
    }

    fn smallestKnownStride(self: *const Scanner) usize {
        var min = self.known_layouts[0].record_stride;
        var i: usize = 1;
        while (i < self.known_layout_count) : (i += 1) {
            min = @min(min, self.known_layouts[i].record_stride);
        }
        return min;
    }

    fn tryStrideBlock(self: *Scanner, data: []const u8, layout: DenseLayout) StrideRun {
        if (data.len < layout.record_stride) return .{};
        const count = consumeStrideBlock(data, layout, self.options.max_line_bytes);
        if (count > 0) {
            self.layout = layout;
            return .{ .count = count, .stride = layout.record_stride };
        }
        return .{};
    }

    fn consumeWithKnownLayouts(self: *Scanner, data: []const u8) StrideRun {
        if (self.known_layout_count == 0) return .{};
        if (data.len < self.smallestKnownStride()) return .{};
        if (data[0] != '@') return .{};

        const header_bytes = headerLineBytesAt(data);
        var tried_stride: ?usize = null;

        if (header_bytes) |h| {
            if (self.layoutByHeader(h)) |layout| {
                tried_stride = layout.record_stride;
                const run = self.tryStrideBlock(data, layout);
                if (run.count > 0) return run;
            }
        }

        if (self.layout) |active| {
            const already_tried = tried_stride != null and active.record_stride == tried_stride.?;
            const header_ok = header_bytes == null or active.header_line_bytes == header_bytes.?;
            if (!already_tried and header_ok) {
                const run = self.tryStrideBlock(data, active);
                if (run.count > 0) return run;
            }
        }

        var i: usize = 0;
        while (i < self.known_layout_count) : (i += 1) {
            const layout = self.known_layouts[i];
            if (tried_stride) |t| {
                if (layout.record_stride == t) continue;
            }
            if (self.layout) |active| {
                if (layout.record_stride == active.record_stride) continue;
            }
            if (header_bytes) |h| {
                if (layout.header_line_bytes != h) continue;
            }
            const run = self.tryStrideBlock(data, layout);
            if (run.count > 0) return run;
        }

        return .{};
    }

    fn layoutByHeader(self: *const Scanner, header_line_bytes: usize) ?DenseLayout {
        if (header_line_bytes >= self.layout_by_header.len) return null;
        const idx = self.layout_by_header[header_line_bytes];
        if (idx == NO_LAYOUT_INDEX or idx >= self.known_layout_count) return null;
        return self.known_layouts[idx];
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

    fn tryFastRecord(self: *Scanner, data: []const u8) parse_error.ReaderError!usize {
        if (data.len == 0) return 0;
        if (data[0] != '@') {
            self.disableFast();
            return 0;
        }

        if (headerLineBytesAt(data)) |h| {
            if (self.layoutByHeader(h)) |layout| {
                if (DenseLayout.validateRecordAt(data, 0, layout, self.options.max_line_bytes)) {
                    self.layout = layout;
                    return layout.record_stride;
                }
            }
            var i: usize = 0;
            while (i < self.known_layout_count) : (i += 1) {
                const layout = self.known_layouts[i];
                if (layout.header_line_bytes != h) continue;
                if (DenseLayout.validateRecordAt(data, 0, layout, self.options.max_line_bytes)) {
                    self.layout = layout;
                    return layout.record_stride;
                }
            }
        }

        var i: usize = 0;
        while (i < self.known_layout_count) : (i += 1) {
            const layout = self.known_layouts[i];
            if (DenseLayout.validateRecordAt(data, 0, layout, self.options.max_line_bytes)) {
                self.layout = layout;
                return layout.record_stride;
            }
        }

        const nl1 = std.mem.findScalar(u8, data, '\n') orelse return 0;
        const hdr_len = lineContentLen(data[0..nl1]);
        if (hdr_len == 0) {
            self.disableFast();
            return 0;
        }
        if (hdr_len > self.options.max_line_bytes) return error.LineTooLong;

        const seq_start = nl1 + 1;
        const nl2_rel = std.mem.indexOfScalar(u8, data[seq_start..], '\n') orelse return 0;
        const seq_len = nl2_rel;
        const after_seq = std.math.add(usize, seq_start, seq_len) catch return 0;
        if (after_seq >= data.len or data[after_seq] != '\n') return 0;
        const sequence_len = lineContentLen(data[seq_start..after_seq]);
        if (sequence_len > self.options.max_line_bytes) return error.LineTooLong;

        const plus_start = std.math.add(usize, after_seq, 1) catch return 0;
        if (plus_start + 1 >= data.len) return 0;
        if (data[plus_start] != '+' or data[plus_start + 1] != '\n') return 0;

        const qual_start = std.math.add(usize, plus_start, 2) catch return 0;
        const after_qual = std.math.add(usize, qual_start, seq_len) catch return 0;
        if (after_qual >= data.len) return 0;
        if (containsNewline(data[qual_start..after_qual])) return 0;
        if (data[after_qual] != '\n') return 0;
        const quality_len = lineContentLen(data[qual_start..after_qual]);
        if (quality_len > self.options.max_line_bytes) return error.LineTooLong;
        if (quality_len != sequence_len) return 0;

        const record_end = after_qual + 1;
        if (DenseLayout.fromRecordEnd(record_end, seq_len)) |layout| {
            self.rememberLayout(layout);
            self.layout = layout;
            self.fast_path_enabled = true;
        }
        return record_end;
    }

    fn rememberLayout(self: *Scanner, layout: DenseLayout) void {
        var i: usize = 0;
        while (i < self.known_layout_count) : (i += 1) {
            if (self.known_layouts[i].record_stride == layout.record_stride) {
                if (layout.header_line_bytes < self.layout_by_header.len) {
                    self.layout_by_header[layout.header_line_bytes] = @intCast(i);
                }
                return;
            }
        }
        if (self.known_layout_count < MAX_KNOWN_LAYOUTS) {
            const idx = self.known_layout_count;
            self.known_layouts[idx] = layout;
            self.known_layout_count += 1;
            if (layout.header_line_bytes < self.layout_by_header.len) {
                self.layout_by_header[layout.header_line_bytes] = idx;
            }
        }
    }

    fn feedSlow(self: *Scanner, data: []const u8) parse_error.ReaderError!usize {
        var pos: usize = 0;

        while (pos < data.len) {
            const rel = std.mem.indexOfScalar(u8, data[pos..], '\n');
            if (rel == null) {
                try self.consumeLineBytes(data[pos..]);
                self.byte_offset += @intCast(data.len - pos);
                return data.len;
            }

            const segment = data[pos .. pos + rel.?];
            try self.consumeLineBytes(segment);
            self.byte_offset += @intCast(segment.len + 1);
            pos += segment.len + 1;
            try self.finishLine();

            if (self.fast_path_enabled and self.machine.expected == .header) {
                return pos;
            }
        }
        return pos;
    }

    fn disableFast(self: *Scanner) void {
        self.fast_path_enabled = false;
        self.layout = null;
        self.known_layout_count = 0;
        self.layout_by_header = [_]u8{NO_LAYOUT_INDEX} ** MAX_HEADER_INDEX_BYTES;
        self.last_header_line_bytes = 0;
        self.current_record_minimal_plus = false;
        self.current_record_dense_eligible = false;
    }

    fn learnDenseLayout(self: *Scanner, seq_len: usize, header_line_bytes: usize) void {
        const doubled = std.math.mul(usize, seq_len, 2) catch return;
        const with_header = std.math.add(usize, header_line_bytes, doubled) catch return;
        const record_bytes = std.math.add(usize, with_header, 4) catch return;
        if (DenseLayout.fromRecordEnd(record_bytes, seq_len)) |layout| {
            self.rememberLayout(layout);
            self.layout = layout;
            self.fast_path_enabled = true;
        }
    }

    fn consumeLineBytes(self: *Scanner, bytes: []const u8) parse_error.ReaderError!void {
        if (bytes.len == 0) return;
        if (self.line_raw_len == 0) self.line_first_byte = bytes[0];
        self.line_raw_len = std.math.add(usize, self.line_raw_len, bytes.len) catch
            return error.LineTooLong;
        self.line_last_byte = bytes[bytes.len - 1];
        if (self.line_raw_len > self.options.max_line_bytes) {
            const excess = self.line_raw_len - self.options.max_line_bytes;
            if (excess > 1 or self.line_last_byte != '\r') return error.LineTooLong;
        }
    }

    fn finishLine(self: *Scanner) parse_error.ReaderError!void {
        const had_cr = self.line_raw_len > 0 and self.line_last_byte == '\r';
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
        const record_ready = self.machine.push(content_len, first_byte) catch |err| {
            return self.structuralError(err, self.line_start_offset);
        };
        if (record_ready) {
            self.record_index += 1;
            if (self.current_record_minimal_plus and
                self.current_record_dense_eligible and
                !had_cr)
            {
                self.learnDenseLayout(self.machine.sequence_len, self.last_header_line_bytes);
            }
        }
        self.line_raw_len = 0;
        self.line_first_byte = 0;
        self.line_last_byte = 0;
        self.line_start_offset = self.byte_offset;
    }

    fn structuralError(
        self: *Scanner,
        err: structure.Error,
        offset: u64,
    ) parse_error.ReaderError {
        const details = structure.diagnostic(err);
        self.storeError(details.code, details.message, details.line, offset);
        return err;
    }

    fn storeError(
        self: *Scanner,
        code: parse_error.LintCode,
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

fn headerLineBytesAt(data: []const u8) ?usize {
    if (data.len == 0 or data[0] != '@') return null;
    const scan_limit = @min(data.len, MAX_HEADER_INDEX_BYTES);
    const rel = std.mem.indexOfScalar(u8, data[0..scan_limit], '\n') orelse return null;
    return rel + 1;
}

pub fn countSlice(
    data: []const u8,
    options: Options,
    scanner: *Scanner,
) parse_error.ReaderError!u64 {
    scanner.* = Scanner.init(options);
    _ = try scanner.feed(data);
    try scanner.finishEof();
    return scanner.record_index;
}
