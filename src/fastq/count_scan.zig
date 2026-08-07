//! Single-pass plain FASTQ record counting.
//!
//! Zero heap allocations; fixed stack buffers only; portable `std.Io` reads.

const std = @import("std");
const limits = @import("../io/limits.zig");
const parse_error = @import("ParseError.zig");

pub const Options = struct {
    max_line_bytes: usize = limits.default_max_line_bytes,
};

const max_known_layouts = 8;
const max_header_index_bytes = 128;
const stride_block_batch = 256;
const no_layout_index = 255;

/// Fixed tail geometry for `+\n` plus lines: `seq\n`, `+\n`, `qual\n`.
const DenseLayout = struct {
    seq_len: usize,
    record_stride: usize,
    header_line_bytes: usize,

    fn fromRecordEnd(record_bytes: usize, seq_len: usize) ?DenseLayout {
        const tail = seq_len * 2 + 4;
        if (record_bytes <= tail) return null;
        return .{
            .seq_len = seq_len,
            .record_stride = record_bytes,
            .header_line_bytes = record_bytes - tail,
        };
    }

    fn validateRecordAt(data: []const u8, off: usize, layout: DenseLayout) bool {
        const end = off + layout.record_stride;
        if (end > data.len) return false;
        const h = layout.header_line_bytes;
        const s = layout.seq_len;
        return data[off] == '@' and
            data[off + h - 1] == '\n' and
            data[off + h + s] == '\n' and
            data[off + h + s + 1] == '+' and
            data[off + h + s + 2] == '\n' and
            data[end - 1] == '\n';
    }
};

const StrideRun = struct {
    count: usize = 0,
    stride: usize = 0,
};

pub const Scanner = struct {
    options: Options,
    state: State = .header,
    seq_len: usize = 0,
    fast_path_enabled: bool = false,
    layout: ?DenseLayout = null,
    known_layouts: [max_known_layouts]DenseLayout = undefined,
    known_layout_count: u8 = 0,
    layout_by_header: [max_header_index_bytes]u8 = [_]u8{no_layout_index} ** max_header_index_bytes,
    last_header_line_bytes: usize = 0,
    current_record_minimal_plus: bool = false,
    record_index: u64 = 0,
    byte_offset: u64 = 0,
    last_error: ?parse_error.ParseError = null,

    const State = enum {
        header,
        sequence,
        plus,
        quality,
    };

    pub fn init(options: Options) Scanner {
        return .{ .options = options };
    }

    pub fn takeLastError(self: *Scanner) ?parse_error.ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    pub fn feed(self: *Scanner, data: []const u8, at_eof: bool) parse_error.ReaderError!usize {
        var pos: usize = 0;
        while (pos < data.len) {
            if (self.fast_path_enabled and self.state == .header) {
                const advanced = try self.feedFast(data[pos..]);
                if (advanced > 0) {
                    self.byte_offset += @intCast(advanced);
                    pos += advanced;
                    continue;
                }
            }

            const advanced = try self.feedSlow(data[pos..], at_eof);
            if (advanced == 0) break;
            pos += advanced;
        }
        return pos;
    }

    pub fn finishEof(self: *Scanner) parse_error.ReaderError!void {
        switch (self.state) {
            .header => {},
            .sequence => {
                self.storeError(.s004_truncated_record, "unexpected end of file in sequence line", 2);
                return error.S004TruncatedRecord;
            },
            .plus => {
                self.storeError(.s004_truncated_record, "unexpected end of file in plus line", 3);
                return error.S004TruncatedRecord;
            },
            .quality => {
                self.storeError(.s004_truncated_record, "unexpected end of file in quality line", 4);
                return error.S004TruncatedRecord;
            },
        }
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
        const count = consumeStrideBlock(data, layout);
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
        if (idx == no_layout_index or idx >= self.known_layout_count) return null;
        return self.known_layouts[idx];
    }

    fn consumeStrideBlock(data: []const u8, layout: DenseLayout) usize {
        if (data.len < layout.record_stride or data[0] != '@') return 0;

        const stride = layout.record_stride;
        const max_records = data.len / stride;
        var accepted: usize = 0;

        while (accepted + stride_block_batch <= max_records) {
            var batch_ok: usize = 0;
            while (batch_ok < stride_block_batch) : (batch_ok += 1) {
                const off = (accepted + batch_ok) * stride;
                if (!DenseLayout.validateRecordAt(data, off, layout)) break;
            }
            if (batch_ok < stride_block_batch) {
                accepted += batch_ok;
                break;
            }
            accepted += stride_block_batch;
        }

        while (accepted < max_records) {
            const off = accepted * stride;
            if (!DenseLayout.validateRecordAt(data, off, layout)) break;
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
                if (DenseLayout.validateRecordAt(data, 0, layout)) {
                    self.state = .header;
                    self.layout = layout;
                    return layout.record_stride;
                }
            }
            var i: usize = 0;
            while (i < self.known_layout_count) : (i += 1) {
                const layout = self.known_layouts[i];
                if (layout.header_line_bytes != h) continue;
                if (DenseLayout.validateRecordAt(data, 0, layout)) {
                    self.state = .header;
                    self.layout = layout;
                    return layout.record_stride;
                }
            }
        }

        var i: usize = 0;
        while (i < self.known_layout_count) : (i += 1) {
            const layout = self.known_layouts[i];
            if (DenseLayout.validateRecordAt(data, 0, layout)) {
                self.state = .header;
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
        const after_seq = seq_start + seq_len;
        if (after_seq >= data.len or data[after_seq] != '\n') return 0;

        const plus_start = after_seq + 1;
        if (plus_start + 1 >= data.len) return 0;
        if (data[plus_start] != '+' or data[plus_start + 1] != '\n') return 0;

        const qual_start = plus_start + 2;
        const after_qual = qual_start + seq_len;
        if (after_qual >= data.len) return 0;
        if (data[after_qual] != '\n') return 0;

        const record_end = after_qual + 1;
        self.state = .header;
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
        if (self.known_layout_count < max_known_layouts) {
            const idx = self.known_layout_count;
            self.known_layouts[idx] = layout;
            self.known_layout_count += 1;
            if (layout.header_line_bytes < self.layout_by_header.len) {
                self.layout_by_header[layout.header_line_bytes] = idx;
            }
        }
    }

    fn feedSlow(self: *Scanner, data: []const u8, final_chunk: bool) parse_error.ReaderError!usize {
        var line_start: usize = 0;
        var pos: usize = 0;

        while (pos < data.len) {
            const rel = std.mem.indexOfScalar(u8, data[pos..], '\n');
            if (rel == null) break;
            pos += rel.?;

            var line_end = pos;
            if (line_end > line_start and data[line_end - 1] == '\r') line_end -= 1;
            try self.ingestLine(data[line_start..line_end]);
            self.byte_offset += @intCast(pos - line_start + 1);
            pos += 1;
            line_start = pos;

            if (self.fast_path_enabled and self.state == .header) {
                return pos;
            }
        }

        if (final_chunk and line_start < data.len) {
            try self.ingestLine(stripTrailingCr(data[line_start..]));
            return data.len;
        }

        return line_start;
    }

    fn disableFast(self: *Scanner) void {
        self.fast_path_enabled = false;
        self.layout = null;
        self.known_layout_count = 0;
        self.layout_by_header = [_]u8{no_layout_index} ** max_header_index_bytes;
        self.last_header_line_bytes = 0;
        self.current_record_minimal_plus = false;
    }

    fn learnDenseLayout(self: *Scanner, seq_len: usize, header_line_bytes: usize) void {
        const record_bytes = header_line_bytes + seq_len * 2 + 4;
        if (DenseLayout.fromRecordEnd(record_bytes, seq_len)) |layout| {
            self.rememberLayout(layout);
            self.layout = layout;
            self.fast_path_enabled = true;
        }
    }

    fn ingestLine(self: *Scanner, line: []const u8) parse_error.ReaderError!void {
        if (line.len > self.options.max_line_bytes) return error.LineTooLong;
        switch (self.state) {
            .header => {
                if (line.len == 0 or line[0] != '@') {
                    self.storeError(.s003_invalid_header, "header line must start with '@'", 1);
                    return error.S003InvalidHeader;
                }
                self.last_header_line_bytes = line.len + 1;
                self.state = .sequence;
            },
            .sequence => {
                self.seq_len = line.len;
                self.state = .plus;
            },
            .plus => {
                self.current_record_minimal_plus = line.len == 1 and line[0] == '+';
                self.state = .quality;
            },
            .quality => {
                if (line.len != self.seq_len) {
                    self.storeError(.s005_length_mismatch, "sequence and quality lengths differ", 4);
                    return error.S005LengthMismatch;
                }
                self.record_index += 1;
                self.state = .header;
                if (self.current_record_minimal_plus) {
                    self.learnDenseLayout(self.seq_len, self.last_header_line_bytes);
                }
            },
        }
    }

    fn storeError(self: *Scanner, code: parse_error.LintCode, message: []const u8, line: u3) void {
        self.last_error = .{
            .code = code,
            .message = message,
            .record_index = self.record_index,
            .byte_offset = self.byte_offset,
            .line_in_record = line,
        };
    }
};

fn stripTrailingCr(line: []const u8) []const u8 {
    if (line.len > 0 and line[line.len - 1] == '\r') return line[0 .. line.len - 1];
    return line;
}

fn lineContentLen(line: []const u8) usize {
    return stripTrailingCr(line).len;
}

fn headerLineBytesAt(data: []const u8) ?usize {
    if (data.len == 0 or data[0] != '@') return null;
    const scan_limit = @min(data.len, max_header_index_bytes);
    const rel = std.mem.indexOfScalar(u8, data[0..scan_limit], '\n') orelse return null;
    return rel + 1;
}

pub fn countPlainFile(
    io: std.Io,
    file: std.Io.File,
    options: Options,
    scanner: *Scanner,
) parse_error.ReaderError!u64 {
    scanner.* = Scanner.init(options);

    var buf: [limits.count_read_buffer_bytes]u8 = undefined;
    var filled: usize = 0;

    while (true) {
        const n = file.readStreaming(io, &.{buf[filled..]}) catch |err| switch (err) {
            error.EndOfStream => 0,
            else => return error.Io,
        };
        const at_eof = n == 0;
        if (filled == 0 and at_eof) break;

        filled += n;
        const consumed = try scanner.feed(buf[0..filled], at_eof);
        const tail = filled - consumed;
        if (tail > 0) std.mem.copyForwards(u8, buf[0..tail], buf[consumed..filled]);
        filled = tail;

        if (at_eof) {
            try scanner.finishEof();
            break;
        }
    }

    return scanner.record_index;
}

/// Count records in an in-memory slice (tests and embedders).
pub fn countSlice(data: []const u8, options: Options, scanner: *Scanner) parse_error.ReaderError!u64 {
    scanner.* = Scanner.init(options);
    _ = try scanner.feed(data, true);
    try scanner.finishEof();
    return scanner.record_index;
}
