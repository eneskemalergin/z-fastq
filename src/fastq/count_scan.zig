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

    fn validateRecord(data: []const u8, layout: DenseLayout) bool {
        if (data.len < layout.record_stride) return false;
        const h = layout.header_line_bytes;
        const s = layout.seq_len;
        const stride = layout.record_stride;
        return data[0] == '@' and
            data[h - 1] == '\n' and
            data[h + s] == '\n' and
            data[h + s + 1] == '+' and
            data[h + s + 2] == '\n' and
            data[stride - 1] == '\n';
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
    uniform_seq_len: ?usize = null,
    layout: ?DenseLayout = null,
    known_layouts: [max_known_layouts]DenseLayout = undefined,
    known_layout_count: u8 = 0,
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
            if (self.uniform_seq_len != null and self.state == .header) {
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
        const seq_len = self.uniform_seq_len.?;
        var cursor: usize = 0;

        while (cursor < data.len) {
            const run = consumeWithKnownLayouts(self, data[cursor..]);
            if (run.count > 0) {
                self.record_index += run.count;
                cursor += run.count * run.stride;
                continue;
            }

            const consumed = try self.tryFastRecord(data[cursor..], seq_len);
            if (consumed == 0) return cursor;
            self.record_index += 1;
            cursor += consumed;
        }
        return cursor;
    }

    fn consumeWithKnownLayouts(self: *Scanner, data: []const u8) StrideRun {
        if (self.layout) |layout| {
            const count = consumeStrideBlock(data, layout);
            if (count > 0) return .{ .count = count, .stride = layout.record_stride };
        }

        var i: usize = 0;
        while (i < self.known_layout_count) : (i += 1) {
            const layout = self.known_layouts[i];
            if (self.layout) |active| {
                if (layout.record_stride == active.record_stride) continue;
            }
            const count = consumeStrideBlock(data, layout);
            if (count > 0) {
                self.layout = layout;
                return .{ .count = count, .stride = layout.record_stride };
            }
        }

        return .{};
    }

    fn consumeStrideBlock(data: []const u8, layout: DenseLayout) usize {
        if (data.len < layout.record_stride or data[0] != '@') return 0;

        const max_records = data.len / layout.record_stride;
        var accepted: usize = 0;
        while (accepted < max_records) {
            const off = accepted * layout.record_stride;
            if (!DenseLayout.validateRecord(data[off..], layout)) break;
            accepted += 1;
        }
        return accepted;
    }

    fn tryFastRecord(self: *Scanner, data: []const u8, seq_len: usize) parse_error.ReaderError!usize {
        if (data.len == 0 or data[0] != '@') {
            self.disableFast();
            return 0;
        }

        if (self.layout) |layout| {
            if (layout.seq_len == seq_len and DenseLayout.validateRecord(data, layout)) {
                self.state = .header;
                return layout.record_stride;
            }
        }

        const nl1 = std.mem.findScalar(u8, data, '\n') orelse return 0;
        const hdr_len = lineContentLen(data[0..nl1]);
        if (hdr_len == 0 or hdr_len > self.options.max_line_bytes) {
            self.disableFast();
            if (hdr_len > self.options.max_line_bytes) return error.LineTooLong;
            return 0;
        }

        const seq_start = nl1 + 1;
        const after_seq = seq_start + seq_len;
        if (after_seq >= data.len or data[after_seq] != '\n') {
            self.disableFast();
            return 0;
        }

        const plus_start = after_seq + 1;
        if (plus_start + 1 >= data.len or data[plus_start] != '+' or data[plus_start + 1] != '\n') {
            self.disableFast();
            return 0;
        }

        const qual_start = plus_start + 2;
        const after_qual = qual_start + seq_len;
        if (after_qual >= data.len or data[after_qual] != '\n') {
            self.disableFast();
            return 0;
        }

        const record_end = after_qual + 1;
        self.state = .header;
        if (DenseLayout.fromRecordEnd(record_end, seq_len)) |layout| {
            self.rememberLayout(layout);
            self.layout = layout;
        }
        return record_end;
    }

    fn rememberLayout(self: *Scanner, layout: DenseLayout) void {
        var i: usize = 0;
        while (i < self.known_layout_count) : (i += 1) {
            if (self.known_layouts[i].record_stride == layout.record_stride) return;
        }
        if (self.known_layout_count < max_known_layouts) {
            self.known_layouts[self.known_layout_count] = layout;
            self.known_layout_count += 1;
        }
    }

    fn feedSlow(self: *Scanner, data: []const u8, final_chunk: bool) parse_error.ReaderError!usize {
        var line_start: usize = 0;
        var pos: usize = 0;

        while (pos < data.len) {
            if (data[pos] != '\n') {
                pos += 1;
                continue;
            }

            var line_end = pos;
            if (line_end > line_start and data[line_end - 1] == '\r') line_end -= 1;
            try self.ingestLine(data[line_start..line_end]);
            self.byte_offset += @intCast(pos - line_start + 1);
            pos += 1;
            line_start = pos;

            if (self.uniform_seq_len != null and self.state == .header) {
                return pos;
            }
        }

        if (final_chunk and line_start < data.len) {
            try self.ingestLine(trimCr(data[line_start..]));
            return data.len;
        }

        return line_start;
    }

    fn disableFast(self: *Scanner) void {
        self.uniform_seq_len = null;
        self.layout = null;
        self.known_layout_count = 0;
        self.last_header_line_bytes = 0;
        self.current_record_minimal_plus = false;
    }

    fn learnDenseLayout(self: *Scanner, seq_len: usize, header_line_bytes: usize) void {
        const record_bytes = header_line_bytes + seq_len * 2 + 4;
        if (DenseLayout.fromRecordEnd(record_bytes, seq_len)) |layout| {
            self.rememberLayout(layout);
            self.layout = layout;
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
                if (self.uniform_seq_len == null) {
                    self.uniform_seq_len = self.seq_len;
                } else if (self.seq_len != self.uniform_seq_len.?) {
                    self.disableFast();
                }
                if (self.uniform_seq_len != null and
                    self.seq_len == self.uniform_seq_len.? and
                    self.current_record_minimal_plus)
                {
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

fn trimCr(line: []const u8) []const u8 {
    if (line.len > 0 and line[line.len - 1] == '\r') return line[0 .. line.len - 1];
    return line;
}

fn lineContentLen(line: []const u8) usize {
    if (line.len > 0 and line[line.len - 1] == '\r') return line.len - 1;
    return line.len;
}

pub fn countPlainFile(
    io: std.Io,
    file: std.Io.File,
    options: Options,
    scanner: *Scanner,
) parse_error.ReaderError!u64 {
    scanner.* = Scanner.init(options);

    var buf: [limits.count_read_buffer_bytes]u8 = undefined;
    const io_buf = buf[0..limits.file_io_buffer_bytes];
    var read_buf = buf[limits.file_io_buffer_bytes..];
    var file_reader = file.reader(io, io_buf);

    var filled: usize = 0;

    while (true) {
        const n = file_reader.interface.readSliceShort(read_buf[filled..]) catch return error.Io;
        const at_eof = n == 0;
        if (filled == 0 and at_eof) break;

        filled += n;
        const consumed = try scanner.feed(read_buf[0..filled], at_eof);
        const tail = filled - consumed;
        if (tail > 0) {
            std.mem.copyForwards(u8, read_buf[0..tail], read_buf[consumed..filled]);
        }
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
