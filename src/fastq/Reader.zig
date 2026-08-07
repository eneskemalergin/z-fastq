//! Streaming FASTQ reader over `ByteSource`.
//!
//! Accepts a plus line without a leading `+` when it matches the sequence line.
//! Borrowed `Record` slices are valid until the next `next()` call.

const std = @import("std");
const record = @import("Record.zig");
const parse_error = @import("ParseError.zig");
const ByteSource = @import("../io/ByteSource.zig").ByteSource;
const limits = @import("../io/limits.zig");

pub const Options = struct {
    max_line_bytes: usize = limits.default_max_line_bytes,
};

const State = enum {
    header,
    sequence,
    plus,
    quality,
};

pub const Reader = struct {
    allocator: std.mem.Allocator,
    source: *const ByteSource,
    buf: []u8,
    fill_end: usize,
    cursor: usize,
    record_index: u64,
    byte_offset: u64,
    options: Options,
    state: State,
    line_accum: ?[]u8,
    line_accum_len: usize,
    last_error: ?parse_error.ParseError,
    header_line: []const u8 = undefined,
    sequence_line: []const u8 = undefined,
    plus_line: []const u8 = undefined,
    quality_line: []const u8 = undefined,

    pub fn init(
        allocator: std.mem.Allocator,
        source: *const ByteSource,
        options: Options,
    ) !Reader {
        const buf = try allocator.alloc(u8, limits.default_reader_buffer_bytes);
        return .{
            .allocator = allocator,
            .source = source,
            .buf = buf,
            .fill_end = 0,
            .cursor = 0,
            .record_index = 0,
            .byte_offset = 0,
            .options = options,
            .state = .header,
            .line_accum = null,
            .line_accum_len = 0,
            .last_error = null,
        };
    }

    pub fn deinit(self: *Reader, allocator: std.mem.Allocator) void {
        if (self.line_accum) |accum| allocator.free(accum);
        allocator.free(self.buf);
        self.* = undefined;
    }

    pub fn recordIndex(self: *const Reader) u64 {
        return self.record_index;
    }

    pub fn byteOffset(self: *const Reader) u64 {
        return self.byte_offset;
    }

    pub fn takeLastError(self: *Reader) ?parse_error.ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    pub fn next(self: *Reader) parse_error.ReaderError! ?record.Record {
        while (true) {
            const line = try self.readLine();
            switch (try self.ingestLine(line)) {
                .eof => return null,
                .continue_ => {},
                .record_ready => return record.Record{
                    .header = self.header_line,
                    .id = record.firstToken(self.header_line),
                    .sequence = self.sequence_line,
                    .plus = self.plus_line,
                    .quality = self.quality_line,
                },
            }
        }
    }

    /// Advance one record without building a `Record`.
    pub fn advance(self: *Reader) parse_error.ReaderError!bool {
        while (true) {
            const line = try self.readLine();
            switch (try self.ingestLine(line)) {
                .eof => return false,
                .continue_ => {},
                .record_ready => return true,
            }
        }
    }

    const IngestResult = enum {
        eof,
        continue_,
        record_ready,
    };

    fn ingestLine(self: *Reader, line: ?[]const u8) parse_error.ReaderError!IngestResult {
        switch (self.state) {
            .header => {
                if (line == null) return .eof;
                const header_line = line.?;
                if (header_line.len == 0 or header_line[0] != '@') {
                    self.storeError(.s003_invalid_header, "header line must start with '@'", 1);
                    return error.S003InvalidHeader;
                }
                self.header_line = header_line[1..];
                self.state = .sequence;
                return .continue_;
            },
            .sequence => {
                const seq_line = line orelse {
                    self.storeError(.s004_truncated_record, "unexpected end of file in sequence line", 2);
                    return error.S004TruncatedRecord;
                };
                self.sequence_line = seq_line;
                self.state = .plus;
                return .continue_;
            },
            .plus => {
                const plus_line = line orelse {
                    self.storeError(.s004_truncated_record, "unexpected end of file in plus line", 3);
                    return error.S004TruncatedRecord;
                };
                if (plus_line.len > 0 and plus_line[0] == '+') {
                    self.plus_line = plus_line[1..];
                } else {
                    self.plus_line = plus_line;
                }
                self.state = .quality;
                return .continue_;
            },
            .quality => {
                const qual_line = line orelse {
                    self.storeError(.s004_truncated_record, "unexpected end of file in quality line", 4);
                    return error.S004TruncatedRecord;
                };
                self.quality_line = qual_line;
                if (self.quality_line.len != self.sequence_line.len) {
                    self.storeError(.s005_length_mismatch, "sequence and quality lengths differ", 4);
                    return error.S005LengthMismatch;
                }
                self.record_index += 1;
                self.state = .header;
                return .record_ready;
            },
        }
    }

    fn storeError(self: *Reader, code: parse_error.LintCode, message: []const u8, line: u3) void {
        self.last_error = .{
            .code = code,
            .message = message,
            .record_index = self.record_index,
            .byte_offset = self.byte_offset,
            .line_in_record = line,
        };
    }

    fn compactIfNeeded(self: *Reader) void {
        if (self.cursor >= self.fill_end) {
            self.cursor = 0;
            self.fill_end = 0;
            return;
        }
        if (self.cursor == 0) return;

        const tail_len = self.fill_end - self.cursor;
        std.mem.copyForwards(u8, self.buf[0..tail_len], self.buf[self.cursor..self.fill_end]);
        self.fill_end = tail_len;
        self.cursor = 0;
    }

    fn refill(self: *Reader) parse_error.ReaderError!bool {
        self.compactIfNeeded();
        const space = self.buf.len - self.fill_end;
        if (space == 0) return error.LineTooLong;
        const n = self.source.read(self.buf[self.fill_end..]) catch return error.Io;
        self.fill_end += n;
        return n > 0;
    }

    fn readLine(self: *Reader) parse_error.ReaderError!?[]const u8 {
        self.line_accum_len = 0;

        while (true) {
            if (self.cursor < self.fill_end) {
                const haystack = self.buf[self.cursor..self.fill_end];
                if (std.mem.indexOfScalar(u8, haystack, '\n')) |rel| {
                    const line_start = self.cursor;
                    const line_end = self.cursor + rel;
                    var content = self.buf[line_start..line_end];
                    self.cursor = line_end + 1;
                    self.byte_offset += @as(u64, @intCast(rel + 1));
                    if (content.len > 0 and content[content.len - 1] == '\r') {
                        content = content[0 .. content.len - 1];
                    }
                    if (self.line_accum_len > 0) {
                        try self.appendToLineAccum(content);
                        const out = self.line_accum.?[0..self.line_accum_len];
                        self.line_accum_len = 0;
                        return out;
                    }
                    if (content.len > self.options.max_line_bytes) return error.LineTooLong;
                    return content;
                }

                const partial = haystack;
                if (partial.len > 0) {
                    try self.appendToLineAccum(partial);
                    self.cursor = self.fill_end;
                }
            }

            const got_data = try self.refill();
            if (!got_data) {
                if (self.line_accum_len > 0) {
                    if (self.line_accum_len > self.options.max_line_bytes) return error.LineTooLong;
                    const out = self.line_accum.?[0..self.line_accum_len];
                    self.line_accum_len = 0;
                    return out;
                }
                return null;
            }
        }
    }

    fn appendToLineAccum(self: *Reader, chunk: []const u8) parse_error.ReaderError!void {
        const new_len = std.math.add(usize, self.line_accum_len, chunk.len) catch return error.LineTooLong;
        if (new_len > self.options.max_line_bytes) return error.LineTooLong;
        try self.ensureLineAccumCapacity(new_len);
        const accum = self.line_accum.?;
        @memcpy(accum[self.line_accum_len..][0..chunk.len], chunk);
        self.line_accum_len = new_len;
    }

    fn ensureLineAccumCapacity(self: *Reader, needed: usize) parse_error.ReaderError!void {
        if (self.line_accum) |accum| {
            if (needed <= accum.len) return;
            const grown_len = @min(
                self.options.max_line_bytes,
                @max(accum.len * 2, needed),
            );
            self.line_accum = self.allocator.realloc(accum, grown_len) catch return error.OutOfMemory;
            return;
        }

        const initial = @min(self.options.max_line_bytes, @max(needed, self.buf.len));
        self.line_accum = self.allocator.alloc(u8, initial) catch return error.OutOfMemory;
    }
};
