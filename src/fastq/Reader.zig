//! Streaming four-line FASTQ reader over `ByteSource`.
//!
//! Borrowed `Record` slices are valid until the next `next()` or `advance()` call.

const std = @import("std");
const record = @import("Record.zig");
const parse_error = @import("ParseError.zig");
const structure = @import("structure.zig");
const ByteSource = @import("../io/ByteSource.zig").ByteSource;
const limits = @import("../io/limits.zig");

pub const Options = struct {
    max_line_bytes: usize = limits.DEFAULT_MAX_LINE_BYTES,
};

const Range = struct {
    start: usize,
    end: usize,

    fn slice(self: Range, bytes: []const u8) []const u8 {
        return bytes[self.start..self.end];
    }
};

const Line = struct {
    range: Range,
    start_offset: u64,
};

pub const Reader = struct {
    allocator: std.mem.Allocator,
    source: ByteSource,
    buf: []u8,
    fill_end: usize,
    cursor: usize,
    record_buf: []u8,
    record_len: usize,
    record_index: u64,
    byte_offset: u64,
    options: Options,
    machine: structure.Machine,
    last_error: ?parse_error.ParseError,
    header_range: Range = undefined,
    sequence_range: Range = undefined,
    plus_range: Range = undefined,
    quality_range: Range = undefined,

    /// The source wrapper is copied; its referenced adapter must outlive the reader.
    pub fn init(
        allocator: std.mem.Allocator,
        source: ByteSource,
        options: Options,
    ) !Reader {
        const buf = try allocator.alloc(u8, limits.DEFAULT_READER_BUFFER_BYTES);
        errdefer allocator.free(buf);
        const record_buf = try allocator.alloc(u8, limits.DEFAULT_READER_BUFFER_BYTES);
        return .{
            .allocator = allocator,
            .source = source,
            .buf = buf,
            .fill_end = 0,
            .cursor = 0,
            .record_buf = record_buf,
            .record_len = 0,
            .record_index = 0,
            .byte_offset = 0,
            .options = options,
            .machine = .{},
            .last_error = null,
        };
    }

    pub fn deinit(self: *Reader) void {
        self.allocator.free(self.record_buf);
        self.allocator.free(self.buf);
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

    pub fn next(self: *Reader) parse_error.ReaderError!?record.Record {
        self.beginRecord();
        while (true) {
            const line = try self.readLine();
            switch (try self.ingestLine(line)) {
                .eof => return null,
                .continue_ => {},
                .record_ready => {
                    const bytes = self.record_buf[0..self.record_len];
                    const header = self.header_range.slice(bytes);
                    return .{
                        .header = header,
                        .id = record.firstToken(header),
                        .sequence = self.sequence_range.slice(bytes),
                        .plus = self.plus_range.slice(bytes),
                        .quality = self.quality_range.slice(bytes),
                    };
                },
            }
        }
    }

    pub fn advance(self: *Reader) parse_error.ReaderError!bool {
        self.beginRecord();
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

    fn beginRecord(self: *Reader) void {
        if (self.machine.expected == .header) self.record_len = 0;
    }

    fn ingestLine(self: *Reader, line: ?Line) parse_error.ReaderError!IngestResult {
        const actual_line = line orelse {
            const missing_line = self.machine.missingLine() orelse return .eof;
            self.storeError(
                .s004_truncated_record,
                structure.truncatedMessage(missing_line),
                missing_line,
                self.byte_offset,
            );
            return error.S004TruncatedRecord;
        };
        const line_kind = self.machine.expected;
        const content = actual_line.range.slice(self.record_buf);
        const first_byte = if (content.len == 0) null else content[0];
        const record_ready = self.machine.push(content.len, first_byte) catch |err| {
            return self.structuralError(err, actual_line.start_offset);
        };

        switch (line_kind) {
            .header => {
                self.header_range = .{
                    .start = actual_line.range.start + 1,
                    .end = actual_line.range.end,
                };
            },
            .sequence => {
                self.sequence_range = actual_line.range;
            },
            .plus => {
                self.plus_range = .{
                    .start = actual_line.range.start + 1,
                    .end = actual_line.range.end,
                };
            },
            .quality => {
                self.quality_range = actual_line.range;
            },
        }
        if (record_ready) {
            self.record_index += 1;
            return .record_ready;
        }
        return .continue_;
    }

    fn structuralError(
        self: *Reader,
        err: structure.Error,
        offset: u64,
    ) parse_error.ReaderError {
        const details = structure.diagnostic(err);
        self.storeError(details.code, details.message, details.line, offset);
        return err;
    }

    fn storeError(
        self: *Reader,
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
        std.debug.assert(space > 0);
        const n = self.source.read(self.buf[self.fill_end..]) catch return error.Io;
        if (n > space) return error.Io;
        self.fill_end += n;
        return n > 0;
    }

    fn readLine(self: *Reader) parse_error.ReaderError!?Line {
        const content_start = self.record_len;
        const line_start_offset = self.byte_offset;

        while (true) {
            if (self.cursor < self.fill_end) {
                const haystack = self.buf[self.cursor..self.fill_end];
                if (std.mem.indexOfScalar(u8, haystack, '\n')) |rel| {
                    try self.appendLineBytes(content_start, haystack[0..rel]);
                    self.cursor += rel + 1;
                    self.byte_offset += @intCast(rel + 1);
                    self.stripLineCr(content_start);
                    return .{
                        .range = .{ .start = content_start, .end = self.record_len },
                        .start_offset = line_start_offset,
                    };
                }

                try self.appendLineBytes(content_start, haystack);
                self.cursor = self.fill_end;
                self.byte_offset += @intCast(haystack.len);
            }

            const got_data = try self.refill();
            if (!got_data) {
                if (self.record_len > content_start) {
                    if (self.record_len - content_start > self.options.max_line_bytes) {
                        return error.LineTooLong;
                    }
                    return .{
                        .range = .{ .start = content_start, .end = self.record_len },
                        .start_offset = line_start_offset,
                    };
                }
                return null;
            }
        }
    }

    fn appendLineBytes(
        self: *Reader,
        content_start: usize,
        chunk: []const u8,
    ) parse_error.ReaderError!void {
        if (chunk.len == 0) return;
        const new_len = std.math.add(usize, self.record_len, chunk.len) catch
            return error.LineTooLong;
        const line_len = new_len - content_start;
        if (line_len > self.options.max_line_bytes) {
            const excess = line_len - self.options.max_line_bytes;
            if (excess > 1 or chunk[chunk.len - 1] != '\r') return error.LineTooLong;
        }
        try self.ensureRecordCapacity(new_len);
        @memcpy(self.record_buf[self.record_len..new_len], chunk);
        self.record_len = new_len;
    }

    fn stripLineCr(self: *Reader, content_start: usize) void {
        if (self.record_len > content_start and self.record_buf[self.record_len - 1] == '\r') {
            self.record_len -= 1;
        }
    }

    fn ensureRecordCapacity(self: *Reader, needed: usize) parse_error.ReaderError!void {
        if (needed <= self.record_buf.len) return;
        const doubled = std.math.add(usize, self.record_buf.len, self.record_buf.len) catch needed;
        const grown_len = @max(doubled, needed);
        self.record_buf = self.allocator.realloc(self.record_buf, grown_len) catch
            return error.OutOfMemory;
    }
};
