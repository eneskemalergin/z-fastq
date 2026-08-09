//! Four-line FASTQ records, streaming reader and writer, and structural errors.

const std = @import("std");
const io_layer = @import("io.zig");

const ByteSource = io_layer.ByteSource;
const ByteSink = io_layer.ByteSink;
const WriteError = io_layer.WriteError;

/// Stable identifiers for structural failures in the supported four-line grammar.
pub const LintCode = enum {
    s001_invalid_plus_line,
    s003_invalid_header,
    s004_truncated_record,
    s005_length_mismatch,
};

/// Structural failure details with zero-based record and byte positions and a one-based line.
pub const ParseError = struct {
    code: LintCode,
    message: []const u8,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u3,
};

/// Returns the stable ASCII tag associated with `code`.
pub fn codeTag(code: LintCode) []const u8 {
    return switch (code) {
        .s001_invalid_plus_line => "S001",
        .s003_invalid_header => "S003",
        .s004_truncated_record => "S004",
        .s005_length_mismatch => "S005",
    };
}

/// Errors returned by the general reader and specialized count scanner.
pub const ReaderError = error{
    S001InvalidPlusLine,
    S003InvalidHeader,
    S004TruncatedRecord,
    S005LengthMismatch,
    LineTooLong,
    OutOfMemory,
    Io,
};

/// Borrowed FASTQ fields with structural prefix bytes and line endings removed.
///
/// Values returned by `Reader.next` remain valid only until the reader advances
/// again or is deinitialized. `id` is the first space- or tab-delimited token of
/// `header` and aliases the same storage.
pub const Record = struct {
    header: []const u8,
    id: []const u8,
    sequence: []const u8,
    plus: []const u8,
    quality: []const u8,
};

/// Independently allocated fields released together by `deinit`.
pub const OwnedRecord = struct {
    allocator: std.mem.Allocator,
    header: []u8,
    id: []u8,
    sequence: []u8,
    plus: []u8,
    quality: []u8,

    pub fn deinit(self: *OwnedRecord) void {
        self.allocator.free(self.header);
        self.allocator.free(self.id);
        self.allocator.free(self.sequence);
        self.allocator.free(self.plus);
        self.allocator.free(self.quality);
        self.* = undefined;
    }
};

fn firstToken(header: []const u8) []const u8 {
    const end = std.mem.indexOfAny(u8, header, "\t ") orelse header.len;
    return header[0..end];
}

/// Duplicates all record fields; the caller must invoke `OwnedRecord.deinit`.
pub fn toOwned(allocator: std.mem.Allocator, record: Record) !OwnedRecord {
    const header = try allocator.dupe(u8, record.header);
    errdefer allocator.free(header);
    const id = try allocator.dupe(u8, record.id);
    errdefer allocator.free(id);
    const sequence = try allocator.dupe(u8, record.sequence);
    errdefer allocator.free(sequence);
    const plus = try allocator.dupe(u8, record.plus);
    errdefer allocator.free(plus);
    const quality = try allocator.dupe(u8, record.quality);
    errdefer allocator.free(quality);
    return .{
        .allocator = allocator,
        .header = header,
        .id = id,
        .sequence = sequence,
        .plus = plus,
        .quality = quality,
    };
}

// --- Reader ---

pub const Options = struct {
    /// Maximum content bytes in one logical line, excluding LF and the optional CR.
    max_line_bytes: usize = io_layer.DEFAULT_MAX_LINE_BYTES,
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

/// Streaming parser constructed with `init`; fields are implementation state.
/// The copied source wrapper's referenced adapter must outlive the reader.
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
    machine: Machine,
    last_error: ?ParseError,
    header_range: Range = undefined,
    sequence_range: Range = undefined,
    plus_range: Range = undefined,
    quality_range: Range = undefined,

    pub fn init(
        allocator: std.mem.Allocator,
        source: ByteSource,
        options: Options,
    ) !Reader {
        const buf = try allocator.alloc(u8, io_layer.DEFAULT_READER_BUFFER_BYTES);
        errdefer allocator.free(buf);
        const record_buf = try allocator.alloc(u8, io_layer.DEFAULT_READER_BUFFER_BYTES);
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

    /// Returns the number of complete records consumed.
    pub fn recordIndex(self: *const Reader) u64 {
        return self.record_index;
    }

    /// Returns the number of input bytes consumed.
    pub fn byteOffset(self: *const Reader) u64 {
        return self.byte_offset;
    }

    /// Returns and clears structural details retained after a parse error.
    pub fn takeLastError(self: *Reader) ?ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    /// Returns the next borrowed record, or null at a clean EOF boundary.
    pub fn next(self: *Reader) ReaderError!?Record {
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
                        .id = firstToken(header),
                        .sequence = self.sequence_range.slice(bytes),
                        .plus = self.plus_range.slice(bytes),
                        .quality = self.quality_range.slice(bytes),
                    };
                },
            }
        }
    }

    /// Consumes one record without returning its fields, or returns false at clean EOF.
    pub fn advance(self: *Reader) ReaderError!bool {
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

    fn ingestLine(self: *Reader, line: ?Line) ReaderError!IngestResult {
        const actual_line = line orelse {
            const missing_line = self.machine.missingLine() orelse return .eof;
            self.storeError(
                .s004_truncated_record,
                truncatedMessage(missing_line),
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
        err: Error,
        offset: u64,
    ) ReaderError {
        const details = diagnostic(err);
        self.storeError(details.code, details.message, details.line, offset);
        return err;
    }

    fn storeError(
        self: *Reader,
        code: LintCode,
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

    fn refill(self: *Reader) ReaderError!bool {
        self.compactIfNeeded();
        const space = self.buf.len - self.fill_end;
        std.debug.assert(space > 0);
        const n = self.source.read(self.buf[self.fill_end..]) catch return error.Io;
        if (n > space) return error.Io;
        self.fill_end += n;
        return n > 0;
    }

    fn readLine(self: *Reader) ReaderError!?Line {
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
    ) ReaderError!void {
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

    fn ensureRecordCapacity(self: *Reader, needed: usize) ReaderError!void {
        if (needed <= self.record_buf.len) return;
        const doubled = std.math.add(usize, self.record_buf.len, self.record_buf.len) catch needed;
        const grown_len = @max(doubled, needed);
        self.record_buf = self.allocator.realloc(self.record_buf, grown_len) catch
            return error.OutOfMemory;
    }
};

// --- Writer ---

pub const WriterError = WriteError || error{InvalidRecord};

/// Streaming writer constructed with `init`; its referenced sink adapter must outlive it.
pub const Writer = struct {
    sink: ByteSink,

    pub fn init(sink: ByteSink) Writer {
        return .{ .sink = sink };
    }

    /// Validates and writes one complete record using LF line endings.
    pub fn writeRecord(self: *Writer, record: Record) WriterError!void {
        if (record.sequence.len != record.quality.len or
            !isWritableField(record.header) or
            !isWritableField(record.sequence) or
            !isWritableField(record.plus) or
            !isWritableField(record.quality))
        {
            return error.InvalidRecord;
        }

        try self.sink.write("@");
        try self.sink.write(record.header);
        try self.sink.write("\n");
        try self.sink.write(record.sequence);
        try self.sink.write("\n+");
        try self.sink.write(record.plus);
        try self.sink.write("\n");
        try self.sink.write(record.quality);
        try self.sink.write("\n");
    }

    /// Flushes the underlying sink when it provides a flush callback.
    pub fn flush(self: *Writer) WriteError!void {
        return self.sink.flush();
    }
};

fn isWritableField(bytes: []const u8) bool {
    return std.mem.indexOfScalar(u8, bytes, '\n') == null and
        (bytes.len == 0 or bytes[bytes.len - 1] != '\r');
}

// --- Structural validation ---

pub const ExpectedLine = enum {
    header,
    sequence,
    plus,
    quality,
};

pub const Error = error{
    S001InvalidPlusLine,
    S003InvalidHeader,
    S005LengthMismatch,
};

pub const Diagnostic = struct {
    code: LintCode,
    message: []const u8,
    line: u3,
};

pub fn diagnostic(err: Error) Diagnostic {
    return switch (err) {
        error.S001InvalidPlusLine => .{
            .code = .s001_invalid_plus_line,
            .message = "plus line must start with '+'",
            .line = 3,
        },
        error.S003InvalidHeader => .{
            .code = .s003_invalid_header,
            .message = "header line must start with '@'",
            .line = 1,
        },
        error.S005LengthMismatch => .{
            .code = .s005_length_mismatch,
            .message = "sequence and quality lengths differ",
            .line = 4,
        },
    };
}

pub fn truncatedMessage(line: u3) []const u8 {
    return switch (line) {
        2 => "unexpected end of file in sequence line",
        3 => "unexpected end of file in plus line",
        4 => "unexpected end of file in quality line",
        else => "unexpected end of file in record",
    };
}

pub const Machine = struct {
    expected: ExpectedLine = .header,
    sequence_len: usize = 0,

    pub fn push(self: *Machine, line_len: usize, first_byte: ?u8) Error!bool {
        switch (self.expected) {
            .header => {
                if (first_byte != '@') return error.S003InvalidHeader;
                self.expected = .sequence;
            },
            .sequence => {
                self.sequence_len = line_len;
                self.expected = .plus;
            },
            .plus => {
                if (first_byte != '+') return error.S001InvalidPlusLine;
                self.expected = .quality;
            },
            .quality => {
                if (line_len != self.sequence_len) return error.S005LengthMismatch;
                self.expected = .header;
                return true;
            },
        }
        return false;
    }

    pub fn missingLine(self: *const Machine) ?u3 {
        return switch (self.expected) {
            .header => null,
            .sequence => 2,
            .plus => 3,
            .quality => 4,
        };
    }
};
