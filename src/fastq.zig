//! Four-line FASTQ records, streaming reader and writer, and record validation.

const std = @import("std");
const io_layer = @import("io.zig");

const ByteSource = io_layer.ByteSource;
const ByteSink = io_layer.ByteSink;
const WriteError = io_layer.WriteError;

pub const LintCode = enum {
    s001_invalid_plus_line,
    s002_invalid_sequence_alphabet,
    s003_invalid_header,
    s004_truncated_record,
    s005_length_mismatch,
    s006_invalid_quality_range,
};
pub const ParseError = struct {
    code: LintCode,
    message: []const u8,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u3,
};

pub fn codeTag(code: LintCode) []const u8 {
    return switch (code) {
        .s001_invalid_plus_line => "S001",
        .s002_invalid_sequence_alphabet => "S002",
        .s003_invalid_header => "S003",
        .s004_truncated_record => "S004",
        .s005_length_mismatch => "S005",
        .s006_invalid_quality_range => "S006",
    };
}

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

pub const Alphabet = enum {
    iupac,
    acgtn,
};

pub const ValidationOptions = struct {
    alphabet: Alphabet = .iupac,
};

pub const SemanticField = enum {
    sequence,
    quality,
};

/// Describes the first semantic failure at a zero-based field-relative byte index.
pub const SemanticError = struct {
    code: LintCode,
    message: []const u8,
    field: SemanticField,
    byte_index: usize,
};

/// Checks only sequence alphabet and Phred+33 bytes without allocating.
/// Structural prefixes, completeness, and equal field lengths remain caller-owned.
pub fn validateRecord(record: Record, options: ValidationOptions) ?SemanticError {
    if (record.sequence.len == record.quality.len) {
        if (std.simd.suggestVectorLength(u8)) |vector_len| {
            return switch (options.alphabet) {
                .iupac => validateRecordVector(vector_len, .iupac, record),
                .acgtn => validateRecordVector(vector_len, .acgtn, record),
            };
        }
    }
    if (firstInvalidSequence(record.sequence, options.alphabet)) |byte_index| {
        return .{
            .code = .s002_invalid_sequence_alphabet,
            .message = "sequence byte is outside the selected alphabet",
            .field = .sequence,
            .byte_index = byte_index,
        };
    }
    if (firstInvalidQuality(record.quality)) |byte_index| {
        return .{
            .code = .s006_invalid_quality_range,
            .message = "quality byte must be ASCII 33 through 126",
            .field = .quality,
            .byte_index = byte_index,
        };
    }
    return null;
}

fn validateRecordVector(
    comptime vector_len: comptime_int,
    comptime alphabet: Alphabet,
    record: Record,
) ?SemanticError {
    const Bytes = @Vector(vector_len, u8);
    const minimum: Bytes = @splat(33);
    const maximum: Bytes = @splat(126);

    var quality_failure: ?usize = null;
    var byte_index: usize = 0;
    while (record.sequence.len - byte_index >= vector_len) : (byte_index += vector_len) {
        const sequence: Bytes = record.sequence[byte_index..][0..vector_len].*;
        const sequence_invalid = invalidSequenceVector(vector_len, alphabet, sequence);
        if (@reduce(.Or, sequence_invalid)) {
            return .{
                .code = .s002_invalid_sequence_alphabet,
                .message = "sequence byte is outside the selected alphabet",
                .field = .sequence,
                .byte_index = byte_index + std.simd.firstTrue(sequence_invalid).?,
            };
        }

        if (quality_failure == null) {
            const quality: Bytes = record.quality[byte_index..][0..vector_len].*;
            const quality_invalid = (quality < minimum) | (quality > maximum);
            if (@reduce(.Or, quality_invalid)) {
                quality_failure = byte_index + std.simd.firstTrue(quality_invalid).?;
            }
        }
    }
    while (byte_index < record.sequence.len) : (byte_index += 1) {
        if (!alphabetAccepts(alphabet, record.sequence[byte_index])) {
            return .{
                .code = .s002_invalid_sequence_alphabet,
                .message = "sequence byte is outside the selected alphabet",
                .field = .sequence,
                .byte_index = byte_index,
            };
        }
        if (quality_failure == null and
            (record.quality[byte_index] < 33 or record.quality[byte_index] > 126))
        {
            quality_failure = byte_index;
        }
    }
    if (quality_failure) |invalid_index| {
        return .{
            .code = .s006_invalid_quality_range,
            .message = "quality byte must be ASCII 33 through 126",
            .field = .quality,
            .byte_index = invalid_index,
        };
    }
    return null;
}

/// Decodes one Phred+33 byte and rejects values outside ASCII 33 through 126.
pub fn decodePhred33(quality_byte: u8) error{InvalidQuality}!u8 {
    if (quality_byte < 33 or quality_byte > 126) return error.InvalidQuality;
    return quality_byte - 33;
}

fn firstInvalidSequence(sequence: []const u8, alphabet: Alphabet) ?usize {
    if (std.simd.suggestVectorLength(u8)) |vector_len| {
        return switch (alphabet) {
            .iupac => firstInvalidSequenceVector(vector_len, .iupac, sequence),
            .acgtn => firstInvalidSequenceVector(vector_len, .acgtn, sequence),
        };
    }
    return firstInvalidSequenceScalar(sequence, alphabet, 0);
}

fn firstInvalidSequenceVector(
    comptime vector_len: comptime_int,
    comptime alphabet: Alphabet,
    sequence: []const u8,
) ?usize {
    const Bytes = @Vector(vector_len, u8);

    var byte_index: usize = 0;
    while (sequence.len - byte_index >= vector_len) : (byte_index += vector_len) {
        const bytes: Bytes = sequence[byte_index..][0..vector_len].*;
        const invalid = invalidSequenceVector(vector_len, alphabet, bytes);
        if (@reduce(.Or, invalid)) {
            return byte_index + std.simd.firstTrue(invalid).?;
        }
    }
    return firstInvalidSequenceScalar(sequence[byte_index..], alphabet, byte_index);
}

fn invalidSequenceVector(
    comptime vector_len: comptime_int,
    comptime alphabet: Alphabet,
    bytes: @Vector(vector_len, u8),
) @Vector(vector_len, bool) {
    const Bytes = @Vector(vector_len, u8);
    const normalized = bytes & @as(Bytes, @splat(0xdf));
    return switch (alphabet) {
        .iupac => (normalized -% @as(Bytes, @splat('A')) > @as(Bytes, @splat(3))) &
            (normalized -% @as(Bytes, @splat('G')) > @as(Bytes, @splat(1))) &
            (normalized != @as(Bytes, @splat('K'))) &
            (normalized -% @as(Bytes, @splat('M')) > @as(Bytes, @splat(1))) &
            (normalized -% @as(Bytes, @splat('R')) > @as(Bytes, @splat(5))) &
            (normalized != @as(Bytes, @splat('Y'))),
        .acgtn => (normalized != @as(Bytes, @splat('A'))) &
            (normalized != @as(Bytes, @splat('C'))) &
            (normalized != @as(Bytes, @splat('G'))) &
            (normalized != @as(Bytes, @splat('T'))) &
            (normalized != @as(Bytes, @splat('N'))),
    };
}

fn firstInvalidSequenceScalar(
    sequence: []const u8,
    alphabet: Alphabet,
    start_index: usize,
) ?usize {
    for (sequence, start_index..) |byte, byte_index| {
        if (!alphabetAccepts(alphabet, byte)) return byte_index;
    }
    return null;
}

fn firstInvalidQuality(quality: []const u8) ?usize {
    if (std.simd.suggestVectorLength(u8)) |vector_len| {
        return firstInvalidQualityVector(vector_len, quality);
    }
    return firstInvalidQualityScalar(quality, 0);
}

fn firstInvalidQualityVector(
    comptime vector_len: comptime_int,
    quality: []const u8,
) ?usize {
    const Bytes = @Vector(vector_len, u8);
    const minimum: Bytes = @splat(33);
    const maximum: Bytes = @splat(126);

    var byte_index: usize = 0;
    while (quality.len - byte_index >= vector_len) : (byte_index += vector_len) {
        const bytes: Bytes = quality[byte_index..][0..vector_len].*;
        const invalid = (bytes < minimum) | (bytes > maximum);
        if (@reduce(.Or, invalid)) {
            return byte_index + std.simd.firstTrue(invalid).?;
        }
    }
    return firstInvalidQualityScalar(quality[byte_index..], byte_index);
}

fn firstInvalidQualityScalar(quality: []const u8, start_index: usize) ?usize {
    for (quality, start_index..) |byte, byte_index| {
        _ = decodePhred33(byte) catch return byte_index;
    }
    return null;
}

fn alphabetAccepts(alphabet: Alphabet, byte: u8) bool {
    return switch (alphabet) {
        .iupac => switch (byte) {
            'A',
            'C',
            'G',
            'T',
            'U',
            'R',
            'Y',
            'S',
            'W',
            'K',
            'M',
            'B',
            'D',
            'H',
            'V',
            'N',
            'a',
            'c',
            'g',
            't',
            'u',
            'r',
            'y',
            's',
            'w',
            'k',
            'm',
            'b',
            'd',
            'h',
            'v',
            'n',
            => true,
            else => false,
        },
        .acgtn => switch (byte) {
            'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n' => true,
            else => false,
        },
    };
}

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
    content: []const u8,
    start_offset: u64,
};

const FallbackField = struct {
    storage: []u8 = &.{},
    len: usize = 0,
};

const BufferedRecord = struct {
    ranges: [4]Range,
    canonical_range: ?Range,
};

const BufferedRecordResult = union(enum) {
    incomplete,
    eof,
    record: BufferedRecord,
};

pub const RecordOffsets = struct {
    header: u64,
    sequence: u64,
    plus: u64,
    quality: u64,
};

/// Streaming parser constructed with `init`; fields are implementation state.
/// The copied source wrapper's referenced adapter must outlive the reader.
pub const Reader = struct {
    allocator: std.mem.Allocator,
    source: ByteSource,
    buf: []u8,
    fill_end: usize,
    cursor: usize,
    fallback_fields: [4]FallbackField,
    record_index: u64,
    byte_offset: u64,
    options: Options,
    machine: Machine,
    last_error: ?ParseError,
    record_offsets: RecordOffsets = undefined,
    current_record_offsets: ?RecordOffsets = null,

    pub fn init(
        allocator: std.mem.Allocator,
        source: ByteSource,
        options: Options,
    ) !Reader {
        const buf = try allocator.alloc(u8, io_layer.DEFAULT_READER_BUFFER_BYTES);
        errdefer allocator.free(buf);
        const fallback_sequence = try allocator.alloc(u8, io_layer.DEFAULT_READER_BUFFER_BYTES);
        return .{
            .allocator = allocator,
            .source = source,
            .buf = buf,
            .fill_end = 0,
            .cursor = 0,
            .fallback_fields = .{
                .{},
                .{ .storage = fallback_sequence },
                .{},
                .{},
            },
            .record_index = 0,
            .byte_offset = 0,
            .options = options,
            .machine = .{},
            .last_error = null,
        };
    }

    pub fn deinit(self: *Reader) void {
        for (self.fallback_fields) |field| {
            if (field.storage.len != 0) self.allocator.free(field.storage);
        }
        self.allocator.free(self.buf);
        self.* = undefined;
    }

    pub fn recordIndex(self: *const Reader) u64 {
        return self.record_index;
    }

    pub fn byteOffset(self: *const Reader) u64 {
        return self.byte_offset;
    }

    /// Returns decoded-stream line-start offsets for the current borrowed record.
    pub fn currentRecordOffsets(self: *const Reader) ?RecordOffsets {
        return self.current_record_offsets;
    }

    /// Returns and clears structural details retained after a parse error.
    pub fn takeLastError(self: *Reader) ?ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    /// Returns the next borrowed record, or null at a clean EOF boundary.
    pub fn next(self: *Reader) ReaderError!?Record {
        var canonical_span: ?[]const u8 = null;
        return self.nextWithCanonicalSpan(&canonical_span, true);
    }

    fn nextWithCanonicalSpan(
        self: *Reader,
        canonical_span: *?[]const u8,
        comptime derive_id: bool,
    ) ReaderError!?Record {
        canonical_span.* = null;
        self.beginRecord();
        switch (try self.readBufferedRecord()) {
            .incomplete => return self.nextFallback(derive_id),
            .eof => return null,
            .record => |buffered| return self.finishBufferedRecord(
                buffered,
                canonical_span,
                derive_id,
            ),
        }
    }

    fn finishBufferedRecord(
        self: *Reader,
        buffered: BufferedRecord,
        canonical_span: *?[]const u8,
        comptime derive_id: bool,
    ) Record {
        self.current_record_offsets = self.record_offsets;
        const ranges = buffered.ranges;
        const header = self.buf[ranges[0].start + 1 .. ranges[0].end];
        canonical_span.* = if (buffered.canonical_range) |range|
            range.slice(self.buf)
        else
            null;
        return .{
            .header = header,
            .id = if (derive_id) firstToken(header) else header[0..0],
            .sequence = ranges[1].slice(self.buf),
            .plus = self.buf[ranges[2].start + 1 .. ranges[2].end],
            .quality = ranges[3].slice(self.buf),
        };
    }

    fn nextFallback(self: *Reader, comptime derive_id: bool) ReaderError!?Record {
        self.beginFallbackRecord();
        while (true) {
            const line = try self.readLine();
            switch (try self.ingestLine(line)) {
                .eof => return null,
                .continue_ => {},
                .record_ready => {
                    self.current_record_offsets = self.record_offsets;
                    const header_field = &self.fallback_fields[0];
                    const plus_field = &self.fallback_fields[2];
                    const header_bytes = header_field.storage[0..header_field.len];
                    const plus_bytes = plus_field.storage[0..plus_field.len];
                    const header = header_bytes[1..];
                    return .{
                        .header = header,
                        .id = if (derive_id) firstToken(header) else header[0..0],
                        .sequence = self.fallback_fields[1].storage[0..self.fallback_fields[1].len],
                        .plus = plus_bytes[1..],
                        .quality = self.fallback_fields[3].storage[0..self.fallback_fields[3].len],
                    };
                },
            }
        }
    }

    /// Consumes one record without returning its fields, or returns false at clean EOF.
    pub fn advance(self: *Reader) ReaderError!bool {
        self.beginRecord();
        switch (try self.readBufferedRecord()) {
            .incomplete => return self.advanceFallback(),
            .eof => return false,
            .record => return true,
        }
    }

    fn advanceFallback(self: *Reader) ReaderError!bool {
        self.beginFallbackRecord();
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
        if (self.machine.expected == .header) {
            self.current_record_offsets = null;
        }
    }

    fn beginFallbackRecord(self: *Reader) void {
        if (self.machine.expected != .header) return;
        for (&self.fallback_fields) |*field| field.len = 0;
    }

    fn readBufferedRecord(self: *Reader) ReaderError!BufferedRecordResult {
        if (self.machine.expected != .header) return .incomplete;
        if (self.cursor == self.fill_end and !try self.refill()) return .eof;

        var starts: [4]usize = undefined;
        var ends: [4]usize = undefined;
        var search_start = self.cursor;
        for (0..4) |line_index| {
            const rel = std.mem.indexOfScalar(
                u8,
                self.buf[search_start..self.fill_end],
                '\n',
            ) orelse return .incomplete;
            starts[line_index] = search_start;
            ends[line_index] = search_start + rel;
            search_start = ends[line_index] + 1;
        }

        var ranges: [4]Range = undefined;
        var canonical = true;
        for (0..4) |line_index| {
            const start = starts[line_index];
            const raw_end = ends[line_index];
            const end = if (raw_end > start and self.buf[raw_end - 1] == '\r')
                raw_end - 1
            else
                raw_end;
            canonical = canonical and end == raw_end;
            if (end - start > self.options.max_line_bytes) return error.LineTooLong;

            const start_offset = self.byte_offset;
            self.cursor = raw_end + 1;
            self.byte_offset += @intCast(raw_end + 1 - start);

            const line_kind = self.machine.expected;
            const content = self.buf[start..end];
            const first_byte = if (content.len == 0) null else content[0];
            const record_ready = self.machine.push(content.len, first_byte) catch |err| {
                return self.structuralError(err, start_offset);
            };
            ranges[line_index] = .{ .start = start, .end = end };
            switch (line_kind) {
                .header => self.record_offsets.header = start_offset,
                .sequence => self.record_offsets.sequence = start_offset,
                .plus => self.record_offsets.plus = start_offset,
                .quality => self.record_offsets.quality = start_offset,
            }
            std.debug.assert(record_ready == (line_index == 3));
        }

        self.record_index += 1;
        return .{ .record = .{
            .ranges = ranges,
            .canonical_range = if (canonical)
                .{ .start = starts[0], .end = search_start }
            else
                null,
        } };
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
        const content = actual_line.content;
        const first_byte = if (content.len == 0) null else content[0];
        const record_ready = self.machine.push(content.len, first_byte) catch |err| {
            return self.structuralError(err, actual_line.start_offset);
        };

        switch (line_kind) {
            .header => {
                self.record_offsets.header = actual_line.start_offset;
            },
            .sequence => {
                self.record_offsets.sequence = actual_line.start_offset;
            },
            .plus => {
                self.record_offsets.plus = actual_line.start_offset;
            },
            .quality => {
                self.record_offsets.quality = actual_line.start_offset;
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
        const field_index = @intFromEnum(self.machine.expected);
        const field = &self.fallback_fields[field_index];
        const content_start = field.len;
        const line_start_offset = self.byte_offset;

        while (true) {
            if (self.cursor < self.fill_end) {
                const haystack = self.buf[self.cursor..self.fill_end];
                if (std.mem.indexOfScalar(u8, haystack, '\n')) |rel| {
                    try self.appendLineBytes(field_index, content_start, haystack[0..rel], true);
                    self.cursor += rel + 1;
                    self.byte_offset += @intCast(rel + 1);
                    self.stripLineCr(field_index, content_start);
                    return .{
                        .content = field.storage[0..field.len],
                        .start_offset = line_start_offset,
                    };
                }

                try self.appendLineBytes(field_index, content_start, haystack, false);
                self.cursor = self.fill_end;
                self.byte_offset += @intCast(haystack.len);
            }

            const got_data = try self.refill();
            if (!got_data) {
                if (field.len > content_start) {
                    if (field.len - content_start > self.options.max_line_bytes) {
                        return error.LineTooLong;
                    }
                    return .{
                        .content = field.storage[0..field.len],
                        .start_offset = line_start_offset,
                    };
                }
                return null;
            }
        }
    }

    fn appendLineBytes(
        self: *Reader,
        field_index: usize,
        content_start: usize,
        chunk: []const u8,
        complete: bool,
    ) ReaderError!void {
        if (chunk.len == 0) return;
        const field = &self.fallback_fields[field_index];
        const new_len = std.math.add(usize, field.len, chunk.len) catch
            return error.LineTooLong;
        const line_len = new_len - content_start;
        if (line_len > self.options.max_line_bytes) {
            const excess = line_len - self.options.max_line_bytes;
            if (excess > 1 or chunk[chunk.len - 1] != '\r') return error.LineTooLong;
        }
        try self.ensureFieldCapacity(field_index, new_len, complete);
        @memcpy(field.storage[field.len..new_len], chunk);
        field.len = new_len;
    }

    fn stripLineCr(self: *Reader, field_index: usize, content_start: usize) void {
        const field = &self.fallback_fields[field_index];
        if (field.len > content_start and field.storage[field.len - 1] == '\r') {
            field.len -= 1;
        }
    }

    fn ensureFieldCapacity(
        self: *Reader,
        field_index: usize,
        needed: usize,
        complete: bool,
    ) ReaderError!void {
        const field = &self.fallback_fields[field_index];
        if (needed <= field.storage.len) return;
        const storage_limit = fieldStorageLimit(self.options.max_line_bytes);
        if (needed > storage_limit) return error.LineTooLong;
        const doubled = std.math.add(usize, field.storage.len, field.storage.len) catch needed;
        var grown_len = if (complete)
            needed
        else
            @min(
                @max(doubled, @max(io_layer.DEFAULT_READER_BUFFER_BYTES, needed)),
                storage_limit,
            );
        if (self.machine.expected == .quality) {
            const quality_capacity = std.math.add(usize, self.machine.sequence_len, 1) catch
                return error.LineTooLong;
            if (quality_capacity >= needed) grown_len = quality_capacity;
        }
        field.storage = if (field.storage.len == 0)
            self.allocator.alloc(u8, grown_len) catch return error.OutOfMemory
        else
            self.allocator.realloc(field.storage, grown_len) catch return error.OutOfMemory;
    }
};

fn fieldStorageLimit(max_line_bytes: usize) usize {
    return std.math.add(usize, max_line_bytes, 1) catch std.math.maxInt(usize);
}

pub fn nextWithoutId(
    reader: *Reader,
    canonical_span: *?[]const u8,
) ReaderError!?Record {
    return reader.nextWithCanonicalSpan(canonical_span, false);
}

pub fn nextBufferedWithoutId(
    reader: *Reader,
    canonical_span: *?[]const u8,
) ReaderError!?Record {
    canonical_span.* = null;
    reader.beginRecord();
    if (reader.cursor == reader.fill_end) return null;
    return switch (try reader.readBufferedRecord()) {
        .incomplete => null,
        .eof => unreachable,
        .record => |buffered| reader.finishBufferedRecord(
            buffered,
            canonical_span,
            false,
        ),
    };
}

// --- Writer ---

pub const WriterError = WriteError || error{InvalidRecord};

/// Streaming writer constructed with `init`; its referenced sink adapter must outlive it.
pub const Writer = struct {
    sink: ByteSink,

    pub fn init(sink: ByteSink) Writer {
        return .{ .sink = sink };
    }

    /// Rejects invalid fields before output; accepted records use LF line endings.
    pub fn writeRecord(self: *Writer, record: Record) WriterError!void {
        if (record.sequence.len != record.quality.len or
            !isWritableField(record.header) or
            !isWritableField(record.sequence) or
            !isWritableField(record.plus) or
            !isWritableField(record.quality))
        {
            return error.InvalidRecord;
        }

        return writeRecordFields(self, record);
    }

    /// Flushes the underlying sink when it provides a flush callback.
    pub fn flush(self: *Writer) WriteError!void {
        return self.sink.flush();
    }
};

pub fn writeValidatedRecord(writer: *Writer, record: Record) WriteError!void {
    return writeRecordFields(writer, record);
}

pub fn writeCanonicalRecordSpan(writer: *Writer, span: []const u8) WriteError!void {
    return writer.sink.write(span);
}

fn writeRecordFields(writer: *Writer, record: Record) WriteError!void {
    try writer.sink.write("@");
    try writer.sink.write(record.header);
    try writer.sink.write("\n");
    try writer.sink.write(record.sequence);
    try writer.sink.write("\n+");
    try writer.sink.write(record.plus);
    try writer.sink.write("\n");
    try writer.sink.write(record.quality);
    try writer.sink.write("\n");
}

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

pub const CheckScannerError = error{
    Format,
    LineTooLong,
    ArithmeticLimit,
};

const STRUCTURAL_BLOCK_BYTES = 64;

fn newlineMask(block: *const [STRUCTURAL_BLOCK_BYTES]u8) u64 {
    const Bytes = @Vector(STRUCTURAL_BLOCK_BYTES, u8);
    const bytes: Bytes = block.*;
    return @bitCast(bytes == @as(Bytes, @splat('\n')));
}

pub const CheckScanner = struct {
    max_line_bytes: usize,
    alphabet: Alphabet,
    machine: Machine = .{},
    line_len: usize = 0,
    first_byte: ?u8 = null,
    pending_cr: bool = false,
    line_start_offset: u64 = 0,
    sequence_start_offset: u64 = 0,
    quality_start_offset: u64 = 0,
    sequence_failure: ?usize = null,
    quality_failure: ?usize = null,
    record_index: u64 = 0,
    byte_offset: u64 = 0,
    last_error: ?ParseError = null,

    pub fn init(options: Options, validation_options: ValidationOptions) CheckScanner {
        return .{
            .max_line_bytes = options.max_line_bytes,
            .alphabet = validation_options.alphabet,
        };
    }

    pub fn feed(self: *CheckScanner, data: []const u8) CheckScannerError!usize {
        var segment_start: usize = 0;
        var block_start: usize = 0;
        while (data.len - block_start >= STRUCTURAL_BLOCK_BYTES) {
            var line_feeds = newlineMask(data[block_start..][0..STRUCTURAL_BLOCK_BYTES]);
            while (line_feeds != 0) {
                const line_end = block_start + @as(usize, @intCast(@ctz(line_feeds)));
                try self.consumeLineBytes(data[segment_start..line_end]);
                try self.advanceOffset(line_end - segment_start + 1);
                try self.finishLine(true);
                segment_start = line_end + 1;
                line_feeds &= line_feeds - 1;
            }
            block_start += STRUCTURAL_BLOCK_BYTES;
        }

        for (data[block_start..], block_start..) |byte, line_end| {
            if (byte != '\n') continue;
            try self.consumeLineBytes(data[segment_start..line_end]);
            try self.advanceOffset(line_end - segment_start + 1);
            try self.finishLine(true);
            segment_start = line_end + 1;
        }
        if (segment_start < data.len) {
            try self.consumeLineBytes(data[segment_start..]);
            try self.advanceOffset(data.len - segment_start);
        }
        return data.len;
    }

    pub fn finishEof(self: *CheckScanner) CheckScannerError!void {
        if (self.line_len != 0) {
            if (self.pending_cr) {
                self.pending_cr = false;
                if (self.line_len > self.max_line_bytes) return error.LineTooLong;
                self.classifyByte('\r', self.line_len - 1);
            }
            try self.finishLine(false);
        }
        const missing_line = self.machine.missingLine() orelse return;
        self.storeError(
            .s004_truncated_record,
            truncatedMessage(missing_line),
            missing_line,
            self.byte_offset,
        );
        return error.Format;
    }

    pub fn takeLastError(self: *CheckScanner) ?ParseError {
        const err = self.last_error;
        self.last_error = null;
        return err;
    }

    fn consumeLineBytes(self: *CheckScanner, bytes: []const u8) CheckScannerError!void {
        if (bytes.len == 0) return;

        if (self.pending_cr) {
            self.pending_cr = false;
            self.classifyByte('\r', self.line_len - 1);
        }

        const previous_len = self.line_len;
        self.line_len = std.math.add(usize, self.line_len, bytes.len) catch
            return error.LineTooLong;
        if (self.first_byte == null) self.first_byte = bytes[0];

        const semantic_end = if (bytes[bytes.len - 1] == '\r') bytes.len - 1 else bytes.len;
        self.pending_cr = semantic_end != bytes.len;
        const content_len = self.line_len - @intFromBool(self.pending_cr);
        if (content_len > self.max_line_bytes) return error.LineTooLong;
        try self.classifyBytes(bytes[0..semantic_end], previous_len);
    }

    fn classifyBytes(
        self: *CheckScanner,
        bytes: []const u8,
        start_index: usize,
    ) CheckScannerError!void {
        switch (self.machine.expected) {
            .sequence => {
                if (self.sequence_failure != null) return;
                const relative = firstInvalidSequence(bytes, self.alphabet) orelse return;
                self.sequence_failure = std.math.add(usize, start_index, relative) catch
                    return error.ArithmeticLimit;
            },
            .quality => {
                if (self.quality_failure != null) return;
                const relative = firstInvalidQuality(bytes) orelse return;
                self.quality_failure = std.math.add(usize, start_index, relative) catch
                    return error.ArithmeticLimit;
            },
            .header, .plus => {},
        }
    }

    fn classifyByte(self: *CheckScanner, byte: u8, byte_index: usize) void {
        switch (self.machine.expected) {
            .sequence => if (self.sequence_failure == null and
                !alphabetAccepts(self.alphabet, byte))
            {
                self.sequence_failure = byte_index;
            },
            .quality => if (self.quality_failure == null and
                (byte < 33 or byte > 126))
            {
                self.quality_failure = byte_index;
            },
            .header, .plus => {},
        }
    }

    fn finishLine(self: *CheckScanner, terminated: bool) CheckScannerError!void {
        const content_len = self.line_len - @intFromBool(terminated and self.pending_cr);
        const content_first = if (content_len == 0) null else self.first_byte;
        const line_kind = self.machine.expected;
        const record_ready = self.machine.push(content_len, content_first) catch |err| {
            const details = diagnostic(err);
            self.storeError(
                details.code,
                details.message,
                details.line,
                self.line_start_offset,
            );
            return error.Format;
        };

        switch (line_kind) {
            .sequence => self.sequence_start_offset = self.line_start_offset,
            .quality => self.quality_start_offset = self.line_start_offset,
            .header, .plus => {},
        }
        if (record_ready) try self.finishRecord();

        self.line_len = 0;
        self.first_byte = null;
        self.pending_cr = false;
        self.line_start_offset = self.byte_offset;
    }

    fn finishRecord(self: *CheckScanner) CheckScannerError!void {
        if (self.sequence_failure) |relative_offset| {
            try self.storeSemanticError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                2,
                self.sequence_start_offset,
                relative_offset,
            );
            return error.Format;
        }
        if (self.quality_failure) |relative_offset| {
            try self.storeSemanticError(
                .s006_invalid_quality_range,
                "quality byte must be ASCII 33 through 126",
                4,
                self.quality_start_offset,
                relative_offset,
            );
            return error.Format;
        }
        self.record_index = std.math.add(u64, self.record_index, 1) catch
            return error.ArithmeticLimit;
        self.sequence_failure = null;
        self.quality_failure = null;
    }

    fn storeSemanticError(
        self: *CheckScanner,
        code: LintCode,
        message: []const u8,
        line: u3,
        field_offset: u64,
        relative_offset: usize,
    ) CheckScannerError!void {
        const relative_u64 = std.math.cast(u64, relative_offset) orelse
            return error.ArithmeticLimit;
        const byte_offset = std.math.add(u64, field_offset, relative_u64) catch
            return error.ArithmeticLimit;
        self.storeError(code, message, line, byte_offset);
    }

    fn storeError(
        self: *CheckScanner,
        code: LintCode,
        message: []const u8,
        line: u3,
        byte_offset: u64,
    ) void {
        self.last_error = .{
            .code = code,
            .message = message,
            .record_index = self.record_index,
            .byte_offset = byte_offset,
            .line_in_record = line,
        };
    }

    fn advanceOffset(self: *CheckScanner, amount: usize) CheckScannerError!void {
        const amount_u64 = std.math.cast(u64, amount) orelse return error.ArithmeticLimit;
        self.byte_offset = std.math.add(u64, self.byte_offset, amount_u64) catch
            return error.ArithmeticLimit;
    }
};

const CheckTestOutcome = union(enum) {
    valid: u64,
    parse_error: ParseError,
    line_too_long,
};

const ReaderSpillFixture = struct {
    fn init(
        allocator: std.mem.Allocator,
        sequence_len: usize,
        quality_len: usize,
        ending: []const u8,
        final_ending: bool,
    ) ![]u8 {
        const ending_count: usize = if (final_ending) 4 else 3;
        const total_len = 4 + sequence_len + 1 + quality_len + ending_count * ending.len;
        const input = try allocator.alloc(u8, total_len);
        var offset: usize = 0;
        @memcpy(input[offset..][0..4], "@abc");
        offset += 4;
        @memcpy(input[offset..][0..ending.len], ending);
        offset += ending.len;
        @memset(input[offset..][0..sequence_len], 'A');
        offset += sequence_len;
        @memcpy(input[offset..][0..ending.len], ending);
        offset += ending.len;
        input[offset] = '+';
        offset += 1;
        @memcpy(input[offset..][0..ending.len], ending);
        offset += ending.len;
        @memset(input[offset..][0..quality_len], 'I');
        offset += quality_len;
        if (final_ending) {
            @memcpy(input[offset..][0..ending.len], ending);
            offset += ending.len;
        }
        std.debug.assert(offset == input.len);
        return input;
    }
};

fn readSpillForAllocationCheck(allocator: std.mem.Allocator, input: []const u8) !void {
    var source = io_layer.SliceSource.init(input);
    var reader = try Reader.init(allocator, source.byteSource(), .{});
    defer reader.deinit();
    _ = try reader.next();
}

test "[failure] - [reader]: partial fallback ownership is released after allocation failure" {
    const input = try ReaderSpillFixture.init(
        std.testing.allocator,
        512 * 1024,
        512 * 1024,
        "\n",
        true,
    );
    defer std.testing.allocator.free(input);

    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        readSpillForAllocationCheck,
        .{input},
    );
}

test "[edge] - [reader]: spill field capacities stop at the line limit" {
    const max_line_bytes = 80 * 1024;
    const storage_limit = fieldStorageLimit(max_line_bytes);
    const fill_bytes = [_]u8{ 'h', 'A', 'p', 'I' };
    const prefix_bytes = [_]?u8{ '@', null, '+', null };

    for ([_]bool{ false, true }) |crlf| {
        const ending_bytes: usize = if (crlf) 2 else 1;
        const input = try std.testing.allocator.alloc(
            u8,
            4 * (max_line_bytes + ending_bytes),
        );
        defer std.testing.allocator.free(input);

        var offset: usize = 0;
        for (fill_bytes, prefix_bytes) |fill, prefix| {
            @memset(input[offset..][0..max_line_bytes], fill);
            if (prefix) |byte| input[offset] = byte;
            offset += max_line_bytes;
            if (crlf) {
                input[offset] = '\r';
                offset += 1;
            }
            input[offset] = '\n';
            offset += 1;
        }
        try std.testing.expectEqual(input.len, offset);

        var source = io_layer.SliceSource.init(input);
        var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        var reader = try Reader.init(
            tracking.allocator(),
            source.byteSource(),
            .{ .max_line_bytes = max_line_bytes },
        );
        defer reader.deinit();

        const record = (try reader.next()).?;
        try std.testing.expectEqual(max_line_bytes - 1, record.header.len);
        try std.testing.expectEqual(max_line_bytes, record.sequence.len);
        try std.testing.expectEqual(max_line_bytes - 1, record.plus.len);
        try std.testing.expectEqual(max_line_bytes, record.quality.len);
        try std.testing.expect(std.mem.allEqual(u8, record.sequence, 'A'));
        try std.testing.expect(std.mem.allEqual(u8, record.quality, 'I'));
        for (&reader.fallback_fields, 0..) |*field, field_index| {
            try std.testing.expectEqual(max_line_bytes, field.len);
            const expected_capacity = if (field_index == 1)
                io_layer.DEFAULT_READER_BUFFER_BYTES
            else if (crlf or field_index == 3)
                storage_limit
            else
                max_line_bytes;
            try std.testing.expectEqual(expected_capacity, field.storage.len);
        }
        try std.testing.expectError(
            ReaderError.LineTooLong,
            reader.ensureFieldCapacity(0, storage_limit + 1, false),
        );
    }
}

test "[edge] - [reader]: field storage limit saturates after checked overflow" {
    const max = std.math.maxInt(usize);

    try std.testing.expectEqual(1024 * 1024 + 1, fieldStorageLimit(1024 * 1024));
    try std.testing.expectEqual(max, fieldStorageLimit(max - 1));
    try std.testing.expectEqual(max, fieldStorageLimit(max));
}

test "[failure] - [reader]: failed spill growth preserves owned storage" {
    var source = io_layer.SliceSource.init("");
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 2,
        .resize_fail_index = 0,
    });
    var reader = try Reader.init(
        failing.allocator(),
        source.byteSource(),
        .{ .max_line_bytes = 80 * 1024 },
    );
    defer reader.deinit();

    try std.testing.expectError(
        ReaderError.OutOfMemory,
        reader.ensureFieldCapacity(0, 1, false),
    );
    try std.testing.expectEqual(@as(usize, 0), reader.fallback_fields[0].storage.len);
}

test "[edge] - [reader]: quality spill reserves only its valid tail" {
    const sequence_len = 131_070;
    const Case = struct {
        quality_len: usize,
        ending: []const u8,
        final_ending: bool,
        expected_error: ?ReaderError,
        expected_quality_capacity: usize,
    };
    const cases = [_]Case{
        .{
            .quality_len = sequence_len,
            .ending = "\n",
            .final_ending = true,
            .expected_error = null,
            .expected_quality_capacity = sequence_len + 1,
        },
        .{
            .quality_len = sequence_len,
            .ending = "\r\n",
            .final_ending = true,
            .expected_error = null,
            .expected_quality_capacity = sequence_len + 1,
        },
        .{
            .quality_len = sequence_len,
            .ending = "\n",
            .final_ending = false,
            .expected_error = null,
            .expected_quality_capacity = sequence_len + 1,
        },
        .{
            .quality_len = sequence_len - 1,
            .ending = "\n",
            .final_ending = true,
            .expected_error = error.S005LengthMismatch,
            .expected_quality_capacity = sequence_len + 1,
        },
        .{
            .quality_len = sequence_len + 1,
            .ending = "\n",
            .final_ending = true,
            .expected_error = error.S005LengthMismatch,
            .expected_quality_capacity = sequence_len + 1,
        },
        .{
            .quality_len = sequence_len + 2,
            .ending = "\n",
            .final_ending = true,
            .expected_error = error.S005LengthMismatch,
            .expected_quality_capacity = sequence_len + 2,
        },
    };

    for (cases) |case| {
        const input = try ReaderSpillFixture.init(
            std.testing.allocator,
            sequence_len,
            case.quality_len,
            case.ending,
            case.final_ending,
        );
        defer std.testing.allocator.free(input);
        var source = io_layer.SliceSource.init(input);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();

        if (case.expected_error) |expected_error| {
            try std.testing.expectError(expected_error, reader.next());
        } else {
            const record = (try reader.next()).?;
            try std.testing.expectEqual(sequence_len, record.sequence.len);
            try std.testing.expectEqual(sequence_len, record.quality.len);
        }
        try std.testing.expectEqual(
            case.expected_quality_capacity,
            reader.fallback_fields[3].storage.len,
        );
    }

    for ([_]usize{ 512 * 1024, 1024 * 1024 }) |long_sequence_len| {
        const input = try ReaderSpillFixture.init(
            std.testing.allocator,
            long_sequence_len,
            long_sequence_len,
            "\n",
            true,
        );
        defer std.testing.allocator.free(input);
        var source = io_layer.SliceSource.init(input);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();

        const record = (try reader.next()).?;
        try std.testing.expectEqual(long_sequence_len, record.sequence.len);
        try std.testing.expectEqual(long_sequence_len, record.quality.len);
        try std.testing.expectEqual(@as(usize, 4), reader.fallback_fields[0].storage.len);
        try std.testing.expectEqual(long_sequence_len, reader.fallback_fields[1].storage.len);
        try std.testing.expectEqual(@as(usize, 1), reader.fallback_fields[2].storage.len);
        try std.testing.expectEqual(
            long_sequence_len + 1,
            reader.fallback_fields[3].storage.len,
        );
    }
}

test "[failure] - [reader]: failed quality-tail reserve preserves owned storage" {
    const sequence_len = 131_070;
    const input = try ReaderSpillFixture.init(
        std.testing.allocator,
        sequence_len,
        sequence_len,
        "\n",
        true,
    );
    defer std.testing.allocator.free(input);

    var source = io_layer.SliceSource.init(input);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 4,
    });
    var reader = try Reader.init(failing.allocator(), source.byteSource(), .{});
    defer reader.deinit();

    try std.testing.expectError(ReaderError.OutOfMemory, reader.next());
    try std.testing.expectEqual(@as(usize, 4), reader.fallback_fields[0].storage.len);
    try std.testing.expectEqual(
        io_layer.DEFAULT_READER_BUFFER_BYTES,
        reader.fallback_fields[1].storage.len,
    );
    try std.testing.expectEqual(@as(usize, 1), reader.fallback_fields[2].storage.len);
    try std.testing.expectEqual(@as(usize, 0), reader.fallback_fields[3].storage.len);
    try std.testing.expect(failing.has_induced_failure);
}

test "[property] - [record delivery]: omits identifiers and preserves buffered canonical spans" {
    const cases = [_]struct {
        input: []const u8,
        expected_output: []const u8,
        has_span: bool,
    }{
        .{
            .input = "@r\rdesc\nACGT\n+note\rx\n!#$%\n",
            .expected_output = "@r\rdesc\nACGT\n+note\rx\n!#$%\n",
            .has_span = true,
        },
        .{
            .input = "@crlf\r\nACGT\r\n+\r\n!!!!\r\n",
            .expected_output = "@crlf\nACGT\n+\n!!!!\n",
            .has_span = false,
        },
        .{
            .input = "@eof\nA\n+\n!",
            .expected_output = "@eof\nA\n+\n!\n",
            .has_span = false,
        },
    };

    for (cases) |case| {
        var source = io_layer.SliceSource.init(case.input);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        var canonical_span: ?[]const u8 = null;
        const record = (try nextWithoutId(&reader, &canonical_span)).?;
        try std.testing.expect(validateRecord(record, .{}) == null);
        try std.testing.expectEqual(@as(usize, 0), record.id.len);
        try std.testing.expectEqual(case.has_span, canonical_span != null);

        var output: [64]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = Writer.init(sink.byteSink());
        if (canonical_span) |span| {
            try writeCanonicalRecordSpan(&writer, span);
        } else {
            try writeValidatedRecord(&writer, record);
        }
        try std.testing.expectEqualStrings(case.expected_output, sink.written());
        try std.testing.expect((try nextWithoutId(&reader, &canonical_span)) == null);
        try std.testing.expect(canonical_span == null);
    }

    {
        const input1 = "@pair/1\nAC\n+left\n!!\n";
        const input2 = "@pair/2\nGT\n+right\n##\n";
        const input = input1 ++ input2;
        var source = io_layer.SliceSource.init(input);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        var canonical_span1: ?[]const u8 = null;
        const record1 = (try nextWithoutId(&reader, &canonical_span1)).?;
        var canonical_span2: ?[]const u8 = null;
        const record2 = (try nextBufferedWithoutId(&reader, &canonical_span2)).?;

        try std.testing.expectEqualStrings("pair/1", record1.header);
        try std.testing.expectEqualStrings("AC", record1.sequence);
        try std.testing.expectEqualStrings("left", record1.plus);
        try std.testing.expectEqualStrings("!!", record1.quality);
        try std.testing.expectEqualStrings("pair/2", record2.header);
        try std.testing.expectEqualStrings(input1, canonical_span1.?);
        try std.testing.expectEqualStrings(input2, canonical_span2.?);
    }

    {
        const input1 = "@pair/1\nAC\n+left\n!!\n";
        const record2_prefix = "@pair/2\nGT\n+right\n";
        const record2_suffix = "##\n";
        const initial_input = input1 ++ record2_prefix;
        const complete_input = initial_input ++ record2_suffix;
        var source = io_layer.SliceSource.init(initial_input);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        var canonical_span1: ?[]const u8 = null;
        const record1 = (try nextWithoutId(&reader, &canonical_span1)).?;
        source.data = complete_input;
        const source_position = source.pos;
        var canonical_span2: ?[]const u8 = null;

        try std.testing.expect(
            (try nextBufferedWithoutId(&reader, &canonical_span2)) == null,
        );
        try std.testing.expectEqual(source_position, source.pos);
        try std.testing.expectEqualStrings("pair/1", record1.header);
        try std.testing.expectEqualStrings(input1, canonical_span1.?);

        const record2 = (try nextWithoutId(&reader, &canonical_span2)).?;
        try std.testing.expectEqualStrings("pair/2", record2.header);
        try std.testing.expectEqualStrings("GT", record2.sequence);
        try std.testing.expectEqualStrings("right", record2.plus);
        try std.testing.expectEqualStrings("##", record2.quality);
        try std.testing.expect(canonical_span2 == null);
    }

    const header_len = io_layer.DEFAULT_READER_BUFFER_BYTES;
    const input = try std.testing.allocator.alloc(u8, header_len + 8);
    defer std.testing.allocator.free(input);
    @memset(input, 'h');
    input[0] = '@';
    input[header_len + 1] = '\n';
    input[header_len + 2] = 'A';
    input[header_len + 3] = '\n';
    input[header_len + 4] = '+';
    input[header_len + 5] = '\n';
    input[header_len + 6] = '!';
    input[header_len + 7] = '\n';

    var source = io_layer.SliceSource.init(input);
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    var canonical_span: ?[]const u8 = null;
    const record = (try nextWithoutId(&reader, &canonical_span)).?;
    try std.testing.expectEqual(header_len, record.header.len);
    try std.testing.expectEqual(@as(usize, 0), record.id.len);
    try std.testing.expect(canonical_span == null);
}

test "[integration] - [writer]: trusted Reader records match checked serialization" {
    const input = "@r\rdesc\r\nACGT\r\n+note\rx\r\n!#$%\r\n";
    const expected = "@r\rdesc\nACGT\n+note\rx\n!#$%\n";
    var source = io_layer.SliceSource.init(input);
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    const record = (try reader.next()).?;
    try std.testing.expect(validateRecord(record, .{}) == null);

    var checked_bytes: [64]u8 = undefined;
    var checked_sink = io_layer.SliceSink.init(&checked_bytes);
    var checked_writer = Writer.init(checked_sink.byteSink());
    try checked_writer.writeRecord(record);

    var trusted_bytes: [64]u8 = undefined;
    var trusted_sink = io_layer.SliceSink.init(&trusted_bytes);
    var trusted_writer = Writer.init(trusted_sink.byteSink());
    try writeValidatedRecord(&trusted_writer, record);

    try std.testing.expectEqualStrings(expected, checked_sink.written());
    try std.testing.expectEqualStrings(checked_sink.written(), trusted_sink.written());
}

test "[unit] - [check scanner]: state remains fixed-size" {
    try std.testing.expect(@sizeOf(CheckScanner) <= 160);
}

test "[property] - [check scanner]: structural masks preserve newline positions" {
    var block: [STRUCTURAL_BLOCK_BYTES]u8 = @splat('A');
    for (0..STRUCTURAL_BLOCK_BYTES) |lane| {
        block[lane] = '\n';
        try std.testing.expectEqual(@as(u64, 1) << @intCast(lane), newlineMask(&block));
        block[lane] = 'A';
    }

    var expected: u64 = 0;
    for ([_]usize{ 0, 1, 15, 16, 31, 32, 47, 48, 62, 63 }) |lane| {
        block[lane] = '\n';
        expected |= @as(u64, 1) << @intCast(lane);
    }
    try std.testing.expectEqual(expected, newlineMask(&block));
}

test "[property] - [check scanner]: fragmented results match independent expectations and Reader" {
    const Case = struct {
        data: []const u8,
        options: Options = .{},
        validation_options: ValidationOptions = .{},
        expected: CheckTestOutcome,
    };
    const cases = [_]Case{
        .{ .data = "", .expected = .{ .valid = 0 } },
        .{ .data = "@r\nA\n+\n!\n", .expected = .{ .valid = 1 } },
        .{ .data = "@r\r\nA\r\n+\r\n!\r\n", .expected = .{ .valid = 1 } },
        .{ .data = "@r\nA\n+\n!", .expected = .{ .valid = 1 } },
        .{ .data = "@r\nR\n+\n!\n", .expected = .{ .valid = 1 } },
        .{
            .data = "@r\nR\n+\n!\n",
            .validation_options = .{ .alphabet = .acgtn },
            .expected = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                0,
                3,
                2,
            ) },
        },
        .{
            .data = "r\nA\n+\n!\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s003_invalid_header,
                "header line must start with '@'",
                0,
                0,
                1,
            ) },
        },
        .{
            .data = "@r\n.\nx\n!\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s001_invalid_plus_line,
                "plus line must start with '+'",
                0,
                5,
                3,
            ) },
        },
        .{
            .data = "@r\n.\n+\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s004_truncated_record,
                "unexpected end of file in quality line",
                0,
                7,
                4,
            ) },
        },
        .{
            .data = "@r\n.\n+\n!!\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s005_length_mismatch,
                "sequence and quality lengths differ",
                0,
                7,
                4,
            ) },
        },
        .{
            .data = "@r\n.\n+\n\x7f\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                0,
                3,
                2,
            ) },
        },
        .{
            .data = "@r\nA\n+\n\r",
            .expected = .{ .parse_error = expectedCheckError(
                .s006_invalid_quality_range,
                "quality byte must be ASCII 33 through 126",
                0,
                7,
                4,
            ) },
        },
        .{
            .data = "@r\nAA\r\n+\r\n!!\r\n",
            .options = .{ .max_line_bytes = 2 },
            .expected = .{ .valid = 1 },
        },
        .{
            .data = "@r\nAAA\n+\n!!!\n",
            .options = .{ .max_line_bytes = 2 },
            .expected = .line_too_long,
        },
        .{
            .data = "@r\nAA\n+\n!!\r",
            .options = .{ .max_line_bytes = 2 },
            .expected = .line_too_long,
        },
        .{
            .data = "@a\nA\n+\n!\n@b\nT\n+\n~\n",
            .expected = .{ .valid = 2 },
        },
    };

    for (cases) |case| {
        const reference = try referenceCheckOutcome(
            case.data,
            case.options,
            case.validation_options,
        );
        try expectCheckOutcome(case.expected, reference);
        for (1..case.data.len + 2) |chunk_len| {
            const direct = directCheckOutcome(
                case.data,
                chunk_len,
                case.options,
                case.validation_options,
            );
            try expectCheckOutcome(case.expected, direct);
        }
    }
}

test "[property] - [check scanner]: structural block boundaries match Reader" {
    for ([_]usize{ 62, 63, 64, 65, 126, 127, 128, 129 }) |field_len| {
        var data: std.ArrayList(u8) = .empty;
        defer data.deinit(std.testing.allocator);
        try data.appendSlice(std.testing.allocator, "@record\n");
        try data.appendNTimes(std.testing.allocator, 'A', field_len);
        try data.appendSlice(std.testing.allocator, "\n+description\n");
        try data.appendNTimes(std.testing.allocator, '!', field_len);
        try data.append(std.testing.allocator, '\n');

        const expected = try referenceCheckOutcome(data.items, .{}, .{});
        try expectCheckOutcome(.{ .valid = 1 }, expected);
        for ([_]usize{ 63, 64, 65, 127, 128, 129 }) |chunk_len| {
            try expectCheckOutcome(
                expected,
                directCheckOutcome(data.items, chunk_len, .{}, .{}),
            );
        }
        try expectCheckOutcome(
            expected,
            directCheckOutcome(data.items, data.items.len, .{}, .{}),
        );
    }
}

test "[property] - [check scanner]: refill-spanning records match Reader" {
    const field_len = io_layer.DEFAULT_READER_BUFFER_BYTES + 17;
    const data_len = 3 + field_len + 3 + field_len + 1;
    const data = try std.testing.allocator.alloc(u8, data_len);
    defer std.testing.allocator.free(data);

    var cursor: usize = 0;
    @memcpy(data[cursor..][0..3], "@r\n");
    cursor += 3;
    @memset(data[cursor..][0..field_len], 'A');
    cursor += field_len;
    @memcpy(data[cursor..][0..3], "\n+\n");
    cursor += 3;
    @memset(data[cursor..][0..field_len], '!');
    cursor += field_len;
    data[cursor] = '\n';

    const expected: CheckTestOutcome = .{ .valid = 1 };
    try expectCheckOutcome(expected, try referenceCheckOutcome(data, .{}, .{}));
    for ([_]usize{
        io_layer.DEFAULT_READER_BUFFER_BYTES - 1,
        io_layer.DEFAULT_READER_BUFFER_BYTES,
        io_layer.DEFAULT_READER_BUFFER_BYTES + 1,
    }) |chunk_len| {
        try expectCheckOutcome(expected, directCheckOutcome(data, chunk_len, .{}, .{}));
    }

    const invalid_index = io_layer.DEFAULT_READER_BUFFER_BYTES;
    data[3 + invalid_index] = '.';
    const invalid: CheckTestOutcome = .{ .parse_error = expectedCheckError(
        .s002_invalid_sequence_alphabet,
        "sequence byte is outside the selected alphabet",
        0,
        3 + invalid_index,
        2,
    ) };
    try expectCheckOutcome(invalid, try referenceCheckOutcome(data, .{}, .{}));
    try expectCheckOutcome(
        invalid,
        directCheckOutcome(data, io_layer.DEFAULT_READER_BUFFER_BYTES, .{}, .{}),
    );
}

test "[property] - [check scanner]: generated semantic mutations retain exact locations" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse 16;
    for ([_]usize{
        1,
        vector_len - 1,
        vector_len,
        vector_len + 1,
        2 * vector_len + 1,
    }) |field_len| {
        var data: std.ArrayList(u8) = .empty;
        defer data.deinit(std.testing.allocator);
        try data.appendSlice(std.testing.allocator, "@ok\nA\n+\n!\n@r\n");
        const sequence_start = data.items.len;
        try data.appendNTimes(std.testing.allocator, 'A', field_len);
        try data.appendSlice(std.testing.allocator, "\n+\n");
        const quality_start = data.items.len;
        try data.appendNTimes(std.testing.allocator, '!', field_len);
        try data.append(std.testing.allocator, '\n');

        for (0..field_len) |invalid_index| {
            data.items[sequence_start + invalid_index] = '.';
            const sequence_error: CheckTestOutcome = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                1,
                sequence_start + invalid_index,
                2,
            ) };
            try expectCheckOutcome(
                sequence_error,
                try referenceCheckOutcome(data.items, .{}, .{}),
            );
            for ([_]usize{ 1, vector_len, vector_len + 1 }) |chunk_len| {
                try expectCheckOutcome(
                    sequence_error,
                    directCheckOutcome(data.items, chunk_len, .{}, .{}),
                );
            }
            data.items[sequence_start + invalid_index] = 'A';

            data.items[quality_start + invalid_index] = 127;
            const quality_error: CheckTestOutcome = .{ .parse_error = expectedCheckError(
                .s006_invalid_quality_range,
                "quality byte must be ASCII 33 through 126",
                1,
                quality_start + invalid_index,
                4,
            ) };
            try expectCheckOutcome(
                quality_error,
                try referenceCheckOutcome(data.items, .{}, .{}),
            );
            for ([_]usize{ 1, vector_len, vector_len + 1 }) |chunk_len| {
                try expectCheckOutcome(
                    quality_error,
                    directCheckOutcome(data.items, chunk_len, .{}, .{}),
                );
            }
            data.items[quality_start + invalid_index] = '!';
        }
    }
}

test "[edge] - [check scanner]: arithmetic limits fail explicitly" {
    var offset = CheckScanner.init(.{}, .{});
    offset.byte_offset = std.math.maxInt(u64);
    try std.testing.expectError(error.ArithmeticLimit, offset.feed("A"));

    var records = CheckScanner.init(.{}, .{});
    records.record_index = std.math.maxInt(u64);
    try std.testing.expectError(error.ArithmeticLimit, records.feed("@r\nA\n+\n!\n"));
}

fn directCheckOutcome(
    data: []const u8,
    chunk_len: usize,
    options: Options,
    validation_options: ValidationOptions,
) CheckTestOutcome {
    var scanner = CheckScanner.init(options, validation_options);
    var cursor: usize = 0;
    while (cursor < data.len) {
        const end = cursor + @min(chunk_len, data.len - cursor);
        _ = scanner.feed(data[cursor..end]) catch |err| return scannerErrorOutcome(
            &scanner,
            err,
        );
        cursor = end;
    }
    scanner.finishEof() catch |err| return scannerErrorOutcome(&scanner, err);
    return .{ .valid = scanner.record_index };
}

fn referenceCheckOutcome(
    data: []const u8,
    options: Options,
    validation_options: ValidationOptions,
) !CheckTestOutcome {
    var source = io_layer.SliceSource.init(data);
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), options);
    defer reader.deinit();

    while (reader.next() catch |err| return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => .{ .parse_error = reader.takeLastError().? },
        error.LineTooLong => .line_too_long,
        error.OutOfMemory, error.Io => return err,
    }) |record| {
        const semantic_error = validateRecord(record, validation_options) orelse continue;
        const offsets = reader.currentRecordOffsets().?;
        const field_offset = switch (semantic_error.field) {
            .sequence => offsets.sequence,
            .quality => offsets.quality,
        };
        return .{ .parse_error = .{
            .code = semantic_error.code,
            .message = semantic_error.message,
            .record_index = reader.recordIndex() - 1,
            .byte_offset = field_offset + semantic_error.byte_index,
            .line_in_record = switch (semantic_error.field) {
                .sequence => 2,
                .quality => 4,
            },
        } };
    }
    return .{ .valid = reader.recordIndex() };
}

fn scannerErrorOutcome(
    scanner: *CheckScanner,
    err: CheckScannerError,
) CheckTestOutcome {
    return switch (err) {
        error.Format => .{ .parse_error = scanner.takeLastError().? },
        error.LineTooLong => .line_too_long,
        error.ArithmeticLimit => unreachable,
    };
}

fn expectedCheckError(
    code: LintCode,
    message: []const u8,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u3,
) ParseError {
    return .{
        .code = code,
        .message = message,
        .record_index = record_index,
        .byte_offset = byte_offset,
        .line_in_record = line_in_record,
    };
}

fn expectCheckOutcome(expected: CheckTestOutcome, actual: CheckTestOutcome) !void {
    try std.testing.expectEqual(std.meta.activeTag(expected), std.meta.activeTag(actual));
    switch (expected) {
        .valid => |count| try std.testing.expectEqual(count, actual.valid),
        .parse_error => |details| try std.testing.expectEqualDeep(details, actual.parse_error),
        .line_too_long => {},
    }
}

test "[property] - [record validation]: vector quality validation matches scalar results" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse 16;
    const max_len = 2 * vector_len + 1;
    const quality = try std.testing.allocator.alloc(u8, max_len);
    defer std.testing.allocator.free(quality);

    @memset(quality, '!');
    for ([_]usize{
        0,
        1,
        vector_len - 1,
        vector_len,
        vector_len + 1,
        2 * vector_len - 1,
        2 * vector_len,
        2 * vector_len + 1,
    }) |length| {
        try std.testing.expectEqual(null, firstInvalidQualityScalar(quality[0..length], 0));
        try std.testing.expectEqual(null, firstInvalidQuality(quality[0..length]));
    }

    for (0..256) |value| {
        @memset(quality, @intCast(value));
        const expected: ?usize = if (value < 33 or value > 126) 0 else null;
        try std.testing.expectEqual(expected, firstInvalidQualityScalar(quality, 0));
        try std.testing.expectEqual(expected, firstInvalidQuality(quality));
    }

    @memset(quality, '!');
    for (0..quality.len) |invalid_index| {
        quality[invalid_index] = if (invalid_index % 2 == 0) 32 else 127;
        try std.testing.expectEqual(
            invalid_index,
            firstInvalidQualityScalar(quality, 0).?,
        );
        try std.testing.expectEqual(invalid_index, firstInvalidQuality(quality).?);
        quality[invalid_index] = '!';
    }

    quality[vector_len] = 127;
    quality[1] = 32;
    try std.testing.expectEqual(@as(usize, 1), firstInvalidQualityScalar(quality, 0).?);
    try std.testing.expectEqual(@as(usize, 1), firstInvalidQuality(quality).?);
}

test "[property] - [record validation]: vector sequence validation matches scalar policies" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse 16;
    const max_len = 2 * vector_len + 1;
    const sequence = try std.testing.allocator.alloc(u8, max_len);
    defer std.testing.allocator.free(sequence);

    for ([_]Alphabet{ .iupac, .acgtn }) |alphabet| {
        @memset(sequence, 'A');
        for ([_]usize{
            0,
            1,
            vector_len - 1,
            vector_len,
            vector_len + 1,
            2 * vector_len - 1,
            2 * vector_len,
            2 * vector_len + 1,
        }) |length| {
            try std.testing.expectEqual(
                null,
                firstInvalidSequenceScalar(sequence[0..length], alphabet, 0),
            );
            try std.testing.expectEqual(null, firstInvalidSequence(sequence[0..length], alphabet));
        }

        for (0..256) |value| {
            @memset(sequence, @intCast(value));
            const expected: ?usize = if (alphabetAccepts(alphabet, @intCast(value))) null else 0;
            try std.testing.expectEqual(
                expected,
                firstInvalidSequenceScalar(sequence, alphabet, 0),
            );
            try std.testing.expectEqual(expected, firstInvalidSequence(sequence, alphabet));
        }

        @memset(sequence, 'A');
        for (0..sequence.len) |invalid_index| {
            sequence[invalid_index] = '.';
            try std.testing.expectEqual(
                invalid_index,
                firstInvalidSequenceScalar(sequence, alphabet, 0).?,
            );
            try std.testing.expectEqual(
                invalid_index,
                firstInvalidSequence(sequence, alphabet).?,
            );
            sequence[invalid_index] = 'A';
        }

        sequence[vector_len] = '.';
        sequence[1] = 0xff;
        try std.testing.expectEqual(
            @as(usize, 1),
            firstInvalidSequenceScalar(sequence, alphabet, 0).?,
        );
        try std.testing.expectEqual(@as(usize, 1), firstInvalidSequence(sequence, alphabet).?);
    }
}

test "[property] - [record validation]: fused vectors preserve field precedence" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse return;
    const length = 2 * vector_len + 1;
    const sequence = try std.testing.allocator.alloc(u8, length);
    defer std.testing.allocator.free(sequence);
    const quality = try std.testing.allocator.alloc(u8, length);
    defer std.testing.allocator.free(quality);
    @memset(sequence, 'A');
    @memset(quality, '!');

    const record = Record{
        .header = "record",
        .id = "record",
        .sequence = sequence,
        .plus = "",
        .quality = quality,
    };
    try std.testing.expect(validateRecord(record, .{}) == null);

    quality[1] = 32;
    sequence[vector_len] = '.';
    const sequence_error = validateRecord(record, .{}).?;
    try std.testing.expectEqual(SemanticField.sequence, sequence_error.field);
    try std.testing.expectEqual(vector_len, sequence_error.byte_index);

    sequence[vector_len] = 'A';
    const quality_error = validateRecord(record, .{}).?;
    try std.testing.expectEqual(SemanticField.quality, quality_error.field);
    try std.testing.expectEqual(@as(usize, 1), quality_error.byte_index);
}
