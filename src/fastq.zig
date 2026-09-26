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
    ArithmeticLimit,
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

pub const RecordPayload = struct {
    sequence: []const u8,
    quality: []const u8,
};

pub const PredictedPayload = struct {
    payload: RecordPayload,
    checkpoint: ?PayloadCheckpoint = null,
};

const PayloadCheckpoint = struct {
    cursor: usize,
    byte_offset: u64,
    record_index: u64,
    machine: Machine,
    record_offsets: RecordOffsets,

    pub fn reread(self: PayloadCheckpoint, reader: *Reader) ReaderError!RecordPayload {
        reader.cursor = self.cursor;
        reader.byte_offset = self.byte_offset;
        reader.record_index = self.record_index;
        reader.machine = self.machine;
        reader.record_offsets = self.record_offsets;
        return (try reader.nextPayload()).?;
    }
};

pub const ValidatedHeader = struct {
    header: []const u8,
    semantic_error: ?SemanticError,
};

pub const ValidatedRecord = struct {
    record: Record,
    canonical_span: ?[]const u8,
    semantic_error: ?SemanticError,
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
    return switch (options.alphabet) {
        .iupac => validateRecordFor(.iupac, record),
        .acgtn => validateRecordFor(.acgtn, record),
    };
}

fn validateRecordFor(comptime alphabet: Alphabet, record: Record) ?SemanticError {
    if (record.sequence.len == record.quality.len) {
        if (std.simd.suggestVectorLength(u8)) |vector_len| {
            return validateRecordVector(vector_len, alphabet, record);
        }
    }
    if (firstInvalidSequence(record.sequence, alphabet)) |byte_index| {
        return semanticSequenceError(byte_index);
    }
    if (firstInvalidQuality(record.quality)) |byte_index| {
        return semanticQualityError(byte_index);
    }
    return null;
}

pub const AdaptiveRecordValidator = struct {
    alphabet: Alphabet,
    use_full_iupac: bool = false,

    pub fn init(options: ValidationOptions) AdaptiveRecordValidator {
        return .{ .alphabet = options.alphabet };
    }

    pub fn validate(self: *AdaptiveRecordValidator, record: Record) ?SemanticError {
        if (self.alphabet == .acgtn) return validateRecordFor(.acgtn, record);
        if (self.use_full_iupac) return validateRecordFor(.iupac, record);

        const narrow_error = validateRecordFor(.acgtn, record) orelse return null;
        if (narrow_error.field == .quality) return narrow_error;

        const full_error = validateRecordFor(.iupac, record);
        if (full_error == null or full_error.?.field == .quality) {
            self.use_full_iupac = true;
        }
        return full_error;
    }
};

fn validateRecordVector(
    comptime vector_len: comptime_int,
    comptime alphabet: Alphabet,
    record: Record,
) ?SemanticError {
    const Bytes = @Vector(vector_len, u8);
    const Mask = @Vector(vector_len, bool);

    if (record.sequence.len < vector_len) {
        if (firstInvalidSequenceScalar(record.sequence, alphabet, 0)) |byte_index| {
            return semanticSequenceError(byte_index);
        }
        if (firstInvalidQualityScalar(record.quality, 0)) |byte_index| {
            return semanticQualityError(byte_index);
        }
        return null;
    }

    var sequence_invalid: Mask = @splat(false);
    var quality_invalid: Mask = @splat(false);
    var byte_index: usize = 0;
    while (record.sequence.len - byte_index >= 4 * vector_len) : (byte_index += 4 * vector_len) {
        const sequence0: Bytes = record.sequence[byte_index..][0..vector_len].*;
        const sequence1: Bytes = record.sequence[byte_index + vector_len ..][0..vector_len].*;
        const sequence2: Bytes = record.sequence[byte_index + 2 * vector_len ..][0..vector_len].*;
        const sequence3: Bytes = record.sequence[byte_index + 3 * vector_len ..][0..vector_len].*;
        sequence_invalid |= invalidSequenceVector(vector_len, alphabet, sequence0) |
            invalidSequenceVector(vector_len, alphabet, sequence1) |
            invalidSequenceVector(vector_len, alphabet, sequence2) |
            invalidSequenceVector(vector_len, alphabet, sequence3);

        const quality0: Bytes = record.quality[byte_index..][0..vector_len].*;
        const quality1: Bytes = record.quality[byte_index + vector_len ..][0..vector_len].*;
        const quality2: Bytes = record.quality[byte_index + 2 * vector_len ..][0..vector_len].*;
        const quality3: Bytes = record.quality[byte_index + 3 * vector_len ..][0..vector_len].*;
        quality_invalid |= invalidQualityVector(vector_len, quality0) |
            invalidQualityVector(vector_len, quality1) |
            invalidQualityVector(vector_len, quality2) |
            invalidQualityVector(vector_len, quality3);
    }
    while (record.sequence.len - byte_index >= 2 * vector_len) : (byte_index += 2 * vector_len) {
        const sequence0: Bytes = record.sequence[byte_index..][0..vector_len].*;
        const sequence1: Bytes = record.sequence[byte_index + vector_len ..][0..vector_len].*;
        sequence_invalid |= invalidSequenceVector(vector_len, alphabet, sequence0) |
            invalidSequenceVector(vector_len, alphabet, sequence1);

        const quality0: Bytes = record.quality[byte_index..][0..vector_len].*;
        const quality1: Bytes = record.quality[byte_index + vector_len ..][0..vector_len].*;
        quality_invalid |= invalidQualityVector(vector_len, quality0) |
            invalidQualityVector(vector_len, quality1);
    }
    if (record.sequence.len - byte_index >= vector_len) {
        const sequence: Bytes = record.sequence[byte_index..][0..vector_len].*;
        sequence_invalid |= invalidSequenceVector(vector_len, alphabet, sequence);
        const quality: Bytes = record.quality[byte_index..][0..vector_len].*;
        quality_invalid |= invalidQualityVector(vector_len, quality);
        byte_index += vector_len;
    }
    if (byte_index < record.sequence.len) {
        const tail_start = record.sequence.len - vector_len;
        const active = uncheckedTailMask(vector_len, byte_index - tail_start);
        const sequence: Bytes = record.sequence[tail_start..][0..vector_len].*;
        sequence_invalid |= invalidSequenceVector(vector_len, alphabet, sequence) & active;
        const quality: Bytes = record.quality[tail_start..][0..vector_len].*;
        quality_invalid |= invalidQualityVector(vector_len, quality) & active;
    }

    if (@reduce(.Or, sequence_invalid)) {
        return semanticSequenceError(
            firstInvalidSequenceScalar(record.sequence, alphabet, 0).?,
        );
    }
    if (@reduce(.Or, quality_invalid)) {
        return semanticQualityError(firstInvalidQualityScalar(record.quality, 0).?);
    }
    return null;
}

fn semanticSequenceError(byte_index: usize) SemanticError {
    return .{
        .code = .s002_invalid_sequence_alphabet,
        .message = "sequence byte is outside the selected alphabet",
        .field = .sequence,
        .byte_index = byte_index,
    };
}

pub fn semanticQualityError(byte_index: usize) SemanticError {
    return .{
        .code = .s006_invalid_quality_range,
        .message = "quality byte must be ASCII 33 through 126",
        .field = .quality,
        .byte_index = byte_index,
    };
}

pub fn semanticParseError(
    semantic_error: SemanticError,
    record_index: u64,
    field_offset: u64,
) error{ArithmeticLimit}!ParseError {
    return .{
        .code = semantic_error.code,
        .message = semantic_error.message,
        .record_index = record_index,
        .byte_offset = try progressAfter(field_offset, semantic_error.byte_index),
        .line_in_record = switch (semantic_error.field) {
            .sequence => 2,
            .quality => 4,
        },
    };
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
    const Mask = @Vector(vector_len, bool);

    if (sequence.len < vector_len) {
        return firstInvalidSequenceScalar(sequence, alphabet, 0);
    }

    var invalid: Mask = @splat(false);
    var byte_index: usize = 0;
    while (sequence.len - byte_index >= 2 * vector_len) : (byte_index += 2 * vector_len) {
        const bytes0: Bytes = sequence[byte_index..][0..vector_len].*;
        const bytes1: Bytes = sequence[byte_index + vector_len ..][0..vector_len].*;
        invalid |= invalidSequenceVector(vector_len, alphabet, bytes0) |
            invalidSequenceVector(vector_len, alphabet, bytes1);
    }
    if (sequence.len - byte_index >= vector_len) {
        const bytes: Bytes = sequence[byte_index..][0..vector_len].*;
        invalid |= invalidSequenceVector(vector_len, alphabet, bytes);
        byte_index += vector_len;
    }
    if (byte_index < sequence.len) {
        const tail_start = sequence.len - vector_len;
        const active = uncheckedTailMask(vector_len, byte_index - tail_start);
        const bytes: Bytes = sequence[tail_start..][0..vector_len].*;
        invalid |= invalidSequenceVector(vector_len, alphabet, bytes) & active;
    }
    if (!@reduce(.Or, invalid)) return null;
    return firstInvalidSequenceScalar(sequence, alphabet, 0);
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
    const Mask = @Vector(vector_len, bool);

    if (quality.len < vector_len) return firstInvalidQualityScalar(quality, 0);

    var invalid: Mask = @splat(false);
    var byte_index: usize = 0;
    while (quality.len - byte_index >= 2 * vector_len) : (byte_index += 2 * vector_len) {
        const bytes0: Bytes = quality[byte_index..][0..vector_len].*;
        const bytes1: Bytes = quality[byte_index + vector_len ..][0..vector_len].*;
        invalid |= invalidQualityVector(vector_len, bytes0) |
            invalidQualityVector(vector_len, bytes1);
    }
    if (quality.len - byte_index >= vector_len) {
        const bytes: Bytes = quality[byte_index..][0..vector_len].*;
        invalid |= invalidQualityVector(vector_len, bytes);
        byte_index += vector_len;
    }
    if (byte_index < quality.len) {
        const tail_start = quality.len - vector_len;
        const active = uncheckedTailMask(vector_len, byte_index - tail_start);
        const bytes: Bytes = quality[tail_start..][0..vector_len].*;
        invalid |= invalidQualityVector(vector_len, bytes) & active;
    }
    if (!@reduce(.Or, invalid)) return null;
    return firstInvalidQualityScalar(quality, 0);
}

fn invalidQualityVector(
    comptime vector_len: comptime_int,
    quality: @Vector(vector_len, u8),
) @Vector(vector_len, bool) {
    const Bytes = @Vector(vector_len, u8);
    return (quality < @as(Bytes, @splat(33))) | (quality > @as(Bytes, @splat(126)));
}

fn uncheckedTailMask(
    comptime vector_len: comptime_int,
    checked_prefix: usize,
) @Vector(vector_len, bool) {
    const Indexes = @Vector(vector_len, usize);
    return std.simd.iota(usize, vector_len) >= @as(Indexes, @splat(checked_prefix));
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
    const end = std.mem.findAny(u8, header, "\t ") orelse header.len;
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
    content_len: usize,
    first_byte: ?u8,
    second_byte: ?u8,
    start_offset: u64,
};

const FallbackValidation = struct {
    alphabet: Alphabet,
    use_full_iupac: bool,
    sequence_failure: ?usize = null,
    quality_failure: ?usize = null,
    pending_cr: bool = false,

    fn init(validator: *const AdaptiveRecordValidator) FallbackValidation {
        return .{
            .alphabet = validator.alphabet,
            .use_full_iupac = validator.use_full_iupac,
        };
    }

    fn consume(
        self: *FallbackValidation,
        line_kind: ExpectedLine,
        bytes: []const u8,
        start_index: usize,
    ) ReaderError!void {
        return switch (line_kind) {
            .sequence => self.consumeField(.sequence, bytes, start_index),
            .quality => self.consumeField(.quality, bytes, start_index),
            .header, .plus => {},
        };
    }

    fn consumeField(
        self: *FallbackValidation,
        comptime field: SemanticField,
        bytes: []const u8,
        start_index: usize,
    ) ReaderError!void {
        if (bytes.len == 0) return;
        if (self.pending_cr) {
            std.debug.assert(start_index > 0);
            classifySemanticByte(
                field,
                self.alphabet,
                &self.use_full_iupac,
                self.failure(field),
                '\r',
                start_index - 1,
            );
            self.pending_cr = false;
        }

        const semantic_end = if (bytes[bytes.len - 1] == '\r') bytes.len - 1 else bytes.len;
        self.pending_cr = semantic_end != bytes.len;
        const semantic_bytes = bytes[0..semantic_end];
        if (semantic_bytes.len == 0) return;
        classifySemanticBytes(
            field,
            self.alphabet,
            &self.use_full_iupac,
            self.failure(field),
            semantic_bytes,
            start_index,
        ) catch return error.ArithmeticLimit;
    }

    fn finishLine(
        self: *FallbackValidation,
        line_kind: ExpectedLine,
        terminated: bool,
        content_len: usize,
    ) void {
        switch (line_kind) {
            .sequence => self.finishField(.sequence, terminated, content_len),
            .quality => self.finishField(.quality, terminated, content_len),
            .header, .plus => {},
        }
    }

    fn finishField(
        self: *FallbackValidation,
        comptime field: SemanticField,
        terminated: bool,
        content_len: usize,
    ) void {
        if (self.pending_cr and !terminated) {
            std.debug.assert(content_len > 0);
            classifySemanticByte(
                field,
                self.alphabet,
                &self.use_full_iupac,
                self.failure(field),
                '\r',
                content_len - 1,
            );
        }
        self.pending_cr = false;
    }

    fn failure(
        self: *FallbackValidation,
        comptime field: SemanticField,
    ) *?usize {
        return switch (field) {
            .sequence => &self.sequence_failure,
            .quality => &self.quality_failure,
        };
    }

    fn result(self: *const FallbackValidation) ?SemanticError {
        if (self.sequence_failure) |byte_index| return semanticSequenceError(byte_index);
        if (self.quality_failure) |byte_index| return semanticQualityError(byte_index);
        return null;
    }

    fn commit(self: *const FallbackValidation, validator: *AdaptiveRecordValidator) void {
        if (self.sequence_failure == null) {
            validator.use_full_iupac = self.use_full_iupac;
        }
    }
};

const FallbackField = struct {
    storage: []u8 = &.{},
    len: usize = 0,
};

pub const RetainedRecordStorage = struct {
    fields: [4]FallbackField = .{ .{}, .{}, .{}, .{} },

    pub fn deinit(self: *RetainedRecordStorage, allocator: std.mem.Allocator) void {
        for (self.fields) |field| {
            if (field.storage.len != 0) allocator.free(field.storage);
        }
        self.* = undefined;
    }
};

const BufferedRecord = struct {
    bytes: []const u8,
    ranges: [4]Range,
    canonical_range: ?Range,
};

const BufferedValidatedRecord = struct {
    record: BufferedRecord,
    semantic_error: ?SemanticError,
};

const BufferedPayload = struct {
    sequence: Range,
    quality: Range,
};

const BufferedProjection = enum { full, payload, predicted_payload, validated_header, validated };

fn BufferedRecordResult(comptime projection: BufferedProjection) type {
    return union(enum) {
        incomplete,
        eof,
        record: switch (projection) {
            .full => BufferedRecord,
            .payload => BufferedPayload,
            .predicted_payload => struct { payload: BufferedPayload, checkpoint: ?PayloadCheckpoint },
            .validated_header => struct { header: Range, semantic_error: ?SemanticError },
            .validated => BufferedValidatedRecord,
        },
    };
}

// Saved positions remain valid only while the buffer's bytes and indexes stay unchanged.
const LineFeedSearch = struct {
    block_start: usize = 0,
    block_end: usize = 0,
    mask: @Int(.unsigned, @max(64, std.simd.suggestVectorLength(u8) orelse 1)) = 0,

    fn find(self: *LineFeedSearch, comptime lanes: usize, bytes: []const u8, start: usize) ?usize {
        var offset = start;
        if (start >= self.block_start and start < self.block_end and self.block_end <= bytes.len) {
            const remaining = self.mask >> @intCast(start - self.block_start);
            if (remaining != 0) return start + @as(usize, @intCast(@ctz(remaining)));
            offset = self.block_end;
        }

        while (bytes.len - offset >= lanes) : (offset += lanes) {
            self.block_start = offset;
            self.block_end = offset + lanes;
            self.mask = newlineMask(lanes, bytes[offset..][0..lanes]);
            if (self.mask != 0) return offset + @as(usize, @intCast(@ctz(self.mask)));
        }
        for (bytes[offset..], offset..) |byte, index| {
            if (byte == '\n') return index;
        }
        return null;
    }

    fn findLines(
        self: *LineFeedSearch,
        comptime line_count: usize,
        bytes: []const u8,
        start: usize,
        ends: *[line_count]usize,
    ) bool {
        var offset = start;
        for (ends) |*end| {
            const line_end = self.find(std.simd.suggestVectorLength(u8) orelse 1, bytes, offset) orelse
                return false;
            end.* = line_end - start;
            offset = line_end + 1;
        }
        return true;
    }
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
    line_feeds: LineFeedSearch = .{},
    fallback_fields: [4]FallbackField,
    record_index: u64,
    byte_offset: u64,
    options: Options,
    machine: Machine,
    last_error: ?ParseError,
    record_offsets: RecordOffsets = undefined,
    current_record_offsets: ?RecordOffsets = null,
    transport_storage: []u8,
    borrowed_gzip: ?*io_layer.GzipSource,
    pending_refill: ?struct { bytes: []u8, cursor: usize, end: usize } = null,

    pub fn init(
        allocator: std.mem.Allocator,
        source: ByteSource,
        options: Options,
    ) !Reader {
        const buf = try allocator.alloc(u8, io_layer.DEFAULT_READER_BUFFER_BYTES);
        return initBuffer(allocator, source, options, buf);
    }

    fn initBuffer(
        allocator: std.mem.Allocator,
        source: ByteSource,
        options: Options,
        buf: []u8,
    ) Reader {
        return .{
            .allocator = allocator,
            .source = source,
            .buf = buf,
            .transport_storage = buf,
            .borrowed_gzip = null,
            .fill_end = 0,
            .cursor = 0,
            .fallback_fields = .{ .{}, .{}, .{}, .{} },
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
        if (self.transport_storage.len != 0) self.allocator.free(self.transport_storage);
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
        var unused_canonical_span: ?[]const u8 = null;
        return self.nextRecord(&unused_canonical_span, true, false);
    }

    fn nextWithCanonicalSpan(
        self: *Reader,
        canonical_span: *?[]const u8,
        comptime derive_id: bool,
    ) ReaderError!?Record {
        return self.nextRecord(canonical_span, derive_id, true);
    }

    fn nextRecord(
        self: *Reader,
        canonical_span: *?[]const u8,
        comptime derive_id: bool,
        comptime include_canonical_span: bool,
    ) ReaderError!?Record {
        canonical_span.* = null;
        self.beginRecord();
        switch (try self.readBufferedRecord(.full, include_canonical_span, null, true)) {
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
        const header = buffered.bytes[ranges[0].start + 1 .. ranges[0].end];
        canonical_span.* = if (buffered.canonical_range) |range|
            range.slice(buffered.bytes)
        else
            null;
        return .{
            .header = header,
            .id = if (derive_id) firstToken(header) else header[0..0],
            .sequence = ranges[1].slice(buffered.bytes),
            .plus = buffered.bytes[ranges[2].start + 1 .. ranges[2].end],
            .quality = ranges[3].slice(buffered.bytes),
        };
    }

    fn nextPayload(self: *Reader) ReaderError!?RecordPayload {
        self.beginRecord();
        switch (try self.readBufferedRecord(.payload, false, null, false)) {
            .incomplete => return self.nextFallbackPayload(),
            .eof => return null,
            .record => |buffered| {
                self.current_record_offsets = self.record_offsets;
                return .{
                    .sequence = buffered.sequence.slice(self.buf),
                    .quality = buffered.quality.slice(self.buf),
                };
            },
        }
    }

    fn nextPredictedPayload(self: *Reader) ReaderError!?PredictedPayload {
        self.beginRecord();
        switch (try self.readBufferedRecord(.predicted_payload, false, null, false)) {
            .incomplete => return if (try self.nextFallbackPayload()) |payload|
                .{ .payload = payload }
            else
                null,
            .eof => return null,
            .record => |buffered| {
                self.current_record_offsets = self.record_offsets;
                return .{
                    .payload = .{
                        .sequence = buffered.payload.sequence.slice(self.buf),
                        .quality = buffered.payload.quality.slice(self.buf),
                    },
                    .checkpoint = buffered.checkpoint,
                };
            },
        }
    }

    fn nextValidatedHeader(
        self: *Reader,
        validator: *AdaptiveRecordValidator,
    ) ReaderError!?ValidatedHeader {
        self.beginRecord();
        switch (try self.readBufferedRecord(.validated_header, false, validator, false)) {
            .incomplete => return self.nextFallbackValidatedHeader(validator),
            .eof => return null,
            .record => |buffered| {
                self.current_record_offsets = self.record_offsets;
                return .{
                    .header = self.buf[buffered.header.start + 1 .. buffered.header.end],
                    .semantic_error = buffered.semantic_error,
                };
            },
        }
    }

    fn nextFallback(self: *Reader, comptime derive_id: bool) ReaderError!?Record {
        if (!try self.readFallbackRecord(0b1111)) return null;
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
    }

    fn nextFallbackPayload(self: *Reader) ReaderError!?RecordPayload {
        if (!try self.readFallbackRecord(0b1010)) return null;
        return .{
            .sequence = self.fallback_fields[1].storage[0..self.fallback_fields[1].len],
            .quality = self.fallback_fields[3].storage[0..self.fallback_fields[3].len],
        };
    }

    fn nextFallbackValidatedHeader(
        self: *Reader,
        validator: *AdaptiveRecordValidator,
    ) ReaderError!?ValidatedHeader {
        self.beginFallbackRecord();
        var validation = FallbackValidation.init(validator);
        while (true) {
            const line = try self.readValidatedFallbackLine(&validation);
            switch (try self.ingestLine(line)) {
                .eof => return null,
                .continue_ => {},
                .record_ready => {
                    self.current_record_offsets = self.record_offsets;
                    validation.commit(validator);
                    const header = self.fallback_fields[0];
                    return .{
                        .header = header.storage[1..header.len],
                        .semantic_error = validation.result(),
                    };
                },
            }
        }
    }

    fn readFallbackRecord(self: *Reader, comptime retained_fields: u4) ReaderError!bool {
        self.beginFallbackRecord();
        while (true) {
            const line = try self.readFallbackLine(retained_fields);
            switch (try self.ingestLine(line)) {
                .eof => return false,
                .continue_ => {},
                .record_ready => {
                    self.current_record_offsets = self.record_offsets;
                    return true;
                },
            }
        }
    }

    /// Consumes one record without returning its fields, or returns false at clean EOF.
    pub fn advance(self: *Reader) ReaderError!bool {
        self.beginRecord();
        switch (try self.readBufferedRecord(.full, false, null, false)) {
            .incomplete => return self.advanceFallback(),
            .eof => return false,
            .record => return true,
        }
    }

    fn advanceFallback(self: *Reader) ReaderError!bool {
        self.beginFallbackRecord();
        while (true) {
            const line = try self.readFallbackLine(0);
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

    fn readBufferedRecord(
        self: *Reader,
        comptime projection: BufferedProjection,
        comptime include_canonical_span: bool,
        validator: ?*AdaptiveRecordValidator,
        comptime allow_refill: bool,
    ) ReaderError!BufferedRecordResult(projection) {
        if (self.machine.expected != .header) return .incomplete;
        if (self.cursor == self.fill_end and !try self.refill()) return .eof;

        var relative_ends: [4]usize = undefined;
        var checkpoint: ?PayloadCheckpoint = null;
        const complete = if (projection == .validated or projection == .validated_header)
            scanCompleteRecord(
                self.buf[0..self.fill_end],
                self.cursor,
                &self.line_feeds,
                self.options.max_line_bytes,
                validator.?.alphabet,
                validator.?.use_full_iupac,
            )
        else
            null;
        if (projection == .predicted_payload) {
            const bytes = self.buf[self.cursor..self.fill_end];
            if (!self.line_feeds.findLines(3, self.buf[0..self.fill_end], self.cursor, relative_ends[0..3])) return .incomplete;
            if (self.predictQualityEnd(&relative_ends)) {
                checkpoint = .{
                    .cursor = self.cursor,
                    .byte_offset = self.byte_offset,
                    .record_index = self.record_index,
                    .machine = self.machine,
                    .record_offsets = self.record_offsets,
                };
            } else {
                const quality_start = relative_ends[2] + 1;
                const quality_lf = std.mem.findScalar(u8, bytes[quality_start..], '\n') orelse
                    return .incomplete;
                relative_ends[3] = quality_start + quality_lf;
            }
        } else if (complete) |checked| {
            relative_ends = checked.line_ends;
        } else if (!self.line_feeds.findLines(4, self.buf[0..self.fill_end], self.cursor, &relative_ends)) {
            if (allow_refill) {
                if (try self.bufferIncompleteRecord(0)) |bytes| {
                    return self.readSavedRecord(bytes, projection, include_canonical_span, validator);
                }
            }
            return .incomplete;
        }

        var ranges: [4]Range = undefined;
        var canonical = true;
        const record_start = self.cursor;
        for (0..4) |line_index| {
            const start = if (line_index == 0)
                record_start
            else
                record_start + relative_ends[line_index - 1] + 1;
            const raw_end = record_start + relative_ends[line_index];
            const end = start + lineContentLen(self.buf[start..raw_end]);
            if (include_canonical_span) canonical = canonical and end == raw_end;
            ranges[line_index] = .{ .start = start, .end = end };
        }

        try self.consumeBufferedRecord(ranges, relative_ends, complete != null);
        const semantic_error = if (projection == .validated or projection == .validated_header) result: {
            if (complete) |checked| {
                validator.?.use_full_iupac = checked.use_full_iupac;
                break :result null;
            }
            break :result validator.?.validate(.{
                .header = "",
                .id = "",
                .sequence = ranges[1].slice(self.buf),
                .plus = "",
                .quality = ranges[3].slice(self.buf),
            });
        } else null;
        if (projection == .full or projection == .validated) {
            const record: BufferedRecord = .{
                .bytes = self.buf,
                .ranges = ranges,
                .canonical_range = if (include_canonical_span and canonical)
                    .{ .start = record_start, .end = record_start + relative_ends[3] + 1 }
                else
                    null,
            };
            return if (projection == .validated)
                .{ .record = .{ .record = record, .semantic_error = semantic_error } }
            else
                .{ .record = record };
        }
        const payload: BufferedPayload = .{ .sequence = ranges[1], .quality = ranges[3] };
        return if (projection == .validated_header)
            .{ .record = .{ .header = ranges[0], .semantic_error = semantic_error } }
        else if (projection == .predicted_payload)
            .{ .record = .{ .payload = payload, .checkpoint = checkpoint } }
        else
            .{ .record = payload };
    }

    fn consumeBufferedRecord(
        self: *Reader,
        ranges: [4]Range,
        relative_ends: [4]usize,
        already_validated: bool,
    ) ReaderError!void {
        const record_start = self.cursor;
        const record_len = relative_ends[3] + 1;
        const sequence_len = ranges[1].end - ranges[1].start;
        const valid = already_validated or valid: {
            const header = ranges[0].slice(self.buf);
            const plus = ranges[2].slice(self.buf);
            break :valid header.len <= self.options.max_line_bytes and
                sequence_len <= self.options.max_line_bytes and
                plus.len <= self.options.max_line_bytes and
                headerPrefixIsValid(
                    if (header.len == 0) null else header[0],
                    if (header.len < 2) null else header[1],
                ) and plus.len != 0 and plus[0] == '+' and
                sequence_len == ranges[3].end - ranges[3].start;
        };
        const next_offset = if (valid) self.offsetAfter(record_len) catch null else null;
        if (next_offset) |end_offset| {
            // The checked end also bounds every intermediate line offset.
            self.record_offsets = .{
                .header = self.byte_offset,
                .sequence = self.byte_offset + relative_ends[0] + 1,
                .plus = self.byte_offset + relative_ends[1] + 1,
                .quality = self.byte_offset + relative_ends[2] + 1,
            };
            self.cursor = record_start + record_len;
            self.byte_offset = end_offset;
            self.machine.sequence_len = sequence_len;
            self.record_index = std.math.add(u64, self.record_index, 1) catch
                return error.ArithmeticLimit;
            return;
        }

        // Replay rejected records in line order to preserve errors and partial progress.
        for (ranges, relative_ends) |range, relative_end| {
            const content = range.slice(self.buf);
            if (content.len > self.options.max_line_bytes) return error.LineTooLong;
            const raw_end = record_start + relative_end;
            const start_offset = self.byte_offset;
            const line_end_offset = try self.offsetAfter(raw_end + 1 - range.start);
            self.cursor = raw_end + 1;
            self.byte_offset = line_end_offset;
            _ = try self.ingestLine(.{
                .content_len = content.len,
                .first_byte = if (content.len == 0) null else content[0],
                .second_byte = if (content.len < 2) null else content[1],
                .start_offset = start_offset,
            });
        }
    }

    fn predictQualityEnd(self: *const Reader, ends: *[4]usize) bool {
        const bytes = self.buf[self.cursor..self.fill_end];
        const sequence_start = ends[0] + 1;
        const sequence_len = lineContentLen(bytes[sequence_start..ends[1]]);
        const quality_start = ends[2] + 1;
        if (sequence_len >= bytes.len - quality_start) return false;
        const quality_end = quality_start + sequence_len;
        const quality_lf = quality_end + @intFromBool(bytes[quality_end] == '\r');
        if (quality_lf == bytes.len or bytes[quality_lf] != '\n') return false;
        // Stripping a guessed CR could reject the length before quality validation can retry.
        if (quality_lf == quality_end and sequence_len != 0 and bytes[quality_end - 1] == '\r')
            return false;

        // A guessed end must not hide an earlier length error with counter overflow.
        if (self.record_index == std.math.maxInt(u64)) return false;
        const consumed = std.math.cast(u64, quality_lf + 1) orelse return false;
        _ = std.math.add(u64, self.byte_offset, consumed) catch return false;
        ends[3] = quality_lf;
        return true;
    }

    const RefillRecord = struct {
        start: usize,
        len: usize = 0,
        stopped: bool = false,
    };

    fn readSavedRecord(
        self: *Reader,
        bytes: []u8,
        comptime projection: BufferedProjection,
        comptime include_canonical_span: bool,
        validator: ?*AdaptiveRecordValidator,
    ) ReaderError!BufferedRecordResult(projection) {
        const buf = self.buf;
        const cursor = self.cursor;
        const fill_end = self.fill_end;
        self.buf = bytes;
        self.line_feeds = .{};
        self.cursor = 0;
        self.fill_end = bytes.len;
        defer {
            self.buf = buf;
            self.line_feeds = .{};
            self.cursor = cursor;
            self.fill_end = fill_end;
        }
        return self.readBufferedRecord(projection, include_canonical_span, validator, false);
    }

    fn canBufferIncompleteRecord(self: *const Reader) bool {
        return self.borrowed_gzip == null and self.pending_refill == null and
            self.fill_end - self.cursor != self.buf.len;
    }

    fn ensureRefillCapacity(self: *Reader) ReaderError!void {
        const field = &self.fallback_fields[1];
        const capacity = io_layer.DEFAULT_READER_BUFFER_BYTES;
        if (field.storage.len >= capacity) return;
        field.storage = self.allocator.realloc(field.storage, capacity) catch
            return error.OutOfMemory;
    }

    fn bufferIncompleteRecord(self: *Reader, prefix_len: usize) ReaderError!?[]u8 {
        if (!self.canBufferIncompleteRecord()) return null;
        try self.ensureRefillCapacity();
        const byte_offset = self.byte_offset;
        const record_index = self.record_index;
        const machine = self.machine;
        const record_offsets = self.record_offsets;
        var record: RefillRecord = .{ .start = self.cursor, .len = prefix_len };
        while (true) {
            const line = try self.readLine(false, null, &record);
            if (record.stopped) break;
            if (try self.ingestLine(line) == .record_ready) break;
        }
        if (!record.stopped) self.copyRefillRecordBytes(&record);
        self.byte_offset = byte_offset;
        self.record_index = record_index;
        self.machine = machine;
        self.record_offsets = record_offsets;
        const bytes = self.fallback_fields[1].storage[0..record.len];
        if (!record.stopped) return bytes[prefix_len..];
        self.cursor = record.start;
        if (record.len != 0) {
            self.pending_refill = .{ .bytes = self.buf, .cursor = self.cursor, .end = self.fill_end };
            self.buf = bytes;
            self.line_feeds = .{};
            self.cursor = prefix_len;
            self.fill_end = bytes.len;
        }
        return null;
    }

    fn copyRefillRecordBytes(self: *Reader, record: *RefillRecord) void {
        const bytes = self.buf[record.start..self.cursor];
        if (bytes.len > io_layer.DEFAULT_READER_BUFFER_BYTES - record.len) {
            record.stopped = true;
            return;
        }
        @memcpy(self.fallback_fields[1].storage[record.len..][0..bytes.len], bytes);
        record.len += bytes.len;
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
        const record_ready = self.machine.push(
            actual_line.content_len,
            actual_line.first_byte,
            actual_line.second_byte,
        ) catch |err| {
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
            self.record_index = std.math.add(u64, self.record_index, 1) catch
                return error.ArithmeticLimit;
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
        if (self.borrowed_gzip != null) {
            std.debug.assert(self.cursor == self.fill_end);
            self.cursor = 0;
            self.fill_end = 0;
            return;
        }
        if (self.cursor >= self.fill_end) {
            self.cursor = 0;
            self.fill_end = 0;
            return;
        }
        if (self.cursor == 0) return;

        const tail_len = self.fill_end - self.cursor;
        @memmove(self.buf[0..tail_len], self.buf[self.cursor..self.fill_end]);
        self.fill_end = tail_len;
        self.cursor = 0;
    }

    fn refill(self: *Reader) ReaderError!bool {
        self.line_feeds = .{};
        if (self.pending_refill) |pending| {
            std.debug.assert(self.cursor == self.fill_end);
            self.buf = pending.bytes;
            self.cursor = pending.cursor;
            self.fill_end = pending.end;
            self.pending_refill = null;
            return self.cursor != self.fill_end;
        }
        self.compactIfNeeded();
        if (self.borrowed_gzip) |source| {
            const decoded = io_layer.readGzipChunk(
                source,
                &source.decompressor_buffer,
            ) catch return error.Io;
            self.buf = if (decoded) |bytes|
                @constCast(bytes)
            else
                source.decompressor_buffer[0..0];
            self.fill_end = self.buf.len;
            return self.fill_end != 0;
        }
        const space = self.buf.len - self.fill_end;
        std.debug.assert(space > 0);
        const n = self.source.read(self.buf[self.fill_end..]) catch return error.Io;
        self.fill_end += n;
        return n > 0;
    }

    fn readFallbackLine(self: *Reader, comptime retained_fields: u4) ReaderError!?Line {
        return switch (self.machine.expected) {
            .header => self.readLine(retained_fields & 0b0001 != 0, null, null),
            .sequence => self.readLine(retained_fields & 0b0010 != 0, null, null),
            .plus => self.readLine(retained_fields & 0b0100 != 0, null, null),
            .quality => self.readLine(retained_fields & 0b1000 != 0, null, null),
        };
    }

    fn readValidatedFallbackLine(
        self: *Reader,
        validation: *FallbackValidation,
    ) ReaderError!?Line {
        return switch (self.machine.expected) {
            .header => self.readLine(true, validation, null),
            .sequence, .plus, .quality => self.readLine(false, validation, null),
        };
    }

    fn readLine(
        self: *Reader,
        comptime retain: bool,
        validation: ?*FallbackValidation,
        refill_record: ?*RefillRecord,
    ) ReaderError!?Line {
        const field_index = @intFromEnum(self.machine.expected);
        const field = &self.fallback_fields[field_index];
        const content_start = if (retain) field.len else 0;
        const line_start_offset = self.byte_offset;
        var discarded_len: usize = 0;
        var first_byte: ?u8 = null;
        var second_byte: ?u8 = null;
        var last_byte: ?u8 = null;

        while (true) {
            if (self.cursor < self.fill_end) {
                const haystack = self.buf[self.cursor..self.fill_end];
                if (std.mem.findScalar(u8, haystack, '\n')) |rel| {
                    if (retain) {
                        try self.appendLineBytes(field_index, content_start, haystack[0..rel], true);
                    } else {
                        const chunk_start = discarded_len;
                        try self.discardLineBytes(
                            &discarded_len,
                            &first_byte,
                            &second_byte,
                            &last_byte,
                            haystack[0..rel],
                        );
                        if (validation) |state| {
                            try state.consume(self.machine.expected, haystack[0..rel], chunk_start);
                        }
                    }
                    const next_offset = try self.offsetAfter(rel + 1);
                    self.cursor += rel + 1;
                    self.byte_offset = next_offset;
                    if (retain) {
                        field.len = content_start + lineContentLen(field.storage[content_start..field.len]);
                        const content = field.storage[0..field.len];
                        return .{
                            .content_len = content.len,
                            .first_byte = if (content.len == 0) null else content[0],
                            .second_byte = if (content.len < 2) null else content[1],
                            .start_offset = line_start_offset,
                        };
                    }
                    if (validation) |state| {
                        state.finishLine(self.machine.expected, true, discarded_len);
                    }
                    if (last_byte == '\r') discarded_len -= 1;
                    return .{
                        .content_len = discarded_len,
                        .first_byte = if (discarded_len == 0) null else first_byte,
                        .second_byte = if (discarded_len < 2) null else second_byte,
                        .start_offset = line_start_offset,
                    };
                }

                if (retain) {
                    try self.appendLineBytes(field_index, content_start, haystack, false);
                } else {
                    const chunk_start = discarded_len;
                    try self.discardLineBytes(
                        &discarded_len,
                        &first_byte,
                        &second_byte,
                        &last_byte,
                        haystack,
                    );
                    if (validation) |state| {
                        try state.consume(self.machine.expected, haystack, chunk_start);
                    }
                }
                const next_offset = try self.offsetAfter(haystack.len);
                self.cursor = self.fill_end;
                self.byte_offset = next_offset;
            }

            if (refill_record) |record| {
                self.copyRefillRecordBytes(record);
                if (record.stopped) return null;
                const got_data = try self.refill();
                record.start = 0;
                if (!got_data) {
                    record.stopped = true;
                    return null;
                }
                continue;
            }

            const got_data = try self.refill();
            if (!got_data) {
                const content_len = if (retain) field.len - content_start else discarded_len;
                if (content_len > 0) {
                    if (content_len > self.options.max_line_bytes) {
                        return error.LineTooLong;
                    }
                    if (validation) |state| {
                        state.finishLine(self.machine.expected, false, content_len);
                    }
                    if (!retain) return .{
                        .content_len = content_len,
                        .first_byte = first_byte,
                        .second_byte = second_byte,
                        .start_offset = line_start_offset,
                    };
                    const content = field.storage[0..field.len];
                    return .{
                        .content_len = content.len,
                        .first_byte = content[0],
                        .second_byte = if (content.len < 2) null else content[1],
                        .start_offset = line_start_offset,
                    };
                }
                return null;
            }
        }
    }

    fn offsetAfter(self: *const Reader, amount: usize) ReaderError!u64 {
        return progressAfter(self.byte_offset, amount);
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
        const line_len = try self.lineLengthAfter(field.len - content_start, chunk);
        const new_len = std.math.add(usize, content_start, line_len) catch
            return error.LineTooLong;
        try self.ensureFieldCapacity(field_index, new_len, complete);
        if (self.buf.ptr == field.storage.ptr) {
            // A replayed record prefix shares the sequence fallback allocation.
            @memmove(field.storage[field.len..new_len], chunk);
        } else @memcpy(field.storage[field.len..new_len], chunk);
        field.len = new_len;
    }

    fn discardLineBytes(
        self: *const Reader,
        line_len: *usize,
        first_byte: *?u8,
        second_byte: *?u8,
        last_byte: *?u8,
        chunk: []const u8,
    ) ReaderError!void {
        if (chunk.len == 0) return;
        const old_len = line_len.*;
        line_len.* = try self.lineLengthAfter(old_len, chunk);
        if (old_len == 0) {
            first_byte.* = chunk[0];
            if (chunk.len > 1) second_byte.* = chunk[1];
        } else if (old_len == 1) {
            second_byte.* = chunk[0];
        }
        last_byte.* = chunk[chunk.len - 1];
    }

    fn lineLengthAfter(
        self: *const Reader,
        current_len: usize,
        chunk: []const u8,
    ) ReaderError!usize {
        const new_len = std.math.add(usize, current_len, chunk.len) catch
            return error.LineTooLong;
        if (new_len > self.options.max_line_bytes) {
            const excess = new_len - self.options.max_line_bytes;
            if (excess > 1 or chunk[chunk.len - 1] != '\r') return error.LineTooLong;
        }
        return new_len;
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
        if (needed > io_layer.DEFAULT_READER_BUFFER_BYTES) {
            var available: ?*FallbackField = null;
            for (&self.fallback_fields) |*donor| {
                if (donor == field or donor.len > field.storage.len or donor.storage.len < grown_len or
                    donor.storage.len > storage_limit) continue;
                if (available == null or donor.len < available.?.len) available = donor;
            }
            if (available) |donor| {
                // No field slices escape before this fallback record is complete.
                const shared = @min(field.len, donor.len);
                for (field.storage[0..shared], donor.storage[0..shared]) |*a, *b| {
                    std.mem.swap(u8, a, b);
                }
                if (field.len > shared) {
                    @memcpy(donor.storage[shared..field.len], field.storage[shared..field.len]);
                } else {
                    @memcpy(field.storage[shared..donor.len], donor.storage[shared..donor.len]);
                }
                std.mem.swap([]u8, &field.storage, &donor.storage);
                return;
            }
        }
        field.storage = if (field.storage.len == 0)
            self.allocator.alloc(u8, grown_len) catch return error.OutOfMemory
        else
            self.allocator.realloc(field.storage, grown_len) catch return error.OutOfMemory;
    }
};

pub fn initBorrowedGzipReader(
    allocator: std.mem.Allocator,
    source: *io_layer.GzipSource,
    options: Options,
) !Reader {
    std.debug.assert(source.decompressor_buffer.len != 0);
    var reader = Reader.initBuffer(allocator, source.byteSource(), options, source.decompressor_buffer[0..0]);
    reader.borrowed_gzip = source;
    return reader;
}

fn fieldStorageLimit(max_line_bytes: usize) usize {
    return std.math.add(usize, max_line_bytes, 1) catch std.math.maxInt(usize);
}

pub fn nextWithoutId(
    reader: *Reader,
    canonical_span: *?[]const u8,
) ReaderError!?Record {
    return reader.nextWithCanonicalSpan(canonical_span, false);
}

pub fn nextValidatedRecord(
    reader: *Reader,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    reader.beginRecord();
    return switch (try reader.readBufferedRecord(.validated, true, validator, true)) {
        .incomplete => nextFallbackValidatedRecord(reader, validator),
        .eof => null,
        .record => |buffered| finishValidatedRecord(reader, buffered),
    };
}

fn finishValidatedRecord(
    reader: *Reader,
    buffered: BufferedValidatedRecord,
) ValidatedRecord {
    var canonical_span: ?[]const u8 = null;
    const record = reader.finishBufferedRecord(buffered.record, &canonical_span, false);
    return .{
        .record = record,
        .canonical_span = canonical_span,
        .semantic_error = buffered.semantic_error,
    };
}

pub fn nextRecordWithoutId(reader: *Reader) ReaderError!?Record {
    var unused_canonical_span: ?[]const u8 = null;
    return reader.nextRecord(&unused_canonical_span, false, false);
}

pub fn nextPayload(reader: *Reader) ReaderError!?RecordPayload {
    return reader.nextPayload();
}

/// Validate predicted quality as Phred+33 before advancing the reader again.
/// On rejection, reread from the checkpoint to recover structural errors first.
pub fn nextPredictedPayload(reader: *Reader) ReaderError!?PredictedPayload {
    return reader.nextPredictedPayload();
}

pub fn nextValidatedHeader(
    reader: *Reader,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedHeader {
    return reader.nextValidatedHeader(validator);
}

pub fn nextBufferedWithoutId(
    reader: *Reader,
    canonical_span: *?[]const u8,
) ReaderError!?Record {
    return nextBufferedRecord(reader, canonical_span, true);
}

pub fn nextBufferedRecordWithoutId(reader: *Reader) ReaderError!?Record {
    var unused_canonical_span: ?[]const u8 = null;
    return nextBufferedRecord(reader, &unused_canonical_span, false);
}

pub fn nextBufferedValidatedRecord(
    reader: *Reader,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    reader.beginRecord();
    if (reader.cursor == reader.fill_end) return null;
    return switch (try reader.readBufferedRecord(.validated, true, validator, false)) {
        .incomplete, .eof => null,
        .record => |buffered| finishValidatedRecord(reader, buffered),
    };
}

/// Updates mate 1's views when moving it; both mates borrow until the next advance.
/// On null, preserve the updated mate 1 before continuing with ordinary reads.
pub fn nextPairedValidatedRecord(
    reader: *Reader,
    first: *ValidatedRecord,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    if (try nextBufferedValidatedRecord(reader, validator)) |second| return second;
    if (!reader.canBufferIncompleteRecord()) return null;
    // Separate fallback fields can overlap the destination before they are copied.
    if (first.record.sequence.ptr == reader.fallback_fields[1].storage.ptr) return null;
    const capacity = io_layer.DEFAULT_READER_BUFFER_BYTES;
    const span = first.canonical_span;
    var size: usize = 6;
    if (span) |bytes| {
        size = bytes.len;
    } else {
        inline for (.{ "header", "sequence", "plus", "quality" }) |field| {
            const len = @field(first.record, field).len;
            if (len > capacity - size) return null;
            size += len;
        }
    }
    if (size >= capacity) return null;
    try reader.ensureRefillCapacity();
    const saved = reader.fallback_fields[1].storage[0..size];
    if (span) |bytes| @memmove(saved, bytes) else saved[0] = '@';
    var offset: usize = 1;
    inline for (.{ "header", "sequence", "plus", "quality" }, 0..) |field, index| {
        const bytes = &@field(first.record, field);
        const start = if (span) |raw| @intFromPtr(bytes.ptr) - @intFromPtr(raw.ptr) else offset;
        const dest = saved[start..][0..bytes.len];
        if (span == null) {
            @memmove(dest, bytes.*);
            offset += bytes.len;
            saved[offset] = '\n';
            offset += 1;
            if (index == 1) {
                saved[offset] = '+';
                offset += 1;
            }
        }
        bytes.* = dest;
    }
    first.record.id = first.record.header[0..first.record.id.len];
    // Staging writes one canonical span, including when the input uses CRLF.
    first.canonical_span = if (recordHasTerminalCr(first.record)) null else saved;
    return nextSavedPairMate(reader, size, validator);
}

/// Keeps only mate 1's updated header alive through the mate-2 read.
/// On null, preserve that header before continuing with ordinary reads.
pub fn nextPairedValidatedHeader(
    reader: *Reader,
    first_header: *[]const u8,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    if (try nextBufferedValidatedRecord(reader, validator)) |second| return second;
    if (!reader.canBufferIncompleteRecord() or
        first_header.len >= io_layer.DEFAULT_READER_BUFFER_BYTES) return null;
    try reader.ensureRefillCapacity();
    const saved = reader.fallback_fields[1].storage[0..first_header.len];
    @memmove(saved, first_header.*);
    first_header.* = saved;
    return nextSavedPairMate(reader, saved.len, validator);
}

fn nextSavedPairMate(
    reader: *Reader,
    prefix_len: usize,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    const bytes = try reader.bufferIncompleteRecord(prefix_len) orelse return null;
    const buffered = try reader.readSavedRecord(bytes, .validated, true, validator);
    std.debug.assert(buffered == .record);
    return finishValidatedRecord(reader, buffered.record);
}

fn nextBufferedRecord(
    reader: *Reader,
    canonical_span: *?[]const u8,
    comptime include_canonical_span: bool,
) ReaderError!?Record {
    canonical_span.* = null;
    reader.beginRecord();
    if (reader.cursor == reader.fill_end) return null;
    return switch (try reader.readBufferedRecord(.full, include_canonical_span, null, false)) {
        .incomplete => null,
        .eof => unreachable,
        .record => |buffered| reader.finishBufferedRecord(
            buffered,
            canonical_span,
            false,
        ),
    };
}

pub fn retainFallbackRecordStorage(
    reader: *Reader,
    retained: *RetainedRecordStorage,
    record: Record,
) bool {
    const fields = &reader.fallback_fields;
    if (fields[0].len == 0 or
        fields[2].len == 0 or
        fields[0].len - 1 != record.header.len or
        fields[1].len != record.sequence.len or
        fields[2].len - 1 != record.plus.len or
        fields[3].len != record.quality.len or
        record.header.ptr != fields[0].storage.ptr + 1 or
        record.sequence.ptr != fields[1].storage.ptr or
        record.plus.ptr != fields[2].storage.ptr + 1 or
        record.quality.ptr != fields[3].storage.ptr)
    {
        return false;
    }
    std.mem.swap([4]FallbackField, fields, &retained.fields);
    return true;
}

pub fn restoreFallbackRecordStorage(
    reader: *Reader,
    retained: *RetainedRecordStorage,
) void {
    std.mem.swap([4]FallbackField, &reader.fallback_fields, &retained.fields);
}

pub fn nextBufferedAfterFallbackTransfer(
    reader: *Reader,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    reader.beginRecord();
    if (reader.borrowed_gzip != null) return null;
    if (reader.cursor == 0 and reader.fill_end == reader.buf.len) {
        return null;
    }
    const got_data = try reader.refill();
    if (!got_data and reader.cursor == reader.fill_end) return null;
    return nextBufferedValidatedRecord(reader, validator);
}

pub fn nextFallbackValidatedRecord(
    reader: *Reader,
    validator: *AdaptiveRecordValidator,
) ReaderError!?ValidatedRecord {
    reader.beginRecord();
    const record = try reader.nextFallback(false) orelse return null;
    return .{
        .record = record,
        .canonical_span = null,
        .semantic_error = validator.validate(record),
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
            !identifierFirstByteIsValid(if (record.header.len == 0) null else record.header[0]) or
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

pub fn recordHasTerminalCr(record: Record) bool {
    inline for (.{ record.header, record.sequence, record.plus, record.quality }) |field| {
        if (std.mem.endsWith(u8, field, "\r")) return true;
    }
    return false;
}

pub fn writeCanonicalRecordSpan(writer: *Writer, span: []const u8) WriteError!void {
    return writer.sink.write(span);
}

// CLI callers check both mates, including terminal CR, before writing either record.
pub fn writeRecordFields(writer: *Writer, record: Record) WriteError!void {
    const fields = [_][]const u8{
        "@",
        record.header,
        "\n",
        record.sequence,
        "\n+",
        record.plus,
        "\n",
        record.quality,
        "\n",
    };
    try writer.sink.writeVec(&fields);
}

fn isWritableField(bytes: []const u8) bool {
    return std.mem.findScalar(u8, bytes, '\n') == null and
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
            .message = "header line must start with '@' and contain a nonempty identifier",
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

    pub fn push(
        self: *Machine,
        line_len: usize,
        first_byte: ?u8,
        second_byte: ?u8,
    ) Error!bool {
        switch (self.expected) {
            .header => {
                if (!headerPrefixIsValid(first_byte, second_byte)) {
                    return error.S003InvalidHeader;
                }
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

/// Excludes one trailing CR; callers check LF termination and handle EOF separately.
pub fn lineContentLen(line: []const u8) usize {
    if (line.len > 0 and line[line.len - 1] == '\r') return line.len - 1;
    return line.len;
}

pub fn progressAfter(current: u64, amount: usize) error{ArithmeticLimit}!u64 {
    const amount_u64 = std.math.cast(u64, amount) orelse return error.ArithmeticLimit;
    return std.math.add(u64, current, amount_u64) catch error.ArithmeticLimit;
}

pub fn headerPrefixIsValid(first_byte: ?u8, identifier_first_byte: ?u8) bool {
    return first_byte == '@' and identifierFirstByteIsValid(identifier_first_byte);
}

fn identifierFirstByteIsValid(byte: ?u8) bool {
    return byte != null and byte != ' ' and byte != '\t';
}

pub const CheckScannerError = error{
    Format,
    LineTooLong,
    ArithmeticLimit,
};

const STRUCTURAL_BLOCK_BYTES = 64;

pub fn newlineMask(comptime lanes: usize, block: *const [lanes]u8) @Int(.unsigned, lanes) {
    const Bytes = @Vector(lanes, u8);
    const bytes: Bytes = block.*;
    return @bitCast(bytes == @as(Bytes, @splat('\n')));
}

const SequenceLineScan = union(enum) {
    line_end: usize,
    invalid_start: usize,
    incomplete,
};

fn firstInvalidCheckSequence(
    sequence: []const u8,
    alphabet: Alphabet,
    use_full_iupac: *bool,
) ?usize {
    return switch (alphabet) {
        .acgtn => firstInvalidSequence(sequence, .acgtn),
        .iupac => if (use_full_iupac.*)
            firstInvalidSequence(sequence, .iupac)
        else
            firstInvalidNarrowIupacSequence(sequence, use_full_iupac),
    };
}

fn firstInvalidNarrowIupacSequence(
    sequence: []const u8,
    use_full_iupac: *bool,
) ?usize {
    return switch (scanSequenceLineFor(.acgtn, sequence, false)) {
        .incomplete => null,
        .line_end => |line_end| line_end,
        .invalid_start => |start| switch (scanSequenceLineFor(
            .iupac,
            sequence[start..],
            false,
        )) {
            .incomplete => result: {
                use_full_iupac.* = true;
                break :result null;
            },
            .line_end => |line_end| start + line_end,
            .invalid_start => firstInvalidSequenceScalar(
                sequence[start..],
                .iupac,
                start,
            ),
        },
    };
}

fn firstValidCheckSequenceLineEnd(
    bytes: []const u8,
    alphabet: Alphabet,
    use_full_iupac: *bool,
) ?usize {
    return switch (alphabet) {
        .acgtn => firstValidSequenceLineEndFor(.acgtn, bytes),
        .iupac => if (use_full_iupac.*)
            firstValidSequenceLineEndFor(.iupac, bytes)
        else
            firstValidNarrowIupacSequenceLineEnd(bytes, use_full_iupac),
    };
}

fn firstValidNarrowIupacSequenceLineEnd(
    bytes: []const u8,
    use_full_iupac: *bool,
) ?usize {
    return switch (scanSequenceLineFor(.acgtn, bytes, true)) {
        .line_end => |line_end| line_end,
        .incomplete => null,
        .invalid_start => |start| switch (scanSequenceLineFor(.iupac, bytes[start..], true)) {
            .line_end => |line_end| result: {
                use_full_iupac.* = true;
                break :result start + line_end;
            },
            .incomplete => result: {
                use_full_iupac.* = true;
                break :result null;
            },
            .invalid_start => null,
        },
    };
}

fn firstValidSequenceLineEndFor(
    comptime alphabet: Alphabet,
    bytes: []const u8,
) ?usize {
    return switch (scanSequenceLineFor(alphabet, bytes, true)) {
        .line_end => |line_end| line_end,
        .invalid_start, .incomplete => null,
    };
}

fn scanSequenceLineFor(
    comptime alphabet: Alphabet,
    bytes: []const u8,
    comptime allow_crlf: bool,
) SequenceLineScan {
    const Bytes = @Vector(STRUCTURAL_BLOCK_BYTES, u8);
    var block_start: usize = 0;
    while (bytes.len - block_start >= STRUCTURAL_BLOCK_BYTES) {
        const block: Bytes = bytes[block_start..][0..STRUCTURAL_BLOCK_BYTES].*;
        const invalid: u64 = @bitCast(invalidSequenceVector(
            STRUCTURAL_BLOCK_BYTES,
            alphabet,
            block,
        ));
        if (invalid != 0) {
            const stop = block_start + @as(usize, @intCast(@ctz(invalid)));
            if (bytes[stop] == '\n') return .{ .line_end = stop };
            if (allow_crlf and bytes[stop] == '\r' and
                stop + 1 < bytes.len and bytes[stop + 1] == '\n')
            {
                return .{ .line_end = stop + 1 };
            }
            return .{ .invalid_start = block_start };
        }
        block_start += STRUCTURAL_BLOCK_BYTES;
    }
    for (bytes[block_start..], block_start..) |byte, byte_index| {
        if (byte == '\n') return .{ .line_end = byte_index };
        if (allow_crlf and byte == '\r' and
            byte_index + 1 < bytes.len and bytes[byte_index + 1] == '\n')
        {
            return .{ .line_end = byte_index + 1 };
        }
        if (!alphabetAccepts(alphabet, byte)) return .{ .invalid_start = byte_index };
    }
    return .incomplete;
}

fn classifySemanticBytes(
    comptime field: SemanticField,
    alphabet: Alphabet,
    use_full_iupac: *bool,
    failure: *?usize,
    bytes: []const u8,
    start_index: usize,
) error{ArithmeticLimit}!void {
    if (failure.* != null) return;
    const relative = switch (field) {
        .sequence => firstInvalidCheckSequence(bytes, alphabet, use_full_iupac),
        .quality => firstInvalidQuality(bytes),
    } orelse return;
    failure.* = std.math.add(usize, start_index, relative) catch
        return error.ArithmeticLimit;
}

fn classifySemanticByte(
    comptime field: SemanticField,
    alphabet: Alphabet,
    use_full_iupac: *bool,
    failure: *?usize,
    byte: u8,
    byte_index: usize,
) void {
    if (failure.* != null) return;
    switch (field) {
        .sequence => {
            if (alphabet == .acgtn or use_full_iupac.*) {
                if (!alphabetAccepts(alphabet, byte)) failure.* = byte_index;
                return;
            }
            if (alphabetAccepts(.acgtn, byte)) return;
            if (alphabetAccepts(.iupac, byte)) {
                use_full_iupac.* = true;
                return;
            }
            failure.* = byte_index;
        },
        .quality => {
            _ = decodePhred33(byte) catch {
                failure.* = byte_index;
            };
        },
    }
}

const CompleteRecord = struct {
    line_ends: [4]usize,
    use_full_iupac: bool,
};

fn scanCompleteRecord(
    buffer: []const u8,
    start: usize,
    line_feeds: *LineFeedSearch,
    max_line_bytes: usize,
    alphabet: Alphabet,
    initial_full_iupac: bool,
) ?CompleteRecord {
    const data = buffer[start..];
    const header_end = (line_feeds.find(STRUCTURAL_BLOCK_BYTES, buffer, start) orelse return null) - start;
    const header_len = lineContentLen(data[0..header_end]);
    const header = data[0..header_len];
    if (header.len > max_line_bytes) return null;
    if (!headerPrefixIsValid(
        if (header.len == 0) null else header[0],
        if (header.len < 2) null else header[1],
    )) return null;

    const sequence_start = header_end + 1;
    var use_full_iupac = initial_full_iupac;
    const sequence_raw_len = firstValidCheckSequenceLineEnd(
        data[sequence_start..],
        alphabet,
        &use_full_iupac,
    ) orelse return null;
    const sequence_len = lineContentLen(data[sequence_start..][0..sequence_raw_len]);
    if (sequence_len > max_line_bytes) return null;

    const plus_start = sequence_start + sequence_raw_len + 1;
    const plus_raw_len = (line_feeds.find(STRUCTURAL_BLOCK_BYTES, buffer, start + plus_start) orelse
        return null) - start - plus_start;
    const plus_len = lineContentLen(data[plus_start..][0..plus_raw_len]);
    const plus = data[plus_start..][0..plus_len];
    if (plus.len > max_line_bytes) return null;
    if (plus.len == 0 or plus[0] != '+') return null;

    const quality_start = plus_start + plus_raw_len + 1;
    if (sequence_len >= data.len - quality_start) return null;
    const quality_end = quality_start + sequence_len;
    const quality_lf = quality_end + @intFromBool(data[quality_end] == '\r');
    if (quality_lf >= data.len or data[quality_lf] != '\n') return null;
    if (firstInvalidQuality(data[quality_start..quality_end]) != null) return null;

    return .{
        .line_ends = .{ header_end, plus_start - 1, quality_start - 1, quality_lf },
        .use_full_iupac = use_full_iupac,
    };
}

pub const CheckScanner = struct {
    max_line_bytes: usize,
    alphabet: Alphabet,
    use_full_iupac: bool = false,
    machine: Machine = .{},
    line_len: usize = 0,
    first_byte: ?u8 = null,
    second_byte: ?u8 = null,
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
        var consumed: usize = 0;
        var boundary_proved = false;
        var line_feeds: LineFeedSearch = .{};
        while (consumed < data.len) {
            if (self.consumeCompleteRecord(data, consumed, &line_feeds, boundary_proved)) |record_len| {
                boundary_proved = true;
                consumed += record_len;
                continue;
            }
            boundary_proved = false;
            consumed += try self.consumeIncrementalRecord(data, consumed, &line_feeds);
        }
        return data.len;
    }

    fn consumeCompleteRecord(
        self: *CheckScanner,
        data: []const u8,
        start: usize,
        line_feeds: *LineFeedSearch,
        boundary_proved: bool,
    ) ?usize {
        if (!boundary_proved and !self.atRecordBoundary()) return null;

        const complete = scanCompleteRecord(
            data,
            start,
            line_feeds,
            self.max_line_bytes,
            self.alphabet,
            self.use_full_iupac,
        ) orelse return null;
        const record_len = complete.line_ends[3] + 1;
        const next_offset = progressAfter(self.byte_offset, record_len) catch return null;
        const next_record_index = std.math.add(u64, self.record_index, 1) catch return null;

        self.use_full_iupac = complete.use_full_iupac;
        self.record_index = next_record_index;
        self.byte_offset = next_offset;
        self.line_start_offset = next_offset;
        return record_len;
    }

    fn atRecordBoundary(self: *const CheckScanner) bool {
        return self.machine.expected == .header and
            self.line_len == 0 and
            self.first_byte == null and
            self.second_byte == null and
            !self.pending_cr and
            self.sequence_failure == null and
            self.quality_failure == null and
            self.line_start_offset == self.byte_offset;
    }

    fn consumeIncrementalRecord(
        self: *CheckScanner,
        data: []const u8,
        start: usize,
        line_feeds: *LineFeedSearch,
    ) CheckScannerError!usize {
        const initial_record_index = self.record_index;
        var segment_start = start;
        while (line_feeds.find(STRUCTURAL_BLOCK_BYTES, data, segment_start)) |line_end| {
            try self.consumeLineBytes(data[segment_start..line_end]);
            try self.advanceOffset(line_end - segment_start + 1);
            try self.finishLine(true);
            segment_start = line_end + 1;
            if (self.record_index != initial_record_index) return segment_start - start;
        }
        if (segment_start < data.len) {
            try self.consumeLineBytes(data[segment_start..]);
            try self.advanceOffset(data.len - segment_start);
        }
        return data.len - start;
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
        if (self.second_byte == null and previous_len < 2 and self.line_len >= 2) {
            self.second_byte = bytes[1 - previous_len];
        }

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
            .sequence => try classifySemanticBytes(
                .sequence,
                self.alphabet,
                &self.use_full_iupac,
                &self.sequence_failure,
                bytes,
                start_index,
            ),
            .quality => try classifySemanticBytes(
                .quality,
                self.alphabet,
                &self.use_full_iupac,
                &self.quality_failure,
                bytes,
                start_index,
            ),
            .header, .plus => {},
        }
    }

    fn classifyByte(self: *CheckScanner, byte: u8, byte_index: usize) void {
        switch (self.machine.expected) {
            .sequence => classifySemanticByte(
                .sequence,
                self.alphabet,
                &self.use_full_iupac,
                &self.sequence_failure,
                byte,
                byte_index,
            ),
            .quality => classifySemanticByte(
                .quality,
                self.alphabet,
                &self.use_full_iupac,
                &self.quality_failure,
                byte,
                byte_index,
            ),
            .header, .plus => {},
        }
    }

    fn finishLine(self: *CheckScanner, terminated: bool) CheckScannerError!void {
        const content_len = self.line_len - @intFromBool(terminated and self.pending_cr);
        const content_first = if (content_len == 0) null else self.first_byte;
        const content_second = if (content_len < 2) null else self.second_byte;
        const line_kind = self.machine.expected;
        const record_ready = self.machine.push(
            content_len,
            content_first,
            content_second,
        ) catch |err| {
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
        self.second_byte = null;
        self.pending_cr = false;
        self.line_start_offset = self.byte_offset;
    }

    fn finishRecord(self: *CheckScanner) CheckScannerError!void {
        if (self.sequence_failure) |relative_offset| {
            self.last_error = try semanticParseError(
                semanticSequenceError(relative_offset),
                self.record_index,
                self.sequence_start_offset,
            );
            return error.Format;
        }
        if (self.quality_failure) |relative_offset| {
            self.last_error = try semanticParseError(
                semanticQualityError(relative_offset),
                self.record_index,
                self.quality_start_offset,
            );
            return error.Format;
        }
        self.record_index = std.math.add(u64, self.record_index, 1) catch
            return error.ArithmeticLimit;
        self.sequence_failure = null;
        self.quality_failure = null;
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
        self.byte_offset = try progressAfter(self.byte_offset, amount);
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

const ProjectionTestSource = struct {
    data: []const u8,
    pos: usize = 0,
    split: usize,
    split_pending: bool,
    fail_at: ?usize,
    chunk_limit: usize = std.math.maxInt(usize),
    read_count: usize = 0,

    fn init(data: []const u8, split: usize, fail_at: ?usize) ProjectionTestSource {
        return .{
            .data = data,
            .split = split,
            .split_pending = split > 0 and split < data.len,
            .fail_at = fail_at,
        };
    }

    fn byteSource(self: *ProjectionTestSource) ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *ProjectionTestSource = @ptrCast(@alignCast(ctx));
        self.read_count += 1;
        if (self.fail_at) |fail_at| {
            if (self.pos >= fail_at) return error.ReadFailed;
        }
        if (self.pos == self.data.len) return 0;

        var end = self.pos + @min(self.chunk_limit, @min(dest.len, self.data.len - self.pos));
        if (self.split_pending) {
            end = @min(end, self.split);
            self.split_pending = end != self.split;
        }
        if (self.fail_at) |fail_at| end = @min(end, fail_at);
        const bytes = self.data[self.pos..end];
        @memcpy(dest[0..bytes.len], bytes);
        self.pos = end;
        return bytes.len;
    }
};

fn expectProjectionErrorEqual(expected: ?ParseError, actual: ?ParseError) !void {
    if (expected) |expected_error| {
        const actual_error = actual orelse return error.TestExpectedEqual;
        try std.testing.expectEqual(expected_error.code, actual_error.code);
        try std.testing.expectEqualStrings(expected_error.message, actual_error.message);
        try std.testing.expectEqual(expected_error.record_index, actual_error.record_index);
        try std.testing.expectEqual(expected_error.byte_offset, actual_error.byte_offset);
        try std.testing.expectEqual(expected_error.line_in_record, actual_error.line_in_record);
    } else {
        try std.testing.expect(actual == null);
    }
}

fn expectSemanticErrorEqual(expected: ?SemanticError, actual: ?SemanticError) !void {
    if (expected) |expected_error| {
        const actual_error = actual orelse return error.TestExpectedEqual;
        try std.testing.expectEqual(expected_error.code, actual_error.code);
        try std.testing.expectEqualStrings(expected_error.message, actual_error.message);
        try std.testing.expectEqual(expected_error.field, actual_error.field);
        try std.testing.expectEqual(expected_error.byte_index, actual_error.byte_index);
    } else {
        try std.testing.expect(actual == null);
    }
}

fn expectPayloadProjection(
    input: []const u8,
    split: usize,
    fail_at: ?usize,
    options: Options,
    validation_options: ValidationOptions,
) !void {
    var full_source = ProjectionTestSource.init(input, split, fail_at);
    var full_reader = try Reader.init(
        std.testing.allocator,
        full_source.byteSource(),
        options,
    );
    defer full_reader.deinit();

    var payload_source = ProjectionTestSource.init(input, split, fail_at);
    var payload_reader = try Reader.init(
        std.testing.allocator,
        payload_source.byteSource(),
        options,
    );
    defer payload_reader.deinit();

    var record_source = ProjectionTestSource.init(input, split, fail_at);
    var record_reader = try Reader.init(
        std.testing.allocator,
        record_source.byteSource(),
        options,
    );
    defer record_reader.deinit();

    var advance_source = ProjectionTestSource.init(input, split, fail_at);
    var advance_reader = try Reader.init(
        std.testing.allocator,
        advance_source.byteSource(),
        options,
    );
    defer advance_reader.deinit();

    var validated_source = ProjectionTestSource.init(input, split, fail_at);
    var validated_reader = try Reader.init(
        std.testing.allocator,
        validated_source.byteSource(),
        options,
    );
    defer validated_reader.deinit();
    var full_validator = AdaptiveRecordValidator.init(validation_options);
    var projected_validator = AdaptiveRecordValidator.init(validation_options);

    while (true) {
        const full_result = full_reader.next();
        const payload_result = nextPayload(&payload_reader);
        const record_result = nextRecordWithoutId(&record_reader);
        const advance_result = advance_reader.advance();
        const validated_result = nextValidatedHeader(&validated_reader, &projected_validator);
        if (full_result) |full_record| {
            const payload_record = try payload_result;
            const projected_record = try record_result;
            const advanced = try advance_result;
            const validated_record = try validated_result;
            try std.testing.expectEqual(full_record == null, payload_record == null);
            try std.testing.expectEqual(full_record == null, projected_record == null);
            try std.testing.expectEqual(full_record == null, validated_record == null);
            try std.testing.expectEqual(full_record != null, advanced);
            if (full_record) |record| {
                const payload = payload_record.?;
                const projected = projected_record.?;
                const validated = validated_record.?;
                try std.testing.expectEqualStrings(record.sequence, payload.sequence);
                try std.testing.expectEqualStrings(record.quality, payload.quality);
                try std.testing.expectEqualStrings(record.header, projected.header);
                try std.testing.expectEqualStrings(record.sequence, projected.sequence);
                try std.testing.expectEqualStrings(record.plus, projected.plus);
                try std.testing.expectEqualStrings(record.quality, projected.quality);
                try std.testing.expectEqual(@as(usize, 0), projected.id.len);
                try std.testing.expectEqualStrings(record.header, validated.header);
                try expectSemanticErrorEqual(
                    full_validator.validate(record),
                    validated.semantic_error,
                );
            }
            try std.testing.expectEqual(full_reader.recordIndex(), payload_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), payload_reader.byteOffset());
            try std.testing.expectEqual(
                full_reader.currentRecordOffsets(),
                payload_reader.currentRecordOffsets(),
            );
            try std.testing.expectEqual(full_reader.recordIndex(), record_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), record_reader.byteOffset());
            try std.testing.expectEqual(
                full_reader.currentRecordOffsets(),
                record_reader.currentRecordOffsets(),
            );
            try std.testing.expectEqual(full_reader.recordIndex(), advance_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), advance_reader.byteOffset());
            try std.testing.expectEqual(full_reader.recordIndex(), validated_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), validated_reader.byteOffset());
            try std.testing.expectEqual(
                full_reader.currentRecordOffsets(),
                validated_reader.currentRecordOffsets(),
            );
            if (full_record == null) return;
        } else |expected_error| {
            try std.testing.expectError(expected_error, payload_result);
            try std.testing.expectError(expected_error, record_result);
            try std.testing.expectError(expected_error, advance_result);
            try std.testing.expectError(expected_error, validated_result);
            try std.testing.expectEqual(full_reader.recordIndex(), payload_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), payload_reader.byteOffset());
            try std.testing.expectEqual(
                full_reader.currentRecordOffsets(),
                payload_reader.currentRecordOffsets(),
            );
            const full_error = full_reader.takeLastError();
            try expectProjectionErrorEqual(full_error, payload_reader.takeLastError());
            try std.testing.expectEqual(full_reader.recordIndex(), record_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), record_reader.byteOffset());
            try std.testing.expectEqual(
                full_reader.currentRecordOffsets(),
                record_reader.currentRecordOffsets(),
            );
            try expectProjectionErrorEqual(full_error, record_reader.takeLastError());
            try std.testing.expectEqual(full_reader.recordIndex(), advance_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), advance_reader.byteOffset());
            try expectProjectionErrorEqual(full_error, advance_reader.takeLastError());
            try std.testing.expectEqual(full_reader.recordIndex(), validated_reader.recordIndex());
            try std.testing.expectEqual(full_reader.byteOffset(), validated_reader.byteOffset());
            try std.testing.expectEqual(
                full_reader.currentRecordOffsets(),
                validated_reader.currentRecordOffsets(),
            );
            try expectProjectionErrorEqual(full_error, validated_reader.takeLastError());
            return;
        }
    }
}

fn expectBufferedLineProgress(
    input: []const u8,
    max_line_bytes: usize,
    progress: struct { records: u64 = 0, bytes: u64 = 0 },
    comptime delivery: enum { record, validated, advance },
) !void {
    errdefer std.debug.print("buffered lines: input {any}, limit {d}, progress {any}, delivery {t}\n", .{
        input, max_line_bytes, progress, delivery,
    });
    var reference_source = ProjectionTestSource.init(input, 0, null);
    var source = reference_source;
    var reference = try Reader.init(std.testing.allocator, reference_source.byteSource(), .{
        .max_line_bytes = max_line_bytes,
    });
    defer reference.deinit();
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 1 });
    var reader = try Reader.init(failing.allocator(), source.byteSource(), .{
        .max_line_bytes = max_line_bytes,
    });
    defer reader.deinit();
    for ([_]*Reader{ &reference, &reader }) |current| {
        current.record_index = progress.records;
        current.byte_offset = progress.bytes;
        current.machine.sequence_len = 123;
        current.record_offsets = .{ .header = 1, .sequence = 2, .plus = 3, .quality = 4 };
        current.current_record_offsets = current.record_offsets;
        current.last_error = .{
            .code = .s005_length_mismatch,
            .message = "earlier error",
            .record_index = 0,
            .byte_offset = 0,
            .line_in_record = 4,
        };
        // Complete records avoid the different progress rules for an incomplete line.
        try std.testing.expect(try current.refill());
        try std.testing.expectEqual(input.len, current.fill_end);
    }
    var reference_validator = AdaptiveRecordValidator.init(.{});
    var validator = AdaptiveRecordValidator.init(.{});
    while (true) {
        const record_start = reader.cursor;
        reference.beginRecord();
        var finished = false;
        if (delivery == .advance) {
            const expected = reference.advanceFallback();
            const actual = reader.advance();
            if (expected) |has_record| {
                try std.testing.expectEqual(has_record, try actual);
                finished = !has_record;
            } else |err| {
                try std.testing.expectError(err, actual);
                finished = true;
            }
        } else {
            const expected = reference.nextFallback(false);
            var span: ?[]const u8 = null;
            var semantic_error: ?SemanticError = null;
            const actual: ReaderError!?Record = if (delivery == .validated) result: {
                const validated = nextValidatedRecord(&reader, &validator) catch |err|
                    break :result err;
                if (validated) |record| {
                    span = record.canonical_span;
                    semantic_error = record.semantic_error;
                    break :result record.record;
                }
                break :result null;
            } else nextWithoutId(&reader, &span);
            if (expected) |record| {
                try std.testing.expectEqualDeep(record, try actual);
                if (record) |fields| {
                    if (delivery == .validated) {
                        try expectSemanticErrorEqual(reference_validator.validate(fields), semantic_error);
                    }
                    const raw = input[record_start..reader.cursor];
                    var canonical = true;
                    var line_start: usize = 0;
                    for (raw, 0..) |byte, index| {
                        if (byte != '\n') continue;
                        if (index > line_start and raw[index - 1] == '\r') canonical = false;
                        line_start = index + 1;
                    }
                    if (canonical) {
                        try std.testing.expectEqualStrings(raw, span.?);
                    } else try std.testing.expect(span == null);
                } else finished = true;
            } else |err| {
                try std.testing.expectError(err, actual);
                try std.testing.expect(span == null);
                finished = true;
            }
        }
        try std.testing.expectEqual(reference.cursor, reader.cursor);
        try std.testing.expectEqual(reference.fill_end, reader.fill_end);
        try std.testing.expectEqual(reference.recordIndex(), reader.recordIndex());
        try std.testing.expectEqual(reference.byteOffset(), reader.byteOffset());
        try std.testing.expectEqual(reference.currentRecordOffsets(), reader.currentRecordOffsets());
        try std.testing.expectEqualDeep(reference.record_offsets, reader.record_offsets);
        try std.testing.expectEqualDeep(reference.machine, reader.machine);
        try std.testing.expectEqualDeep(reference_validator, validator);
        try expectProjectionErrorEqual(reference.last_error, reader.last_error);
        try std.testing.expectEqual(reference_source.pos, source.pos);
        try std.testing.expectEqual(reference_source.read_count, source.read_count);
        try std.testing.expect(!failing.has_induced_failure);
        if (finished) return;
    }
}

test "[property] - [reader]: complete buffered records preserve per-line results and progress" {
    const cases = [_][4][]const u8{
        .{ "@r note", "AC", "+note", "!~" },
        .{ "@r", "R", "+", "~" },
        .{ "@r", "", "+", "" },
        .{ "@r", "A\rA", "+text\rinside", "!\r!" },
        .{ "@\x00", "?\xff", "+", " \x7f" },
        .{ "", "AC", "+", "!!" },
        .{ "@", "AC", "+", "!!" },
        .{ "@ r", "AC", "+", "!!" },
        .{ "@\tr", "AC", "+", "!!" },
        .{ "bad", "AC", "bad", "!" },
        .{ "@r", "AC", "", "!!" },
        .{ "@r", "AC", "bad", "!" },
        .{ "@r", "AC", "+", "!" },
        .{ "@r", "A\r", "+", "!!" },
        .{ "@r", "AA", "+", "!\r" },
    };
    const prefix = "@previous\nN\n+\nI\n";
    var storage: [128]u8 = undefined;
    @memcpy(storage[0..prefix.len], prefix);
    for (cases) |fields| {
        for (0..16) |mask| {
            var end: usize = prefix.len;
            for (fields, 0..) |field, index| {
                @memcpy(storage[end..][0..field.len], field);
                end += field.len;
                if (mask & (@as(usize, 1) << @intCast(index)) != 0) {
                    storage[end] = '\r';
                    end += 1;
                }
                storage[end] = '\n';
                end += 1;
            }
            const input = storage[prefix.len..end];
            inline for (.{ .record, .validated, .advance }) |delivery| {
                for ([_]usize{ 0, 1, 2, 3, 4, 7, std.math.maxInt(usize) }) |limit| {
                    try expectBufferedLineProgress(input, limit, .{}, delivery);
                }
                try expectBufferedLineProgress(storage[0..end], 128, .{}, delivery);
                try expectBufferedLineProgress(input, 128, .{ .bytes = (1 << 32) - 3 }, delivery);
                for ([_]u64{ (1 << 32) - 1, std.math.maxInt(u64) - 1, std.math.maxInt(u64) }) |count| {
                    try expectBufferedLineProgress(input, 128, .{ .records = count }, delivery);
                }
                if (mask == 0 or mask == 15) {
                    for (0..input.len + 1) |remaining| {
                        for ([_]usize{ 2, 128 }) |limit| {
                            try expectBufferedLineProgress(input, limit, .{
                                .bytes = std.math.maxInt(u64) - remaining,
                            }, delivery);
                        }
                    }
                }
            }
        }
    }
}

fn readPredictedPayloadForTest(reader: *Reader) ReaderError!?RecordPayload {
    const predicted = try nextPredictedPayload(reader) orelse return null;
    for (predicted.payload.quality) |byte| {
        if (byte < 33 or byte > 126) {
            if (predicted.checkpoint) |checkpoint| return try checkpoint.reread(reader);
            break;
        }
    }
    return predicted.payload;
}

fn expectPredictedPayloads(
    input: []const u8,
    split: usize,
    chunk_limit: usize,
    fail_at: ?usize,
    options: Options,
    progress: struct { records: u64 = 0, bytes: u64 = 0 },
) !void {
    var reference_source = ProjectionTestSource.init(input, split, fail_at);
    reference_source.chunk_limit = chunk_limit;
    var reference = try Reader.init(std.testing.allocator, reference_source.byteSource(), options);
    defer reference.deinit();
    var source = ProjectionTestSource.init(input, split, fail_at);
    source.chunk_limit = chunk_limit;
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), options);
    defer reader.deinit();
    for ([_]*Reader{ &reference, &reader }) |current| {
        current.record_index = progress.records;
        current.byte_offset = progress.bytes;
        current.record_offsets = .{ .header = 1, .sequence = 2, .plus = 3, .quality = 4 };
        current.last_error = .{
            .code = .s005_length_mismatch,
            .message = "earlier error",
            .record_index = 0,
            .byte_offset = 0,
            .line_in_record = 4,
        };
    }

    while (true) {
        const expected_result = nextPayload(&reference);
        const actual_result = readPredictedPayloadForTest(&reader);
        var finished = false;
        if (expected_result) |expected| {
            try std.testing.expectEqualDeep(expected, try actual_result);
            finished = expected == null;
        } else |err| {
            try std.testing.expectError(err, actual_result);
            finished = true;
        }
        try std.testing.expectEqual(reference.cursor, reader.cursor);
        try std.testing.expectEqual(reference.fill_end, reader.fill_end);
        try std.testing.expectEqual(reference.recordIndex(), reader.recordIndex());
        try std.testing.expectEqual(reference.byteOffset(), reader.byteOffset());
        try std.testing.expectEqual(reference.currentRecordOffsets(), reader.currentRecordOffsets());
        try std.testing.expectEqualDeep(reference.record_offsets, reader.record_offsets);
        try std.testing.expectEqualDeep(reference.machine, reader.machine);
        try expectProjectionErrorEqual(reference.last_error, reader.last_error);
        try std.testing.expectEqual(reference_source.pos, source.pos);
        try std.testing.expectEqual(reference_source.read_count, source.read_count);
        if (finished) return;
    }
}

test "[property] - [reader]: predicted payloads preserve parsing and source progress" {
    const first = "@r\nAC\n+\n!~\n";
    const unlimited = std.math.maxInt(usize);
    for ([_][]const u8{
        first,
        "@empty\n\n+\n\n",
        "@r\r\nac?\r\n+note\r\n!~I\r\n",
        "@r\n\x00\x7f\r\n+\r\n!~\n",
        "@r\nAC\n+\n!!",
        "@r\nAC\n+\n \n\n",
        "@r\nACG\n+\n!\n!\n",
        "@r\nACG\n+\n!\n\r\n",
        "@r\nACG\n+\n!\n!\r\n",
        "@r\nAC\n+\n!\x7f\r\n",
        "@r\nAC\n+\n!\r\r\n",
        "@r\nAC\n+\n!\r\n",
        "@r\nAC\n+\n!!!\n",
        "bad\nAC\n+\n!~\n",
        "@r\nAC\nbad\n!~\n",
        first ++ "@r\nACG\n+\n!\n!\n",
        first ++ "@r\r\nAC\r\n+\r\n!\x7f\r\n",
    }) |input| {
        for (0..input.len + 1) |position| {
            try expectPredictedPayloads(input, position, unlimited, null, .{}, .{});
            try expectPredictedPayloads(input[0..position], 0, unlimited, null, .{}, .{});
            try expectPredictedPayloads(input, 0, unlimited, position, .{}, .{});
        }
        for (1..input.len + 1) |chunk_limit| {
            try expectPredictedPayloads(input, 0, chunk_limit, null, .{}, .{});
        }
        for (0..9) |limit| {
            try expectPredictedPayloads(input, 0, unlimited, null, .{ .max_line_bytes = limit }, .{});
        }
        for ([_]u64{ (1 << 32) - 1, std.math.maxInt(u64) - 1, std.math.maxInt(u64) }) |counter| {
            try expectPredictedPayloads(input, 0, unlimited, null, .{}, .{ .records = counter });
        }
        for (0..input.len + 1) |remaining| {
            try expectPredictedPayloads(input, 0, unlimited, null, .{}, .{
                .bytes = std.math.maxInt(u64) - remaining,
            });
        }
        try expectPredictedPayloads(input, 0, unlimited, null, .{}, .{ .bytes = (1 << 32) - 1 });
    }
}

test "[property] - [reader]: predicted payloads retain refill and oversized fallbacks" {
    const capacity = io_layer.DEFAULT_READER_BUFFER_BYTES;
    for ([_]usize{ capacity / 2 - 8, capacity / 2, capacity + 1 }) |len| {
        for ([_][]const u8{ "\n", "\r\n" }) |ending| {
            const input = try ReaderSpillFixture.init(std.testing.allocator, len, len, ending, true);
            defer std.testing.allocator.free(input);
            for ([_]usize{ 0, 5, len + 5, capacity - 1, capacity, input.len - 1 }) |split| {
                try expectPredictedPayloads(input, split, capacity, null, .{}, .{});
            }
            try expectPredictedPayloads(input, 0, capacity - 1, capacity, .{}, .{});
            try expectPredictedPayloads(input, 0, capacity, null, .{ .max_line_bytes = len - 1 }, .{});
            try expectPredictedPayloads(input[0 .. input.len - ending.len], 0, capacity, null, .{}, .{});
        }
    }
}

test "[unit] - [reader]: prediction and retry use the current buffer without allocating" {
    for ([_][]const u8{ "\n", "\r\n" }) |ending| {
        for ([_]usize{ 0, 1, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129 }) |len| {
            const input = try ReaderSpillFixture.init(std.testing.allocator, len, len, ending, true);
            defer std.testing.allocator.free(input);
            var source = ProjectionTestSource.init(input, 0, null);
            var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 1 });
            var reader = try Reader.init(failing.allocator(), source.byteSource(), .{});
            defer reader.deinit();
            const predicted = (try nextPredictedPayload(&reader)).?;
            try std.testing.expect(predicted.checkpoint != null);
            try std.testing.expectEqual(len, predicted.payload.sequence.len);
            try std.testing.expectEqual(len, predicted.payload.quality.len);
            const buffer = reader.buf.ptr;
            const offset = reader.byteOffset();
            const actual = try predicted.checkpoint.?.reread(&reader);
            try std.testing.expectEqualDeep(predicted.payload, actual);
            try std.testing.expectEqual(buffer, reader.buf.ptr);
            try std.testing.expectEqual(offset, reader.byteOffset());
            try std.testing.expectEqual(@as(usize, 1), source.read_count);
            try std.testing.expect(!failing.has_induced_failure);
        }
    }

    const input = "@r\nACG\n+\n!\n!\n";
    var source = ProjectionTestSource.init(input, 0, null);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 1 });
    var reader = try Reader.init(failing.allocator(), source.byteSource(), .{});
    defer reader.deinit();
    const predicted = (try nextPredictedPayload(&reader)).?;
    try std.testing.expectEqualStrings("!\n!", predicted.payload.quality);
    try std.testing.expectError(error.S005LengthMismatch, predicted.checkpoint.?.reread(&reader));
    try std.testing.expectEqual(@as(u64, 9), reader.takeLastError().?.byte_offset);
    try std.testing.expectEqual(@as(u64, 11), reader.byteOffset());
    try std.testing.expectEqual(@as(usize, 1), source.read_count);
    try std.testing.expect(!failing.has_induced_failure);
}

fn expectValidatedProjection(
    input: []const u8,
    split: usize,
    fail_at: ?usize,
    options: Options,
    validation_options: ValidationOptions,
    progress: struct { records: u64 = 0, bytes: u64 = 0 },
) !void {
    var reference_source = ProjectionTestSource.init(input, split, fail_at);
    var reference = try Reader.init(std.testing.allocator, reference_source.byteSource(), options);
    defer reference.deinit();
    var source = ProjectionTestSource.init(input, split, fail_at);
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), options);
    defer reader.deinit();
    reference.record_index = progress.records;
    reference.byte_offset = progress.bytes;
    reader.record_index = progress.records;
    reader.byte_offset = progress.bytes;
    var reference_validator = AdaptiveRecordValidator.init(validation_options);
    var validator = AdaptiveRecordValidator.init(validation_options);

    while (true) {
        var expected_span: ?[]const u8 = null;
        const expected_result = nextWithoutId(&reference, &expected_span);
        const actual_result = nextValidatedRecord(&reader, &validator);
        var finished = false;
        if (expected_result) |expected_record| {
            const actual = try actual_result;
            try std.testing.expectEqual(expected_record == null, actual == null);
            if (expected_record) |record| {
                try std.testing.expectEqualDeep(record, actual.?.record);
                try std.testing.expectEqual(expected_span == null, actual.?.canonical_span == null);
                if (expected_span) |span| {
                    try std.testing.expectEqualStrings(span, actual.?.canonical_span.?);
                }
                try expectSemanticErrorEqual(reference_validator.validate(record), actual.?.semantic_error);
            } else finished = true;
        } else |err| {
            try std.testing.expectError(err, actual_result);
            finished = true;
        }
        try std.testing.expectEqualDeep(reference_validator, validator);
        try std.testing.expectEqual(reference.recordIndex(), reader.recordIndex());
        try std.testing.expectEqual(reference.byteOffset(), reader.byteOffset());
        try std.testing.expectEqual(reference.currentRecordOffsets(), reader.currentRecordOffsets());
        try std.testing.expectEqualDeep(reference.machine, reader.machine);
        try expectProjectionErrorEqual(reference.takeLastError(), reader.takeLastError());
        try std.testing.expectEqual(reference_source.pos, source.pos);
        if (finished) return;
    }
}

fn readSpillForAllocationCheck(allocator: std.mem.Allocator, input: []const u8, chunk_limit: usize) !void {
    var source = ProjectionTestSource.init(input, 0, null);
    source.chunk_limit = chunk_limit;
    var reader = try Reader.init(allocator, source.byteSource(), .{});
    defer reader.deinit();
    _ = try reader.next();
}

fn expectRefillDelivery(
    input: []const u8,
    split: usize,
    chunk_limit: usize,
    fail_at: ?usize,
    options: Options,
    progress: struct { records: u64 = 0, bytes: u64 = 0 },
) !void {
    errdefer std.debug.print("refill comparison: bytes {d}, split {d}, chunk {d}, failure {?d}, limit {d}, progress {any}\n", .{
        input.len, split, chunk_limit, fail_at, options.max_line_bytes, progress,
    });
    var reference_source = ProjectionTestSource.init(input, split, fail_at);
    reference_source.chunk_limit = chunk_limit;
    var source = reference_source;
    var reference = try Reader.init(std.testing.allocator, reference_source.byteSource(), options);
    defer reference.deinit();
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), options);
    defer reader.deinit();
    reference.record_index = progress.records;
    reader.record_index = progress.records;
    reference.byte_offset = progress.bytes;
    reader.byte_offset = progress.bytes;
    var reference_validator = AdaptiveRecordValidator.init(.{});
    var validator = AdaptiveRecordValidator.init(.{});

    while (true) {
        const record_start = reader.byteOffset() - progress.bytes;
        var reference_span: ?[]const u8 = null;
        reference.beginRecord();
        // Keep the pre-refill parser as a reference for fields and failures.
        const expected: ReaderError!?Record = result: {
            const buffered = reference.readBufferedRecord(.full, true, null, false) catch |err|
                break :result err;
            break :result switch (buffered) {
                .incomplete => reference.nextFallback(false),
                .eof => null,
                .record => |record| reference.finishBufferedRecord(record, &reference_span, false),
            };
        };
        const actual = nextValidatedRecord(&reader, &validator);
        var finished = false;
        if (expected) |expected_record| {
            const delivered = try actual;
            try std.testing.expectEqual(expected_record == null, delivered == null);
            if (expected_record) |record| {
                try std.testing.expectEqualDeep(record, delivered.?.record);
                try expectSemanticErrorEqual(reference_validator.validate(record), delivered.?.semantic_error);
                if (reference_span) |span| try std.testing.expectEqualStrings(span, delivered.?.canonical_span.?);
                if (delivered.?.canonical_span) |span| {
                    try std.testing.expectEqualStrings(input[@intCast(record_start)..@intCast(reader.byteOffset() - progress.bytes)], span);
                }
            } else finished = true;
        } else |err| {
            try std.testing.expectError(err, actual);
            finished = true;
        }
        try std.testing.expectEqual(reference.recordIndex(), reader.recordIndex());
        try std.testing.expectEqual(reference.byteOffset(), reader.byteOffset());
        try std.testing.expectEqual(reference.currentRecordOffsets(), reader.currentRecordOffsets());
        try std.testing.expectEqualDeep(reference.machine, reader.machine);
        try std.testing.expectEqualDeep(reference_validator, validator);
        try expectProjectionErrorEqual(reference.takeLastError(), reader.takeLastError());
        try std.testing.expectEqual(reference_source.pos, source.pos);
        try std.testing.expectEqual(reference_source.read_count, source.read_count);
        if (finished) return;
    }
}

fn expectPairedRefill(
    comptime header_only: bool,
    input: []const u8,
    split: usize,
    chunk_limit: usize,
    fail_at: ?usize,
    options: Options,
    progress: struct { records: u64 = 0, bytes: u64 = 0 },
) !void {
    errdefer std.debug.print("pair refill: header {}, bytes {d}, split {d}, chunk {d}, failure {?d}, progress {any}\n", .{
        header_only, input.len, split, chunk_limit, fail_at, progress,
    });
    var reference_source = ProjectionTestSource.init(input, split, fail_at);
    reference_source.chunk_limit = chunk_limit;
    var source = reference_source;
    var reference = try Reader.init(std.testing.allocator, reference_source.byteSource(), options);
    defer reference.deinit();
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), options);
    defer reader.deinit();
    reference.record_index = progress.records;
    reader.record_index = progress.records;
    reference.byte_offset = progress.bytes;
    reader.byte_offset = progress.bytes;
    var reference_validator = AdaptiveRecordValidator.init(.{});
    var validator = AdaptiveRecordValidator.init(.{});

    while (true) {
        var finished = false;
        pair: {
            const expected_first = nextValidatedRecord(&reference, &reference_validator);
            const actual_first = nextValidatedRecord(&reader, &validator);
            const first = expected_first catch |err| {
                try std.testing.expectError(err, actual_first);
                finished = true;
                break :pair;
            };
            var delivered = try actual_first;
            try std.testing.expectEqual(first == null, delivered == null);
            if (first == null) {
                finished = true;
                break :pair;
            }
            // The reference owns mate 1 before an ordinary read invalidates it.
            var saved = try toOwned(std.testing.allocator, first.?.record);
            defer saved.deinit();
            const expected_second = nextValidatedRecord(&reference, &reference_validator);
            const paired_second = if (header_only)
                nextPairedValidatedHeader(&reader, &delivered.?.record.header, &validator)
            else
                nextPairedValidatedRecord(&reader, &delivered.?, &validator);
            const actual_second: ReaderError!?ValidatedRecord = second: {
                const buffered = paired_second catch |err| break :second err;
                try std.testing.expectEqualStrings(saved.header, delivered.?.record.header);
                if (!header_only) {
                    inline for (.{ "id", "sequence", "plus", "quality" }) |field| {
                        try std.testing.expectEqualStrings(@field(saved, field), @field(delivered.?.record, field));
                    }
                    if (delivered.?.canonical_span) |span| {
                        const output = try std.mem.concat(std.testing.allocator, u8, &.{
                            "@", saved.header, "\n", saved.sequence, "\n+", saved.plus, "\n", saved.quality, "\n",
                        });
                        defer std.testing.allocator.free(output);
                        try std.testing.expectEqualStrings(output, span);
                    }
                }
                break :second buffered orelse nextValidatedRecord(&reader, &validator);
            };
            if (expected_second) |second| {
                const actual = try actual_second;
                try std.testing.expectEqual(second == null, actual == null);
                if (second) |expected| {
                    try std.testing.expectEqualDeep(expected.record, actual.?.record);
                    try expectSemanticErrorEqual(expected.semantic_error, actual.?.semantic_error);
                } else finished = true;
            } else |err| {
                try std.testing.expectError(err, actual_second);
                finished = true;
            }
        }
        try std.testing.expectEqual(reference.recordIndex(), reader.recordIndex());
        try std.testing.expectEqual(reference.byteOffset(), reader.byteOffset());
        try std.testing.expectEqual(reference.currentRecordOffsets(), reader.currentRecordOffsets());
        try std.testing.expectEqualDeep(reference.machine, reader.machine);
        try std.testing.expectEqualDeep(reference_validator, validator);
        try expectProjectionErrorEqual(reference.takeLastError(), reader.takeLastError());
        try std.testing.expectEqual(reference_source.pos, source.pos);
        try std.testing.expectEqual(reference_source.read_count, source.read_count);
        if (finished) return;
    }
}

test "[property] - [reader]: paired refills preserve fields, progress, and failures" {
    const first = "@r/1\nAR\n+left\n!~\n";
    const first_crlf = "@r/1\r\nAR\r\n+left\r\n!~\r\n";
    inline for (.{ false, true }) |header_only| {
        for ([_][]const u8{
            first ++ "@r/2\nC\n+right\n#\n" ++ first,
            first_crlf ++ "@r/2\r\nC\r\n+right\r\n#\r\n" ++ first,
            first ++ "@r/2\n\n+\n\n",
            first_crlf ++ "@r/2\nR\n+right\n \n",
            first ++ "@r/2\n?\n-\n \n",
            first ++ "@r/2\nAC\n+\n!\n",
            first ++ "@r/2\r\r\nR\rA\n+\r\r\n!!!\n",
            "@r/1\r\r\nAR\n+left\r\r\n!~\n" ++ first,
        }) |input| {
            for (0..input.len + 1) |split| {
                try expectPairedRefill(header_only, input, split, std.math.maxInt(usize), null, .{}, .{});
                try expectPairedRefill(header_only, input[0..split], 0, 3, null, .{}, .{});
                try expectPairedRefill(header_only, input, 0, 3, split, .{}, .{});
            }
            for (1..input.len + 1) |chunk_limit| {
                try expectPairedRefill(header_only, input, 0, chunk_limit, null, .{ .max_line_bytes = 8 }, .{});
            }
        }
        for ([_]u64{ (1 << 32) - 3, std.math.maxInt(u64) - first.len, std.math.maxInt(u64) - 1 }) |counter| {
            try expectPairedRefill(header_only, first ** 2, first.len, 3, null, .{}, .{ .bytes = counter });
            try expectPairedRefill(header_only, first ** 2, first.len, 3, null, .{}, .{ .records = counter });
        }
    }
}

test "[property] - [reader]: paired refills retain oversized and EOF fallbacks" {
    const capacity = io_layer.DEFAULT_READER_BUFFER_BYTES;
    const first = "@r/1\nAR\n+left\n!~\n";
    const last = "@last\nN\n+\n!\n";
    for ([_]usize{ capacity - first.len - 1, capacity - first.len, capacity - first.len + 1, capacity - 3, capacity - 2, capacity + 1, 2 * capacity + 11 }) |record_len| {
        const field_len = (record_len - 9) / 2;
        const second = try ReaderSpillFixture.init(std.testing.allocator, field_len, field_len, "\n", true);
        defer std.testing.allocator.free(second);
        const padding = if (second.len == record_len) "" else "x";
        const input = try std.mem.concat(std.testing.allocator, u8, &.{ first, second[0..4], padding, second[4..], last });
        defer std.testing.allocator.free(input);
        inline for (.{ false, true }) |header_only| {
            if (record_len + (if (header_only) @as(usize, 3) else first.len) <= capacity) {
                var source = ProjectionTestSource.init(input, first.len, null);
                source.chunk_limit = 17;
                var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
                var reader = try Reader.init(failing.allocator(), source.byteSource(), .{});
                defer reader.deinit();
                var validator = AdaptiveRecordValidator.init(.{});
                var mate1 = (try nextValidatedRecord(&reader, &validator)).?;
                const mate2 = (try if (header_only)
                    nextPairedValidatedHeader(&reader, &mate1.record.header, &validator)
                else
                    nextPairedValidatedRecord(&reader, &mate1, &validator)).?;
                try std.testing.expectEqualStrings("r/1", mate1.record.header);
                try std.testing.expectEqual(field_len, mate2.record.sequence.len);
                try std.testing.expectEqualStrings(input[first.len .. first.len + record_len], mate2.canonical_span.?);
                try std.testing.expect(!failing.has_induced_failure);
            }
            for ([_]usize{ first.len, capacity - 1, capacity, capacity + 1 }) |split| {
                try expectPairedRefill(header_only, input, split, std.math.maxInt(usize), null, .{}, .{});
            }
            try expectPairedRefill(header_only, input, first.len, 17, null, .{}, .{});
            try expectPairedRefill(header_only, input[0 .. first.len + record_len - 1], first.len, 37, null, .{}, .{});
            try expectPairedRefill(header_only, input, first.len, 37, null, .{ .max_line_bytes = field_len - 1 }, .{});
        }
    }
}

test "[failure] - [reader]: lazy refill allocation preserves progress and buffered mates" {
    const first = "@r/1\nAR\n+left\n!~\n";
    const second = "@r/2\nC\n+right\n#\n";
    for (0..3) |delivery| {
        var source = ProjectionTestSource.init(first ++ second ** 2, first.len + 5, null);
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 1 });
        var reader = try Reader.init(failing.allocator(), source.byteSource(), .{ .max_line_bytes = 8 });
        defer reader.deinit();
        var validator = AdaptiveRecordValidator.init(.{});
        var mate1 = (try nextValidatedRecord(&reader, &validator)).?;
        source.chunk_limit = 3;
        const position = source.pos;
        const reads = source.read_count;

        try std.testing.expect((try nextBufferedValidatedRecord(&reader, &validator)) == null);
        try std.testing.expect(!failing.has_induced_failure);
        const result = switch (delivery) {
            0 => nextValidatedRecord(&reader, &validator),
            1 => nextPairedValidatedRecord(&reader, &mate1, &validator),
            2 => nextPairedValidatedHeader(&reader, &mate1.record.header, &validator),
            else => unreachable,
        };
        try std.testing.expectError(error.OutOfMemory, result);
        try std.testing.expect(failing.has_induced_failure);
        try std.testing.expectEqual(@as(usize, 0), reader.fallback_fields[1].storage.len);
        try std.testing.expectEqual(position, source.pos);
        try std.testing.expectEqual(reads, source.read_count);
        try std.testing.expectEqual(@as(u64, first.len), reader.byteOffset());
        try std.testing.expectEqual(@as(u64, 1), reader.recordIndex());
        try std.testing.expect(reader.takeLastError() == null);
        try std.testing.expectEqualStrings("r/1", mate1.record.header);
        try std.testing.expectEqualStrings("AR", mate1.record.sequence);
        try std.testing.expectEqualStrings("left", mate1.record.plus);
        try std.testing.expectEqualStrings("!~", mate1.record.quality);

        failing.fail_index = std.math.maxInt(usize);
        const mate2 = (try switch (delivery) {
            0 => nextValidatedRecord(&reader, &validator),
            1 => nextPairedValidatedRecord(&reader, &mate1, &validator),
            2 => nextPairedValidatedHeader(&reader, &mate1.record.header, &validator),
            else => unreachable,
        }).?;
        try std.testing.expectEqualStrings(second, mate2.canonical_span.?);
        try std.testing.expect(mate2.semantic_error == null);
        if (delivery != 0) try std.testing.expectEqualStrings("r/1", mate1.record.header);
        if (delivery == 1) {
            try std.testing.expectEqualStrings(first, mate1.canonical_span.?);
            try std.testing.expectEqualStrings("AR", mate1.record.sequence);
            try std.testing.expectEqualStrings("!~", mate1.record.quality);
        }
        try std.testing.expectEqual(@as(usize, 2), failing.allocations);
        try std.testing.expectEqual(2 * io_layer.DEFAULT_READER_BUFFER_BYTES, failing.allocated_bytes);
        const reserve = reader.fallback_fields[1].storage.ptr;
        try std.testing.expectEqualStrings(second, (try nextValidatedRecord(&reader, &validator)).?.canonical_span.?);
        try std.testing.expect((try nextValidatedRecord(&reader, &validator)) == null);
        try std.testing.expectEqual(reserve, reader.fallback_fields[1].storage.ptr);
        try std.testing.expectEqual(@as(usize, 2), failing.allocations);
    }
}

test "[edge] - [reader]: refill reserve grows from projected fields and retains earlier capacity" {
    const old = "@old\nAC\n+\n!!\n";
    const first = "@r/1\nAR\n+left\n!~\n";
    const second = "@r/2\nC\n+right\n#\n";
    var source = ProjectionTestSource.init(old ++ first ++ second, old.len - 1, null);
    source.chunk_limit = old.len - 1;
    var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{ .resize_fail_index = 0 });
    var reader = try Reader.init(tracking.allocator(), source.byteSource(), .{ .max_line_bytes = 8 });
    defer reader.deinit();
    const payload = (try reader.nextPayload()).?;
    try std.testing.expectEqualStrings("AC", payload.sequence);
    try std.testing.expectEqualStrings("!!", payload.quality);
    try std.testing.expectEqual(@as(usize, 2), reader.fallback_fields[1].storage.len);
    try std.testing.expectEqual(@as(usize, 3), reader.fallback_fields[3].storage.len);
    const sequence = reader.fallback_fields[1].storage.ptr;
    const quality = reader.fallback_fields[3].storage.ptr;
    const position = source.pos;

    tracking.fail_index = tracking.allocations;
    var validator = AdaptiveRecordValidator.init(.{});
    try std.testing.expectError(error.OutOfMemory, nextValidatedRecord(&reader, &validator));
    try std.testing.expect(tracking.has_induced_failure);
    try std.testing.expectEqual(sequence, reader.fallback_fields[1].storage.ptr);
    try std.testing.expectEqualStrings("AC", payload.sequence);
    try std.testing.expectEqualStrings("!!", payload.quality);
    try std.testing.expectEqual(position, source.pos);
    try std.testing.expectEqual(@as(u64, old.len), reader.byteOffset());
    try std.testing.expectEqual(@as(u64, 1), reader.recordIndex());

    tracking.fail_index = std.math.maxInt(usize);
    source.chunk_limit = 3;
    var mate1 = (try nextValidatedRecord(&reader, &validator)).?;
    const allocations = tracking.allocations;
    const mate2 = (try nextPairedValidatedRecord(&reader, &mate1, &validator)).?;
    try std.testing.expectEqualStrings(first, mate1.canonical_span.?);
    try std.testing.expectEqualStrings(second, mate2.canonical_span.?);
    try std.testing.expect(mate1.semantic_error == null and mate2.semantic_error == null);
    try std.testing.expectEqual(allocations, tracking.allocations);
    try std.testing.expectEqual(quality, reader.fallback_fields[3].storage.ptr);
    try std.testing.expectEqual(@as(usize, 3), reader.fallback_fields[3].storage.len);
    try std.testing.expectEqual(2 * io_layer.DEFAULT_READER_BUFFER_BYTES + 3, tracking.allocated_bytes - tracking.freed_bytes);
    try std.testing.expect((try nextValidatedRecord(&reader, &validator)) == null);
}

test "[property] - [reader]: structural masks match scalar line boundaries" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse 1;
    var storage: [5 * vector_len + 8]u8 = undefined;

    for (0..vector_len) |alignment| {
        for (0..vector_len) |first_len| {
            @memset(&storage, 'x');
            const bytes = storage[alignment..];
            const line_lengths = [4]usize{ first_len, 0, vector_len - 1, vector_len };
            var input_len: usize = 0;
            for (line_lengths) |line_len| {
                input_len += line_len;
                bytes[input_len] = '\n';
                input_len += 1;
            }
            const input = bytes[0..input_len];

            var expected: [4]usize = undefined;
            var search_start: usize = 0;
            for (&expected) |*line_end| {
                const relative = std.mem.findScalar(
                    u8,
                    input[search_start..],
                    '\n',
                ).?;
                line_end.* = search_start + relative;
                search_start = line_end.* + 1;
            }

            var actual: [4]usize = undefined;
            var line_feeds: LineFeedSearch = .{};
            try std.testing.expect(line_feeds.findLines(4, input, 0, &actual));
            try std.testing.expectEqual(expected, actual);
        }
    }

    var incomplete_ends: [4]usize = undefined;
    var line_feeds: LineFeedSearch = .{};
    try std.testing.expect(!line_feeds.findLines(4, "a\nb\nc\n", 0, &incomplete_ends));
}

test "[property] - [reader]: saved LF masks preserve searches after skips and rewinds" {
    var storage: [193]u8 = undefined;
    inline for (.{ std.simd.suggestVectorLength(u8) orelse 1, STRUCTURAL_BLOCK_BYTES }) |lanes| {
        for (0..storage.len) |newline| {
            @memset(&storage, 'x');
            storage[newline] = '\n';
            storage[storage.len - 1] = '\n';
            var line_feeds: LineFeedSearch = .{};
            for (0..storage.len + 1) |start| {
                const expected = if (std.mem.findScalar(u8, storage[start..], '\n')) |end|
                    start + end
                else
                    null;
                try std.testing.expectEqual(expected, line_feeds.find(lanes, &storage, start));
                try std.testing.expectEqual(expected, line_feeds.find(1, &storage, start));
                try std.testing.expectEqual(@as(?usize, newline), line_feeds.find(lanes, &storage, 0));
            }
        }
        @memset(&storage, '\n');
        var line_feeds: LineFeedSearch = .{};
        for (0..lanes) |start| {
            try std.testing.expectEqual(@as(?usize, start), line_feeds.find(lanes, &storage, start));
            try std.testing.expectEqual(@as(usize, lanes), line_feeds.block_end);
        }
    }
}

test "[property] - [reader]: short reads borrow complete records without field allocations" {
    const cases = [_][]const u8{
        "@r\nAR\n+note\n!~\n@next\nC\n+\n#\n",
        "@r\r\nAR\r\n+note\r\n!~\r\n@next\r\nC\r\n+\r\n#\r\n",
        "@r\r\nAR\n+note\r\n!~\n@next\nC\r\n+\n#\r\n",
    };
    for (cases) |data| {
        for (1..data.len + 1) |chunk_limit| {
            var source = ProjectionTestSource.init(data, 0, null);
            source.chunk_limit = chunk_limit;
            var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
            var reader = try Reader.init(failing.allocator(), source.byteSource(), .{});
            defer reader.deinit();
            var validator = AdaptiveRecordValidator.init(.{});

            const first = (try nextValidatedRecord(&reader, &validator)).?;
            try std.testing.expectEqualStrings("r", first.record.header);
            try std.testing.expectEqualStrings("AR", first.record.sequence);
            try std.testing.expectEqualStrings("note", first.record.plus);
            try std.testing.expectEqualStrings("!~", first.record.quality);
            try std.testing.expect(first.semantic_error == null);
            try std.testing.expect(validator.use_full_iupac);
            if (std.mem.findScalar(u8, data, '\r') == null) {
                try std.testing.expectEqualStrings("@r\nAR\n+note\n!~\n", first.canonical_span.?);
            } else try std.testing.expect(first.canonical_span == null);
            const second = (try nextValidatedRecord(&reader, &validator)).?;
            try std.testing.expectEqualStrings("next", second.record.header);
            try std.testing.expectEqualStrings("C", second.record.sequence);
            try std.testing.expectEqualStrings("#", second.record.quality);
            try std.testing.expect((try nextValidatedRecord(&reader, &validator)) == null);
            try std.testing.expectEqual(@as(u64, data.len), reader.byteOffset());
            try std.testing.expectEqual(@as(u64, 2), reader.recordIndex());
            try std.testing.expect(!failing.has_induced_failure);
        }
    }
}

test "[property] - [reader]: refill retries preserve records, progress, and failure order" {
    const unlimited = std.math.maxInt(usize);
    for ([_][]const u8{
        "",                      "@r\nR\n+\n!\n",      "@r\r\nRY\r\n+note\r\n!~\r\n", "@r\nA\r\n+\n!\r\n",
        "@r\n\n+\n\n",           "r\nA\n+\n!\n",       "@ \nA\n+\n!\n",               "@r\nA\n-\n!\n",
        "@r\nR\n+\n \n",         "@r\n?\n-\n \n",      "@r\nAAA\n+\n!\n!\n",          "@r\nA\n+\n!\r",
        "@r\r\r\nA\n+\r\r\n!\n", "@r\nR\rA\n+\n!!!\n",
    }) |input| {
        for ([_]Options{ .{}, .{ .max_line_bytes = 2 } }) |options| {
            for (0..input.len + 1) |split| {
                try expectRefillDelivery(input, split, unlimited, null, options, .{});
                try expectRefillDelivery(input[0..split], 0, 1, null, options, .{});
                try expectRefillDelivery(input, 0, unlimited, split, options, .{});
            }
            for (1..input.len + 1) |chunk_limit| {
                try expectRefillDelivery(input, 0, chunk_limit, null, options, .{});
            }
        }
    }

    const input = "@r\nR\n+\n!\n" ** 2;
    for ([_]u64{ (1 << 32) - 3, std.math.maxInt(u64) - 8, std.math.maxInt(u64) - 1, std.math.maxInt(u64) }) |counter| {
        for (0..input.len + 1) |split| {
            try expectRefillDelivery(input, split, unlimited, null, .{}, .{ .bytes = counter });
            try expectRefillDelivery(input, split, unlimited, null, .{}, .{ .records = counter });
        }
    }
}

test "[property] - [reader]: refill retries cover the transport capacity and oversized fallback" {
    const window = io_layer.DEFAULT_READER_BUFFER_BYTES;
    const prefix = "@first\nA\n+\n!\n";
    const suffix = "@last\nC\n+\n#\n";
    for ([_]usize{ window - 1, window, window + 1, 2 * window + 11 }) |record_len| {
        const sequence_len = (record_len - 9) / 2;
        const base = try ReaderSpillFixture.init(std.testing.allocator, sequence_len, sequence_len, "\n", true);
        defer std.testing.allocator.free(base);
        const extra_header = if (base.len == record_len) "" else "x";
        const record = try std.mem.concat(std.testing.allocator, u8, &.{ base[0..4], extra_header, base[4..] });
        defer std.testing.allocator.free(record);
        try std.testing.expectEqual(record_len, record.len);
        const input = try std.mem.concat(std.testing.allocator, u8, &.{ prefix, record, suffix });
        defer std.testing.allocator.free(input);
        for ([_]usize{ window - 1, window, window + 1 }) |split| {
            try expectRefillDelivery(input, split, std.math.maxInt(usize), null, .{}, .{ .bytes = (1 << 32) - window });
            try expectRefillDelivery(input, split, std.math.maxInt(usize), null, .{ .max_line_bytes = sequence_len - 1 }, .{});
        }
        for ([_]usize{ 5 + extra_header.len, 6 + extra_header.len + sequence_len, record.len - 1 }) |changed| {
            const saved = input[prefix.len + changed];
            input[prefix.len + changed] = if (changed == record.len - 1) 'x' else '?';
            try expectRefillDelivery(input, window - 1, std.math.maxInt(usize), null, .{}, .{});
            input[prefix.len + changed] = saved;
        }
        if (record_len <= window) {
            var source = ProjectionTestSource.init(input, window - 1, null);
            source.chunk_limit = 17;
            var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
            var reader = try Reader.init(failing.allocator(), source.byteSource(), .{});
            defer reader.deinit();
            try std.testing.expect(try reader.advance());
            var span: ?[]const u8 = null;
            const parsed = (try nextWithoutId(&reader, &span)).?;
            try std.testing.expectEqual(sequence_len, parsed.sequence.len);
            try std.testing.expectEqual(sequence_len, parsed.quality.len);
            try std.testing.expectEqualStrings(record, span.?);
            try std.testing.expect(try reader.advance());
            try std.testing.expect(!try reader.advance());
            try std.testing.expect(!failing.has_induced_failure);
        }
    }
}

test "[property] - [reader]: field projections preserve delivery and failures" {
    const valid =
        "@first description\nACGT\n+annotated\n!!!!\n" ++
        "@opaque\tmetadata /1\r\nacgtn\r\n+different annotation\r\n!#~AB\r\n" ++
        "@empty\r\n\r\n+\r\n\r\n" ++
        "@tail\nN\n+\n#";
    for (0..valid.len + 1) |split| {
        try expectPayloadProjection(valid, split, null, .{}, .{});
    }

    const semantic_cases = [_][]const u8{
        "@bad-sequence\nAC?T\n+\n!!!!\n",
        "@bad-quality\nACGT\n+\n!!\x1f!\n",
        "@wider-iupac\nARYN\n+\n!!!!\n@next\nACGT\n+\n!!!!\n",
        "@embedded-cr\nAC\rGT\n+\n!!!!!\n",
        "@quality-cr\nACGT\n+\n!!\r!!\n",
        "@valid\nAC\n+annotation\n!!\n@bad-both\nA?\n+other\n!\x1f\n",
        "@valid\nAC\n+\n!!\n@bad-structure\nA?\n-\n!\x1f\n",
    };
    for (semantic_cases) |input| {
        for (0..input.len + 1) |split| {
            try expectPayloadProjection(input, split, null, .{}, .{});
        }
    }

    const acgtn_input = "@acgtn-policy\nARYN\n+\n!!!!\n";
    for (0..acgtn_input.len + 1) |split| {
        try expectPayloadProjection(
            acgtn_input,
            split,
            null,
            .{},
            .{ .alphabet = .acgtn },
        );
    }

    const exact_limit_crlf = "@r\r\nA\r\n+\r\n!\r\n";
    for (0..exact_limit_crlf.len + 1) |split| {
        try expectPayloadProjection(
            exact_limit_crlf,
            split,
            null,
            .{ .max_line_bytes = 2 },
            .{},
        );
    }

    const malformed = [_]struct {
        input: []const u8,
        split: usize,
        options: Options = .{},
    }{
        .{ .input = "@r\nAC\n+\n", .split = 5 },
        .{ .input = "r\nA\n+\n!\n", .split = 3 },
        .{ .input = "@r\nA\n-\n!\n", .split = 6 },
        .{ .input = "@r\nAA\n+\n!\n", .split = 7 },
        .{
            .input = "@long\nA\n+\n!\n",
            .split = 4,
            .options = .{ .max_line_bytes = 4 },
        },
        .{
            .input = "@r\nA\n+\n!\r",
            .split = 10,
            .options = .{ .max_line_bytes = 1 },
        },
    };
    for (malformed) |case| {
        try expectPayloadProjection(case.input, 0, null, case.options, .{});
        try expectPayloadProjection(case.input, case.split, null, case.options, .{});
    }

    try expectPayloadProjection(valid, 4, 9, .{}, .{});
}

test "[property] - [reader]: validated records match separate parsing and validation" {
    const prefix = "@first\nAC\n+note\n!~\n";
    for ([_][]const u8{
        "",                            "@r\n\n+\n\n",         "@r\nR\n+\n!\n",  "@r\nR\n+\n \n", "@r\nR?\n+\n!!\n",
        "@r\r\nRY\r\n+note\r\n!~\r\n", "@r\nAA\r\n+\n!!\r\n", "@r\nA\n+\n!",    "@r\nA\n+\n\r",  "@r\r\r\nA\n+\n!\n",
        "@r\nA\rA\n+\n!!!\n",          "r\nA\n+\n!\n",        "@ \nA\n+\n!\n",  "@r\n?\n-\n \n", "@r\nAA\n+\n!\n",
        "@r\nAAA\n+\n!\n!\n",          "@r\nR\n+\n",          "@r\nR\n+\n!!\n",
    }) |suffix| {
        var storage: [128]u8 = undefined;
        const input = try std.fmt.bufPrint(&storage, "{s}{s}{s}", .{ prefix, suffix, prefix });
        for (0..input.len + 1) |split| {
            try expectValidatedProjection(input, split, null, .{}, .{}, .{});
        }
        for (0..suffix.len + 1) |split| {
            try expectValidatedProjection(suffix, split, null, .{}, .{}, .{});
        }
        try expectValidatedProjection(suffix, 0, null, .{}, .{ .alphabet = .acgtn }, .{});
        try expectValidatedProjection(suffix, 0, null, .{ .max_line_bytes = 2 }, .{}, .{});
        for (0..suffix.len + 1) |fail_at| {
            try expectValidatedProjection(suffix, 0, fail_at, .{}, .{}, .{});
        }
    }

    const vector_len = std.simd.suggestVectorLength(u8) orelse 16;
    for ([_]usize{ 15, 16, 17, vector_len - 1, vector_len, vector_len + 1, 63, 64, 65 }) |len| {
        const input = try ReaderSpillFixture.init(std.testing.allocator, len, len, "\n", true);
        defer std.testing.allocator.free(input);
        try expectValidatedProjection(input, 0, null, .{}, .{}, .{});
        const sequence_start = std.mem.findScalar(u8, input, '\n').? + 1;
        const quality_start = std.mem.findScalarPos(u8, input, sequence_start + len + 1, '\n').? + 1;
        for ([_]usize{ sequence_start, quality_start }) |start| {
            for ([_]usize{ 0, len - 1 }) |index| {
                const saved = input[start + index];
                defer input[start + index] = saved;
                for ([_]u8{ 'R', '?', 0, '\r', '\n', 127 }) |byte| {
                    input[start + index] = byte;
                    for ([_]usize{ 0, start + index, start + index + 1 }) |split| {
                        try expectValidatedProjection(input, split, null, .{}, .{}, .{});
                    }
                }
            }
        }
    }

    for ([_]u64{ std.math.maxInt(u32) - 1, std.math.maxInt(u64) - 1, std.math.maxInt(u64) }) |counter| {
        try expectValidatedProjection(prefix ++ prefix, 0, null, .{}, .{}, .{ .records = counter });
        try expectValidatedProjection(prefix ++ prefix, 0, null, .{}, .{}, .{ .bytes = counter });
    }
}

test "[edge] - [reader]: buffered validation preserves borrowed mates and validator state" {
    const first = "@r/1\nA\n+\n!\n";
    const second = "@r/2\nR\n+\n~\n";
    for ([_]struct { split: usize, second_buffered: bool }{
        .{ .split = 1, .second_buffered = true },
        .{ .split = first.len, .second_buffered = false },
        .{ .split = first.len + 8, .second_buffered = false },
        .{ .split = first.len + second.len, .second_buffered = true },
    }) |case| {
        var source = ProjectionTestSource.init(first ++ second, case.split, null);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        var validator = AdaptiveRecordValidator.init(.{});
        const mate1 = (try nextValidatedRecord(&reader, &validator)).?;
        try std.testing.expectEqualStrings(first, mate1.canonical_span.?);
        try std.testing.expect(mate1.semantic_error == null);
        const old_pos = source.pos;
        const old_cursor = reader.cursor;
        const buffered = try nextBufferedValidatedRecord(&reader, &validator);
        try std.testing.expectEqual(old_pos, source.pos);
        try std.testing.expectEqualStrings(first, mate1.canonical_span.?);
        if (!case.second_buffered) {
            try std.testing.expect(buffered == null);
            try std.testing.expectEqual(old_cursor, reader.cursor);
            try std.testing.expectEqual(@as(u64, 1), reader.recordIndex());
            try std.testing.expectEqual(@as(u64, first.len), reader.byteOffset());
            try std.testing.expect(!validator.use_full_iupac);
        }
        const mate2 = buffered orelse (try nextValidatedRecord(&reader, &validator)).?;
        try std.testing.expectEqualStrings("r/2", mate2.record.header);
        try std.testing.expectEqualStrings("R", mate2.record.sequence);
        try std.testing.expectEqualStrings("", mate2.record.plus);
        try std.testing.expectEqualStrings("~", mate2.record.quality);
        try std.testing.expect(mate2.semantic_error == null);
        try std.testing.expect(validator.use_full_iupac);
        try std.testing.expectEqualDeep(RecordOffsets{
            .header = first.len,
            .sequence = first.len + 5,
            .plus = first.len + 7,
            .quality = first.len + 9,
        }, reader.currentRecordOffsets().?);
    }

    var source = io_layer.SliceSource.init("@r\nR\n+\n \n");
    var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    var validator = AdaptiveRecordValidator.init(.{});
    const invalid = (try nextValidatedRecord(&reader, &validator)).?;
    try expectSemanticErrorEqual(.{
        .code = .s006_invalid_quality_range,
        .message = "quality byte must be ASCII 33 through 126",
        .field = .quality,
        .byte_index = 0,
    }, invalid.semantic_error);
    try std.testing.expect(validator.use_full_iupac);
}

test "[edge] - [reader]: fallback projections retain only requested fields" {
    const field_len = io_layer.DEFAULT_READER_BUFFER_BYTES + 17;
    const input = try ReaderSpillFixture.init(
        std.testing.allocator,
        field_len,
        field_len,
        "\n",
        true,
    );
    defer std.testing.allocator.free(input);

    var advance_source = io_layer.SliceSource.init(input);
    var advance_tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    var advance_reader = try Reader.init(
        advance_tracking.allocator(),
        advance_source.byteSource(),
        .{},
    );
    defer advance_reader.deinit();
    const initial_allocations = advance_tracking.allocations;

    try std.testing.expect(try advance_reader.advance());
    try std.testing.expectEqual(initial_allocations, advance_tracking.allocations);
    try std.testing.expectEqual(@as(u64, 1), advance_reader.recordIndex());
    for (advance_reader.fallback_fields) |field| {
        try std.testing.expectEqual(@as(usize, 0), field.len);
    }

    var payload_source = io_layer.SliceSource.init(input);
    var payload_reader = try Reader.init(
        std.testing.allocator,
        payload_source.byteSource(),
        .{},
    );
    defer payload_reader.deinit();

    const payload = (try payload_reader.nextPayload()).?;
    try std.testing.expectEqual(field_len, payload.sequence.len);
    try std.testing.expectEqual(field_len, payload.quality.len);
    try std.testing.expectEqual(@as(usize, 0), payload_reader.fallback_fields[0].len);
    try std.testing.expectEqual(@as(usize, 0), payload_reader.fallback_fields[2].len);
    try std.testing.expectEqual(@as(usize, 0), payload_reader.fallback_fields[0].storage.len);
    try std.testing.expectEqual(@as(usize, 0), payload_reader.fallback_fields[2].storage.len);

    var validated_source = io_layer.SliceSource.init(input);
    var validated_reader = try Reader.init(
        std.testing.allocator,
        validated_source.byteSource(),
        .{},
    );
    defer validated_reader.deinit();
    var validator = AdaptiveRecordValidator.init(.{});

    const validated = (try validated_reader.nextValidatedHeader(&validator)).?;
    try std.testing.expectEqualStrings("abc", validated.header);
    try std.testing.expect(validated.semantic_error == null);
    try std.testing.expectEqual(@as(usize, 0), validated_reader.fallback_fields[1].storage.len);
    try std.testing.expectEqual(@as(usize, 0), validated_reader.fallback_fields[1].len);
    try std.testing.expectEqual(@as(usize, 0), validated_reader.fallback_fields[2].storage.len);
    try std.testing.expectEqual(@as(usize, 0), validated_reader.fallback_fields[3].storage.len);
}

test "[edge] - [reader]: counts and locations cross 32 bits through buffered and fallback reads" {
    const crossing: u64 = 1 << 32;
    const data = "@r\nA\n+\n!\n" ** 3 ++ "@bad\nA\n?\n!\n";
    for ([_]usize{ 1, 28, data.len }) |split| {
        for ([_]bool{ false, true }) |advance| {
            errdefer std.debug.print("reader 32-bit crossing: split {d}, advance {}\n", .{
                split, advance,
            });
            var source = ProjectionTestSource.init(data, split, null);
            var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
            defer reader.deinit();
            reader.byte_offset = crossing - 3;
            reader.record_index = crossing - 1;

            for (0..3) |index| {
                const start = crossing - 3 + 9 * index;
                if (advance) {
                    try std.testing.expect(try reader.advance());
                } else {
                    try std.testing.expectEqualDeep(Record{
                        .header = "r",
                        .id = "r",
                        .sequence = "A",
                        .plus = "",
                        .quality = "!",
                    }, (try reader.next()).?);
                    try std.testing.expectEqualDeep(RecordOffsets{
                        .header = start,
                        .sequence = start + 3,
                        .plus = start + 5,
                        .quality = start + 7,
                    }, reader.currentRecordOffsets().?);
                }
                try std.testing.expectEqual(crossing + index, reader.recordIndex());
                try std.testing.expectEqual(start + 9, reader.byteOffset());
            }
            if (advance) {
                try std.testing.expectError(error.S001InvalidPlusLine, reader.advance());
            } else {
                try std.testing.expectError(error.S001InvalidPlusLine, reader.next());
            }
            try std.testing.expectEqual(crossing + 2, reader.recordIndex());
            try std.testing.expectEqual(crossing + 33, reader.byteOffset());
            try std.testing.expectEqualDeep(ParseError{
                .code = .s001_invalid_plus_line,
                .message = "plus line must start with '+'",
                .record_index = crossing + 2,
                .byte_offset = crossing + 31,
                .line_in_record = 3,
            }, reader.takeLastError().?);
            try std.testing.expect(reader.takeLastError() == null);
        }
    }
}

test "[edge] - [reader]: buffered progress rejects maximum offset and record count" {
    const data = "@r\nA\n+\n!\n";

    var offset_source = io_layer.SliceSource.init(data);
    var offset_reader = try Reader.init(
        std.testing.allocator,
        offset_source.byteSource(),
        .{},
    );
    defer offset_reader.deinit();
    offset_reader.byte_offset = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, offset_reader.next());
    try std.testing.expectEqual(std.math.maxInt(u64), offset_reader.byte_offset);
    try std.testing.expectEqual(@as(u64, 0), offset_reader.record_index);

    var count_source = io_layer.SliceSource.init(data);
    var count_reader = try Reader.init(
        std.testing.allocator,
        count_source.byteSource(),
        .{},
    );
    defer count_reader.deinit();
    count_reader.record_index = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, count_reader.next());
    try std.testing.expectEqual(std.math.maxInt(u64), count_reader.record_index);
    try std.testing.expectEqual(@as(u64, data.len), count_reader.byte_offset);
}

test "[edge] - [reader]: fallback progress rejects maximum offset and record count" {
    const data = "@r\nA\n+\n!\n";

    var offset_source = ProjectionTestSource.init(data, 1, null);
    var offset_reader = try Reader.init(
        std.testing.allocator,
        offset_source.byteSource(),
        .{},
    );
    defer offset_reader.deinit();
    offset_reader.byte_offset = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, offset_reader.next());
    try std.testing.expectEqual(std.math.maxInt(u64), offset_reader.byte_offset);
    try std.testing.expectEqual(@as(u64, 0), offset_reader.record_index);

    var count_source = ProjectionTestSource.init(data, 1, null);
    var count_reader = try Reader.init(
        std.testing.allocator,
        count_source.byteSource(),
        .{},
    );
    defer count_reader.deinit();
    count_reader.record_index = std.math.maxInt(u64);

    try std.testing.expectError(error.ArithmeticLimit, count_reader.next());
    try std.testing.expectEqual(std.math.maxInt(u64), count_reader.record_index);
    try std.testing.expectEqual(@as(u64, data.len), count_reader.byte_offset);
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

    for ([_][]const u8{ input, "@r\nA\n+\n!" }) |bytes| {
        for ([_]usize{ std.math.maxInt(usize), 37 }) |chunk_limit| {
            try std.testing.checkAllAllocationFailures(
                std.testing.allocator,
                readSpillForAllocationCheck,
                .{ bytes, chunk_limit },
            );
        }
    }
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

        var repeated: std.ArrayList(u8) = .empty;
        defer repeated.deinit(std.testing.allocator);
        for (0..4) |_| try repeated.appendSlice(std.testing.allocator, input);
        var source = io_layer.SliceSource.init(repeated.items);
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
            const expected_capacity = if (crlf or field_index == 3)
                storage_limit
            else
                max_line_bytes;
            try std.testing.expectEqual(expected_capacity, field.storage.len);
        }
        try std.testing.expectError(
            ReaderError.LineTooLong,
            reader.ensureFieldCapacity(0, storage_limit + 1, false),
        );
        var warm_allocated = tracking.allocated_bytes;
        var count: usize = 1;
        while (try reader.next()) |next| {
            try std.testing.expectEqual(max_line_bytes - 1, next.header.len);
            try std.testing.expectEqual(max_line_bytes - 1, next.plus.len);
            try std.testing.expectEqual(max_line_bytes, next.sequence.len);
            try std.testing.expectEqual(max_line_bytes, next.quality.len);
            try std.testing.expect(std.mem.allEqual(u8, next.sequence, 'A'));
            try std.testing.expect(std.mem.allEqual(u8, next.quality, 'I'));
            for (reader.fallback_fields, 0..) |field, field_index| {
                try std.testing.expect(field.storage.len <= if (field_index == 1)
                    @max(io_layer.DEFAULT_READER_BUFFER_BYTES, storage_limit)
                else
                    storage_limit);
            }
            try std.testing.expect(tracking.allocated_bytes - tracking.freed_bytes <=
                2 * io_layer.DEFAULT_READER_BUFFER_BYTES + 3 * storage_limit);
            // The second record first attempts whole-record refill assembly.
            try std.testing.expectEqual(io_layer.DEFAULT_READER_BUFFER_BYTES, reader.fallback_fields[1].storage.len);
            if (crlf and count > 1) try std.testing.expectEqual(warm_allocated, tracking.allocated_bytes);
            warm_allocated = tracking.allocated_bytes;
            count += 1;
        }
        try std.testing.expectEqual(@as(usize, 4), count);
    }
}

test "[property] - [reader]: repeated initialization releases both source modes at boundary line limits" {
    const Exercise = struct {
        fn run(allocator: std.mem.Allocator, limit: usize, borrowed: bool) !void {
            var empty = std.Io.Reader.fixed("");
            var gzip = io_layer.GzipSource.init(&empty);
            var source = io_layer.SliceSource.init("");
            var reader = if (borrowed)
                try initBorrowedGzipReader(allocator, &gzip, .{ .max_line_bytes = limit })
            else
                try Reader.init(allocator, source.byteSource(), .{ .max_line_bytes = limit });
            defer reader.deinit();

            for (reader.fallback_fields) |field| {
                try std.testing.expectEqual(@as(usize, 0), field.storage.len);
            }
            try std.testing.expectEqual(
                @as(usize, if (borrowed) 0 else io_layer.DEFAULT_READER_BUFFER_BYTES),
                reader.transport_storage.len,
            );
        }
    };
    for ([_]usize{ 0, 8192, 80 * 1024, io_layer.DEFAULT_READER_BUFFER_BYTES - 1, io_layer.DEFAULT_READER_BUFFER_BYTES, std.math.maxInt(usize) }) |limit| {
        for ([_]bool{ false, true }) |borrowed| {
            if (borrowed and @sizeOf(@FieldType(io_layer.GzipSource, "decompressor_buffer")) == 0) continue;
            try std.testing.checkAllAllocationFailures(std.testing.allocator, Exercise.run, .{ limit, borrowed });
            var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
            for (0..3) |iteration| {
                try Exercise.run(tracking.allocator(), limit, borrowed);
                try std.testing.expectEqual(@as(usize, if (borrowed) 0 else iteration + 1), tracking.allocations);
                try std.testing.expectEqual(
                    @as(usize, if (borrowed) 0 else (iteration + 1) * io_layer.DEFAULT_READER_BUFFER_BYTES),
                    tracking.allocated_bytes,
                );
                try std.testing.expectEqual(tracking.allocated_bytes, tracking.freed_bytes);
            }
        }
    }
}

test "[edge] - [reader]: field storage limit saturates after checked overflow" {
    const max = std.math.maxInt(usize);

    try std.testing.expectEqual(1024 * 1024 + 1, fieldStorageLimit(1024 * 1024));
    try std.testing.expectEqual(max, fieldStorageLimit(max - 1));
    try std.testing.expectEqual(max, fieldStorageLimit(max));
}

test "[property] - [reader]: borrowed gzip allocates only for retained fallback fields" {
    // The separate parser test executable always uses native gzip.
    if (@sizeOf(@FieldType(io_layer.GzipSource, "decompressor_buffer")) == 0) return;
    const Exercise = struct {
        const bytes = "@r\r\nAR\r\n+left\r\n!~\r\n@next\nC\n+\n#\n";

        fn run(allocator: std.mem.Allocator, split: usize) !void {
            var compressed: [128]u8 = undefined;
            var output = std.Io.Writer.fixed(&compressed);
            for ([_][]const u8{ bytes[0..split], bytes[split..] }) |part| {
                const len: u16 = @intCast(part.len);
                try output.writeAll(&.{ 0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255, 1 });
                try output.writeInt(u16, len, .little);
                try output.writeInt(u16, ~len, .little);
                try output.writeAll(part);
                try output.writeInt(u32, std.hash.Crc32.hash(part), .little);
                try output.writeInt(u32, len, .little);
            }
            var input = std.Io.Reader.fixed(output.buffered());
            var gzip = io_layer.GzipSource.init(&input);
            var reader = try initBorrowedGzipReader(allocator, &gzip, .{ .max_line_bytes = 8 });
            defer reader.deinit();
            const first = (try reader.next()).?;
            try std.testing.expectEqualStrings("r", first.header);
            try std.testing.expectEqualStrings("r", first.id);
            try std.testing.expectEqualStrings("AR", first.sequence);
            try std.testing.expectEqualStrings("left", first.plus);
            try std.testing.expectEqualStrings("!~", first.quality);
            const second = (try reader.next()).?;
            try std.testing.expectEqualStrings("next", second.header);
            try std.testing.expectEqualStrings("C", second.sequence);
            try std.testing.expectEqualStrings("#", second.quality);
            try std.testing.expect((try reader.next()) == null);
            try std.testing.expectEqual(@as(u64, bytes.len), reader.byteOffset());
            for (reader.fallback_fields) |field| try std.testing.expect(field.storage.len <= 9);
        }
    };
    for ([_]usize{ 0, 1, 5, 14, Exercise.bytes.len - 1 }) |split| {
        try std.testing.checkAllAllocationFailures(std.testing.allocator, Exercise.run, .{split});
        var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        try Exercise.run(tracking.allocator(), split);
        try std.testing.expectEqual(tracking.allocated_bytes, tracking.freed_bytes);
        if (split == 0) {
            try std.testing.expectEqual(@as(usize, 0), tracking.allocations);
        } else try std.testing.expect(tracking.allocations > 0);
    }
}

test "[unit] - [reader]: field storage exchanges preserve both live prefixes" {
    const window = io_layer.DEFAULT_READER_BUFFER_BYTES;
    var source = io_layer.SliceSource.init("");
    var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    var reader = try Reader.init(tracking.allocator(), source.byteSource(), .{});
    defer reader.deinit();
    reader.machine.expected = .sequence;
    try reader.ensureFieldCapacity(1, window, false);
    reader.fallback_fields[1].len = window;
    @memset(reader.fallback_fields[1].storage, 'A');
    reader.fallback_fields[2].storage = try tracking.allocator().alloc(u8, 2 * window);
    const donor_pointer = reader.fallback_fields[2].storage.ptr;
    const old_pointer = reader.fallback_fields[1].storage.ptr;
    const allocations = tracking.allocations;

    try reader.ensureFieldCapacity(1, window + 1, false);

    try std.testing.expectEqual(allocations, tracking.allocations);
    try std.testing.expectEqual(donor_pointer, reader.fallback_fields[1].storage.ptr);
    try std.testing.expectEqual(old_pointer, reader.fallback_fields[2].storage.ptr);
    try std.testing.expectEqual(window, reader.fallback_fields[1].len);
    try std.testing.expect(std.mem.allEqual(u8, reader.fallback_fields[1].storage[0..window], 'A'));
    try std.testing.expectEqual(@as(usize, 0), reader.fallback_fields[2].len);

    reader.fallback_fields[3].storage = try tracking.allocator().alloc(u8, 4 * window);
    reader.fallback_fields[3].len = 1;
    reader.fallback_fields[3].storage[0] = '!';
    const live_pointer = reader.fallback_fields[3].storage.ptr;
    try reader.ensureFieldCapacity(1, 2 * window + 1, false);
    try std.testing.expectEqual(live_pointer, reader.fallback_fields[1].storage.ptr);
    try std.testing.expectEqual(@as(u8, '!'), reader.fallback_fields[3].storage[0]);
    try std.testing.expect(std.mem.allEqual(u8, reader.fallback_fields[1].storage[0..window], 'A'));
}

test "[failure] - [reader]: failed spill growth preserves owned storage" {
    var source = io_layer.SliceSource.init("");
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 1,
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

test "[property] - [reader]: storage exchange preserves unequal prefixes" {
    const window = io_layer.DEFAULT_READER_BUFFER_BYTES;
    for ([_]usize{ 0, 1, window }) |live| {
        for ([_]usize{ 0, 1, window }) |other_live| {
            var source = io_layer.SliceSource.init("");
            var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
            defer reader.deinit();
            try reader.ensureFieldCapacity(1, window, false);
            const field = &reader.fallback_fields[1];
            const donor = &reader.fallback_fields[0];
            field.len = live;
            @memset(field.storage[0..live], 'A');
            donor.storage = try std.testing.allocator.alloc(u8, 2 * window);
            donor.len = other_live;
            @memset(donor.storage[0..other_live], 'h');
            const pointer = donor.storage.ptr;

            try reader.ensureFieldCapacity(1, window + 1, false);

            try std.testing.expectEqual(pointer, field.storage.ptr);
            try std.testing.expectEqual(live, field.len);
            try std.testing.expectEqual(other_live, donor.len);
            try std.testing.expect(std.mem.allEqual(u8, field.storage[0..live], 'A'));
            try std.testing.expect(std.mem.allEqual(u8, donor.storage[0..other_live], 'h'));
        }
    }
}

test "[failure] - [reader]: heterogeneous storage reuse cleans up every allocation failure" {
    const large = io_layer.DEFAULT_READER_BUFFER_BYTES * 2 + 17;
    const geometry = [_][4]usize{
        .{ large, 151, large, 151 },
        .{ large, 151, large, 151 },
        .{ 8, large, 8, large },
        .{ 8, large, 8, large },
        .{ large, large, large, large },
        .{ large, large, large, large },
    };
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(std.testing.allocator);
    for (geometry, 0..) |lengths, index| {
        for (lengths, "hApI", 0..) |length, byte, field| {
            if (field == 0 or field == 2) {
                try input.append(std.testing.allocator, if (field == 0) '@' else '+');
            }
            try input.appendNTimes(std.testing.allocator, byte, length);
            if (index % 2 == 1) try input.append(std.testing.allocator, '\r');
            try input.append(std.testing.allocator, '\n');
        }
    }
    const Exercise = struct {
        fn check(record: Record, lengths: [4]usize) !void {
            for ([_][]const u8{ record.header, record.sequence, record.plus, record.quality }, lengths, "hApI") |bytes, length, fill| {
                try std.testing.expectEqual(length, bytes.len);
                try std.testing.expect(std.mem.allEqual(u8, bytes, fill));
            }
        }

        fn run(allocator: std.mem.Allocator, bytes: []const u8, retain: bool) !void {
            var source = io_layer.SliceSource.init(bytes);
            var reader = try Reader.init(allocator, source.byteSource(), .{});
            defer reader.deinit();
            var retained: RetainedRecordStorage = .{};
            defer retained.deinit(allocator);
            var index: usize = 0;
            while (try reader.next()) |first| {
                try check(first, geometry[index]);
                index += 1;
                if (retain) {
                    try std.testing.expect(retainFallbackRecordStorage(&reader, &retained, first));
                    const second = (try reader.next()).?;
                    try check(first, geometry[index - 1]);
                    try check(second, geometry[index]);
                    restoreFallbackRecordStorage(&reader, &retained);
                    try check(second, geometry[index]);
                    index += 1;
                }
            }
            try std.testing.expectEqual(geometry.len, index);
        }
    };
    for ([_]bool{ false, true }) |retain| {
        try std.testing.checkAllAllocationFailures(std.testing.allocator, Exercise.run, .{ input.items, retain });
    }
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
        sequence_len,
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
            try std.testing.expect(!recordHasTerminalCr(record));
            try writeRecordFields(&writer, record);
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
        try std.testing.expectEqualStrings(record2_prefix ++ record2_suffix, canonical_span2.?);
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

test "[integration] - [record delivery]: retained fallback storage survives a refill" {
    for ([_]usize{ 2, io_layer.DEFAULT_READER_BUFFER_BYTES }) |mate2_len| {
        const field_len = io_layer.DEFAULT_READER_BUFFER_BYTES - 7;
        var input: std.ArrayList(u8) = .empty;
        defer input.deinit(std.testing.allocator);
        try input.appendSlice(std.testing.allocator, "@pair/1\n");
        try input.appendNTimes(std.testing.allocator, 'A', field_len);
        try input.appendSlice(std.testing.allocator, "\n+\n");
        try input.appendNTimes(std.testing.allocator, '!', field_len);
        try input.appendSlice(std.testing.allocator, "\n@pair/2\n");
        try input.appendNTimes(std.testing.allocator, 'T', mate2_len);
        try input.appendSlice(std.testing.allocator, "\n+right\n");
        try input.appendNTimes(std.testing.allocator, '#', mate2_len);
        try input.append(std.testing.allocator, '\n');

        var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
        {
            var source = io_layer.SliceSource.init(input.items);
            var reader = try Reader.init(tracking.allocator(), source.byteSource(), .{});
            defer reader.deinit();
            var retained: RetainedRecordStorage = .{};
            defer retained.deinit(tracking.allocator());

            var canonical_span1: ?[]const u8 = null;
            const record1 = (try nextWithoutId(&reader, &canonical_span1)).?;
            try std.testing.expect(canonical_span1 == null);
            try std.testing.expectEqual(field_len, record1.sequence.len);
            try std.testing.expectEqual(field_len, record1.quality.len);

            var canonical_span2: ?[]const u8 = null;
            try std.testing.expect(
                (try nextBufferedWithoutId(&reader, &canonical_span2)) == null,
            );
            const sequence_ptr = record1.sequence.ptr;
            try std.testing.expect(retainFallbackRecordStorage(&reader, &retained, record1));
            const allocations_before_refill = tracking.allocations;

            var validator = AdaptiveRecordValidator.init(.{});
            const buffered2 = try nextBufferedAfterFallbackTransfer(
                &reader,
                &validator,
            );
            const validated2 = buffered2 orelse (try nextFallbackValidatedRecord(&reader, &validator)).?;
            const record2 = validated2.record;
            try std.testing.expect(validated2.semantic_error == null);
            if (mate2_len == 2) try std.testing.expectEqual(allocations_before_refill, tracking.allocations);
            try std.testing.expectEqualStrings("pair/1", record1.header);
            try std.testing.expect(std.mem.allEqual(u8, record1.sequence, 'A'));
            try std.testing.expect(std.mem.allEqual(u8, record1.quality, '!'));
            try std.testing.expectEqualStrings("pair/2", record2.header);
            try std.testing.expectEqual(mate2_len, record2.sequence.len);
            try std.testing.expect(std.mem.allEqual(u8, record2.sequence, 'T'));
            try std.testing.expectEqualStrings("right", record2.plus);
            try std.testing.expect(std.mem.allEqual(u8, record2.quality, '#'));

            restoreFallbackRecordStorage(&reader, &retained);
            try std.testing.expectEqual(sequence_ptr, reader.fallback_fields[1].storage.ptr);
            try std.testing.expectEqualStrings("pair/2", record2.header);
            try std.testing.expect(std.mem.allEqual(u8, record2.sequence, 'T'));
            try std.testing.expect(std.mem.allEqual(u8, record2.quality, '#'));
            try std.testing.expect((try nextWithoutId(&reader, &canonical_span2)) == null);
        }
        try std.testing.expectEqual(tracking.allocated_bytes, tracking.freed_bytes);
    }
}

test "[integration] - [writer]: prechecked Reader records match checked serialization" {
    const cases = [_]struct { input: []const u8, expected: ?[]const u8 }{
        .{
            .input = "@r\rdesc\r\nACGT\r\n+note\rx\r\n!#$%\r\n",
            .expected = "@r\rdesc\nACGT\n+note\rx\n!#$%\n",
        },
        .{ .input = "@r\r\r\nA\n+\n!\n", .expected = null },
        .{ .input = "@r\nA\n+\r\r\n!\n", .expected = null },
        .{ .input = "@\r\r\nA\n+\n!\n", .expected = null },
    };
    for (cases) |case| {
        var source = io_layer.SliceSource.init(case.input);
        var reader = try Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        const record = (try reader.next()).?;
        try std.testing.expect(validateRecord(record, .{}) == null);

        var checked_bytes: [64]u8 = undefined;
        var checked_sink = io_layer.SliceSink.init(&checked_bytes);
        var checked_writer = Writer.init(checked_sink.byteSink());

        var trusted_bytes: [64]u8 = undefined;
        var trusted_sink = io_layer.SliceSink.init(&trusted_bytes);
        var trusted_writer = Writer.init(trusted_sink.byteSink());
        if (case.expected) |expected| {
            try checked_writer.writeRecord(record);
            try std.testing.expect(!recordHasTerminalCr(record));
            try writeRecordFields(&trusted_writer, record);
            try std.testing.expectEqualStrings(expected, checked_sink.written());
        } else {
            try std.testing.expectError(error.InvalidRecord, checked_writer.writeRecord(record));
            try std.testing.expect(recordHasTerminalCr(record));
            try std.testing.expectEqualStrings("", trusted_sink.written());
        }
        try std.testing.expectEqualStrings(checked_sink.written(), trusted_sink.written());
    }
}

test "[unit] - [check scanner]: state remains fixed-size" {
    try std.testing.expect(@sizeOf(CheckScanner) <= 160);
}

test "[property] - [LF mask]: preserves newline positions across block widths" {
    inline for (.{ 1, 16, STRUCTURAL_BLOCK_BYTES, std.simd.suggestVectorLength(u8) orelse 1 }) |lanes| {
        const Mask = @Int(.unsigned, lanes);
        var block: [lanes]u8 = @splat('A');
        try std.testing.expectEqual(@as(Mask, 0), newlineMask(lanes, &block));
        for (0..lanes) |lane| {
            block[lane] = '\n';
            try std.testing.expectEqual(@as(Mask, 1) << @intCast(lane), newlineMask(lanes, &block));
            block[lane] = 'A';
        }

        var expected: Mask = 0;
        for ([_]usize{ 0, 1, 15, 16, 31, 32, 47, 48, 62, 63 }) |lane| {
            if (lane >= lanes) continue;
            block[lane] = '\n';
            expected |= @as(Mask, 1) << @intCast(lane);
        }
        try std.testing.expectEqual(expected, newlineMask(lanes, &block));
    }
}

fn firstValidSequenceLineEnd(bytes: []const u8, alphabet: Alphabet) ?usize {
    return switch (alphabet) {
        .iupac => firstValidSequenceLineEndFor(.iupac, bytes),
        .acgtn => firstValidSequenceLineEndFor(.acgtn, bytes),
    };
}

test "[property] - [check scanner]: fused sequence scan preserves delimiter boundaries" {
    var bytes: [2 * STRUCTURAL_BLOCK_BYTES + 2]u8 = @splat('A');
    for (0..bytes.len) |line_end| {
        @memset(&bytes, 'A');
        bytes[line_end] = '\n';
        try std.testing.expectEqual(
            line_end,
            firstValidSequenceLineEnd(&bytes, .iupac).?,
        );
        try std.testing.expectEqual(
            line_end,
            firstValidSequenceLineEnd(&bytes, .acgtn).?,
        );

        if (line_end != 0) {
            bytes[line_end - 1] = '.';
            try std.testing.expect(firstValidSequenceLineEnd(&bytes, .iupac) == null);
            bytes[line_end - 1] = 'A';
        }
        if (line_end + 1 < bytes.len) {
            bytes[line_end + 1] = '.';
            try std.testing.expectEqual(
                line_end,
                firstValidSequenceLineEnd(&bytes, .iupac).?,
            );
        }
    }

    @memset(&bytes, 'A');
    bytes[STRUCTURAL_BLOCK_BYTES] = 'R';
    bytes[STRUCTURAL_BLOCK_BYTES + 1] = '\n';
    try std.testing.expectEqual(
        STRUCTURAL_BLOCK_BYTES + 1,
        firstValidSequenceLineEnd(&bytes, .iupac).?,
    );
    try std.testing.expect(firstValidSequenceLineEnd(&bytes, .acgtn) == null);

    for (0..bytes.len - 1) |cr_index| {
        @memset(&bytes, 'A');
        bytes[cr_index] = '\r';
        bytes[cr_index + 1] = '\n';
        for ([_]Alphabet{ .iupac, .acgtn }) |alphabet| {
            try std.testing.expectEqual(cr_index + 1, firstValidSequenceLineEnd(&bytes, alphabet).?);
            try std.testing.expect(firstValidSequenceLineEnd(bytes[0 .. cr_index + 1], alphabet) == null);
        }
        var use_full_iupac = false;
        try std.testing.expectEqual(
            cr_index + 1,
            firstValidCheckSequenceLineEnd(&bytes, .iupac, &use_full_iupac).?,
        );
        try std.testing.expect(!use_full_iupac);
        try std.testing.expectEqual(
            cr_index,
            firstInvalidCheckSequence(bytes[0 .. cr_index + 2], .iupac, &use_full_iupac).?,
        );
        try std.testing.expect(!use_full_iupac);
        if (cr_index != 0) {
            bytes[0] = 'R';
            try std.testing.expectEqual(
                cr_index + 1,
                firstValidCheckSequenceLineEnd(&bytes, .iupac, &use_full_iupac).?,
            );
            try std.testing.expect(use_full_iupac);
            bytes[0] = '.';
            try std.testing.expect(firstValidSequenceLineEnd(&bytes, .iupac) == null);
        }
    }
}

test "[property] - [check scanner]: adaptive IUPAC validation preserves byte policy" {
    for (0..256) |value| {
        const sequence = [_]u8{@intCast(value)};
        var use_full_iupac = false;
        try std.testing.expectEqual(
            firstInvalidSequenceScalar(&sequence, .iupac, 0),
            firstInvalidCheckSequence(&sequence, .iupac, &use_full_iupac),
        );
        try std.testing.expectEqual(
            alphabetAccepts(.iupac, sequence[0]) and !alphabetAccepts(.acgtn, sequence[0]),
            use_full_iupac,
        );

        use_full_iupac = false;
        try std.testing.expectEqual(
            firstInvalidSequenceScalar(&sequence, .acgtn, 0),
            firstInvalidCheckSequence(&sequence, .acgtn, &use_full_iupac),
        );
        try std.testing.expect(!use_full_iupac);
    }

    var block: [STRUCTURAL_BLOCK_BYTES]u8 = @splat('A');
    block[1] = '\n';
    block[2] = 'R';
    var use_full_iupac = false;
    try std.testing.expectEqual(
        @as(usize, 1),
        firstValidCheckSequenceLineEnd(&block, .iupac, &use_full_iupac).?,
    );
    try std.testing.expect(!use_full_iupac);

    block[0] = 'R';
    use_full_iupac = false;
    try std.testing.expectEqual(
        @as(usize, 1),
        firstValidCheckSequenceLineEnd(&block, .iupac, &use_full_iupac).?,
    );
    try std.testing.expect(use_full_iupac);

    var sequence: [2 * STRUCTURAL_BLOCK_BYTES]u8 = @splat('A');
    sequence[31] = 'R';
    sequence[70] = '.';
    use_full_iupac = false;
    try std.testing.expectEqual(
        @as(usize, 70),
        firstInvalidCheckSequence(&sequence, .iupac, &use_full_iupac).?,
    );
    try std.testing.expect(!use_full_iupac);

    sequence[5] = '.';
    try std.testing.expectEqual(
        @as(usize, 5),
        firstInvalidCheckSequence(&sequence, .iupac, &use_full_iupac).?,
    );
    sequence[5] = 'A';
    sequence[70] = 'A';
    try std.testing.expect(firstInvalidCheckSequence(&sequence, .iupac, &use_full_iupac) == null);
    try std.testing.expect(use_full_iupac);
}

test "[unit] - [check scanner]: complete record path commits only proved records" {
    const record1 = "@r one\nACGTN\n+note\n!!!!!\n";
    const record2 = "@s\n\n+\n\n";
    const input = record1 ++ record2;
    var scanner = CheckScanner.init(.{}, .{});
    var line_feeds: LineFeedSearch = .{};

    try std.testing.expectEqual(record1.len, scanner.consumeCompleteRecord(input, 0, &line_feeds, false).?);
    try std.testing.expectEqual(@as(u64, 1), scanner.record_index);
    try std.testing.expectEqual(@as(u64, record1.len), scanner.byte_offset);
    try std.testing.expectEqual(
        record2.len,
        scanner.consumeCompleteRecord(input, record1.len, &line_feeds, false).?,
    );
    try std.testing.expectEqual(@as(u64, 2), scanner.record_index);
    try std.testing.expectEqual(@as(u64, input.len), scanner.byte_offset);

    for ([_][]const u8{
        "@r\r\nA\r\n+\r\n!\r",
        "@r\nA\n+\n",
        "@r\nR\n+\n",
        "r\nA\n+\n!\n",
        "@r\n.\n+\n!\n",
        "@r\nA\n+\n\x7f\n",
        "@r\nAA\n+\n!\n",
    }) |data| {
        var fallback = CheckScanner.init(.{}, .{});
        const before = fallback;
        line_feeds = .{};
        try std.testing.expect(fallback.consumeCompleteRecord(data, 0, &line_feeds, false) == null);
        try std.testing.expectEqualDeep(before, fallback);
    }

    var limited = CheckScanner.init(.{ .max_line_bytes = 1 }, .{});
    const before = limited;
    line_feeds = .{};
    try std.testing.expect(limited.consumeCompleteRecord("@r\nA\n+\n!\n", 0, &line_feeds, false) == null);
    try std.testing.expectEqualDeep(before, limited);

    for ([_][]const u8{ "@r\nR\n+\n!\n", "@r\r\nR\r\n+\r\n!\r\n" }) |data| {
        var wide = CheckScanner.init(.{}, .{});
        line_feeds = .{};
        try std.testing.expectEqual(data.len, wide.consumeCompleteRecord(data, 0, &line_feeds, false).?);
        try std.testing.expect(wide.use_full_iupac);
    }
}

test "[unit] - [check scanner]: adaptive IUPAC state survives a chunk seam" {
    var scanner = CheckScanner.init(.{}, .{});
    _ = try scanner.feed("@r\nR");
    try std.testing.expect(scanner.use_full_iupac);
    _ = try scanner.feed("\n+\n!\n");
    try scanner.finishEof();
    try std.testing.expectEqual(@as(u64, 1), scanner.record_index);
}

test "[property] - [check scanner]: complete records accept mixed LF and CRLF endings" {
    for (0..16) |ending_mask| {
        var endings: [4][]const u8 = undefined;
        for (&endings, 0..) |*ending, index| {
            ending.* = if (ending_mask & (@as(usize, 1) << @intCast(index)) == 0)
                "\n"
            else
                "\r\n";
        }
        var storage: [32]u8 = undefined;
        const data = try std.fmt.bufPrint(&storage, "@r{s}AC{s}+{s}!~{s}", .{
            endings[0], endings[1], endings[2], endings[3],
        });
        const sequence_start = 2 + endings[0].len;
        const plus_start = sequence_start + 2 + endings[1].len;
        const quality_start = plus_start + 1 + endings[2].len;
        var line_feeds: LineFeedSearch = .{};
        const complete = scanCompleteRecord(data, 0, &line_feeds, 2, .iupac, false);
        try std.testing.expect(complete != null);
        try std.testing.expectEqualDeep([4]usize{
            sequence_start - 1, plus_start - 1, quality_start - 1, data.len - 1,
        }, complete.?.line_ends);
        try std.testing.expect(!complete.?.use_full_iupac);

        var scanner = CheckScanner.init(.{ .max_line_bytes = 2 }, .{});
        try std.testing.expectEqual(data.len, scanner.consumeCompleteRecord(data, 0, &line_feeds, false).?);
        try std.testing.expectEqual(@as(u64, 1), scanner.record_index);
        try std.testing.expectEqual(@as(u64, data.len), scanner.byte_offset);
        try std.testing.expect(scanner.atRecordBoundary());
        for (0..data.len) |cut| {
            line_feeds = .{};
            try std.testing.expect(scanCompleteRecord(data[0..cut], 0, &line_feeds, 2, .iupac, false) == null);
        }
        for (1..data.len + 2) |chunk_len| {
            try expectCheckOutcome(
                .{ .valid = 1 },
                directCheckOutcome(data, chunk_len, .{ .max_line_bytes = 2 }, .{}),
            );
        }
        for (0..data.len + 1) |split| {
            try expectValidatedProjection(data, split, null, .{ .max_line_bytes = 2 }, .{}, .{});
        }
    }
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
        .{ .data = "@r\r\n\r\n+\r\n\r\n", .expected = .{ .valid = 1 } },
        .{ .data = "@r\rb\r\nA\r\n+note\rb\r\n!\r\n", .expected = .{ .valid = 1 } },
        .{ .data = "@\r\r\nA\r\n+\r\r\n!\r\n", .expected = .{ .valid = 1 } },
        .{ .data = "@r\nA\n+\n!", .expected = .{ .valid = 1 } },
        .{ .data = "@r\nR\n+\n!\n", .expected = .{ .valid = 1 } },
        .{
            .data = "@r\r\nURYSWKMBDHVuryswkmbdhv\r\n+\r\n!!!!!!!!!!!!!!!!!!!!!!\r\n",
            .expected = .{ .valid = 1 },
        },
        .{
            .data = "@r\n.R\n+\n!!\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                0,
                3,
                2,
            ) },
        },
        .{
            .data = "@r\r\nA\rA\r\n+\r\n!!!\r\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                0,
                5,
                2,
            ) },
        },
        .{
            .data = "@r\r\nA\r\r\n+\r\n!!\r\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                0,
                5,
                2,
            ) },
        },
        .{
            .data = "@r\r\nAAA\r\n+\r\n!\r!\r\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s006_invalid_quality_range,
                "quality byte must be ASCII 33 through 126",
                0,
                13,
                4,
            ) },
        },
        .{
            .data = "@r\r\nAAA\r\n+\r\n!\n!\r\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s005_length_mismatch,
                "sequence and quality lengths differ",
                0,
                12,
                4,
            ) },
        },
        .{
            .data = "@\r\nA\r\n+\r\n!\r\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s003_invalid_header,
                "header line must start with '@' and contain a nonempty identifier",
                0,
                0,
                1,
            ) },
        },
        .{
            .data = "@r\r",
            .expected = .{ .parse_error = expectedCheckError(
                .s004_truncated_record,
                "unexpected end of file in sequence line",
                0,
                3,
                2,
            ) },
        },
        .{
            .data = "@r\r\nA\r",
            .expected = .{ .parse_error = expectedCheckError(
                .s004_truncated_record,
                "unexpected end of file in plus line",
                0,
                6,
                3,
            ) },
        },
        .{
            .data = "@r\r\nA\r\n+\r",
            .expected = .{ .parse_error = expectedCheckError(
                .s004_truncated_record,
                "unexpected end of file in quality line",
                0,
                9,
                4,
            ) },
        },
        .{
            .data = "@r\nR.\n+\n!!\n",
            .expected = .{ .parse_error = expectedCheckError(
                .s002_invalid_sequence_alphabet,
                "sequence byte is outside the selected alphabet",
                0,
                4,
                2,
            ) },
        },
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
                "header line must start with '@' and contain a nonempty identifier",
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
    for (1..2 * vector_len + 2) |field_len| {
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
            for (1..2 * vector_len + 2) |chunk_len| {
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
            for (1..2 * vector_len + 2) |chunk_len| {
                try expectCheckOutcome(
                    quality_error,
                    directCheckOutcome(data.items, chunk_len, .{}, .{}),
                );
            }
            data.items[quality_start + invalid_index] = '!';
        }
    }
}

test "[edge] - [check scanner]: counts and error locations cross 32 bits" {
    const crossing: u64 = 1 << 32;
    const valid = "@r\nA\n+\n!\n" ** 3;
    const cases = [_]struct {
        data: []const u8,
        code: LintCode,
        message: []const u8,
        offset: u64,
        line: u3,
        at_eof: bool = false,
    }{
        .{
            .data = "@bad\nA.\n+\n!!\n",
            .code = .s002_invalid_sequence_alphabet,
            .message = "sequence byte is outside the selected alphabet",
            .offset = crossing + 30,
            .line = 2,
        },
        .{
            .data = "@bad\nAA\n+\n! \n",
            .code = .s006_invalid_quality_range,
            .message = "quality byte must be ASCII 33 through 126",
            .offset = crossing + 35,
            .line = 4,
        },
        .{
            .data = "@bad\nA\n+\n",
            .code = .s004_truncated_record,
            .message = "unexpected end of file in quality line",
            .offset = crossing + 33,
            .line = 4,
            .at_eof = true,
        },
    };
    for ([_]usize{ 1, valid.len }) |chunk_len| {
        for (cases) |case| {
            errdefer std.debug.print("check 32-bit crossing: chunk {d}, code {s}\n", .{
                chunk_len, codeTag(case.code),
            });
            var scanner = CheckScanner.init(.{}, .{});
            scanner.byte_offset = crossing - 3;
            scanner.line_start_offset = scanner.byte_offset;
            scanner.record_index = crossing - 1;
            var cursor: usize = 0;
            while (cursor < valid.len) {
                const end = @min(cursor + chunk_len, valid.len);
                try std.testing.expectEqual(end - cursor, try scanner.feed(valid[cursor..end]));
                cursor = end;
            }
            try scanner.finishEof();
            try std.testing.expectEqual(crossing + 2, scanner.record_index);
            try std.testing.expectEqual(crossing + 24, scanner.byte_offset);

            if (case.at_eof) {
                try std.testing.expectEqual(case.data.len, try scanner.feed(case.data));
                try std.testing.expectError(error.Format, scanner.finishEof());
            } else {
                try std.testing.expectError(error.Format, scanner.feed(case.data));
            }
            try std.testing.expectEqual(crossing + 2, scanner.record_index);
            try std.testing.expectEqualDeep(ParseError{
                .code = case.code,
                .message = case.message,
                .record_index = crossing + 2,
                .byte_offset = case.offset,
                .line_in_record = case.line,
            }, scanner.takeLastError().?);
            try std.testing.expect(scanner.takeLastError() == null);
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
        error.ArithmeticLimit, error.OutOfMemory, error.Io => return err,
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

test "[property] - [record validation]: quality checks match scalar results" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse 16;
    const max_len = 4 * vector_len - 1;
    const quality = try std.testing.allocator.alloc(u8, max_len);
    defer std.testing.allocator.free(quality);

    @memset(quality, '!');
    for (0..max_len + 1) |length| {
        try std.testing.expectEqual(null, firstInvalidQualityScalar(quality[0..length], 0));
        try std.testing.expectEqual(null, firstInvalidQuality(quality[0..length]));
    }

    for (0..256) |value| {
        @memset(quality, @intCast(value));
        const expected: ?usize = if (value < 33 or value > 126) 0 else null;
        try std.testing.expectEqual(expected, firstInvalidQualityScalar(quality, 0));
        try std.testing.expectEqual(expected, firstInvalidQuality(quality));
        var failure: ?usize = null;
        var use_full_iupac = false;
        classifySemanticByte(.quality, .iupac, &use_full_iupac, &failure, @intCast(value), 0);
        try std.testing.expectEqual(expected, failure);
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

    for (1..vector_len) |remainder| {
        const length = 2 * vector_len + remainder;
        const overlap_index = length - vector_len;
        quality[overlap_index] = 32;
        try std.testing.expectEqual(
            overlap_index,
            firstInvalidQualityScalar(quality[0..length], 0).?,
        );
        try std.testing.expectEqual(overlap_index, firstInvalidQuality(quality[0..length]).?);
        quality[overlap_index] = '!';

        quality[length - 1] = 127;
        try std.testing.expectEqual(
            length - 1,
            firstInvalidQualityScalar(quality[0..length], 0).?,
        );
        try std.testing.expectEqual(length - 1, firstInvalidQuality(quality[0..length]).?);
        quality[length - 1] = '!';
    }

    quality[0] = 32;
    try std.testing.expectEqual(@as(usize, 0), firstInvalidQuality(quality[0..1]).?);
    quality[0] = '!';

    quality[vector_len] = 127;
    quality[1] = 32;
    try std.testing.expectEqual(@as(usize, 1), firstInvalidQualityScalar(quality, 0).?);
    try std.testing.expectEqual(@as(usize, 1), firstInvalidQuality(quality).?);
}

test "[property] - [record validation]: vector sequence validation matches scalar policies" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse 16;
    const max_len = 4 * vector_len - 1;
    const sequence = try std.testing.allocator.alloc(u8, max_len);
    defer std.testing.allocator.free(sequence);

    for ([_]Alphabet{ .iupac, .acgtn }) |alphabet| {
        @memset(sequence, 'A');
        for (0..max_len + 1) |length| {
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

        for (1..vector_len) |remainder| {
            const length = 2 * vector_len + remainder;
            const overlap_index = length - vector_len;
            sequence[overlap_index] = '.';
            try std.testing.expectEqual(
                overlap_index,
                firstInvalidSequenceScalar(sequence[0..length], alphabet, 0).?,
            );
            try std.testing.expectEqual(
                overlap_index,
                firstInvalidSequence(sequence[0..length], alphabet).?,
            );
            sequence[overlap_index] = 'A';

            sequence[length - 1] = 0xff;
            try std.testing.expectEqual(
                length - 1,
                firstInvalidSequenceScalar(sequence[0..length], alphabet, 0).?,
            );
            try std.testing.expectEqual(
                length - 1,
                firstInvalidSequence(sequence[0..length], alphabet).?,
            );
            sequence[length - 1] = 'A';
        }

        sequence[0] = '.';
        try std.testing.expectEqual(
            @as(usize, 0),
            firstInvalidSequence(sequence[0..1], alphabet).?,
        );
        sequence[0] = 'A';

        sequence[vector_len] = '.';
        sequence[1] = 0xff;
        try std.testing.expectEqual(
            @as(usize, 1),
            firstInvalidSequenceScalar(sequence, alphabet, 0).?,
        );
        try std.testing.expectEqual(@as(usize, 1), firstInvalidSequence(sequence, alphabet).?);
    }
}

test "[property] - [record validation]: adaptive IUPAC state preserves byte policy" {
    for (0..256) |value| {
        const sequence = [_]u8{@intCast(value)};
        const record: Record = .{
            .header = "r",
            .id = "r",
            .sequence = &sequence,
            .plus = "",
            .quality = "!",
        };
        var validator = AdaptiveRecordValidator.init(.{});
        try std.testing.expectEqualDeep(
            validateRecord(record, .{}),
            validator.validate(record),
        );
        try std.testing.expectEqual(
            alphabetAccepts(.iupac, sequence[0]) and
                !alphabetAccepts(.acgtn, sequence[0]),
            validator.use_full_iupac,
        );
    }

    var validator = AdaptiveRecordValidator.init(.{});
    const bad_quality: Record = .{
        .header = "r",
        .id = "r",
        .sequence = "A",
        .plus = "",
        .quality = " ",
    };
    try std.testing.expectEqualDeep(
        validateRecord(bad_quality, .{}),
        validator.validate(bad_quality),
    );
    try std.testing.expect(!validator.use_full_iupac);

    const wider: Record = .{
        .header = "r",
        .id = "r",
        .sequence = "R",
        .plus = "",
        .quality = "!",
    };
    try std.testing.expect(validator.validate(wider) == null);
    try std.testing.expect(validator.use_full_iupac);
    const invalid: Record = .{
        .header = "s",
        .id = "s",
        .sequence = ".",
        .plus = "",
        .quality = "!",
    };
    try std.testing.expectEqualDeep(
        validateRecord(invalid, .{}),
        validator.validate(invalid),
    );

    var independent = AdaptiveRecordValidator.init(.{});
    try std.testing.expect(!independent.use_full_iupac);
    const wider_bad_quality: Record = .{
        .header = "t",
        .id = "t",
        .sequence = "R",
        .plus = "",
        .quality = " ",
    };
    try std.testing.expectEqualDeep(
        validateRecord(wider_bad_quality, .{}),
        independent.validate(wider_bad_quality),
    );
    try std.testing.expect(independent.use_full_iupac);

    var acgtn = AdaptiveRecordValidator.init(.{ .alphabet = .acgtn });
    try std.testing.expectEqualDeep(
        validateRecord(wider, .{ .alphabet = .acgtn }),
        acgtn.validate(wider),
    );
    try std.testing.expect(!acgtn.use_full_iupac);
}

test "[property] - [record validation]: fused vectors preserve field precedence" {
    const vector_len = std.simd.suggestVectorLength(u8) orelse return;
    const length = 8 * vector_len - 1;
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
    for ([_]Alphabet{ .iupac, .acgtn }) |alphabet| {
        for (0..length + 1) |field_len| {
            const current = Record{
                .header = record.header,
                .id = record.id,
                .sequence = record.sequence[0..field_len],
                .plus = record.plus,
                .quality = record.quality[0..field_len],
            };
            const options: ValidationOptions = .{ .alphabet = alphabet };
            try std.testing.expect(validateRecord(current, options) == null);
            if (field_len == 0) continue;

            quality[0] = 32;
            sequence[field_len - 1] = '.';
            const sequence_error = validateRecord(current, options).?;
            try std.testing.expectEqual(SemanticField.sequence, sequence_error.field);
            try std.testing.expectEqual(field_len - 1, sequence_error.byte_index);

            sequence[field_len - 1] = 'A';
            const quality_error = validateRecord(current, options).?;
            try std.testing.expectEqual(SemanticField.quality, quality_error.field);
            try std.testing.expectEqual(@as(usize, 0), quality_error.byte_index);
            quality[0] = '!';
        }
    }
}
