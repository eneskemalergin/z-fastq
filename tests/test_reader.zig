//! Library exports, byte sources, Reader and count scanner behavior, and record validation.

const std = @import("std");
const zfastq = @import("z-fastq");

const GZIP_EMPTY = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x01, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00,
};

const GZIP_OPTIONAL_X = [_]u8{
    0x1f, 0x8b, 0x08, 0x1e, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x02, 0x00, 'x',  'y',  'n',  0x00, 'c',  0x00, 0xca, 0x4e,
    0x01, 0x01, 0x00, 0xfe, 0xff, 'x',  0x83, 0x16, 0xdc, 0x8c,
    0x01, 0x00, 0x00, 0x00,
};

const GZIP_FIXED = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x03,
    0xf3, 0x48, 0xcd, 0xc9, 0xc9, 0x57, 0x28, 0xcf, 0x2f, 0xca,
    0x49, 0xe1, 0x02, 0x00, 0xd5, 0xe0, 0x39, 0xb7, 0x0c, 0x00,
    0x00, 0x00,
};

const GZIP_DYNAMIC = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
    0x3d, 0xc6, 0x39, 0x11, 0x00, 0x00, 0x0c, 0x02, 0x30, 0x2b,
    0xb5, 0x52, 0x1e, 0xff, 0x96, 0x38, 0x16, 0x96, 0x5c, 0x1e,
    0x94, 0xcb, 0x6d, 0x01, 0x17, 0x1c, 0x39, 0xb4, 0x13, 0x00,
    0x00, 0x00,
};

const GZIP_CASES = [_]struct {
    name: []const u8,
    compressed: []const u8,
    decoded: []const u8,
    prefix: ?struct { compressed_len: usize, decoded_len: usize } = null,
}{
    .{ .name = "empty member", .compressed = &GZIP_EMPTY, .decoded = "" },
    .{ .name = "stored with optional header", .compressed = &GZIP_OPTIONAL_X, .decoded = "x" },
    .{ .name = "fixed", .compressed = &GZIP_FIXED, .decoded = "Hello world\n" },
    .{ .name = "dynamic", .compressed = &GZIP_DYNAMIC, .decoded = "ABCDEABCD ABCDEABCD" },
    .{
        .name = "stored and dynamic members",
        .compressed = &(GZIP_OPTIONAL_X ++ GZIP_DYNAMIC),
        .decoded = "xABCDEABCD ABCDEABCD",
        .prefix = .{ .compressed_len = GZIP_OPTIONAL_X.len, .decoded_len = 1 },
    },
};

const GZIP_CHUNKS = [_]struct { input: usize, output: usize }{
    .{ .input = 1, .output = 1 },
    .{ .input = 31, .output = 17 },
};

const GZIP_HEADER_FIELD_OFFSETS = [_]usize{ 0, 2, 3, 4, 8, 9, 10, 12, 14, 16, 18 };

const GZIP_FLAG_FCOMMENT: u8 = 0x10;

const GZIP_FIELD_LIMIT = 64 * 1024;

test "[unit] - [root]: every exported declaration is analyzable" {
    std.testing.refAllDecls(zfastq);
}

test "[unit] - [root]: limit namespaces and the count buffer name remain aliases" {
    try std.testing.expect(zfastq.limits == zfastq.io.limits);
    try std.testing.expectEqual(
        zfastq.limits.DEFAULT_READER_BUFFER_BYTES,
        zfastq.io.limits.COUNT_READ_BUFFER_BYTES,
    );
}

test "[unit] - [root]: library version matches the manifest and retains its sentinel" {
    const expected = @import("test_options").package_version;
    try std.testing.expectEqualStrings(expected, zfastq.VERSION);
    try std.testing.expect(@TypeOf(zfastq.VERSION) == *const [expected.len:0]u8);
}

test "[property] - [gzip source]: optional member chains decode at boundary input chunk sizes" {
    const gzip_xx = GZIP_OPTIONAL_X ++ GZIP_OPTIONAL_X;

    for ([_]usize{ 0, 2, 16 }) |buffer_len| {
        for ([_]usize{ 1, 2, gzip_xx.len }) |chunk_len| {
            var input_buffer: [16]u8 = undefined;
            var input = FragmentReader.init(&gzip_xx, chunk_len, input_buffer[0..buffer_len]);
            if (buffer_len >= 2) {
                try std.testing.expectEqualSlices(u8, gzip_xx[0..2], try input.interface.peek(2));
            }
            var gzip = zfastq.io.gzip.ReaderSource.init(&input.interface);
            const source = gzip.byteSource();
            var output: [2]u8 = undefined;

            try std.testing.expectEqual(@as(usize, 2), try source.read(&output));
            try std.testing.expectEqualStrings("xx", &output);
            try std.testing.expectEqual(@as(usize, 0), try source.read(&output));
        }
    }
}

test "[property] - [gzip source]: truncations accept only complete member prefixes" {
    var output: [4096]u8 = undefined;
    for (GZIP_CASES) |case| {
        const first_end = if (case.prefix) |prefix| prefix.compressed_len else case.compressed.len;
        const cuts = [_]usize{
            0,
            1,
            first_end - 1,
            first_end,
            @min(first_end + 1, case.compressed.len),
            case.compressed.len - 1,
            case.compressed.len,
        };
        for (cuts) |end| {
            for (GZIP_CHUNKS) |chunks| {
                errdefer std.debug.print("gzip {s}, cut {d}, chunks {d}/{d}\n", .{
                    case.name, end, chunks.input, chunks.output,
                });
                var expected: ?[]const u8 = null;
                if (end == case.compressed.len) expected = case.decoded;
                if (case.prefix) |prefix| {
                    if (end == prefix.compressed_len) {
                        expected = case.decoded[0..prefix.decoded_len];
                    }
                }

                const result = decodeGzip(case.compressed[0..end], chunks.input, chunks.output, &output);
                if (expected) |bytes| {
                    try std.testing.expectEqualStrings(bytes, try result);
                } else {
                    try std.testing.expectError(error.ReadFailed, result);
                }
            }
        }
    }
}

test "[property] - [gzip source]: one mutation per header field fails or preserves frozen output" {
    var mutated: [GZIP_OPTIONAL_X.len + GZIP_DYNAMIC.len]u8 = undefined;
    var output: [4096]u8 = undefined;

    for (GZIP_CASES) |case| {
        const compressed = mutated[0..case.compressed.len];
        @memcpy(compressed, case.compressed);
        for (GZIP_HEADER_FIELD_OFFSETS) |offset| {
            compressed[offset] = case.compressed[offset] ^ 0xff;
            for (GZIP_CHUNKS) |chunks| {
                errdefer std.debug.print("gzip {s}, byte {d}, chunks {d}/{d}\n", .{
                    case.name, offset, chunks.input, chunks.output,
                });
                const decoded = decodeGzip(compressed, chunks.input, chunks.output, &output) catch |err| switch (err) {
                    error.ReadFailed => continue,
                    else => return err,
                };
                try std.testing.expectEqualStrings(case.decoded, decoded);
            }
            compressed[offset] = case.compressed[offset];
        }
    }
}

test "[failure] - [gzip source]: header CRC, CRC32, and ISIZE are checked in every member" {
    const valid = GZIP_OPTIONAL_X ++ GZIP_OPTIONAL_X;
    var corrupted = valid;
    var output: [16]u8 = undefined;
    const checked_bytes = [_]usize{ 18, 19, 26, 27, 28, 29, 30, 31, 32, 33 };

    for ([_]usize{ 0, GZIP_OPTIONAL_X.len }) |member_start| {
        for (checked_bytes) |byte| {
            const offset = member_start + byte;
            for (0..8) |bit| {
                corrupted[offset] = valid[offset] ^ (@as(u8, 1) << @intCast(bit));
                for (GZIP_CHUNKS) |chunks| {
                    errdefer std.debug.print("gzip checksum byte {d}, bit {d}, chunks {d}/{d}\n", .{
                        offset, bit, chunks.input, chunks.output,
                    });
                    try std.testing.expectError(
                        error.ReadFailed,
                        decodeGzip(&corrupted, chunks.input, chunks.output, &output),
                    );
                }
            }
            corrupted[offset] = valid[offset];
        }
    }
}

test "[edge] - [gzip source]: history refills preserve bytes and validate the final trailer" {
    const compressed = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0xff,
        0xed, 0xc3, 0x31, 0x09, 0x00, 0x00, 0x0c, 0x03, 0x30, 0x6d,
        0xa5, 0xc7, 0x0c, 0xd4, 0xbf, 0x96, 0x19, 0x49, 0x20, 0xe9,
        0x2d, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
        0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
        0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
        0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
        0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
        0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa,
        0xaa, 0xaa, 0xaa, 0xaa, 0xda, 0xdb, 0x03, 0x57, 0xe1, 0x9b,
        0x36, 0x04, 0x00, 0x01, 0x00,
    };
    const expected = "ACGT" ** 16385;
    var output: [expected.len + 1]u8 = undefined;
    const trailer_offsets = [_]usize{ compressed.len - 8, compressed.len - 4 };

    for ([_]usize{ 17, 257, 32769 }) |output_chunk_len| {
        try std.testing.expectEqualStrings(
            expected,
            try decodeGzip(&compressed, 31, output_chunk_len, &output),
        );
        for (trailer_offsets) |offset| {
            var corrupted = compressed;
            corrupted[offset] ^= 1;
            try std.testing.expectError(
                error.ReadFailed,
                decodeGzip(&corrupted, 31, output_chunk_len, &output),
            );
        }
        try std.testing.expectError(
            error.ReadFailed,
            decodeGzip(compressed[0 .. compressed.len - 1], 31, output_chunk_len, &output),
        );
    }
}

test "[failure] - [gzip source]: a read failure after decoded output propagates" {
    var input_buffer: [10]u8 = undefined;
    var input = FragmentReader.init(&GZIP_OPTIONAL_X, 1, &input_buffer);
    input.fail_at = GZIP_OPTIONAL_X.len - 4;
    var gzip = zfastq.io.gzip.ReaderSource.init(&input.interface);
    const source = gzip.byteSource();
    var output: [1]u8 = undefined;

    try std.testing.expectEqual(@as(usize, 1), try source.read(&output));
    try std.testing.expectEqualStrings("x", &output);
    try std.testing.expectError(error.ReadFailed, source.read(&output));
}

test "[failure] - [gzip source]: bytes after the last member fail" {
    var output: [1]u8 = undefined;
    const trailing_text = GZIP_EMPTY ++ [_]u8{'x'};
    try std.testing.expectError(error.ReadFailed, decodeGzip(&trailing_text, 31, 1, &output));
    const zero_padding = GZIP_EMPTY ++ [_]u8{0};
    try std.testing.expectError(error.ReadFailed, decodeGzip(&zero_padding, 31, 1, &output));
}

test "[edge] - [gzip source]: a comment field with its terminator at the limit is accepted and one byte more fails" {
    var output: [1]u8 = undefined;
    for ([_]usize{ GZIP_FIELD_LIMIT - 1, GZIP_FIELD_LIMIT }) |comment_len| {
        const compressed = try std.testing.allocator.alloc(u8, comment_len + GZIP_EMPTY.len + 1);
        defer std.testing.allocator.free(compressed);
        @memcpy(compressed[0..10], GZIP_EMPTY[0..10]);
        compressed[3] = GZIP_FLAG_FCOMMENT;
        @memset(compressed[10..][0..comment_len], 'c');
        compressed[10 + comment_len] = 0;
        @memcpy(compressed[11 + comment_len ..], GZIP_EMPTY[10..]);

        var input = std.Io.Reader.fixed(compressed);
        var gzip = zfastq.io.gzip.ReaderSource.init(&input);
        const source = gzip.byteSource();
        if (comment_len + 1 <= GZIP_FIELD_LIMIT) {
            try std.testing.expectEqual(@as(usize, 0), try source.read(&output));
        } else {
            try std.testing.expectError(error.ReadFailed, source.read(&output));
        }
    }
}

test "[property] - [plain reader source]: sniffed prefixes survive short reads and EOF" {
    const cases = [_][]const u8{ "", "x", "@r", "\x1f", "\x1f\x8a", "@r\nA\n+\n!\n" };
    for (cases) |data| {
        for ([_]usize{ 1, 3 }) |input_chunk_len| {
            for ([_]usize{ 1, 3, 16 }) |output_chunk_len| {
                errdefer std.debug.print("plain prefix length {d}, chunks {d}/{d}\n", .{
                    data.len, input_chunk_len, output_chunk_len,
                });
                var input_buffer: [2]u8 = undefined;
                var input = FragmentReader.init(data, input_chunk_len, &input_buffer);
                if (data.len >= 2) {
                    try std.testing.expectEqualStrings(data[0..2], try input.interface.peek(2));
                } else {
                    _ = input.interface.peek(2) catch {};
                }

                var plain = zfastq.io.plain.ReaderSource.init(&input.interface);
                const source = plain.byteSource();
                var output: [16]u8 = undefined;
                var offset: usize = 0;
                for (0..data.len + 1) |_| {
                    const n = try source.read(output[0..output_chunk_len]);
                    try std.testing.expect(n <= data.len - offset);
                    try std.testing.expectEqualStrings(data[offset..][0..n], output[0..n]);
                    offset += n;
                    if (n == 0) break;
                }
                try std.testing.expectEqual(data.len, offset);
                try std.testing.expectEqual(@as(usize, 0), try source.read(&output));
            }
        }
    }
}

test "[unit] - [lint code]: every implemented code has its stable tag" {
    const cases = [_]struct {
        code: zfastq.LintCode,
        tag: []const u8,
    }{
        .{ .code = .s001_invalid_plus_line, .tag = "S001" },
        .{ .code = .s002_invalid_sequence_alphabet, .tag = "S002" },
        .{ .code = .s003_invalid_header, .tag = "S003" },
        .{ .code = .s004_truncated_record, .tag = "S004" },
        .{ .code = .s005_length_mismatch, .tag = "S005" },
        .{ .code = .s006_invalid_quality_range, .tag = "S006" },
    };

    for (cases) |case| {
        try std.testing.expectEqualStrings(case.tag, zfastq.codeTag(case.code));
    }
}

test "[property] - [record validation]: alphabet policies classify every byte" {
    const policies = [_]struct {
        alphabet: zfastq.Alphabet,
        accepted: []const u8,
    }{
        .{
            .alphabet = .iupac,
            .accepted = "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
        },
        .{
            .alphabet = .acgtn,
            .accepted = "ACGTNacgtn",
        },
    };

    for (policies) |policy| {
        try std.testing.expect(validateSequenceQuality("", "", .{ .alphabet = policy.alphabet }) == null);

        for (0..256) |value| {
            const sequence = [1]u8{@intCast(value)};
            const details = validateSequenceQuality(&sequence, "!", .{ .alphabet = policy.alphabet });
            const accepted = std.mem.findScalar(u8, policy.accepted, sequence[0]) != null;
            try std.testing.expectEqual(accepted, details == null);
            if (details) |semantic_error| {
                try std.testing.expectEqual(
                    zfastq.LintCode.s002_invalid_sequence_alphabet,
                    semantic_error.code,
                );
                try std.testing.expectEqual(zfastq.SemanticField.sequence, semantic_error.field);
                try std.testing.expectEqual(@as(usize, 0), semantic_error.byte_index);
                try std.testing.expectEqualStrings(
                    "sequence byte is outside the selected alphabet",
                    semantic_error.message,
                );
            }
        }
    }
}

test "[edge] - [record validation]: each invalid sequence position is reported exactly" {
    var sequence = [_]u8{ 'A', 'C', 'G', 'T' };
    const quality = [_]u8{ '!', '!', '!', '!' };
    for (sequence, 0..) |original, invalid_index| {
        sequence[invalid_index] = '.';
        const details = validateSequenceQuality(&sequence, &quality, .{}).?;
        try std.testing.expectEqual(zfastq.LintCode.s002_invalid_sequence_alphabet, details.code);
        try std.testing.expectEqual(zfastq.SemanticField.sequence, details.field);
        try std.testing.expectEqual(invalid_index, details.byte_index);
        sequence[invalid_index] = original;
    }
}

test "[edge] - [record validation]: each invalid quality position is reported exactly" {
    const sequence = [_]u8{ 'A', 'C', 'G', 'T' };
    var quality = [_]u8{ '!', '!', '!', '!' };
    for (quality, 0..) |original, invalid_index| {
        quality[invalid_index] = if (invalid_index % 2 == 0) ' ' else 0x7f;
        const details = validateSequenceQuality(&sequence, &quality, .{}).?;
        try std.testing.expectEqual(zfastq.LintCode.s006_invalid_quality_range, details.code);
        try std.testing.expectEqual(zfastq.SemanticField.quality, details.field);
        try std.testing.expectEqual(invalid_index, details.byte_index);
        try std.testing.expectEqualStrings(
            "quality byte must be ASCII 33 through 126",
            details.message,
        );
        quality[invalid_index] = original;
    }
}

test "[unit] - [record validation]: a sequence error takes precedence over a quality error" {
    const details = validateSequenceQuality(".", "\x7f", .{}).?;

    try std.testing.expectEqual(zfastq.LintCode.s002_invalid_sequence_alphabet, details.code);
    try std.testing.expectEqual(zfastq.SemanticField.sequence, details.field);
    try std.testing.expectEqual(@as(usize, 0), details.byte_index);
}

test "[property] - [record validation]: quality range classifies every byte" {
    try std.testing.expect(validateSequenceQuality("AC", "!~", .{}) == null);

    for (0..256) |value| {
        const quality = [1]u8{@intCast(value)};
        const details = validateSequenceQuality("A", &quality, .{});
        try std.testing.expectEqual(value >= '!' and value <= '~', details == null);
        if (details) |semantic_error| {
            try std.testing.expectEqual(
                zfastq.LintCode.s006_invalid_quality_range,
                semantic_error.code,
            );
            try std.testing.expectEqual(zfastq.SemanticField.quality, semantic_error.field);
            try std.testing.expectEqual(@as(usize, 0), semantic_error.byte_index);
            try std.testing.expectEqualStrings(
                "quality byte must be ASCII 33 through 126",
                semantic_error.message,
            );
        }
    }
}

test "[unit] - [reader]: current record offsets identify every line start" {
    const data = "@one\r\nA\r\n+\r\n!\r\n@two\nCC\n+x\n##";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    try std.testing.expect(reader.currentRecordOffsets() == null);

    const first = (try reader.next()).?;
    try std.testing.expectEqualStrings("one", first.id);
    try std.testing.expectEqual(
        zfastq.RecordOffsets{ .header = 0, .sequence = 6, .plus = 9, .quality = 12 },
        reader.currentRecordOffsets().?,
    );

    const second = (try reader.next()).?;
    try std.testing.expectEqualStrings("two", second.id);
    try std.testing.expectEqual(
        zfastq.RecordOffsets{ .header = 15, .sequence = 20, .plus = 23, .quality = 26 },
        reader.currentRecordOffsets().?,
    );

    try std.testing.expect((try reader.next()) == null);
    try std.testing.expect(reader.currentRecordOffsets() == null);

    var advance_source = zfastq.io.plain.SliceSource.init("@r\nA\n+\n!\n");
    var advance_reader = try zfastq.Reader.init(
        std.testing.allocator,
        advance_source.byteSource(),
        .{},
    );
    defer advance_reader.deinit();
    try std.testing.expect(try advance_reader.advance());
    try std.testing.expect(advance_reader.currentRecordOffsets() == null);
}

test "[unit] - [reader]: supported record forms parse exactly" {
    const cases = [_]struct {
        data: []const u8,
        header: []const u8,
        id: []const u8,
        sequence: []const u8,
        plus: []const u8,
        quality: []const u8,
    }{
        .{
            .data = "@read1\nACGT\n+\n!!!!\n",
            .header = "read1",
            .id = "read1",
            .sequence = "ACGT",
            .plus = "",
            .quality = "!!!!",
        },
        .{
            .data = "@r1\r\nACGT\r\n+\r\n!!!!\r\n",
            .header = "r1",
            .id = "r1",
            .sequence = "ACGT",
            .plus = "",
            .quality = "!!!!",
        },
        .{
            .data = "@read3\nAAAA\n+repeat_id\n!!!!\n",
            .header = "read3",
            .id = "read3",
            .sequence = "AAAA",
            .plus = "repeat_id",
            .quality = "!!!!",
        },
        .{
            .data = "@final\nACGT\n+\nIIII",
            .header = "final",
            .id = "final",
            .sequence = "ACGT",
            .plus = "",
            .quality = "IIII",
        },
        .{
            .data = "@read1\tcomment\r\nACGT\n+name\r\n!!!!\n",
            .header = "read1\tcomment",
            .id = "read1",
            .sequence = "ACGT",
            .plus = "name",
            .quality = "!!!!",
        },
        .{
            .data = "@empty\n\n+description\n\n",
            .header = "empty",
            .id = "empty",
            .sequence = "",
            .plus = "description",
            .quality = "",
        },
        .{
            .data = "@r\rcomment\nA\rC\n+x\ry\n!\r!\n",
            .header = "r\rcomment",
            .id = "r\rcomment",
            .sequence = "A\rC",
            .plus = "x\ry",
            .quality = "!\r!",
        },
    };

    for (cases) |case| {
        var source = zfastq.io.plain.SliceSource.init(case.data);
        var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();

        const record = (try reader.next()).?;

        try std.testing.expectEqualStrings(case.header, record.header);
        try std.testing.expectEqualStrings(case.id, record.id);
        try std.testing.expectEqualStrings(case.sequence, record.sequence);
        try std.testing.expectEqualStrings(case.plus, record.plus);
        try std.testing.expectEqualStrings(case.quality, record.quality);
        try std.testing.expect((try reader.next()) == null);
    }
}

const ExpectedRecord = struct {
    header: []const u8,
    id: []const u8,
    sequence: []const u8,
    plus: []const u8,
    quality: []const u8,
    offsets: ?zfastq.RecordOffsets = null,
};

test "[property] - [reader]: fields and offsets are identical at every refill boundary" {
    const cases = [_]struct {
        data: []const u8,
        records: []const ExpectedRecord,
    }{
        .{
            .data = "@one comment\r\nACGT\r\n+same\r\n!!!!\r\n" ++
                "@empty\n\n+description\n\n" ++
                "@final\nNN\n+\n##",
            .records = &.{
                .{
                    .header = "one comment",
                    .id = "one",
                    .sequence = "ACGT",
                    .plus = "same",
                    .quality = "!!!!",
                    .offsets = .{ .header = 0, .sequence = 14, .plus = 20, .quality = 27 },
                },
                .{
                    .header = "empty",
                    .id = "empty",
                    .sequence = "",
                    .plus = "description",
                    .quality = "",
                    .offsets = .{ .header = 33, .sequence = 40, .plus = 41, .quality = 54 },
                },
                .{
                    .header = "final",
                    .id = "final",
                    .sequence = "NN",
                    .plus = "",
                    .quality = "##",
                    .offsets = .{ .header = 55, .sequence = 62, .plus = 65, .quality = 67 },
                },
            },
        },
        .{
            .data = "@read1 comment\nACGT\n+repeat\n!!!!\n",
            .records = &.{.{
                .header = "read1 comment",
                .id = "read1",
                .sequence = "ACGT",
                .plus = "repeat",
                .quality = "!!!!",
            }},
        },
        .{
            .data = "@r\r\r\nA\r\r\n+\r\r\n!\r\r\n",
            .records = &.{.{
                .header = "r\r",
                .id = "r\r",
                .sequence = "A\r",
                .plus = "\r",
                .quality = "!\r",
            }},
        },
        .{
            .data = "@r\nAA\n+\n!\r",
            .records = &.{.{
                .header = "r",
                .id = "r",
                .sequence = "AA",
                .plus = "",
                .quality = "!\r",
            }},
        },
    };
    const chunk_sizes = [_]usize{
        1,
        2,
        3,
        15,
        16,
        17,
        31,
        32,
        33,
        63,
        64,
        65,
        4 * 1024,
        zfastq.limits.DEFAULT_READER_BUFFER_BYTES,
    };

    for (cases) |case| {
        for (chunk_sizes) |chunk_len| {
            var chunked = ChunkSource.init(case.data, chunk_len);
            var reader = try zfastq.Reader.init(
                std.testing.allocator,
                chunked.byteSource(),
                .{},
            );
            defer reader.deinit();

            for (case.records) |expected| {
                const parsed = (try reader.next()).?;
                try std.testing.expectEqualStrings(expected.header, parsed.header);
                try std.testing.expectEqualStrings(expected.id, parsed.id);
                try std.testing.expectEqualStrings(expected.sequence, parsed.sequence);
                try std.testing.expectEqualStrings(expected.plus, parsed.plus);
                try std.testing.expectEqualStrings(expected.quality, parsed.quality);
                if (expected.offsets) |offsets| {
                    try std.testing.expectEqual(offsets, reader.currentRecordOffsets().?);
                }
            }
            try std.testing.expect((try reader.next()) == null);
            try std.testing.expectEqual(@as(u64, case.data.len), reader.byteOffset());
        }
    }
}

test "[property] - [parser]: every input partition preserves the format result" {
    const valid = "@r\nA\n+\n!\n";
    const masks = [_]usize{ 0, 0xff, 1 << 2, 1 << 4, 1 << 6 };
    for (masks) |mask| {
        var partitioned = ChunkSource.init(valid, valid.len);
        partitioned.split_mask = mask;
        const source = partitioned.byteSource();
        var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
        defer reader.deinit();
        const parsed = (try reader.next()).?;
        try std.testing.expectEqualStrings("r", parsed.id);
        try std.testing.expectEqualStrings("A", parsed.sequence);
        try std.testing.expectEqualStrings("!", parsed.quality);

        var scan = zfastq.count_scan.Scanner.init(.{});
        var pos: usize = 0;
        while (pos < valid.len) {
            const end = partitionEnd(valid.len, pos, mask);
            _ = try scan.feed(valid[pos..end]);
            pos = end;
        }
        try scan.finishEof();
        try std.testing.expectEqual(@as(u64, 1), scan.record_index);
    }

    const invalid = "@r\nA\nx\n!\n";
    for (masks) |mask| {
        var partitioned = ChunkSource.init(invalid, invalid.len);
        partitioned.split_mask = mask;
        const source = partitioned.byteSource();
        var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
        defer reader.deinit();
        try std.testing.expectError(zfastq.ReaderError.S001InvalidPlusLine, reader.next());
        try std.testing.expectEqual(@as(u64, 5), reader.takeLastError().?.byte_offset);

        var scan = zfastq.count_scan.Scanner.init(.{});
        var pos: usize = 0;
        var found_error = false;
        while (pos < invalid.len) {
            const end = partitionEnd(invalid.len, pos, mask);
            if (scan.feed(invalid[pos..end])) |_| {
                pos = end;
            } else |err| {
                try std.testing.expectEqual(zfastq.ReaderError.S001InvalidPlusLine, err);
                found_error = true;
                break;
            }
        }
        try std.testing.expect(found_error);
        try std.testing.expectEqual(@as(u64, 5), scan.takeLastError().?.byte_offset);
    }
}

test "[edge] - [reader]: fields of a line one byte past the read buffer survive the next read" {
    const line_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES + 1;
    const next_record = "\n@next comment\nAC\n+note\n!!\n";
    const total_len = 3 + line_len + 3 + line_len + next_record.len;
    const data = try std.testing.allocator.alloc(u8, total_len);
    defer std.testing.allocator.free(data);

    var pos: usize = 0;
    @memcpy(data[pos..][0..4], "@r\nA");
    pos += 4;
    @memset(data[pos..][0 .. line_len - 1], 'A');
    pos += line_len - 1;
    @memcpy(data[pos..][0..3], "\n+\n");
    pos += 3;
    @memset(data[pos..][0..line_len], 'I');
    pos += line_len;
    @memcpy(data[pos..][0..next_record.len], next_record);

    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    const parsed = (try reader.next()).?;
    try std.testing.expectEqual(line_len, parsed.sequence.len);
    try std.testing.expectEqual(line_len, parsed.quality.len);
    try std.testing.expect(std.mem.allEqual(u8, parsed.sequence, 'A'));
    try std.testing.expect(std.mem.allEqual(u8, parsed.quality, 'I'));

    const next = (try reader.next()).?;
    try std.testing.expectEqualStrings("next comment", next.header);
}

test "[unit] - [reader]: complete buffered records allocate only the transport buffer" {
    const data = "@r\nACGT\n+\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 1 });
    var reader = try zfastq.Reader.init(failing.allocator(), source.byteSource(), .{});
    defer reader.deinit();

    const record = (try reader.next()).?;
    try std.testing.expectEqualStrings("ACGT", record.sequence);
    try std.testing.expect((try reader.next()) == null);
    try std.testing.expectEqual(@as(usize, 1), failing.allocations);
    try std.testing.expectEqual(zfastq.limits.DEFAULT_READER_BUFFER_BYTES, failing.allocated_bytes);
    try std.testing.expect(!failing.has_induced_failure);
}

test "[unit] - [owned record]: copied fields remain stable after Reader advances" {
    const data = "@one\nAAAA\n+\n!!!!\n@two\nCCCC\n+\n####\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    var owned = try zfastq.toOwned(std.testing.allocator, (try reader.next()).?);
    defer owned.deinit();
    _ = try reader.next();

    try std.testing.expectEqualStrings("one", owned.header);
    try std.testing.expectEqualStrings("AAAA", owned.sequence);
    try std.testing.expectEqualStrings("!!!!", owned.quality);
}

test "[failure] - [reader]: initialization releases partial allocations" {
    var source = zfastq.io.plain.SliceSource.init("@r\nA\n+\n!\n");
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    try std.testing.expectError(
        error.OutOfMemory,
        zfastq.Reader.init(failing.allocator(), source.byteSource(), .{}),
    );
}

test "[failure] - [owned record]: conversion releases partial allocations" {
    const record = zfastq.Record{
        .header = "read comment",
        .id = "read",
        .sequence = "ACGT",
        .plus = "read comment",
        .quality = "!!!!",
    };
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 4 });
    try std.testing.expectError(error.OutOfMemory, zfastq.toOwned(failing.allocator(), record));
}

test "[failure] - [reader]: read failures at every byte retain complete earlier records" {
    const record = "@r\nA\n+\n!\n";
    const data = record ++ record;
    for (0..data.len + 1) |fail_at| {
        errdefer std.debug.print("read failure after {d} bytes\n", .{fail_at});
        var failing_source = ChunkSource.init(data, data.len);
        failing_source.fail_at = fail_at;
        var reader = try zfastq.Reader.init(
            std.testing.allocator,
            failing_source.byteSource(),
            .{},
        );
        defer reader.deinit();

        for (0..fail_at / record.len) |_| {
            const parsed = (try reader.next()).?;
            try std.testing.expectEqualStrings("r", parsed.id);
            try std.testing.expectEqualStrings("A", parsed.sequence);
            try std.testing.expectEqualStrings("!", parsed.quality);
        }
        try std.testing.expectError(zfastq.ReaderError.Io, reader.next());
        try std.testing.expectEqual(@as(u64, fail_at / record.len), reader.recordIndex());
        try std.testing.expect(reader.takeLastError() == null);
    }
}

test "[failure] - [byte source]: a read count beyond its destination is rejected at each layer" {
    var direct = InvalidCountSource{};
    var output: [1]u8 = undefined;
    try std.testing.expectError(error.ReadFailed, direct.byteSource().read(&output));
    try std.testing.expectEqual(@as(usize, 1), direct.read_count);

    var invalid_source = InvalidCountSource{};
    var reader = try zfastq.Reader.init(
        std.testing.allocator,
        invalid_source.byteSource(),
        .{},
    );
    defer reader.deinit();

    try std.testing.expectError(zfastq.ReaderError.Io, reader.next());
}

test "[unit] - [reader]: the ByteSource wrapper is copied by value" {
    var original = zfastq.io.plain.SliceSource.init("@original\nA\n+\n!\n");
    var replacement = zfastq.io.plain.SliceSource.init("@replacement\nC\n+\n#\n");
    var source = original.byteSource();
    var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
    defer reader.deinit();

    source = replacement.byteSource();
    const parsed = (try reader.next()).?;

    try std.testing.expectEqualStrings("original", parsed.id);
}

test "[edge] - [byte source]: an empty destination does not consume input" {
    var adapter = zfastq.io.plain.SliceSource.init("abc");
    var source = adapter.byteSource();
    var empty: [0]u8 = .{};
    var first: [1]u8 = undefined;

    try std.testing.expectEqual(@as(usize, 0), try source.read(&empty));
    try std.testing.expectEqual(@as(usize, 1), try source.read(&first));
    try std.testing.expectEqual(@as(u8, 'a'), first[0]);
}

test "[failure] - [reader]: allocation failure while growing storage propagates" {
    const line_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES + 1;
    const data = try std.testing.allocator.alloc(u8, line_len + 4);
    defer std.testing.allocator.free(data);
    @memcpy(data[0..4], "@r\nA");
    @memset(data[4..], 'A');

    var source = zfastq.io.plain.SliceSource.init(data);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 2,
        .resize_fail_index = 0,
    });
    var reader = try zfastq.Reader.init(failing.allocator(), source.byteSource(), .{});
    defer reader.deinit();

    try std.testing.expectError(zfastq.ReaderError.OutOfMemory, reader.next());
}

test "[edge] - [count slice]: the configured line limit is exact on every record line" {
    const LIMIT = 4;
    const cases = [_]struct { data: []const u8, accepted: bool }{
        .{ .data = "@abc\nA\n+\n!\n", .accepted = true },
        .{ .data = "@abcd\nA\n+\n!\n", .accepted = false },
        .{ .data = "@r\nACGT\n+\n!!!!\n", .accepted = true },
        .{ .data = "@r\nACGTA\n+\n!!!!!\n", .accepted = false },
        .{ .data = "@r\nA\n+abc\n!\n", .accepted = true },
        .{ .data = "@r\nA\n+abcd\n!\n", .accepted = false },
        .{ .data = "@r\nAAAA\n+\n!!!!\n", .accepted = true },
        .{ .data = "@r\nAAAA\n+\n!!!!!\n", .accepted = false },
        .{ .data = "@r\nACGT\n+\n!!!!\r\n", .accepted = true },
        .{ .data = "@r\nACGT\n+\n!!!!\r", .accepted = false },
    };

    for (cases) |case| {
        var source = zfastq.io.plain.SliceSource.init(case.data);
        var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{ .max_line_bytes = LIMIT });
        defer reader.deinit();
        var scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = LIMIT });
        if (case.accepted) {
            _ = try reader.next();
            try std.testing.expectEqual(
                @as(u64, 1),
                try zfastq.count_scan.countSlice(case.data, .{ .max_line_bytes = LIMIT }, &scan),
            );
        } else {
            try std.testing.expectError(zfastq.ReaderError.LineTooLong, reader.next());
            try std.testing.expectError(
                zfastq.ReaderError.LineTooLong,
                zfastq.count_scan.countSlice(case.data, .{ .max_line_bytes = LIMIT }, &scan),
            );
        }
    }

    const lone_cr = "@r\nACGT\n+\n!!!!\r";
    var eof_source = zfastq.io.plain.SliceSource.init(lone_cr);
    var eof_reader = try zfastq.Reader.init(std.testing.allocator, eof_source.byteSource(), .{ .max_line_bytes = LIMIT });
    defer eof_reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.LineTooLong, eof_reader.next());
    try std.testing.expectEqual(@as(u64, lone_cr.len), eof_reader.byteOffset());

    var eof_scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = LIMIT });
    try std.testing.expectError(
        zfastq.ReaderError.LineTooLong,
        zfastq.count_scan.countSlice(lone_cr, .{ .max_line_bytes = LIMIT }, &eof_scan),
    );
    try std.testing.expectEqual(@as(u64, lone_cr.len), eof_scan.byte_offset);

    const dense = "@a\nAC\n+\n!!\n@b\nACGT\n+\n!!!!\n";
    var dense_scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 3 });
    try std.testing.expectError(
        zfastq.ReaderError.LineTooLong,
        zfastq.count_scan.countSlice(dense, .{ .max_line_bytes = 3 }, &dense_scan),
    );
    try std.testing.expectEqual(@as(u64, 1), dense_scan.record_index);

    const span = "@r\nAAAAAAAA\n+\n!!!!!!!!";
    var byte_scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 8 });
    for (span) |byte| _ = try byte_scan.feed(&.{byte});
    try byte_scan.finishEof();
    try std.testing.expectEqual(@as(u64, 1), byte_scan.record_index);

    var long_scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 7 });
    try std.testing.expectError(zfastq.ReaderError.LineTooLong, long_scan.feed(span));
}

test "[property] - [reader]: advance, next, and the count scanner agree on record count" {
    const data =
        \\@read1
        \\ACGT
        \\+
        \\!!!!
        \\@read2
        \\TGCA
        \\+
        \\JJJJ
        \\
    ;

    var scan = zfastq.count_scan.Scanner.init(.{});
    _ = try scan.feed(data);
    try scan.finishEof();
    try std.testing.expectEqual(@as(u64, 2), scan.record_index);

    var advance_source = zfastq.io.plain.SliceSource.init(data);
    var advance_reader = try zfastq.Reader.init(std.testing.allocator, advance_source.byteSource(), .{});
    defer advance_reader.deinit();
    var advance_count: u64 = 0;
    while (try advance_reader.advance()) {
        advance_count += 1;
    }
    try std.testing.expectEqual(@as(u64, 2), advance_count);

    var next_source = zfastq.io.plain.SliceSource.init(data);
    var next_reader = try zfastq.Reader.init(std.testing.allocator, next_source.byteSource(), .{});
    defer next_reader.deinit();
    var next_count: u64 = 0;
    while ((try next_reader.next()) != null) {
        next_count += 1;
    }
    try std.testing.expectEqual(@as(u64, 2), next_count);
}

test "[unit] - [count slice]: reuse replaces progress, options, and old diagnostics" {
    var scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 0 });
    try std.testing.expectError(
        error.S001InvalidPlusLine,
        zfastq.count_scan.countSlice("@r\nA\nx\n!\n", .{}, &scan),
    );

    const valid = "@r\nA\n+\n!\n";
    const options: zfastq.count_scan.Options = .{ .max_line_bytes = 2 };
    try std.testing.expectEqual(
        @as(u64, 1),
        try zfastq.count_scan.countSlice(valid, options, &scan),
    );
    try std.testing.expectEqual(@as(u64, valid.len), scan.byte_offset);
    try std.testing.expect(scan.takeLastError() == null);
    try std.testing.expectEqual(
        @as(u64, 0),
        try zfastq.count_scan.countSlice("", .{ .max_line_bytes = 0 }, &scan),
    );
    try std.testing.expectEqual(@as(u64, 0), scan.byte_offset);
    try std.testing.expect(scan.takeLastError() == null);
    try std.testing.expectError(error.LineTooLong, scan.feed(valid));
}

test "[property] - [parser]: Reader and scanner agree on mixed-ending records and fixed mutations" {
    const valid = [_]struct { data: []const u8, count: u64 }{
        .{ .data = "@r1\r\nACGT\n+x\r\n!!!!\n@r2\nN\n+\n#\n", .count = 2 },
        .{ .data = "@a\nA\n+\n!\n@b\n\n+\n\n", .count = 2 },
        .{ .data = "@r\nACGT\n+\nIIII", .count = 1 },
    };
    for (valid) |case| {
        for ([_]usize{ 1, 2, case.data.len }) |chunk_len| {
            const reader_valid = try readerOutcome(case.data, chunk_len);
            try std.testing.expectEqual(case.count, reader_valid.count);
            try std.testing.expectEqual(@as(?zfastq.ReaderError, null), reader_valid.err);
            try expectSameOutcome(reader_valid, scannerOutcome(case.data, chunk_len));
        }
    }

    const base = "@r1\r\nACGT\n+x\r\n!!!!\n@r2\nN\n+\n#\n";
    const mutations = [_]struct { offset: usize, byte: u8 }{
        .{ .offset = 0, .byte = '\n' },
        .{ .offset = 4, .byte = '\r' },
        .{ .offset = 7, .byte = 0 },
        .{ .offset = 10, .byte = '@' },
        .{ .offset = 14, .byte = '\n' },
        .{ .offset = 19, .byte = '\r' },
        .{ .offset = 26, .byte = 'A' },
        .{ .offset = 27, .byte = '\r' },
    };
    for (mutations) |mutation| {
        var mutated: [base.len]u8 = undefined;
        @memcpy(&mutated, base);
        mutated[mutation.offset] = mutation.byte;
        for ([_]usize{ 1, 2, mutated.len }) |chunk_len| {
            errdefer std.debug.print("mutation at {d} to 0x{x}, chunk {d}\n", .{
                mutation.offset, mutation.byte, chunk_len,
            });
            const reader_mutated = try readerOutcome(&mutated, chunk_len);
            try expectSameOutcome(reader_mutated, scannerOutcome(&mutated, chunk_len));
        }
    }
}

test "[fuzz] - [parser]: arbitrary bytes keep Reader and scanner in agreement" {
    const corpus = [_][]const u8{
        "\x09\x00\x00\x00@r\nA\n+\n!\n",
        "\x08\x00\x00\x00@\nA\n+\n!\n",
        "\x0a\x00\x00\x00@ r\nA\n+\n!\n",
        "\x04\x00\x00\x00\x00\xff\n@",
    };
    try std.testing.fuzz({}, fuzzReaderScannerAgreement, .{ .corpus = &corpus });
}

test "[property] - [parser]: structural details agree across chunk sizes" {
    const valid = "@r1\nAAAA\n+\n!!!!\n@r2\nCC\n+note\n##\n";
    for ([_]usize{ 1, 2, 3, valid.len }) |chunk_len| {
        var scan = zfastq.count_scan.Scanner.init(.{});
        var pos: usize = 0;
        while (pos < valid.len) {
            const end = @min(valid.len, pos + chunk_len);
            _ = try scan.feed(valid[pos..end]);
            pos = end;
        }
        try scan.finishEof();
        try std.testing.expectEqual(@as(u64, 2), scan.record_index);
        try std.testing.expectEqual(@as(u64, valid.len), scan.byte_offset);
    }

    const cases = [_]struct {
        data: []const u8,
        expected_error: zfastq.ReaderError,
        code: zfastq.LintCode,
        record_index: u64 = 0,
        line: u3,
        offset: u64,
    }{
        .{
            .data = "bad\nA\n+\n!\n",
            .expected_error = error.S003InvalidHeader,
            .code = .s003_invalid_header,
            .line = 1,
            .offset = 0,
        },
        .{
            .data = "not_a_header\nACGT\n+\n!!!!\n",
            .expected_error = error.S003InvalidHeader,
            .code = .s003_invalid_header,
            .line = 1,
            .offset = 0,
        },
        .{
            .data = "@bad\nACGT\nnot-plus\n!!!!\n",
            .expected_error = error.S001InvalidPlusLine,
            .code = .s001_invalid_plus_line,
            .line = 3,
            .offset = 10,
        },
        .{
            .data = "@r1\nAAAA\n+\n!!!!\n@r2\nA\nAA\n+\n!!!!\n",
            .expected_error = error.S001InvalidPlusLine,
            .code = .s001_invalid_plus_line,
            .record_index = 1,
            .line = 3,
            .offset = 22,
        },
        .{
            .data = "@truncated\nACGT\n+",
            .expected_error = error.S004TruncatedRecord,
            .code = .s004_truncated_record,
            .line = 4,
            .offset = 17,
        },
        .{
            .data = "@bad\nACGT\n+\nIII\n",
            .expected_error = error.S005LengthMismatch,
            .code = .s005_length_mismatch,
            .line = 4,
            .offset = 12,
        },
        .{
            .data = "@\nA\n+\n!\n",
            .expected_error = error.S003InvalidHeader,
            .code = .s003_invalid_header,
            .line = 1,
            .offset = 0,
        },
        .{
            .data = "@ description\nA\n+\n!\n",
            .expected_error = error.S003InvalidHeader,
            .code = .s003_invalid_header,
            .line = 1,
            .offset = 0,
        },
        .{
            .data = "@r\nA\nx\n!\n",
            .expected_error = error.S001InvalidPlusLine,
            .code = .s001_invalid_plus_line,
            .line = 3,
            .offset = 5,
        },
        .{
            .data = "@r\nA\n+\n",
            .expected_error = error.S004TruncatedRecord,
            .code = .s004_truncated_record,
            .line = 4,
            .offset = 7,
        },
        .{
            .data = "@r\nAA\n+\n!\n",
            .expected_error = error.S005LengthMismatch,
            .code = .s005_length_mismatch,
            .line = 4,
            .offset = 8,
        },
        .{
            .data = "@r\nA\n+\n!\r",
            .expected_error = error.S005LengthMismatch,
            .code = .s005_length_mismatch,
            .line = 4,
            .offset = 7,
        },
    };

    for (cases) |case| {
        for ([_]usize{ 1, 2, 3, case.data.len }) |chunk_len| {
            var chunked = ChunkSource.init(case.data, chunk_len);
            const source = chunked.byteSource();
            var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
            defer reader.deinit();
            for (0..case.record_index) |_| _ = try reader.next();
            if (reader.next()) |_| {
                return error.ExpectedParseError;
            } else |err| {
                try std.testing.expectEqual(case.expected_error, err);
            }
            const reader_details = reader.takeLastError().?;
            try expectDetails(reader_details, case.code, case.record_index, case.line, case.offset);
            try std.testing.expect(reader.takeLastError() == null);

            var scan = zfastq.count_scan.Scanner.init(.{});
            var pos: usize = 0;
            var scan_error: ?zfastq.ReaderError = null;
            while (pos < case.data.len) {
                const end = @min(case.data.len, pos + chunk_len);
                if (scan.feed(case.data[pos..end])) |_| {
                    pos = end;
                } else |err| {
                    scan_error = err;
                    break;
                }
            }
            if (scan_error == null) {
                scan.finishEof() catch |err| {
                    scan_error = err;
                };
            }
            try std.testing.expectEqual(case.expected_error, scan_error.?);
            const scan_details = scan.takeLastError().?;
            try expectDetails(scan_details, case.code, case.record_index, case.line, case.offset);
            try std.testing.expect(scan.takeLastError() == null);
        }
    }
}

test "[property] - [count scanner]: uniform and varying record strides count every record" {
    const hdr_a = "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/1\n";
    const hdr_b = "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/12\n";
    const cases = [_]struct {
        headers: []const []const u8,
        sequence_lens: []const usize,
        count: u64,
    }{
        .{ .headers = &.{ hdr_a, hdr_a, hdr_a }, .sequence_lens = &.{ 100, 100, 100 }, .count = 3 },
        .{ .headers = &.{ hdr_a, hdr_b }, .sequence_lens = &.{ 100, 100 }, .count = 2 },
        .{ .headers = &.{ "@ont_read\n", "@ont_read\n" }, .sequence_lens = &.{ 4, 6 }, .count = 2 },
    };

    var buf: [3 * 256]u8 = undefined;
    for (cases) |case| {
        var pos: usize = 0;
        for (case.headers, case.sequence_lens) |header, sequence_len| {
            pos = appendSizedRecord(&buf, pos, header, sequence_len);
        }
        const data = buf[0..pos];

        var scan = zfastq.count_scan.Scanner.init(.{});
        try std.testing.expectEqual(case.count, try zfastq.count_scan.countSlice(data, .{}, &scan));
        try std.testing.expectEqual(case.count, (try readerOutcome(data, data.len)).count);
    }
}

test "[property] - [count scanner]: dense exits match Reader at tail boundaries and splits" {
    const prefix = "@r1\nAAAA\n+\n!!!!\n" ** 3;
    const tails = [_][]const u8{
        "@longer\nAAAA\n+\n!!!!\n",
        "@r2\nAA\n+\n!!\n",
        "@r2\nAAAA\n-\n!!!!\n",
        "@r2\nAAA\r\n+\n!!!!\n",
    };
    for (tails) |tail| {
        var storage: [256]u8 = undefined;
        var output = std.Io.Writer.fixed(&storage);
        try output.writeAll(prefix);
        try output.writeAll(tail);
        try output.writeAll(tail);
        try output.writeAll(prefix);
        const complete = output.buffered();

        const seams = lineBoundaries(complete);
        for (seams.slice()) |length| {
            if (length < prefix.len) continue;
            const data = complete[0..length];
            try expectScanAtEverySplit(data, seams.slice(), try readerOutcome(data, 1));
        }
    }
}

test "[property] - [count scanner]: dense headers retain an identifier after CRLF normalization" {
    const prefix = "@a\nA\n+\n!\n" ** 2;
    const cases = [_]struct { header: []const u8, valid: bool }{
        .{ .header = "@\r", .valid = false },
        .{ .header = "@", .valid = false },
        .{ .header = "@ ", .valid = false },
        .{ .header = "@\t", .valid = false },
        .{ .header = "@b", .valid = true },
        .{ .header = "@b\r", .valid = true },
        .{ .header = "@\r\r", .valid = true },
        .{ .header = "@\x00", .valid = true },
    };

    for (cases) |case| {
        var storage: [32]u8 = undefined;
        var output = std.Io.Writer.fixed(&storage);
        try output.writeAll(prefix);
        try output.writeAll(case.header);
        try output.writeAll("\nA\n+\n!\n");
        const data = output.buffered();
        const expected = try readerOutcome(data, 1);

        try std.testing.expectEqual(@as(u64, if (case.valid) 3 else 2), expected.count);
        try std.testing.expectEqual(
            @as(?zfastq.ReaderError, if (case.valid) null else error.S003InvalidHeader),
            expected.err,
        );
        if (!case.valid) try expectDetails(expected.details.?, .s003_invalid_header, 2, 1, prefix.len);
        for ([_]usize{ 1, 2, data.len }) |chunk_len| {
            try expectSameOutcome(expected, scannerOutcome(data, chunk_len));
        }
    }
}

test "[property] - [count scanner]: mixed line endings match Reader at line seams" {
    const cases = [_][4][]const u8{
        .{ "@r", "AC", "+note", "!!" },
        .{ "@r", "", "+", "" },
        .{ "@r", "A\r", "+", "!\r" },
        .{ "@r", "A\r", "+", "!" },
        .{ "@", "A", "+", "!" },
        .{ "@ r", "A", "+", "!" },
        .{ "@r", "A", "-", "!" },
        .{ "@r", "AA", "+", "!" },
    };
    const prefix = "@a\r\nA\r\n+\r\n!\r\n";
    const ending_masks = [_]usize{ 0, 15 };
    for (cases) |fields| {
        for (ending_masks) |endings| {
            var storage: [128]u8 = undefined;
            var output = std.Io.Writer.fixed(&storage);
            try output.writeAll(prefix);
            for (fields, 0..) |field, index| {
                try output.writeAll(field);
                if (endings & (@as(usize, 1) << @intCast(index)) != 0) try output.writeByte('\r');
                try output.writeByte('\n');
            }
            try output.writeAll("@b\nA\n+\n!\n@c\nA\n+\n!\n@\nA\n+\n!\n");
            const complete = output.buffered();
            const seams = lineBoundaries(complete);
            for (seams.slice()) |length| {
                const data = complete[0..length];
                try expectScanAtEverySplit(data, seams.slice(), try readerOutcome(data, 1));
            }
        }
    }
}

test "[property] - [parser]: header lookup boundaries are chunk invariant" {
    var data: [420]u8 = undefined;
    var pos: usize = 0;
    for ([_]usize{ 127, 128, 129 }) |header_line_bytes| {
        data[pos] = '@';
        @memset(data[pos + 1 .. pos + header_line_bytes - 1], 'H');
        data[pos + header_line_bytes - 1] = '\n';
        pos += header_line_bytes;
        @memcpy(data[pos .. pos + 12], "AAAA\n+\nIIII\n");
        pos += 12;
    }

    for ([_]usize{ 1, data.len }) |chunk_len| {
        var scan = zfastq.count_scan.Scanner.init(.{});
        var scan_pos: usize = 0;
        while (scan_pos < data.len) {
            const end = @min(data.len, scan_pos + chunk_len);
            _ = try scan.feed(data[scan_pos..end]);
            scan_pos = end;
        }
        try scan.finishEof();

        var chunked = ChunkSource.init(&data, chunk_len);
        var reader = try zfastq.Reader.init(
            std.testing.allocator,
            chunked.byteSource(),
            .{},
        );
        defer reader.deinit();
        var reader_count: u64 = 0;
        while (try reader.advance()) reader_count += 1;

        try std.testing.expectEqual(reader_count, scan.record_index);
        try std.testing.expectEqual(@as(u64, 3), scan.record_index);
    }
}

fn appendSizedRecord(buf: []u8, start: usize, header: []const u8, seq_len: usize) usize {
    var pos = start;
    @memcpy(buf[pos..][0..header.len], header);
    pos += header.len;
    @memset(buf[pos..][0..seq_len], 'A');
    pos += seq_len;
    @memcpy(buf[pos..][0..3], "\n+\n");
    pos += 3;
    @memset(buf[pos..][0..seq_len], 'I');
    pos += seq_len;
    buf[pos] = '\n';
    return pos + 1;
}

const ChunkSource = struct {
    data: []const u8,
    pos: usize = 0,
    chunk_len: usize,
    split_mask: ?usize = null,
    fail_at: ?usize = null,

    fn init(data: []const u8, chunk_len: usize) ChunkSource {
        return .{ .data = data, .chunk_len = chunk_len };
    }

    fn byteSource(self: *ChunkSource) zfastq.io.ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *ChunkSource = @ptrCast(@alignCast(ctx));
        if (self.fail_at) |fail_at| {
            if (self.pos >= fail_at) return error.ReadFailed;
        }
        if (self.pos == self.data.len) return 0;
        const limit = if (self.fail_at) |fail_at| @min(fail_at, self.data.len) else self.data.len;
        const end = if (self.split_mask) |mask|
            partitionEnd(self.data.len, self.pos, mask)
        else
            self.pos + self.chunk_len;
        const n = @min(@min(end, limit) - self.pos, dest.len);
        @memcpy(dest[0..n], self.data[self.pos..][0..n]);
        self.pos += n;
        return n;
    }
};

fn decodeGzip(
    compressed: []const u8,
    input_chunk_len: usize,
    output_chunk_len: usize,
    output: []u8,
) ![]const u8 {
    var input_buffer: [64]u8 = undefined;
    var input = FragmentReader.init(compressed, input_chunk_len, &input_buffer);
    var gzip = zfastq.io.gzip.ReaderSource.init(&input.interface);
    const source = gzip.byteSource();
    var written: usize = 0;
    while (written < output.len) {
        const end = @min(output.len, written + output_chunk_len);
        const n = try source.read(output[written..end]);
        if (n == 0) return output[0..written];
        written += n;
    }
    return error.TestOutputLimit;
}

const FragmentReader = struct {
    interface: std.Io.Reader,
    data: []const u8,
    pos: usize = 0,
    chunk_len: usize,
    fail_at: ?usize = null,

    fn init(
        data: []const u8,
        chunk_len: usize,
        buffer: []u8,
    ) FragmentReader {
        return .{
            .interface = .{
                .vtable = &.{ .stream = stream },
                .buffer = buffer,
                .seek = 0,
                .end = 0,
            },
            .data = data,
            .chunk_len = chunk_len,
        };
    }

    fn stream(
        reader: *std.Io.Reader,
        writer: *std.Io.Writer,
        limit: std.Io.Limit,
    ) std.Io.Reader.StreamError!usize {
        const self: *FragmentReader = @alignCast(@fieldParentPtr("interface", reader));
        if (self.pos == self.data.len) return error.EndOfStream;
        const before_failure = if (self.fail_at) |fail_at| b: {
            if (self.pos >= fail_at) return error.ReadFailed;
            break :b fail_at - self.pos;
        } else self.data.len - self.pos;
        const n = @min(
            @intFromEnum(limit),
            @min(self.chunk_len, @min(before_failure, self.data.len - self.pos)),
        );
        const written = try writer.write(self.data[self.pos..][0..n]);
        self.pos += written;
        return written;
    }
};

const InvalidCountSource = struct {
    read_count: usize = 0,

    fn byteSource(self: *InvalidCountSource) zfastq.io.ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *InvalidCountSource = @ptrCast(@alignCast(ctx));
        self.read_count += 1;
        return dest.len + 1;
    }
};

fn partitionEnd(len: usize, start: usize, mask: usize) usize {
    var end = start + 1;
    while (end < len and mask & (@as(usize, 1) << @intCast(end - 1)) == 0) : (end += 1) {}
    return end;
}

const MAX_LINE_FEEDS = 32;

const LineBoundaries = struct {
    items: [3 * MAX_LINE_FEEDS + 1]usize = undefined,
    len: usize = 0,

    fn slice(self: *const LineBoundaries) []const usize {
        return self.items[0..self.len];
    }
};

fn lineBoundaries(data: []const u8) LineBoundaries {
    var result: LineBoundaries = .{};
    for (data, 0..) |byte, index| {
        if (byte != '\n') continue;
        for ([_]usize{ index, index + 1, index + 2 }) |position| {
            result.items[result.len] = @min(position, data.len);
            result.len += 1;
        }
    }
    result.items[result.len] = data.len;
    result.len += 1;
    return result;
}

fn validateSequenceQuality(
    sequence: []const u8,
    quality: []const u8,
    options: zfastq.ValidationOptions,
) ?zfastq.SemanticError {
    return zfastq.validateRecord(.{
        .header = "r",
        .id = "r",
        .sequence = sequence,
        .plus = "",
        .quality = quality,
    }, options);
}

fn expectDetails(
    details: zfastq.ParseError,
    code: zfastq.LintCode,
    record_index: u64,
    line: u3,
    offset: u64,
) !void {
    try std.testing.expectEqual(code, details.code);
    try std.testing.expectEqual(record_index, details.record_index);
    try std.testing.expectEqual(line, details.line_in_record);
    try std.testing.expectEqual(offset, details.byte_offset);
}

const ParseOutcome = struct {
    count: u64,
    err: ?zfastq.ReaderError = null,
    details: ?zfastq.ParseError = null,
};

fn readerOutcome(data: []const u8, chunk_len: usize) !ParseOutcome {
    var chunked = ChunkSource.init(data, chunk_len);
    var reader = try zfastq.Reader.init(
        std.testing.allocator,
        chunked.byteSource(),
        .{},
    );
    defer reader.deinit();

    while (true) {
        const has_record = reader.advance() catch |err| {
            return .{
                .count = reader.recordIndex(),
                .err = err,
                .details = reader.takeLastError(),
            };
        };
        if (!has_record) return .{ .count = reader.recordIndex() };
    }
}

fn scannerOutcome(data: []const u8, chunk_len: usize) ParseOutcome {
    var scan = zfastq.count_scan.Scanner.init(.{});
    var pos: usize = 0;
    while (pos < data.len) {
        const end = @min(data.len, pos + chunk_len);
        _ = scan.feed(data[pos..end]) catch |err| {
            return .{
                .count = scan.record_index,
                .err = err,
                .details = scan.takeLastError(),
            };
        };
        pos = end;
    }
    scan.finishEof() catch |err| {
        return .{
            .count = scan.record_index,
            .err = err,
            .details = scan.takeLastError(),
        };
    };
    return .{ .count = scan.record_index };
}

fn expectScanAtEverySplit(data: []const u8, splits: []const usize, expected: ParseOutcome) !void {
    for (splits) |split| {
        if (split > data.len) continue;
        var scan = zfastq.count_scan.Scanner.init(.{});
        var failure: ?zfastq.ReaderError = null;
        for ([_][]const u8{ data[0..split], data[split..] }) |chunk| {
            const consumed = scan.feed(chunk) catch |err| {
                failure = err;
                break;
            };
            try std.testing.expectEqual(chunk.len, consumed);
        }
        if (failure == null) scan.finishEof() catch |err| {
            failure = err;
        };
        try expectSameOutcome(expected, .{
            .count = scan.record_index,
            .err = failure,
            .details = scan.takeLastError(),
        });
        if (failure == null) {
            try std.testing.expectEqual(@as(u64, data.len), scan.byte_offset);
        } else if (expected.details) |details| {
            const line_start: usize = @intCast(details.byte_offset);
            const line_end = std.mem.findScalarPos(u8, data, line_start, '\n');
            const consumed = if (line_end) |end| end + 1 else data.len;
            try std.testing.expectEqual(@as(u64, consumed), scan.byte_offset);
        }
    }
}

fn fuzzReaderScannerAgreement(_: void, smith: *std.testing.Smith) !void {
    var storage: [512]u8 = undefined;
    const data = storage[0..smith.slice(&storage)];
    const reader = try readerOutcome(data, 1);
    try expectOutcomeBounds(data, reader);

    const scanner_one = scannerOutcome(data, 1);
    const scanner_block = scannerOutcome(data, 64);
    try expectOutcomeBounds(data, scanner_one);
    try expectOutcomeBounds(data, scanner_block);
    try expectSameOutcome(reader, scanner_one);
    try expectSameOutcome(reader, scanner_block);
}

fn expectOutcomeBounds(data: []const u8, outcome: ParseOutcome) !void {
    const line_feeds = std.mem.count(u8, data, "\n");
    try std.testing.expect(outcome.count <= (line_feeds + 1) / 4);
    try std.testing.expectEqual(outcome.err != null, outcome.details != null);
    if (outcome.details) |details| {
        try std.testing.expectEqual(outcome.count, details.record_index);
        try std.testing.expect(details.line_in_record >= 1);
        try std.testing.expect(details.line_in_record <= 4);
        try std.testing.expect(details.byte_offset <= data.len);
    }
}

fn expectSameOutcome(expected: ParseOutcome, actual: ParseOutcome) !void {
    try std.testing.expectEqual(expected.count, actual.count);
    try std.testing.expectEqual(expected.err, actual.err);
    try std.testing.expectEqual(expected.details == null, actual.details == null);
    if (expected.details) |expected_details| {
        const actual_details = actual.details.?;
        try expectDetails(
            actual_details,
            expected_details.code,
            expected_details.record_index,
            expected_details.line_in_record,
            expected_details.byte_offset,
        );
        try std.testing.expectEqualStrings(expected_details.message, actual_details.message);
    }
}
