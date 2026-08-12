//! Library root, reader, scanner, and validation contracts.

const std = @import("std");
const zfastq = @import("z-fastq");

const GZIP_OPTIONAL_X = [_]u8{
    0x1f, 0x8b, 0x08, 0x1e, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x02, 0x00, 'x',  'y',  'n',  0x00, 'c',  0x00, 0xca, 0x4e,
    0x01, 0x01, 0x00, 0xfe, 0xff, 'x',  0x83, 0x16, 0xdc, 0x8c,
    0x01, 0x00, 0x00, 0x00,
};

test "[unit] - [root]: every exported declaration is analyzable" {
    std.testing.refAllDecls(zfastq);
}

test "[unit] - [root]: version exposes the current internal checkpoint" {
    try std.testing.expectEqualStrings("0.0.8", zfastq.VERSION);
}

test "[property] - [gzip source]: optional member chains decode at every input chunk size" {
    const gzip_xx = GZIP_OPTIONAL_X ++ GZIP_OPTIONAL_X;

    for (1..gzip_xx.len + 1) |chunk_len| {
        var input_buffer: [10]u8 = undefined;
        var input = FragmentReader.init(&gzip_xx, chunk_len, &input_buffer);
        try std.testing.expectEqualSlices(u8, gzip_xx[0..2], try input.interface.peek(2));
        var gzip = zfastq.io.gzip.ReaderSource.init(&input.interface);
        const source = gzip.byteSource();
        var output: [2]u8 = undefined;

        try std.testing.expectEqual(@as(usize, 2), try source.read(&output));
        try std.testing.expectEqualStrings("xx", &output);
        try std.testing.expectEqual(@as(usize, 0), try source.read(&output));
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

test "[property] - [input detection]: a one-byte first read replays plain prefix bytes" {
    const data = "@r\nA\n+\n!\n";
    var input_buffer: [2]u8 = undefined;
    var input = FragmentReader.init(data, 1, &input_buffer);

    try std.testing.expectEqualSlices(u8, data[0..2], try input.interface.peek(2));
    var plain = zfastq.io.plain.ReaderSource.init(&input.interface);
    const source = plain.byteSource();
    var output: [data.len]u8 = undefined;

    try std.testing.expectEqual(data.len, try source.read(&output));
    try std.testing.expectEqualStrings(data, &output);
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

test "[unit] - [record validation]: both alphabets accept their complete symbol sets" {
    const cases = [_]struct {
        options: zfastq.ValidationOptions,
        sequence: []const u8,
    }{
        .{
            .options = .{},
            .sequence = "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
        },
        .{
            .options = .{ .alphabet = .acgtn },
            .sequence = "ACGTNacgtn",
        },
        .{ .options = .{ .alphabet = .iupac }, .sequence = "" },
        .{ .options = .{ .alphabet = .acgtn }, .sequence = "" },
    };

    for (cases) |case| {
        var quality_buffer: [32]u8 = undefined;
        const quality = quality_buffer[0..case.sequence.len];
        @memset(quality, '!');
        const record = zfastq.Record{
            .header = "r",
            .id = "r",
            .sequence = case.sequence,
            .plus = "",
            .quality = quality,
        };

        try std.testing.expect(zfastq.validateRecord(record, case.options) == null);
    }
}

test "[property] - [record validation]: ACGTN rejects every wider IUPAC symbol" {
    for ("URYSWKMBDHVuryswkmbdhv") |byte| {
        const sequence = [1]u8{byte};
        const details = zfastq.validateRecord(.{
            .header = "r",
            .id = "r",
            .sequence = &sequence,
            .plus = "",
            .quality = "!",
        }, .{ .alphabet = .acgtn }).?;

        try std.testing.expectEqual(zfastq.LintCode.s002_invalid_sequence_alphabet, details.code);
        try std.testing.expectEqual(zfastq.SemanticField.sequence, details.field);
        try std.testing.expectEqual(@as(usize, 0), details.byte_index);
        try std.testing.expectEqualStrings(
            "sequence byte is outside the selected alphabet",
            details.message,
        );
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
        for (0..256) |value| {
            const sequence = [1]u8{@intCast(value)};
            const details = zfastq.validateRecord(.{
                .header = "r",
                .id = "r",
                .sequence = &sequence,
                .plus = "",
                .quality = "!",
            }, .{ .alphabet = policy.alphabet });
            const accepted = std.mem.indexOfScalar(u8, policy.accepted, sequence[0]) != null;
            try std.testing.expectEqual(accepted, details == null);
            if (details) |semantic_error| {
                try std.testing.expectEqual(
                    zfastq.LintCode.s002_invalid_sequence_alphabet,
                    semantic_error.code,
                );
                try std.testing.expectEqual(@as(usize, 0), semantic_error.byte_index);
            }
        }
    }
}

test "[property] - [record validation]: every semantic byte position reports exactly" {
    var sequence = [_]u8{ 'A', 'C', 'G', 'T' };
    const quality = [_]u8{ '!', '!', '!', '!' };
    for (sequence, 0..) |original, invalid_index| {
        sequence[invalid_index] = '.';
        const details = zfastq.validateRecord(.{
            .header = "r",
            .id = "r",
            .sequence = &sequence,
            .plus = "",
            .quality = &quality,
        }, .{}).?;
        try std.testing.expectEqual(zfastq.LintCode.s002_invalid_sequence_alphabet, details.code);
        try std.testing.expectEqual(zfastq.SemanticField.sequence, details.field);
        try std.testing.expectEqual(invalid_index, details.byte_index);
        sequence[invalid_index] = original;
    }

    var mutable_quality = quality;
    for (mutable_quality, 0..) |original, invalid_index| {
        mutable_quality[invalid_index] = if (invalid_index % 2 == 0) 32 else 127;
        const details = zfastq.validateRecord(.{
            .header = "r",
            .id = "r",
            .sequence = &sequence,
            .plus = "",
            .quality = &mutable_quality,
        }, .{}).?;
        try std.testing.expectEqual(zfastq.LintCode.s006_invalid_quality_range, details.code);
        try std.testing.expectEqual(zfastq.SemanticField.quality, details.field);
        try std.testing.expectEqual(invalid_index, details.byte_index);
        try std.testing.expectEqualStrings(
            "quality byte must be ASCII 33 through 126",
            details.message,
        );
        mutable_quality[invalid_index] = original;
    }
}

test "[unit] - [record validation]: sequence errors precede quality errors" {
    const details = zfastq.validateRecord(.{
        .header = "r",
        .id = "r",
        .sequence = ".",
        .plus = "",
        .quality = "\x7f",
    }, .{}).?;

    try std.testing.expectEqual(zfastq.LintCode.s002_invalid_sequence_alphabet, details.code);
}

test "[property] - [record validation]: quality range classifies every byte" {
    try std.testing.expect(zfastq.validateRecord(.{
        .header = "r",
        .id = "r",
        .sequence = "AC",
        .plus = "",
        .quality = "!~",
    }, .{}) == null);

    for (0..256) |value| {
        const quality = [1]u8{@intCast(value)};
        const details = zfastq.validateRecord(.{
            .header = "r",
            .id = "r",
            .sequence = "A",
            .plus = "",
            .quality = &quality,
        }, .{});
        try std.testing.expectEqual(value >= 33 and value <= 126, details == null);
        if (details) |semantic_error| {
            try std.testing.expectEqual(
                zfastq.LintCode.s006_invalid_quality_range,
                semantic_error.code,
            );
            try std.testing.expectEqual(@as(usize, 0), semantic_error.byte_index);
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

test "[failure] - [reader]: unequal sequence and quality lengths report S005" {
    const data =
        \\@bad
        \\ACGT
        \\+
        \\III
        \\
    ;
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.S005LengthMismatch, reader.next());
}

test "[failure] - [reader]: a header without at-sign reports S003" {
    const data = "not_a_header\nACGT\n+\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.S003InvalidHeader, reader.next());
}

test "[failure] - [reader]: a truncated record at EOF reports S004" {
    const data = "@truncated\nACGT\n+";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.S004TruncatedRecord, reader.next());
}

test "[failure] - [reader]: an invalid plus line reports S001 at its start" {
    const data = "@bad\nACGT\nnot-plus\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    try std.testing.expectError(zfastq.ReaderError.S001InvalidPlusLine, reader.next());
    const details = reader.takeLastError().?;
    try std.testing.expectEqual(zfastq.LintCode.s001_invalid_plus_line, details.code);
    try std.testing.expectEqual(@as(u64, 0), details.record_index);
    try std.testing.expectEqual(@as(u3, 3), details.line_in_record);
    try std.testing.expectEqual(@as(u64, 10), details.byte_offset);
    try std.testing.expect(reader.takeLastError() == null);
}

test "[property] - [reader]: borrowed fields survive every fixed short-read size" {
    const data = "@read1 comment\nACGT\n+repeat\n!!!!\n";
    for (1..data.len + 1) |chunk_len| {
        var chunked = ChunkSource.init(data, chunk_len);
        const source = chunked.byteSource();
        var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
        defer reader.deinit();

        const parsed = (try reader.next()).?;
        try std.testing.expectEqualStrings("read1 comment", parsed.header);
        try std.testing.expectEqualStrings("read1", parsed.id);
        try std.testing.expectEqualStrings("ACGT", parsed.sequence);
        try std.testing.expectEqualStrings("repeat", parsed.plus);
        try std.testing.expectEqualStrings("!!!!", parsed.quality);
        try std.testing.expectEqual(@as(u64, data.len), reader.byteOffset());
    }
}

test "[property] - [reader]: complete and refill-spanning records are identical" {
    const data =
        "@one comment\r\nACGT\r\n+same\r\n!!!!\r\n" ++
        "@empty\n\n+description\n\n" ++
        "@final\nNN\n+\n##";
    const cases = [_]struct {
        header: []const u8,
        id: []const u8,
        sequence: []const u8,
        plus: []const u8,
        quality: []const u8,
        offsets: zfastq.RecordOffsets,
    }{
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
    };
    const vector_len = @max(std.simd.suggestVectorLength(u8) orelse 16, 2);
    const chunk_sizes = [_]usize{
        1,
        2,
        3,
        vector_len - 1,
        vector_len,
        vector_len + 1,
        4 * 1024,
        zfastq.limits.DEFAULT_READER_BUFFER_BYTES,
    };

    for (chunk_sizes) |chunk_len| {
        var chunked = ChunkSource.init(data, chunk_len);
        var reader = try zfastq.Reader.init(
            std.testing.allocator,
            chunked.byteSource(),
            .{},
        );
        defer reader.deinit();

        for (cases) |case| {
            const parsed = (try reader.next()).?;
            try std.testing.expectEqualStrings(case.header, parsed.header);
            try std.testing.expectEqualStrings(case.id, parsed.id);
            try std.testing.expectEqualStrings(case.sequence, parsed.sequence);
            try std.testing.expectEqualStrings(case.plus, parsed.plus);
            try std.testing.expectEqualStrings(case.quality, parsed.quality);
            try std.testing.expectEqual(case.offsets, reader.currentRecordOffsets().?);
        }
        try std.testing.expect((try reader.next()) == null);
        try std.testing.expectEqual(@as(u64, data.len), reader.byteOffset());
    }
}

test "[property] - [parser]: every input partition preserves the format result" {
    const valid = "@r\nA\n+\n!\n";
    const partition_count = @as(usize, 1) << (valid.len - 1);
    for (0..partition_count) |mask| {
        var partitioned = PartitionSource.init(valid, mask);
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
    for (0..partition_count) |mask| {
        var partitioned = PartitionSource.init(invalid, mask);
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

test "[edge] - [reader]: borrowed fields survive record-buffer growth" {
    const line_len = 1024 * 1024;
    const total_len = 3 + line_len + 3 + line_len;
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

    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    const parsed = (try reader.next()).?;
    try std.testing.expectEqual(line_len, parsed.sequence.len);
    try std.testing.expectEqual(line_len, parsed.quality.len);
    try std.testing.expect(std.mem.allEqual(u8, parsed.sequence, 'A'));
    try std.testing.expect(std.mem.allEqual(u8, parsed.quality, 'I'));
}

test "[unit] - [reader]: a common record allocates only during initialization" {
    const data = "@r\nACGT\n+\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
    var reader = try zfastq.Reader.init(failing.allocator(), source.byteSource(), .{});
    defer reader.deinit();

    _ = try reader.next();
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
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        initReaderForAllocationCheck,
        .{"@r\nA\n+\n!\n"},
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
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        ownRecordForAllocationCheck,
        .{record},
    );
}

test "[failure] - [reader]: source read failure propagates" {
    var failing_source = FailingSource{};
    const source = failing_source.byteSource();
    var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
    defer reader.deinit();

    try std.testing.expectError(zfastq.ReaderError.Io, reader.next());
}

test "[failure] - [reader]: a source count beyond its destination is rejected" {
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

test "[edge] - [reader]: the configured line limit is exact" {
    const at_limit = "@r\nACGT\n+\n!!!!\n";
    var valid_source = zfastq.io.plain.SliceSource.init(at_limit);
    var valid_reader = try zfastq.Reader.init(
        std.testing.allocator,
        valid_source.byteSource(),
        .{ .max_line_bytes = 4 },
    );
    defer valid_reader.deinit();
    _ = try valid_reader.next();

    const above_limit = "@r\nACGTA\n+\n!!!!!\n";
    var invalid_source = zfastq.io.plain.SliceSource.init(above_limit);
    var invalid_reader = try zfastq.Reader.init(
        std.testing.allocator,
        invalid_source.byteSource(),
        .{ .max_line_bytes = 4 },
    );
    defer invalid_reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.LineTooLong, invalid_reader.next());

    const lone_cr_at_eof = "@r\nA\n+\n!\r";
    var eof_source = zfastq.io.plain.SliceSource.init(lone_cr_at_eof);
    var eof_reader = try zfastq.Reader.init(
        std.testing.allocator,
        eof_source.byteSource(),
        .{ .max_line_bytes = 1 },
    );
    defer eof_reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.LineTooLong, eof_reader.next());

    var scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 1 });
    try std.testing.expectError(
        zfastq.ReaderError.LineTooLong,
        zfastq.count_scan.countSlice(lone_cr_at_eof, .{ .max_line_bytes = 1 }, &scan),
    );
}

test "[property] - [count scanner]: inline counts match Reader advance" {
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

    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    var advance_count: u64 = 0;
    while (try reader.advance()) {
        advance_count += 1;
    }
    try std.testing.expectEqual(scan.record_index, advance_count);
}

test "[failure] - [count scanner]: a length mismatch reports S005" {
    const data =
        \\@bad
        \\ACGT
        \\+
        \\III
        \\
    ;
    var scan = zfastq.count_scan.Scanner.init(.{});
    try std.testing.expectError(zfastq.ReaderError.S005LengthMismatch, scan.feed(data));
}

test "[property] - [count scanner]: result and details are chunk invariant" {
    const valid = "@r1\nAAAA\n+\n!!!!\n@r2\nCC\n+note\n##\n";
    for (1..valid.len + 1) |chunk_len| {
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

    const malformed = "@r1\nAAAA\n+\n!!!!\n@r2\nA\nAA\n+\n!!!!\n";
    for (1..malformed.len + 1) |chunk_len| {
        var scan = zfastq.count_scan.Scanner.init(.{});
        var pos: usize = 0;
        var found_error = false;
        while (pos < malformed.len) {
            const end = @min(malformed.len, pos + chunk_len);
            if (scan.feed(malformed[pos..end])) |_| {
                pos = end;
            } else |err| {
                try std.testing.expectEqual(zfastq.ReaderError.S001InvalidPlusLine, err);
                found_error = true;
                break;
            }
        }
        try std.testing.expect(found_error);
        const details = scan.takeLastError().?;
        try std.testing.expectEqual(zfastq.LintCode.s001_invalid_plus_line, details.code);
        try std.testing.expectEqual(@as(u64, 1), details.record_index);
        try std.testing.expectEqual(@as(u64, 22), details.byte_offset);
        try std.testing.expectEqual(@as(u3, 3), details.line_in_record);
    }
}

test "[fuzz] - [parser]: generated mutations keep Reader and scanner in agreement" {
    var prng = std.Random.DefaultPrng.init(0x6a09e667f3bcc909);
    const random = prng.random();
    const mutations = [_]u8{ '\n', '\r', '@', '+', 0, 'A', '!' };

    for (0..96) |_| {
        var storage: [8192]u8 = undefined;
        var output = std.Io.Writer.fixed(&storage);
        const record_count = random.intRangeAtMost(usize, 1, 16);
        for (0..record_count) |record_index| {
            try output.writeAll("@r");
            const header_len = random.intRangeAtMost(usize, 1, 24);
            for (0..header_len) |_| try output.writeByte(randomLetter(random));
            try writeRandomLineEnding(&output, random);

            const sequence_len = random.intRangeAtMost(usize, 0, 64);
            for (0..sequence_len) |_| try output.writeByte(randomBase(random));
            try writeRandomLineEnding(&output, random);

            try output.writeByte('+');
            const plus_len = random.intRangeAtMost(usize, 0, 16);
            for (0..plus_len) |_| try output.writeByte(randomLetter(random));
            try writeRandomLineEnding(&output, random);

            for (0..sequence_len) |_| try output.writeByte(randomQuality(random));
            if (record_index + 1 < record_count or sequence_len == 0 or random.boolean()) {
                try writeRandomLineEnding(&output, random);
            }
        }
        const data = output.buffered();
        const reader_chunk = random.intRangeAtMost(usize, 1, 67);
        const scanner_chunk = random.intRangeAtMost(usize, 1, 67);

        const reader_valid = try readerOutcome(data, reader_chunk);
        const scanner_valid = scannerOutcome(data, scanner_chunk);

        try std.testing.expectEqual(@as(u64, @intCast(record_count)), reader_valid.count);
        try std.testing.expectEqual(@as(?zfastq.ReaderError, null), reader_valid.err);
        try expectSameOutcome(reader_valid, scanner_valid);

        var mutated: [8192]u8 = undefined;
        @memcpy(mutated[0..data.len], data);
        const mutation_count = random.intRangeAtMost(usize, 1, 3);
        for (0..mutation_count) |_| {
            const index = random.uintLessThan(usize, data.len);
            mutated[index] = mutations[random.uintLessThan(usize, mutations.len)];
        }

        const reader_mutated = try readerOutcome(mutated[0..data.len], reader_chunk);
        const scanner_mutated = scannerOutcome(mutated[0..data.len], scanner_chunk);

        try expectSameOutcome(reader_mutated, scanner_mutated);
    }
}

test "[property] - [parser]: structural details agree across chunk sizes" {
    const cases = [_]struct {
        data: []const u8,
        expected_error: zfastq.ReaderError,
        code: zfastq.LintCode,
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
        for (1..case.data.len + 1) |chunk_len| {
            var chunked = ChunkSource.init(case.data, chunk_len);
            const source = chunked.byteSource();
            var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
            defer reader.deinit();
            if (reader.next()) |_| {
                return error.ExpectedParseError;
            } else |err| {
                try std.testing.expectEqual(case.expected_error, err);
            }
            const reader_details = reader.takeLastError().?;
            try expectDetails(reader_details, case.code, 0, case.line, case.offset);
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
            try expectDetails(scan_details, case.code, 0, case.line, case.offset);
            try std.testing.expect(scan.takeLastError() == null);
        }
    }
}

test "[property] - [count scanner]: dense validation checks every field region" {
    const valid = "@r1\nAAAA\n+\n!!!!\n";
    const cases = [_]struct {
        malformed_record: []const u8,
        expected_error: zfastq.ReaderError,
        code: zfastq.LintCode,
        line: u3,
        offset: u64,
    }{
        .{
            .malformed_record = "@\nx\nAAAA\n+\n!!!!\n",
            .expected_error = error.S001InvalidPlusLine,
            .code = .s001_invalid_plus_line,
            .line = 3,
            .offset = 20,
        },
        .{
            .malformed_record = "@r2\nA\nAA\n+\n!!!!\n",
            .expected_error = error.S001InvalidPlusLine,
            .code = .s001_invalid_plus_line,
            .line = 3,
            .offset = 22,
        },
        .{
            .malformed_record = "@r2\nAAAA\n+\n!\n!!\n",
            .expected_error = error.S005LengthMismatch,
            .code = .s005_length_mismatch,
            .line = 4,
            .offset = 27,
        },
    };

    for (cases) |case| {
        var data: [32]u8 = undefined;
        @memcpy(data[0..valid.len], valid);
        @memcpy(data[valid.len..], case.malformed_record);

        for (1..data.len + 1) |chunk_len| {
            var scan = zfastq.count_scan.Scanner.init(.{});
            var pos: usize = 0;
            var found_error: ?zfastq.ReaderError = null;
            while (pos < data.len) {
                const end = @min(data.len, pos + chunk_len);
                if (scan.feed(data[pos..end])) |_| {
                    pos = end;
                } else |err| {
                    found_error = err;
                    break;
                }
            }
            try std.testing.expectEqual(case.expected_error, found_error.?);
            const details = scan.takeLastError().?;
            try expectDetails(details, case.code, 1, case.line, case.offset);
        }
    }
}

test "[edge] - [count scanner]: lines span chunks while retaining the limit" {
    const data = "@r\nAAAAAAAA\n+\n!!!!!!!!";
    var scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 8 });
    for (data) |byte| {
        _ = try scan.feed((&[_]u8{byte})[0..]);
    }
    try scan.finishEof();
    try std.testing.expectEqual(@as(u64, 1), scan.record_index);

    var too_long = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 7 });
    try std.testing.expectError(
        zfastq.ReaderError.LineTooLong,
        too_long.feed(data),
    );
}

test "[edge] - [count scanner]: dense layout learning retains the line limit" {
    const data = "@a\nAC\n+\n!!\n@b\nACGT\n+\n!!!!\n";
    var scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 3 });

    try std.testing.expectError(
        zfastq.ReaderError.LineTooLong,
        zfastq.count_scan.countSlice(data, .{ .max_line_bytes = 3 }, &scan),
    );
    try std.testing.expectEqual(@as(u64, 1), scan.record_index);
}

test "[property] - [count scanner]: dense layout learning preserves length semantics" {
    const data = "@a\nACGT\n+\n!!!!\n@b\nACG\r\n+\n!!!!\n";
    var scan = zfastq.count_scan.Scanner.init(.{});

    try std.testing.expectError(
        zfastq.ReaderError.S005LengthMismatch,
        zfastq.count_scan.countSlice(data, .{}, &scan),
    );
    try expectDetails(
        scan.takeLastError().?,
        .s005_length_mismatch,
        1,
        4,
        25,
    );
}

test "[unit] - [count scanner]: dense stride blocks count uniform records" {
    var buf: [243 * 3]u8 = undefined;
    const pos = appendUniformRecords(&buf, 3, "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/1\n", 100);
    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 3), n);
}
test "[property] - [count scanner]: known stride layouts may alternate" {
    const hdr_a = "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/1\n";
    const hdr_b = "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/12\n";

    var buf: [243 + 244]u8 = undefined;
    var pos: usize = 0;
    pos = appendRecord(&buf, pos, hdr_a, 100);
    pos = appendRecord(&buf, pos, hdr_b, 100);

    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 2), n);
}

test "[property] - [count scanner]: sequence lengths may alternate under one header" {
    const hdr = "@ont_read\n";

    var buf: [128]u8 = undefined;
    var pos: usize = 0;
    pos = appendRecord(&buf, pos, hdr, 4);
    pos = appendRecord(&buf, pos, hdr, 6);

    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 2), n);
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

    for (1..data.len + 1) |chunk_len| {
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

test "[edge] - [count scanner]: a non-minimal plus record does not disable fast scanning" {
    const data =
        \\@read3
        \\AAAA
        \\+repeat_id
        \\!!!!
        \\@dense
        \\CCCC
        \\+
        \\####
        \\
    ;
    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(data, .{}, &scan);
    try std.testing.expectEqual(@as(u64, 2), n);
}

test "[edge] - [count scanner]: a non-minimal plus line uses structural fallback" {
    const data =
        \\@read3
        \\AAAA
        \\+repeat_id
        \\!!!!
        \\
    ;
    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(data, .{}, &scan);
    try std.testing.expectEqual(@as(u64, 1), n);
}

test "[property] - [reader]: advance and next produce the same record count" {
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
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    var advance_count: usize = 0;
    while (try reader.advance()) {
        advance_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 2), advance_count);

    var source2 = zfastq.io.plain.SliceSource.init(data);
    var reader2 = try zfastq.Reader.init(std.testing.allocator, source2.byteSource(), .{});
    defer reader2.deinit();

    var next_count: usize = 0;
    while ((try reader2.next()) != null) {
        next_count += 1;
    }
    try std.testing.expectEqual(advance_count, next_count);
}

fn appendUniformRecords(buf: []u8, count: usize, header: []const u8, seq_len: usize) usize {
    var pos: usize = 0;
    for (0..count) |_| {
        pos = appendRecord(buf, pos, header, seq_len);
    }
    return pos;
}

fn appendRecord(buf: []u8, start: usize, header: []const u8, seq_len: usize) usize {
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

    fn init(data: []const u8, chunk_len: usize) ChunkSource {
        return .{ .data = data, .chunk_len = chunk_len };
    }

    fn byteSource(self: *ChunkSource) zfastq.io.ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *ChunkSource = @ptrCast(@alignCast(ctx));
        if (self.pos == self.data.len) return 0;
        const n = @min(self.chunk_len, @min(dest.len, self.data.len - self.pos));
        @memcpy(dest[0..n], self.data[self.pos..][0..n]);
        self.pos += n;
        return n;
    }
};

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

const FailingSource = struct {
    fn byteSource(self: *FailingSource) zfastq.io.ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        _ = ctx;
        _ = dest;
        return error.ReadFailed;
    }
};

const InvalidCountSource = struct {
    fn byteSource(self: *InvalidCountSource) zfastq.io.ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        _ = ctx;
        return dest.len + 1;
    }
};

const PartitionSource = struct {
    data: []const u8,
    pos: usize = 0,
    mask: usize,

    fn init(data: []const u8, mask: usize) PartitionSource {
        return .{ .data = data, .mask = mask };
    }

    fn byteSource(self: *PartitionSource) zfastq.io.ByteSource {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *PartitionSource = @ptrCast(@alignCast(ctx));
        if (self.pos == self.data.len) return 0;
        const end = @min(partitionEnd(self.data.len, self.pos, self.mask), self.pos + dest.len);
        const n = end - self.pos;
        @memcpy(dest[0..n], self.data[self.pos..end]);
        self.pos = end;
        return n;
    }
};

fn partitionEnd(len: usize, start: usize, mask: usize) usize {
    var end = start + 1;
    while (end < len and mask & (@as(usize, 1) << @intCast(end - 1)) == 0) : (end += 1) {}
    return end;
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

fn initReaderForAllocationCheck(allocator: std.mem.Allocator, data: []const u8) !void {
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(allocator, source.byteSource(), .{});
    defer reader.deinit();
}

fn ownRecordForAllocationCheck(allocator: std.mem.Allocator, record: zfastq.Record) !void {
    var owned = try zfastq.toOwned(allocator, record);
    defer owned.deinit();
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

fn writeRandomLineEnding(output: *std.Io.Writer, random: std.Random) !void {
    if (random.boolean()) {
        try output.writeAll("\r\n");
    } else {
        try output.writeByte('\n');
    }
}

fn randomBase(random: std.Random) u8 {
    return "ACGTN"[random.uintLessThan(usize, 5)];
}

fn randomLetter(random: std.Random) u8 {
    return 'a' + random.uintLessThan(u8, 26);
}

fn randomQuality(random: std.Random) u8 {
    return '!' + random.uintLessThan(u8, 41);
}
