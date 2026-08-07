//! Unit tests for `fastq.Reader`.

const std = @import("std");
const zfastq = @import("z-fastq");

test "reader: parses basic four-line record" {
    const data =
        \\@read1
        \\ACGT
        \\+
        \\!!!!
        \\
    ;
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    const record = (try reader.next()).?;
    try std.testing.expectEqualStrings("read1", record.header);
    try std.testing.expectEqualStrings("read1", record.id);
    try std.testing.expectEqualStrings("ACGT", record.sequence);
    try std.testing.expectEqualStrings("", record.plus);
    try std.testing.expectEqualStrings("!!!!", record.quality);
    try std.testing.expect((try reader.next()) == null);
}

test "reader: accepts CRLF line endings" {
    const data = "@r1\r\nACGT\r\n+\r\n!!!!\r\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    _ = try reader.next();
    try std.testing.expect((try reader.next()) == null);
}

test "reader: accepts plus line description" {
    const data =
        \\@read3
        \\AAAA
        \\+repeat_id
        \\!!!!
        \\
    ;
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    const record = (try reader.next()).?;
    try std.testing.expectEqualStrings("repeat_id", record.plus);
    try std.testing.expectEqualStrings("!!!!", record.quality);
}

test "reader: accepts final quality line without trailing newline" {
    const data = "@truncated\nACGT\n+\nIIII";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    _ = try reader.next();
    try std.testing.expect((try reader.next()) == null);
}

test "reader: rejects S005 when sequence and quality lengths differ" {
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

test "reader: rejects S003 when header does not start with @" {
    const data = "not_a_header\nACGT\n+\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.S003InvalidHeader, reader.next());
}

test "reader: rejects S004 on truncated record at EOF" {
    const data = "@truncated\nACGT\n+";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    try std.testing.expectError(zfastq.ReaderError.S004TruncatedRecord, reader.next());
}

test "reader: rejects S001 with line-start details" {
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
}

test "reader: returned fields survive every fixed short-read size" {
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

test "reader and count_scan: every partition preserves the format result" {
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

test "reader: fields survive record-buffer growth" {
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

test "reader: common record performs no allocation after init" {
    const data = "@r\nACGT\n+\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
    var reader = try zfastq.Reader.init(failing.allocator(), source.byteSource(), .{});
    defer reader.deinit();

    _ = try reader.next();
    try std.testing.expect(!failing.has_induced_failure);
}

test "reader: owned record remains stable after next" {
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

test "reader: propagates source read failure" {
    var failing_source = FailingSource{};
    const source = failing_source.byteSource();
    var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
    defer reader.deinit();

    try std.testing.expectError(zfastq.ReaderError.Io, reader.next());
}

test "reader: rejects an invalid source byte count" {
    var invalid_source = InvalidCountSource{};
    var reader = try zfastq.Reader.init(
        std.testing.allocator,
        invalid_source.byteSource(),
        .{},
    );
    defer reader.deinit();

    try std.testing.expectError(zfastq.ReaderError.Io, reader.next());
}

test "reader: owns the ByteSource wrapper value" {
    var original = zfastq.io.plain.SliceSource.init("@original\nA\n+\n!\n");
    var replacement = zfastq.io.plain.SliceSource.init("@replacement\nC\n+\n#\n");
    var source = original.byteSource();
    var reader = try zfastq.Reader.init(std.testing.allocator, source, .{});
    defer reader.deinit();

    source = replacement.byteSource();
    const parsed = (try reader.next()).?;

    try std.testing.expectEqualStrings("original", parsed.id);
}

test "reader: reports allocation failure while growing record storage" {
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

test "reader: handles mixed line endings and tab-delimited id" {
    const data = "@read1\tcomment\r\nACGT\n+name\r\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    const parsed = (try reader.next()).?;
    try std.testing.expectEqualStrings("read1", parsed.id);
    try std.testing.expectEqualStrings("name", parsed.plus);
    try std.testing.expectEqual(@as(u64, data.len), reader.byteOffset());
}

test "reader: enforces exact line limit" {
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
}

test "count_scan: matches reader advance on inline data" {
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

test "count_scan: rejects S005 on length mismatch" {
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

test "count_scan: result and details are chunk invariant" {
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

test "reader and count_scan: structural error details agree across chunk sizes" {
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
        }
    }
}

test "count_scan: dense validation rejects embedded newlines in every field region" {
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

test "count_scan: carries lines across chunks and enforces limit" {
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

test "count_scan: enforces line limit after dense layout learning" {
    const data = "@a\nAC\n+\n!!\n@b\nACGT\n+\n!!!!\n";
    var scan = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = 3 });

    try std.testing.expectError(
        zfastq.ReaderError.LineTooLong,
        zfastq.count_scan.countSlice(data, .{ .max_line_bytes = 3 }, &scan),
    );
    try std.testing.expectEqual(@as(u64, 1), scan.record_index);
}

test "count_scan: preserves length semantics after dense layout learning" {
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

test "count_scan: dense stride block counts uniform records" {
    var buf: [243 * 3]u8 = undefined;
    const pos = appendUniformRecords(&buf, 3, "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/1\n", 100);
    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 3), n);
}
test "count_scan: alternates between known stride layouts" {
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

test "count_scan: alternates between sequence lengths with same header" {
    const hdr = "@ont_read\n";

    var buf: [128]u8 = undefined;
    var pos: usize = 0;
    pos = appendRecord(&buf, pos, hdr, 4);
    pos = appendRecord(&buf, pos, hdr, 6);

    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 2), n);
}

test "count_scan: retains fast path after non-minimal plus record" {
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

test "count_scan: falls back when plus line is not minimal" {
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

test "reader: advance matches next record count" {
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
