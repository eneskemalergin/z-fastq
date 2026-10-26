//! Writer output bytes, record round trips, sink failures, and plain output adapters.

const std = @import("std");
const zfastq = @import("z-fastq");

const PAGE_BYTES = 4096;
const PATTERN_MULTIPLIER = 31;
const PATTERN_BLOCK = 251;
const STEP_FACTOR = 7;
const STEP_LIMIT = 5003;
const VEC_EVERY = 3;

fn expectRoundTrip(record: zfastq.Record, expected: []const u8) !void {
    var buf: [512]u8 = undefined;
    var sink = zfastq.io.plain.SliceSink.init(&buf);
    var writer = zfastq.Writer.init(sink.byteSink());

    try writer.writeRecord(record);
    try writer.flush();

    try std.testing.expectEqualStrings(expected, sink.written());

    var source = zfastq.io.plain.SliceSource.init(sink.written());
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    const parsed = (try reader.next()).?;
    try std.testing.expectEqualDeep(record, parsed);
    try std.testing.expect((try reader.next()) == null);
}

test "[integration] - [writer]: supported records serialize and round-trip exactly" {
    const cases = [_]struct {
        record: zfastq.Record,
        expected: []const u8,
    }{
        .{
            .record = .{
                .header = "read1 run=1",
                .id = "read1",
                .sequence = "ACGT",
                .plus = "",
                .quality = "!!!!",
            },
            .expected = "@read1 run=1\nACGT\n+\n!!!!\n",
        },
        .{
            .record = .{
                .header = "read3 lane=1",
                .id = "read3",
                .sequence = "AAAA",
                .plus = "repeat_id",
                .quality = "!!!!",
            },
            .expected = "@read3 lane=1\nAAAA\n+repeat_id\n!!!!\n",
        },
        .{
            .record = .{
                .header = "empty",
                .id = "empty",
                .sequence = "",
                .plus = "description",
                .quality = "",
            },
            .expected = "@empty\n\n+description\n\n",
        },
        .{
            .record = .{
                .header = "r\rcomment",
                .id = "r\rcomment",
                .sequence = "A\rC",
                .plus = "x\ry",
                .quality = "!\r!",
            },
            .expected = "@r\rcomment\nA\rC\n+x\ry\n!\r!\n",
        },
        .{
            .record = .{
                .header = "r\x00\xff opaque",
                .id = "r\x00\xff",
                .sequence = "A\x00\xff",
                .plus = "note\x00\xff",
                .quality = "!\x00\xff",
            },
            .expected = "@r\x00\xff opaque\nA\x00\xff\n+note\x00\xff\n!\x00\xff\n",
        },
    };

    for (cases) |case| {
        try expectRoundTrip(case.record, case.expected);
    }
}

test "[failure] - [writer]: invalid fields are rejected before output" {
    const cases = [_]zfastq.Record{
        .{ .header = "", .id = "", .sequence = "A", .plus = "", .quality = "!" },
        .{ .header = " description", .id = "", .sequence = "A", .plus = "", .quality = "!" },
        .{ .header = "\tdescription", .id = "", .sequence = "A", .plus = "", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "AA", .plus = "", .quality = "!" },
        .{ .header = "r\nnext", .id = "r", .sequence = "A", .plus = "", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "A\n", .plus = "", .quality = "!!" },
        .{ .header = "r", .id = "r", .sequence = "A", .plus = "x\n", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "AA", .plus = "", .quality = "!\n" },
        .{ .header = "r\r", .id = "r", .sequence = "A", .plus = "", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "A\r", .plus = "", .quality = "!x" },
        .{ .header = "r", .id = "r", .sequence = "A", .plus = "x\r", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "Ax", .plus = "", .quality = "!\r" },
    };

    for (cases) |record| {
        var buf: [64]u8 = undefined;
        var sink = zfastq.io.plain.SliceSink.init(&buf);
        var writer = zfastq.Writer.init(sink.byteSink());

        try std.testing.expectError(error.InvalidRecord, writer.writeRecord(record));
        try std.testing.expectEqual(@as(usize, 0), sink.written().len);
    }
}

test "[edge] - [writer]: empty, interior-CR, and long fields serialize exactly and round-trip" {
    const cases = [_]struct {
        record: zfastq.Record,
        expected: []const u8,
    }{
        .{
            .record = .{
                .header = "r",
                .id = "r",
                .sequence = "",
                .plus = "",
                .quality = "",
            },
            .expected = "@r\n\n+\n\n",
        },
        .{
            .record = .{
                .header = "r\rh",
                .id = "r\rh",
                .sequence = "A\rC",
                .plus = "p\rq",
                .quality = "!\r~",
            },
            .expected = "@r\rh\nA\rC\n+p\rq\n!\r~\n",
        },
        .{
            .record = .{
                .header = "h" ** 48,
                .id = "h" ** 48,
                .sequence = "A\rCG" ** 24,
                .plus = "p" ** 32,
                .quality = "!\r~" ** 32,
            },
            .expected = "@" ++ ("h" ** 48) ++ "\n" ++ ("A\rCG" ** 24) ++ "\n+" ++ ("p" ** 32) ++ "\n" ++ ("!\r~" ** 32) ++ "\n",
        },
    };

    for (cases) |case| {
        try expectRoundTrip(case.record, case.expected);
    }
}

test "[integration] - [file sink]: writes and flushes after init returns" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(io, "record.fastq", .{ .read = true });
    defer file.close(io);

    var write_buf: [32]u8 = undefined;
    var file_sink = zfastq.io.plain.FileSink.init(io, file, &write_buf);
    var writer = zfastq.Writer.init(file_sink.byteSink());
    try writer.writeRecord(.{
        .header = "read1",
        .id = "read1",
        .sequence = "ACGT",
        .plus = "",
        .quality = "!!!!",
    });
    try writer.flush();

    var actual: [64]u8 = undefined;
    const n = try file.readPositionalAll(io, &actual, 0);
    try std.testing.expectEqualStrings("@read1\nACGT\n+\n!!!!\n", actual[0..n]);
}

test "[property] - [file sink]: page-sized writes keep every byte in order" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var pattern: [3 * PAGE_BYTES + 777]u8 = undefined;
    for (&pattern, 0..) |*byte, index| {
        byte.* = @truncate(index *% PATTERN_MULTIPLIER +% index / PATTERN_BLOCK);
    }

    for ([_]usize{ 100, PAGE_BYTES + 100, 2 * PAGE_BYTES }) |buffer_len| {
        const file = try tmp.dir.createFile(io, "sink", .{ .read = true, .truncate = true });
        defer file.close(io);
        var write_buf: [2 * PAGE_BYTES]u8 = undefined;
        var file_sink = zfastq.io.plain.FileSink.init(io, file, write_buf[0..buffer_len]);
        const sink = file_sink.byteSink();

        var expected: std.ArrayList(u8) = .empty;
        defer expected.deinit(std.testing.allocator);
        var offset: usize = 0;
        var step: usize = 1;
        while (offset < pattern.len) : (step = step * STEP_FACTOR % STEP_LIMIT + 1) {
            const piece = pattern[offset..@min(pattern.len, offset + step)];
            if (step % VEC_EVERY == 0) {
                try sink.writeVec(&.{ piece[0 .. piece.len / 2], piece[piece.len / 2 ..] });
            } else {
                try sink.write(piece);
            }
            try expected.appendSlice(std.testing.allocator, piece);
            offset += piece.len;
        }
        try sink.flush();

        const actual = try tmp.dir.readFileAlloc(io, "sink", std.testing.allocator, .unlimited);
        defer std.testing.allocator.free(actual);
        try std.testing.expectEqualSlices(u8, expected.items, actual);
    }
}

test "[integration] - [plain adapters]: reader, file source, and writer carry records exactly" {
    const input = "@r\nA\n+\n!\n";
    var fixed_reader = std.Io.Reader.fixed(input);
    var reader_source = zfastq.io.plain.ReaderSource.init(&fixed_reader);
    var reader = try zfastq.Reader.init(std.testing.allocator, reader_source.byteSource(), .{});
    defer reader.deinit();
    const parsed = (try reader.next()).?;

    var output: [32]u8 = undefined;
    var fixed_writer = std.Io.Writer.fixed(&output);
    var writer_sink = zfastq.io.plain.WriterSink.init(&fixed_writer);
    var writer = zfastq.Writer.init(writer_sink.byteSink());
    try writer.writeRecord(parsed);
    try writer.flush();
    try std.testing.expectEqualStrings(input, fixed_writer.buffered());

    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(io, "input.fastq", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, "@file\nAC\n+\n!!\n", 0);

    var read_buf: [3]u8 = undefined;
    var file_source = zfastq.io.plain.FileSource.init(io, file, &read_buf);
    var file_reader = try zfastq.Reader.init(std.testing.allocator, file_source.byteSource(), .{});
    defer file_reader.deinit();
    const file_record = (try file_reader.next()).?;
    try std.testing.expectEqualDeep(zfastq.Record{
        .header = "file",
        .id = "file",
        .sequence = "AC",
        .plus = "",
        .quality = "!!",
    }, file_record);
    try std.testing.expect((try file_reader.next()) == null);
}

test "[integration] - [writer]: capacity, write, and flush failures propagate" {
    const parsed = zfastq.Record{
        .header = "r",
        .id = "r",
        .sequence = "A",
        .plus = "",
        .quality = "!",
    };
    const complete = "@r\nA\n+\n!\n";
    const cases = [_]struct { capacity: usize, expected: []const u8 }{
        .{ .capacity = 0, .expected = "" },
        .{ .capacity = 3, .expected = "@r\n" },
        .{ .capacity = complete.len - 1, .expected = "@r\nA\n+\n!" },
        .{ .capacity = complete.len, .expected = complete },
    };
    for (cases) |case| {
        var output: [complete.len]u8 = undefined;
        var slice_sink = zfastq.io.plain.SliceSink.init(output[0..case.capacity]);
        var slice_writer = zfastq.Writer.init(slice_sink.byteSink());
        if (case.capacity < complete.len) {
            try std.testing.expectError(error.WriteFailed, slice_writer.writeRecord(parsed));
        } else {
            try slice_writer.writeRecord(parsed);
        }
        try slice_writer.flush();
        try std.testing.expectEqualStrings(case.expected, slice_sink.written());
    }

    var pending: [1]u8 = undefined;
    var standard_writer: std.Io.Writer = .{
        .vtable = &.{ .drain = std.Io.Writer.failingDrain },
        .buffer = &pending,
    };
    var adapter = zfastq.io.plain.WriterSink.init(&standard_writer);
    var failing_writer = zfastq.Writer.init(adapter.byteSink());
    try std.testing.expectError(error.WriteFailed, failing_writer.writeRecord(parsed));
    try std.testing.expectEqualStrings("@", standard_writer.buffered());
    try std.testing.expectError(error.WriteFailed, failing_writer.flush());
}

test "[integration] - [plain adapters]: access errors preserve borrowed file handles" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const write_only = try tmp.dir.createFile(io, "input", .{});
    defer write_only.close(io);
    var read_buffer: [8]u8 = undefined;
    var file_source = zfastq.io.plain.FileSource.init(io, write_only, &read_buffer);
    var output: [1]u8 = undefined;
    try std.testing.expectError(error.ReadFailed, file_source.byteSource().read(&output));
    try write_only.writeStreamingAll(io, "kept");

    var failing_reader: std.Io.Reader = .failing;
    var reader_source = zfastq.io.plain.ReaderSource.init(&failing_reader);
    try std.testing.expectError(error.ReadFailed, reader_source.byteSource().read(&output));

    const read_only = try tmp.dir.openFile(io, "input", .{});
    defer read_only.close(io);
    for ([_]usize{ 0, 8 }) |buffer_len| {
        var write_buffer: [8]u8 = undefined;
        var file_sink = zfastq.io.plain.FileSink.init(io, read_only, write_buffer[0..buffer_len]);
        const sink = file_sink.byteSink();
        if (buffer_len == 0) {
            try std.testing.expectError(error.WriteFailed, sink.write("x"));
        } else {
            try sink.write("x");
            try std.testing.expectError(error.WriteFailed, sink.flush());
        }
        var contents: [8]u8 = undefined;
        const n = try read_only.readPositionalAll(io, &contents, 0);
        try std.testing.expectEqualStrings("kept", contents[0..n]);
    }
}

test "[unit] - [writer]: the ByteSink wrapper is copied by value" {
    var original_buf: [32]u8 = undefined;
    var replacement_buf: [32]u8 = undefined;
    var original = zfastq.io.plain.SliceSink.init(&original_buf);
    var replacement = zfastq.io.plain.SliceSink.init(&replacement_buf);
    var sink = original.byteSink();
    var writer = zfastq.Writer.init(sink);

    sink = replacement.byteSink();
    try writer.writeRecord(.{
        .header = "r",
        .id = "r",
        .sequence = "A",
        .plus = "",
        .quality = "!",
    });

    try std.testing.expectEqualStrings("@r\nA\n+\n!\n", original.written());
    try std.testing.expectEqual(@as(usize, 0), replacement.written().len);
}
