//! Writer and plain output-adapter contracts.

const std = @import("std");
const zfastq = @import("z-fastq");

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
        var buf: [256]u8 = undefined;
        var sink = zfastq.io.plain.SliceSink.init(&buf);
        var writer = zfastq.Writer.init(sink.byteSink());

        try writer.writeRecord(case.record);
        try writer.flush();

        try std.testing.expectEqualStrings(case.expected, sink.written());

        var source = zfastq.io.plain.SliceSource.init(sink.written());
        var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        const parsed = (try reader.next()).?;
        try std.testing.expectEqualDeep(case.record, parsed);
        try std.testing.expect((try reader.next()) == null);
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

test "[integration] - [writer]: generated fields serialize exactly and round-trip" {
    var prng = std.Random.DefaultPrng.init(0xbb67ae8584caa73b);
    const random = prng.random();

    for (0..96) |case_index| {
        errdefer std.debug.print(
            "writer seed=0xbb67ae8584caa73b runner_seed=0x{x} case={d}\n",
            .{ std.testing.random_seed, case_index },
        );
        var header_storage: [48]u8 = undefined;
        const header_len = random.intRangeAtMost(usize, 1, header_storage.len);
        for (header_storage[0..header_len]) |*byte| {
            byte.* = 'a' + random.uintLessThan(u8, 26);
        }
        if (header_len > 2 and random.boolean()) {
            header_storage[random.intRangeLessThan(usize, 1, header_len - 1)] = '\r';
        }
        const header = header_storage[0..header_len];

        var sequence_storage: [96]u8 = undefined;
        const sequence_len = random.intRangeAtMost(usize, 0, sequence_storage.len);
        for (sequence_storage[0..sequence_len]) |*byte| {
            byte.* = "ACGTN"[random.uintLessThan(usize, 5)];
        }
        if (sequence_len > 2 and random.boolean()) {
            sequence_storage[random.intRangeLessThan(usize, 1, sequence_len - 1)] = '\r';
        }
        const sequence = sequence_storage[0..sequence_len];

        var plus_storage: [32]u8 = undefined;
        const plus_len = random.intRangeAtMost(usize, 0, plus_storage.len);
        for (plus_storage[0..plus_len]) |*byte| {
            byte.* = 'a' + random.uintLessThan(u8, 26);
        }
        if (plus_len > 2 and random.boolean()) {
            plus_storage[random.intRangeLessThan(usize, 1, plus_len - 1)] = '\r';
        }
        const plus = plus_storage[0..plus_len];

        var quality_storage: [96]u8 = undefined;
        for (quality_storage[0..sequence_len]) |*byte| {
            byte.* = '!' + random.uintLessThan(u8, 41);
        }
        if (sequence_len > 2 and random.boolean()) {
            quality_storage[random.intRangeLessThan(usize, 1, sequence_len - 1)] = '\r';
        }
        const quality = quality_storage[0..sequence_len];

        var output: [512]u8 = undefined;
        var sink = zfastq.io.plain.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        try writer.writeRecord(.{
            .header = header,
            .id = header,
            .sequence = sequence,
            .plus = plus,
            .quality = quality,
        });
        try writer.flush();

        var expected_buffer: [512]u8 = undefined;
        const expected = try std.fmt.bufPrint(
            &expected_buffer,
            "@{s}\n{s}\n+{s}\n{s}\n",
            .{ header, sequence, plus, quality },
        );
        try std.testing.expectEqualStrings(expected, sink.written());

        var source = zfastq.io.plain.SliceSource.init(sink.written());
        var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        const round = (try reader.next()).?;

        try std.testing.expectEqualStrings(header, round.header);
        try std.testing.expectEqualStrings(sequence, round.sequence);
        try std.testing.expectEqualStrings(plus, round.plus);
        try std.testing.expectEqualStrings(quality, round.quality);
        try std.testing.expect((try reader.next()) == null);
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

test "[integration] - [plain adapters]: standard reader and writer preserve bytes" {
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
}

test "[integration] - [file source]: supplies records to Reader" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(io, "input.fastq", .{ .read = true });
    defer file.close(io);
    try file.writePositionalAll(io, "@file\nAC\n+\n!!\n", 0);

    var read_buf: [3]u8 = undefined;
    var file_source = zfastq.io.plain.FileSource.init(io, file, &read_buf);
    var reader = try zfastq.Reader.init(std.testing.allocator, file_source.byteSource(), .{});
    defer reader.deinit();
    const parsed = (try reader.next()).?;
    try std.testing.expectEqualStrings("file", parsed.id);
    try std.testing.expectEqualStrings("AC", parsed.sequence);
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

test "[unit] - [byte sink]: flush succeeds when the sink has no flush callback" {
    var sink_state = NoFlushSink{};
    const sink = sink_state.byteSink();

    try sink.write("abc");
    try sink.flush();

    try std.testing.expectEqual(@as(usize, 3), sink_state.bytes_written);
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

const NoFlushSink = struct {
    bytes_written: usize = 0,

    fn byteSink(self: *NoFlushSink) zfastq.io.ByteSink {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSink.VTable{ .write = write };

    fn write(ctx: *anyopaque, data: []const u8) error{WriteFailed}!void {
        const self: *NoFlushSink = @ptrCast(@alignCast(ctx));
        self.bytes_written += data.len;
    }
};
