//! Writer and plain output-adapter contracts.

const std = @import("std");
const zfastq = @import("z-fastq");

test "[unit] - [writer]: supported records use exact LF serialization" {
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
    };

    for (cases) |case| {
        var buf: [256]u8 = undefined;
        var sink = zfastq.io.plain.SliceSink.init(&buf);
        var writer = zfastq.Writer.init(sink.byteSink());

        try writer.writeRecord(case.record);

        try std.testing.expectEqualStrings(case.expected, sink.written());
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

test "[integration] - [writer]: serialized fields round-trip through Reader" {
    const input =
        \\@read3 lane=1
        \\AAAA
        \\+repeat_id
        \\!!!!
        \\
    ;

    var source = zfastq.io.plain.SliceSource.init(input);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    const parsed = (try reader.next()).?;

    var out_buf: [256]u8 = undefined;
    var sink = zfastq.io.plain.SliceSink.init(&out_buf);
    var writer = zfastq.Writer.init(sink.byteSink());
    try writer.writeRecord(parsed);

    var source2 = zfastq.io.plain.SliceSource.init(sink.written());
    var reader2 = try zfastq.Reader.init(std.testing.allocator, source2.byteSource(), .{});
    defer reader2.deinit();
    const round = (try reader2.next()).?;
    try std.testing.expectEqualStrings(parsed.header, round.header);
    try std.testing.expectEqualStrings(parsed.sequence, round.sequence);
    try std.testing.expectEqualStrings(parsed.plus, round.plus);
    try std.testing.expectEqualStrings(parsed.quality, round.quality);
}

test "[property] - [writer]: generated valid fields round-trip through Reader" {
    var prng = std.Random.DefaultPrng.init(0xbb67ae8584caa73b);
    const random = prng.random();

    for (0..96) |_| {
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

test "[failure] - [writer]: capacity and sink failures propagate" {
    var tiny: [3]u8 = undefined;
    var slice_sink = zfastq.io.plain.SliceSink.init(&tiny);
    var slice_writer = zfastq.Writer.init(slice_sink.byteSink());
    const parsed = zfastq.Record{
        .header = "r",
        .id = "r",
        .sequence = "A",
        .plus = "",
        .quality = "!",
    };
    try std.testing.expectError(error.WriteFailed, slice_writer.writeRecord(parsed));

    var failing = FailingSink{};
    var failing_writer = zfastq.Writer.init(failing.byteSink());
    try std.testing.expectError(error.WriteFailed, failing_writer.writeRecord(parsed));
    try std.testing.expectError(error.WriteFailed, failing_writer.flush());
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

const FailingSink = struct {
    fn byteSink(self: *FailingSink) zfastq.io.ByteSink {
        return .{ .vtable = &vtable, .ctx = self };
    }

    const vtable = zfastq.io.ByteSink.VTable{
        .write = write,
        .flush = flush,
    };

    fn write(ctx: *anyopaque, data: []const u8) error{WriteFailed}!void {
        _ = ctx;
        _ = data;
        return error.WriteFailed;
    }

    fn flush(ctx: *anyopaque) error{WriteFailed}!void {
        _ = ctx;
        return error.WriteFailed;
    }
};

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
