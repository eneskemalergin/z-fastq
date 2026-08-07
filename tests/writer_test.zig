//! Unit tests for `fastq.Writer`.

const std = @import("std");
const zfastq = @import("z-fastq");

test "writer: emits LF FASTQ lines" {
    var buf: [256]u8 = undefined;
    var sink = zfastq.io.plain.SliceSink.init(&buf);
    var writer = zfastq.Writer.init(sink.byteSink());
    const record = zfastq.Record{
        .header = "read1 run=1",
        .id = "read1",
        .sequence = "ACGT",
        .plus = "",
        .quality = "!!!!",
    };
    try writer.writeRecord(record);
    try std.testing.expectEqualStrings("@read1 run=1\nACGT\n+\n!!!!\n", sink.written());
}

test "writer: rejects structurally invalid records before writing" {
    const cases = [_]zfastq.Record{
        .{ .header = "r", .id = "r", .sequence = "AA", .plus = "", .quality = "!" },
        .{ .header = "r\nnext", .id = "r", .sequence = "A", .plus = "", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "A\n", .plus = "", .quality = "!!" },
        .{ .header = "r", .id = "r", .sequence = "A", .plus = "x\n", .quality = "!" },
        .{ .header = "r", .id = "r", .sequence = "AA", .plus = "", .quality = "!\n" },
    };

    for (cases) |record| {
        var buf: [64]u8 = undefined;
        var sink = zfastq.io.plain.SliceSink.init(&buf);
        var writer = zfastq.Writer.init(sink.byteSink());

        try std.testing.expectError(error.InvalidRecord, writer.writeRecord(record));
        try std.testing.expectEqual(@as(usize, 0), sink.written().len);
    }
}

test "writer: round-trip preserves record fields" {
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

test "plain adapters: file sink writes and flushes after init return" {
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

test "plain adapters: reader and writer wrappers preserve bytes" {
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

test "plain adapters: file source reads through Reader" {
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

test "writer: propagates capacity and sink failures" {
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

test "writer: owns the ByteSink wrapper value" {
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
