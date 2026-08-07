//! Unit tests for `fastq.Writer`.

const std = @import("std");
const zfastq = @import("z-fastq");

test "writer: emits LF FASTQ lines" {
    var buf: [256]u8 = undefined;
    var sink = zfastq.io.plain.SliceSink.init(&buf);
    var writer = zfastq.Writer.init(&sink.byteSink());
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

test "writer: round-trip preserves record fields" {
    const input =
        \\@read3 lane=1
        \\AAAA
        \\repeat_id
        \\!!!!
        \\
    ;

    var source = zfastq.io.plain.SliceSource.init(input);
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
    const parsed = (try reader.next()).?;

    var out_buf: [256]u8 = undefined;
    var sink = zfastq.io.plain.SliceSink.init(&out_buf);
    var writer = zfastq.Writer.init(&sink.byteSink());
    try writer.writeRecord(parsed);

    var source2 = zfastq.io.plain.SliceSource.init(sink.written());
    var reader2 = try zfastq.Reader.init(std.testing.allocator, &source2.byteSource(), .{});
    defer reader2.deinit(std.testing.allocator);
    const round = (try reader2.next()).?;
    try std.testing.expectEqualStrings(parsed.header, round.header);
    try std.testing.expectEqualStrings(parsed.sequence, round.sequence);
    try std.testing.expectEqualStrings(parsed.plus, round.plus);
    try std.testing.expectEqualStrings(parsed.quality, round.quality);
}
