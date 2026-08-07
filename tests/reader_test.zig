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
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);

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
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
    _ = try reader.next();
    try std.testing.expect((try reader.next()) == null);
}

test "reader: accepts plus line without leading plus sign" {
    const data =
        \\@read3
        \\AAAA
        \\repeat_id
        \\!!!!
        \\
    ;
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
    const record = (try reader.next()).?;
    try std.testing.expectEqualStrings("repeat_id", record.plus);
    try std.testing.expectEqualStrings("!!!!", record.quality);
}

test "reader: accepts final quality line without trailing newline" {
    const data = "@truncated\nACGT\n+\nIIII";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
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
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
    try std.testing.expectError(zfastq.ReaderError.S005LengthMismatch, reader.next());
}

test "reader: rejects S003 when header does not start with @" {
    const data = "not_a_header\nACGT\n+\n!!!!\n";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
    try std.testing.expectError(zfastq.ReaderError.S003InvalidHeader, reader.next());
}

test "reader: rejects S004 on truncated record at EOF" {
    const data = "@truncated\nACGT\n+";
    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);
    try std.testing.expectError(zfastq.ReaderError.S004TruncatedRecord, reader.next());
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
    _ = try scan.feed(data, true);
    try scan.finishEof();
    try std.testing.expectEqual(@as(u64, 2), scan.record_index);

    var source = zfastq.io.plain.SliceSource.init(data);
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);

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
    try std.testing.expectError(zfastq.ReaderError.S005LengthMismatch, scan.feed(data, true));
}

test "count_scan: dense stride block counts uniform records" {
    var buf: [243 * 3]u8 = undefined;
    const pos = appendUniformRecords(&buf, 3, "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/1\n", 100);
    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 3), n);
}
test "count_scan: alternates between known stride layouts" {
    const hdr_a = "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/1\n"; // 39 bytes
    const hdr_b = "@HWI-ST180_0186:3:1:1484:1936#GGCTAC/12\n"; // 40 bytes

    var buf: [243 + 244]u8 = undefined;
    var pos: usize = 0;

    inline for (.{ hdr_a, hdr_b }) |hdr| {
        @memcpy(buf[pos..][0..hdr.len], hdr);
        pos += hdr.len;
        for (0..100) |_| {
            buf[pos] = 'A';
            pos += 1;
        }
        buf[pos] = '\n';
        pos += 1;
        buf[pos] = '+';
        pos += 1;
        buf[pos] = '\n';
        pos += 1;
        for (0..100) |_| {
            buf[pos] = 'I';
            pos += 1;
        }
        buf[pos] = '\n';
        pos += 1;
    }

    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 2), n);
}

test "count_scan: alternates between sequence lengths with same header" {
    const hdr = "@ont_read\n";

    var buf: [128]u8 = undefined;
    var pos: usize = 0;

    // seq len 4
    @memcpy(buf[pos..][0..hdr.len], hdr);
    pos += hdr.len;
    @memcpy(buf[pos..][0..4], "ACGT");
    pos += 4;
    buf[pos] = '\n';
    pos += 1;
    buf[pos] = '+';
    pos += 1;
    buf[pos] = '\n';
    pos += 1;
    @memcpy(buf[pos..][0..4], "!!!!");
    pos += 4;
    buf[pos] = '\n';
    pos += 1;

    // seq len 6
    @memcpy(buf[pos..][0..hdr.len], hdr);
    pos += hdr.len;
    @memcpy(buf[pos..][0..6], "ACGTAC");
    pos += 6;
    buf[pos] = '\n';
    pos += 1;
    buf[pos] = '+';
    pos += 1;
    buf[pos] = '\n';
    pos += 1;
    @memcpy(buf[pos..][0..6], "!!!!!!");
    pos += 6;
    buf[pos] = '\n';
    pos += 1;

    var scan = zfastq.count_scan.Scanner.init(.{});
    const n = try zfastq.count_scan.countSlice(buf[0..pos], .{}, &scan);
    try std.testing.expectEqual(@as(u64, 2), n);
}

test "count_scan: retains fast path after non-minimal plus record" {
    const data =
        \\@read3
        \\AAAA
        \\repeat_id
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
        \\repeat_id
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
    var reader = try zfastq.Reader.init(std.testing.allocator, &source.byteSource(), .{});
    defer reader.deinit(std.testing.allocator);

    var advance_count: usize = 0;
    while (try reader.advance()) {
        advance_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 2), advance_count);

    var source2 = zfastq.io.plain.SliceSource.init(data);
    var reader2 = try zfastq.Reader.init(std.testing.allocator, &source2.byteSource(), .{});
    defer reader2.deinit(std.testing.allocator);

    var next_count: usize = 0;
    while ((try reader2.next()) != null) {
        next_count += 1;
    }
    try std.testing.expectEqual(advance_count, next_count);
}

fn appendUniformRecords(buf: []u8, count: usize, header: []const u8, seq_len: usize) usize {
    var pos: usize = 0;
    for (0..count) |_| {
        @memcpy(buf[pos..][0..header.len], header);
        pos += header.len;
        for (0..seq_len) |_| {
            buf[pos] = 'A';
            pos += 1;
        }
        buf[pos] = '\n';
        pos += 1;
        buf[pos] = '+';
        pos += 1;
        buf[pos] = '\n';
        pos += 1;
        for (0..seq_len) |_| {
            buf[pos] = 'I';
            pos += 1;
        }
        buf[pos] = '\n';
        pos += 1;
    }
    return pos;
}
