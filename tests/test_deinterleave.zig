//! Installed `z-fastq deinterleave` behavior: pairing, two outputs, validation, and exit status.

const std = @import("std");
const cli = @import("utilities.zig");

const EXACT_INPUT =
    "@same left\nA\n+one\n!\n@same right\nT\n+two\n#\n" ++
    "@empty one\n\n+empty-one\n\n@empty two\n\n+empty-two\n\n";
const EXACT_R1 =
    "@same left\nA\n+one\n!\n" ++
    "@empty one\n\n+empty-one\n\n";
const EXACT_R2 =
    "@same right\nT\n+two\n#\n" ++
    "@empty two\n\n+empty-two\n\n";

const READ_BUFFER_BYTES = 256 * 1024;
const OUTPUT_BUFFER_BYTES = 64 * 1024;
const GZIP_MEMBER_INPUT_BYTES = 60 * 1024;

test "[cli] - [deinterleave]: mates are written in input order with canonical LF output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");

    const input =
        "@cluster 1:N:0:index-a\r\nAC\r\n+first annotation\r\n!~\r\n" ++
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n" ++
        "@legacy/2 other\r\nA\r\n+mate\r\n!";
    const expected1 =
        "@cluster 1:N:0:index-a\nAC\n+first annotation\n!~\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n";
    const expected2 =
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/2 other\nA\n+mate\n!\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = input });

    try cli.expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--out1", out1_path, "--out2", out2_path, input_path }),
        0,
        "",
        "",
    );
    try expectFile(allocator, out1_path, expected1);
    try expectFile(allocator, out2_path, expected2);
}

test "[cli] - [deinterleave]: interleaving the outputs and splitting again reproduces them" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");

    const expected1 =
        "@cluster 1:N:0:index-a\nAC\n+first annotation\n!~\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n";
    const expected2 =
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/2 other\nA\n+mate\n!\n";
    const expected_interleaved =
        "@cluster 1:N:0:index-a\nAC\n+first annotation\n!~\n" ++
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n" ++
        "@legacy/2 other\nA\n+mate\n!\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "out1.fastq", .data = expected1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "out2.fastq", .data = expected2 });

    const joined = try cli.run(allocator, &.{ "interleave", out1_path, out2_path });
    try cli.expectResult(joined, 0, expected_interleaved, "");

    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = joined.stdout });
    try tmp.dir.deleteFile(io, "out1.fastq");
    try tmp.dir.deleteFile(io, "out2.fastq");
    try cli.expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--out1", out1_path, "--out2", out2_path, input_path }),
        0,
        "",
        "",
    );
    try expectFile(allocator, out1_path, expected1);
    try expectFile(allocator, out2_path, expected2);
}

test "[cli] - [deinterleave]: terminal CR fields leave both outputs empty" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for (0..2) |bad_mate| {
        var tmp = std.testing.tmpDir(.{});
        defer tmp.cleanup();
        const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
        const out1 = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
        const out2 = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
        const r1 = if (bad_mate == 0) "@pair/1 x\r\r\nA\n+\n!\n" else "@pair/1\nA\n+\n!\n";
        const r2 = if (bad_mate == 1) "@pair/2\nT\n+\r\r\n#\n" else "@pair/2\nT\n+\n#\n";
        const input = try std.mem.concat(allocator, u8, &.{ r1, r2 });
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = input });
        const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: record fields ending in CR cannot be written with LF endings\n", .{input_path});
        try cli.expectResult(try cli.run(allocator, &.{ "deinterleave", "--out1", out1, "--out2", out2, input_path }), 1, "", diagnostic);
        try expectFile(allocator, out1, "");
        try expectFile(allocator, out2, "");
    }
}

test "[cli] - [deinterleave]: exact policy preserves fields from plain, gzip, and stdin input" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const plain_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const gzip_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq.gz");
    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, EXACT_INPUT, .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = EXACT_INPUT });
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq.gz", .data = gzip.items });

    const transports = [_]struct { name: []const u8, input: []const u8, stdin: ?[]const u8 }{
        .{ .name = "plain", .input = plain_path, .stdin = null },
        .{ .name = "gzip", .input = gzip_path, .stdin = null },
        .{ .name = "stdin", .input = "-", .stdin = EXACT_INPUT },
    };
    for (transports) |transport| {
        const out1 = try cli.tempPath(allocator, &tmp.sub_path, try std.fmt.allocPrint(allocator, "{s}-r1.fastq", .{transport.name}));
        const out2 = try cli.tempPath(allocator, &tmp.sub_path, try std.fmt.allocPrint(allocator, "{s}-r2.fastq", .{transport.name}));
        const args = [_][]const u8{
            "deinterleave", "--pair-names", "exact",  "--alphabet", "acgtn",
            "--out1",       out1,           "--out2", out2,         transport.input,
        };
        const result = if (transport.stdin) |stdin_data|
            try cli.runWithStdin(allocator, &args, stdin_data, 1)
        else
            try cli.run(allocator, &args);
        try cli.expectResult(result, 0, "", "");
        try expectFile(allocator, out1, EXACT_R1);
        try expectFile(allocator, out2, EXACT_R2);
    }
}

test "[cli] - [deinterleave]: empty input creates two empty outputs" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const empty_out1 = try cli.tempPath(allocator, &tmp.sub_path, "empty-r1.fastq");
    const empty_out2 = try cli.tempPath(allocator, &tmp.sub_path, "empty-r2.fastq");
    try cli.expectResult(
        try cli.runWithStdin(allocator, &.{ "deinterleave", "--out1", empty_out1, "--out2", empty_out2, "-" }, "", 1),
        0,
        "",
        "",
    );
    try expectEmptyFile(empty_out1);
    try expectEmptyFile(empty_out2);
}

test "[cli] - [deinterleave]: exact policy compares complete first tokens" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const paths = try outputPaths(allocator, &tmp.sub_path, "exact-mismatch");
    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@cluster/1\nA\n+\n!\n@cluster/2\nT\n+\n#\n",
    });
    const expected = try std.fmt.allocPrint(
        allocator,
        "error: {s}: P001: paired identifiers or mate markers do not match (pair 0)\n",
        .{input_path},
    );

    const result = try cli.run(allocator, &.{
        "deinterleave",
        "--pair-names",
        "exact",
        "--out1",
        paths[0],
        "--out2",
        paths[1],
        input_path,
    });
    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqualStrings("", result.stdout);
    try cli.expectFirstLine(expected, result.stderr);
    try expectEmptyFile(paths[0]);
    try expectEmptyFile(paths[1]);
}

test "[cli] - [deinterleave]: the earliest failure decides the diagnostic and what the outputs keep" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const cases = [_]struct {
        name: []const u8,
        input: []const u8,
        message: []const u8,
        first_line: bool,
        out1: []const u8 = "",
        out2: []const u8 = "",
    }{
        .{
            .name = "structure",
            .input = "@bad/1\n.\n+\n!\n@bad/2\nA\nx\n!\n",
            .message = "S001: plus line must start with '+' (record 1, line 3, offset 22)\n",
            .first_line = false,
        },
        .{
            .name = "iupac-quality",
            .input = "@bad/1\nR\n+\n \n@bad/2\nA\nx\n!\n",
            .message = "S001: plus line must start with '+' (record 1, line 3, offset 22)\n",
            .first_line = false,
        },
        .{
            .name = "mismatch",
            .input = "@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n",
            .message = "P001: paired identifiers or mate markers do not match (pair 0)\n",
            .first_line = true,
        },
        .{
            .name = "odd",
            .input = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@odd/1\nA\n+\n!\n",
            .message = "P002: paired input is missing a mate " ++
                "(pair 1, remaining R1, last R1 record 2, last R2 record 1)\n",
            .first_line = true,
            .out1 = "@ok/1\nA\n+\n!\n",
            .out2 = "@ok/2\nT\n+\n#\n",
        },
        .{
            .name = "odd-semantic",
            .input = "@odd/1\n.\n+\n!\n",
            .message = "P002: paired input is missing a mate " ++
                "(pair 0, remaining R1, last R1 record 0, last R2 record none)\n",
            .first_line = true,
        },
        .{
            .name = "odd-structure",
            .input = "@odd/1\nA\n+\n",
            .message = "S004: unexpected end of file in quality line (record 0, line 4, offset 11)\n",
            .first_line = false,
        },
    };
    for (cases) |case| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = case.input });
        const paths = try outputPaths(allocator, &tmp.sub_path, case.name);
        const message = try std.fmt.allocPrint(allocator, "error: {s}: {s}", .{ input_path, case.message });
        const result = try runDeinterleave(allocator, input_path, paths);
        try std.testing.expectEqual(@as(u8, 1), result.exit_code);
        try std.testing.expectEqualStrings("", result.stdout);
        if (case.first_line) {
            try cli.expectFirstLine(message, result.stderr);
        } else {
            try std.testing.expectEqualStrings(message, result.stderr);
        }
        try expectFile(allocator, paths[0], case.out1);
        try expectFile(allocator, paths[1], case.out2);
    }
}

test "[cli] - [deinterleave]: late gzip failure retains the same complete prefix in both outputs" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq.gz");
    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    try appendFastqRecord(allocator, &r1, "ok/1", 'A', '!', 151);
    try appendFastqRecord(allocator, &r2, "ok/2", 'T', '#', 51);
    var input: std.ArrayList(u8) = .empty;
    for (0..2000) |_| {
        try input.appendSlice(allocator, r1.items);
        try input.appendSlice(allocator, r2.items);
    }
    var gzip: std.ArrayList(u8) = .empty;
    try appendGzipChunks(allocator, &gzip, input.items);
    gzip.items[cli.gzipCrcOffset(gzip.items.len)] ^= 1;
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq.gz", .data = gzip.items });
    const paths = try outputPaths(allocator, &tmp.sub_path, "corrupt-prefix");
    const expected_stderr = try std.fmt.allocPrint(allocator, "error: {s}: I/O error\n", .{input_path});

    try cli.expectResult(try runDeinterleave(allocator, input_path, paths), 3, "", expected_stderr);
    const output1 = try tmp.dir.readFileAlloc(io, "corrupt-prefix-r1.fastq", allocator, .limited(1024 * 1024));
    const output2 = try tmp.dir.readFileAlloc(io, "corrupt-prefix-r2.fastq", allocator, .limited(1024 * 1024));
    try std.testing.expect(output1.len > 64 * 1024);
    try std.testing.expect(output2.len > 64 * 1024);
    try std.testing.expectEqual(@as(usize, 0), output1.len % r1.items.len);
    try std.testing.expectEqual(@as(usize, 0), output2.len % r2.items.len);
    const pairs = output1.len / r1.items.len;
    try std.testing.expect(pairs <= 2000);
    try std.testing.expectEqual(pairs, output2.len / r2.items.len);
    for (0..pairs) |index| {
        try std.testing.expectEqualSlices(u8, r1.items, output1[index * r1.items.len ..][0..r1.items.len]);
        try std.testing.expectEqualSlices(u8, r2.items, output2[index * r2.items.len ..][0..r2.items.len]);
    }
}

test "[cli] - [deinterleave]: validation failures in the first pair leave empty outputs" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");

    const cases = [_]struct {
        prefix: []const u8,
        input: []const u8,
        exit_code: u8,
        error_tail: []const u8,
    }{
        .{
            .prefix = "mate1-structure",
            .input = "@bad/1\nA\nx\n!\n@bad/2\nA\n+\n!\n",
            .exit_code = 1,
            .error_tail = "S001: plus line must start with '+' " ++
                "(record 0, line 3, offset 9)\n",
        },
        .{
            .prefix = "mate1-semantic",
            .input = "@bad/1\n.\n+\n!\n@bad/2\nA\n+\n!\n",
            .exit_code = 1,
            .error_tail = "S002: sequence byte is outside the selected alphabet " ++
                "(record 0, line 2, offset 7)\n",
        },
        .{
            .prefix = "mate2-semantic",
            .input = "@bad/1\nA\n+\n!\n@bad/2\n.\n+\n!\n",
            .exit_code = 1,
            .error_tail = "S002: sequence byte is outside the selected alphabet " ++
                "(record 1, line 2, offset 20)\n",
        },
        .{
            .prefix = "both-semantic",
            .input = "@bad/1\n.\n+\n!\n@bad/2\n.\n+\n!\n",
            .exit_code = 1,
            .error_tail = "S002: sequence byte is outside the selected alphabet " ++
                "(record 0, line 2, offset 7)\n",
        },
    };
    for (cases) |case| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = case.input });
        const paths = try outputPaths(allocator, &tmp.sub_path, case.prefix);
        const expected = try std.fmt.allocPrint(allocator, "error: {s}: {s}", .{ input_path, case.error_tail });
        try cli.expectResult(
            try runDeinterleave(allocator, input_path, paths),
            case.exit_code,
            "",
            expected,
        );
        try expectEmptyFile(paths[0]);
        try expectEmptyFile(paths[1]);
    }
}

test "[cli] - [deinterleave]: a later failure keeps the complete pair beyond both output buffers" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const mate_len = OUTPUT_BUFFER_BYTES + 64;
    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    try appendFastqRecord(allocator, &r1, "ok/1", 'A', '!', mate_len);
    try appendFastqRecord(allocator, &r2, "ok/2", 'T', '#', mate_len);
    var input: std.ArrayList(u8) = .empty;
    try input.appendSlice(allocator, r1.items);
    try input.appendSlice(allocator, r2.items);
    const prefix_len = input.items.len;
    const cases = [_]struct {
        tail: []const u8,
        message: []const u8,
        record: usize,
        line: u8,
        offset: usize,
        first_line: bool = false,
    }{
        .{
            .tail = "@bad/1\nA\nx\n!\n@bad/2\nT\n+\n#\n",
            .message = "S001: plus line must start with '+'",
            .record = 2,
            .line = 3,
            .offset = 9,
        },
        .{
            .tail = "@bad/1\nA\n+\n!\n@bad/2\n.\n+\n#\n",
            .message = "S002: sequence byte is outside the selected alphabet",
            .record = 3,
            .line = 2,
            .offset = 20,
        },
        .{
            .tail = "@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n",
            .message = "P001: paired identifiers or mate markers do not match (pair 1)",
            .record = 0,
            .line = 0,
            .offset = 0,
            .first_line = true,
        },
        .{
            .tail = "@odd/1\nA\n+\n!\n",
            .message = "P002: paired input is missing a mate " ++
                "(pair 1, remaining R1, last R1 record 2, last R2 record 1)",
            .record = 0,
            .line = 0,
            .offset = 0,
            .first_line = true,
        },
    };
    for (cases, 0..) |case, index| {
        input.shrinkRetainingCapacity(prefix_len);
        try input.appendSlice(allocator, case.tail);
        const expected_stderr = if (case.first_line)
            try std.fmt.allocPrint(allocator, "error: {s}: {s}\n", .{ input_path, case.message })
        else
            try std.fmt.allocPrint(
                allocator,
                "error: {s}: {s} (record {d}, line {d}, offset {d})\n",
                .{ input_path, case.message, case.record, case.line, prefix_len + case.offset },
            );
        const transports: []const bool = if (index == 0) &.{ false, true } else &.{false};
        for (transports) |gzip_input| {
            var gzip: std.ArrayList(u8) = .empty;
            if (gzip_input) try appendGzipChunks(allocator, &gzip, input.items);
            try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = if (gzip_input) gzip.items else input.items });
            const paths = try outputPaths(allocator, &tmp.sub_path, "prefix");
            const result = try runDeinterleave(allocator, input_path, paths);
            try std.testing.expectEqual(@as(u8, 1), result.exit_code);
            try std.testing.expectEqualStrings("", result.stdout);
            if (case.first_line) {
                try cli.expectFirstLine(expected_stderr, result.stderr);
            } else {
                try std.testing.expectEqualStrings(expected_stderr, result.stderr);
            }
            try expectFile(allocator, paths[0], r1.items);
            try expectFile(allocator, paths[1], r2.items);
            try tmp.dir.deleteFile(io, "prefix-r1.fastq");
            try tmp.dir.deleteFile(io, "prefix-r2.fastq");
        }
    }
}

test "[cli] - [deinterleave]: IUPAC alphabet accepts ambiguous bases in both mates" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");

    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@iupac/1\nR\n+\n!\n@iupac/2\nY\n+\n#\n",
    });
    const iupac_paths = try outputPaths(allocator, &tmp.sub_path, "iupac");
    try cli.expectResult(try runDeinterleave(allocator, input_path, iupac_paths), 0, "", "");
    try expectFile(allocator, iupac_paths[0], "@iupac/1\nR\n+\n!\n");
    try expectFile(allocator, iupac_paths[1], "@iupac/2\nY\n+\n#\n");
}

test "[cli] - [deinterleave]: ACGTN alphabet rejection leaves empty outputs" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");

    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@iupac/1\nR\n+\n!\n@iupac/2\nY\n+\n#\n",
    });

    const acgtn_paths = try outputPaths(allocator, &tmp.sub_path, "acgtn");
    const acgtn_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 9)\n",
        .{input_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--alphabet",
            "acgtn",
            "--out1",
            acgtn_paths[0],
            "--out2",
            acgtn_paths[1],
            input_path,
        }),
        1,
        "",
        acgtn_error,
    );
    try expectEmptyFile(acgtn_paths[0]);
    try expectEmptyFile(acgtn_paths[1]);
}

test "[cli] - [deinterleave]: line limit failures leave empty outputs" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");

    const line_cases = [_]struct {
        prefix: []const u8,
        input: []const u8,
    }{
        .{
            .prefix = "mate1-line",
            .input = "@ok/1\nAAAAAAAAA\n+\n!!!!!!!!!\n@ok/2\nT\n+\n#\n",
        },
        .{
            .prefix = "mate2-line",
            .input = "@ok/1\nA\n+\n!\n@ok/2\nAAAAAAAAA\n+\n!!!!!!!!!\n",
        },
    };
    const line_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: line length limit exceeded\n",
        .{input_path},
    );
    for (line_cases) |case| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = case.input });
        const line_paths = try outputPaths(allocator, &tmp.sub_path, case.prefix);
        try cli.expectResult(
            try cli.run(allocator, &.{
                "deinterleave",
                "--max-line-bytes",
                "8",
                "--out1",
                line_paths[0],
                "--out2",
                line_paths[1],
                input_path,
            }),
            4,
            "",
            line_error,
        );
        try expectEmptyFile(line_paths[0]);
        try expectEmptyFile(line_paths[1]);
    }
}

test "[cli] - [deinterleave]: an existing output is kept and the other output is not created" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const names = [_][]const u8{ "out1.fastq", "out2.fastq" };
    const kept = [_][]const u8{ "keep-one", "keep-two" };
    for (0..2) |existing| {
        var tmp = std.testing.tmpDir(.{});
        defer tmp.cleanup();
        const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
        const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
        const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");
        const paths = [_][]const u8{ out1_path, out2_path };
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n" });
        try tmp.dir.writeFile(io, .{ .sub_path = names[existing], .data = kept[existing] });

        const expected = try std.fmt.allocPrint(allocator, "error: {s}: output path already exists\n", .{paths[existing]});
        try cli.expectResult(try runDeinterleave(allocator, input_path, paths), 3, "", expected);
        try expectFile(allocator, paths[existing], kept[existing]);
        try std.testing.expect(!try pathExists(io, paths[1 - existing]));
    }
}

test "[cli] - [deinterleave]: an output path equal to the other output or to the input is refused" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");
    const payload = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });

    const alias_path = try cli.tempPath(allocator, &tmp.sub_path, "./out1.fastq");
    const alias_error = try std.fmt.allocPrint(allocator, "error: {s}: output path already exists\n", .{alias_path});
    try cli.expectResult(try runDeinterleave(allocator, input_path, .{ out1_path, alias_path }), 3, "", alias_error);
    try std.testing.expect(!try pathExists(io, out1_path));

    const input_error = try std.fmt.allocPrint(allocator, "error: {s}: output path already exists\n", .{input_path});
    try cli.expectResult(try runDeinterleave(allocator, input_path, .{ input_path, out2_path }), 3, "", input_error);
    try std.testing.expect(!try pathExists(io, out2_path));
    try cli.expectResult(try runDeinterleave(allocator, input_path, .{ out1_path, input_path }), 3, "", input_error);
    try std.testing.expect(!try pathExists(io, out1_path));
    try expectFile(allocator, input_path, payload);
}

test "[cli] - [deinterleave]: a symlink at an output path is refused and not followed" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "sentinel", .data = "keep-target" });
    try tmp.dir.symLink(io, "sentinel", "out1.fastq", .{});

    const expected = try std.fmt.allocPrint(allocator, "error: {s}: output path already exists\n", .{out1_path});
    try cli.expectResult(try runDeinterleave(allocator, input_path, .{ out1_path, out2_path }), 3, "", expected);
    const link_stat = try tmp.dir.statFile(io, "out1.fastq", .{ .follow_symlinks = false });
    try std.testing.expectEqual(std.Io.File.Kind.sym_link, link_stat.kind);
    try expectFile(allocator, try cli.tempPath(allocator, &tmp.sub_path, "sentinel"), "keep-target");
    try std.testing.expect(!try pathExists(io, out2_path));
}

test "[cli] - [deinterleave]: an output in a missing directory fails to create output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "missing/out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n" });

    const expected = try std.fmt.allocPrint(allocator, "error: {s}: failed to create output\n", .{out1_path});
    try cli.expectResult(try runDeinterleave(allocator, input_path, .{ out1_path, out2_path }), 3, "", expected);
    try std.testing.expect(!try pathExists(io, out2_path));
}

test "[edge] - [deinterleave]: a line exactly at the line limit is accepted" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");
    const r1 = "@p/1\nAAAAAAAA\n+\n!!!!!!!!\n";
    const r2 = "@p/2\nTTTTTTTT\n+\n########\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = r1 ++ r2 });
    try cli.expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--max-line-bytes", "8", "--out1", out1_path, "--out2", out2_path, input_path }),
        0,
        "",
        "",
    );
    try expectFile(allocator, out1_path, r1);
    try expectFile(allocator, out2_path, r2);
}

test "[cli] - [deinterleave]: failure does not delete replaced or created outputs" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const out1_path = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq");

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const child_io = threaded.io();
    var process = try std.process.spawn(child_io, .{
        .argv = &.{
            cli.ZFASTQ_BIN,
            "deinterleave",
            "--out1",
            out1_path,
            "--out2",
            out2_path,
            "-",
        },
        .stdin = .pipe,
        .stdout = .pipe,
        .stderr = .pipe,
    });
    defer process.kill(child_io);

    var process_deadline = try cli.ProcessDeadline.init(child_io, &process, cli.PROCESS_TIMEOUT);
    defer process_deadline.deinit();

    while (true) {
        const out1_exists = try pathExists(child_io, out1_path);
        const out2_exists = try pathExists(child_io, out2_path);
        if (out1_exists and out2_exists) break;
        if (process_deadline.expired()) return process_deadline.failTimeout();
        try std.Thread.yield();
    }

    try std.Io.Dir.cwd().deleteFile(child_io, out1_path);
    try std.Io.Dir.cwd().writeFile(child_io, .{ .sub_path = out1_path, .data = "replacement" });
    try std.Io.File.writeStreamingAll(
        process.stdin.?,
        child_io,
        "@bad/1\n.\n+\n!\n@bad/2\nA\n+\n!\n",
    );
    process.stdin.?.close(child_io);
    process.stdin = null;

    const result = try cli.finishSpawned(allocator, child_io, &process, &process_deadline);
    const expected_stderr = "error: -: S002: sequence byte is outside the selected alphabet " ++
        "(record 0, line 2, offset 7)\n";

    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(expected_stderr, result.stderr);
    try expectFile(allocator, out1_path, "replacement");
    try expectEmptyFile(out2_path);
}

test "[cli] - [deinterleave]: invalid arguments fail before output creation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct {
        args: []const []const u8,
        error_message: []const u8,
    }{
        .{
            .args = &.{ "deinterleave", "--out1", "one", "input" },
            .error_message = "error: deinterleave requires --out2\n",
        },
        .{
            .args = &.{ "deinterleave", "--out1", "one", "--out2" },
            .error_message = "error: --out2 requires a value\n",
        },
        .{
            .args = &.{ "deinterleave", "--out1", "-", "--out2", "two", "input" },
            .error_message = "error: deinterleave output paths cannot be standard output\n",
        },
        .{
            .args = &.{ "deinterleave", "--out1", "one", "--out2", "-", "input" },
            .error_message = "error: deinterleave output paths cannot be standard output\n",
        },
        .{
            .args = &.{ "deinterleave", "--out1", "one", "--out2", "two" },
            .error_message = "error: deinterleave requires exactly one input\n",
        },
    };
    for (cases) |case| {
        try cli.expectResult(try cli.run(allocator, case.args), 2, "", case.error_message);
    }
}

test "[cli] - [deinterleave]: missing input cannot be created as an output" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const missing = try cli.tempPath(allocator, &tmp.sub_path, "missing.fastq");
    const other = try cli.tempPath(allocator, &tmp.sub_path, "other.fastq");
    const missing_error = try std.fmt.allocPrint(allocator, "error: {s}: file not found\n", .{missing});
    try cli.expectResult(
        try runDeinterleave(allocator, missing, .{ missing, other }),
        3,
        "",
        missing_error,
    );
    try std.testing.expect(!try pathExists(std.testing.io, missing));
    try std.testing.expect(!try pathExists(std.testing.io, other));
}

test "[cli] - [deinterleave]: double dash treats option bytes as the input path" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const paths = try outputPaths(allocator, &tmp.sub_path, "double-dash");

    try cli.expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            paths[0],
            "--out2",
            paths[1],
            "--",
            "--alphabet",
        }),
        3,
        "",
        "error: --alphabet: file not found\n",
    );
    try std.testing.expect(!try pathExists(std.testing.io, paths[0]));
    try std.testing.expect(!try pathExists(std.testing.io, paths[1]));
}

test "[edge] - [deinterleave]: records spanning read-buffer refills on both sides remain exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const paths = try outputPaths(allocator, &tmp.sub_path, "large");
    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    var interleaved: std.ArrayList(u8) = .empty;
    const lengths = [_][2]usize{
        .{ 1, 1 },
        .{ READ_BUFFER_BYTES * 2, 1 },
        .{ READ_BUFFER_BYTES - 7, READ_BUFFER_BYTES },
    };
    for (lengths) |pair_lengths| {
        const start1 = r1.items.len;
        const start2 = r2.items.len;
        try appendFastqRecord(allocator, &r1, "large/1", 'A', '!', pair_lengths[0]);
        try appendFastqRecord(allocator, &r2, "large/2", 'T', '#', pair_lengths[1]);
        try interleaved.appendSlice(allocator, r1.items[start1..]);
        try interleaved.appendSlice(allocator, r2.items[start2..]);
    }
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = interleaved.items });

    try cli.expectResult(try runDeinterleave(allocator, input_path, paths), 0, "", "");
    try expectFile(allocator, paths[0], r1.items);
    try expectFile(allocator, paths[1], r2.items);

    var gzip: std.ArrayList(u8) = .empty;
    try appendGzipChunks(allocator, &gzip, interleaved.items);
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = gzip.items });
    const gzip_paths = try outputPaths(allocator, &tmp.sub_path, "gzip-large");
    try cli.expectResult(try runDeinterleave(allocator, input_path, gzip_paths), 0, "", "");
    try expectFile(allocator, gzip_paths[0], r1.items);
    try expectFile(allocator, gzip_paths[1], r2.items);
}

fn runDeinterleave(
    allocator: std.mem.Allocator,
    input_path: []const u8,
    paths: [2][]const u8,
) !cli.CommandResult {
    return cli.run(allocator, &.{ "deinterleave", "--out1", paths[0], "--out2", paths[1], input_path });
}

fn outputPaths(
    allocator: std.mem.Allocator,
    sub_path: []const u8,
    prefix: []const u8,
) ![2][]const u8 {
    return .{
        try cli.tempPath(allocator, sub_path, try std.fmt.allocPrint(allocator, "{s}-r1.fastq", .{prefix})),
        try cli.tempPath(allocator, sub_path, try std.fmt.allocPrint(allocator, "{s}-r2.fastq", .{prefix})),
    };
}

fn appendGzipChunks(allocator: std.mem.Allocator, gzip: *std.ArrayList(u8), bytes: []const u8) !void {
    var offset: usize = 0;
    while (offset < bytes.len) {
        const end = @min(offset + GZIP_MEMBER_INPUT_BYTES, bytes.len);
        try cli.appendGzipMember(allocator, gzip, bytes[offset..end], .{});
        offset = end;
    }
}

fn expectEmptyFile(path: []const u8) !void {
    const stat = try std.Io.Dir.cwd().statFile(std.testing.io, path, .{});
    try std.testing.expectEqual(std.Io.File.Kind.file, stat.kind);
    try std.testing.expectEqual(@as(u64, 0), stat.size);
}

fn expectFile(
    allocator: std.mem.Allocator,
    path: []const u8,
    expected: []const u8,
) !void {
    const actual = try std.Io.Dir.cwd().readFileAlloc(
        std.testing.io,
        path,
        allocator,
        .limited(4 * 1024 * 1024),
    );
    try std.testing.expectEqual(expected.len, actual.len);
    try std.testing.expectEqualSlices(u8, expected, actual);
}

fn pathExists(io: std.Io, path: []const u8) !bool {
    std.Io.Dir.cwd().access(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return false,
        else => return err,
    };
    return true;
}

fn appendFastqRecord(
    allocator: std.mem.Allocator,
    output: *std.ArrayList(u8),
    header: []const u8,
    base: u8,
    quality: u8,
    field_len: usize,
) !void {
    try output.append(allocator, '@');
    try output.appendSlice(allocator, header);
    try output.append(allocator, '\n');
    try output.appendNTimes(allocator, base, field_len);
    try output.appendSlice(allocator, "\n+\n");
    try output.appendNTimes(allocator, quality, field_len);
    try output.append(allocator, '\n');
}
