//! Installed CLI contracts for `z-fastq interleave`.

const std = @import("std");
const cli = @import("utilities.zig");

const FIRST_PAIR =
    "@ok/1\nA\n+\n!\n" ++
    "@ok/2\nT\n+\n#\n";

test "[cli] - [interleave]: fields, order, line endings, and name policies are exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    const r1 =
        "@cluster 1:N:0:index-a\r\nAC\r\n+first annotation\r\n!~\r\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n";
    const r2 =
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/2 other\r\nA\r\n+mate\r\n!";
    const expected =
        "@cluster 1:N:0:index-a\nAC\n+first annotation\n!~\n" ++
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n" ++
        "@legacy/2 other\nA\n+mate\n!\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });

    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        expected,
        "",
    );

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@same left\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@same right\nT\n+\n#\n" });
    try cli.expectResult(
        try cli.run(allocator, &.{
            "interleave",
            "--pair-names",
            "exact",
            r1_path,
            r2_path,
        }),
        0,
        "@same left\nA\n+\n!\n@same right\nT\n+\n#\n",
        "",
    );

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "" });
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        "",
        "",
    );
}

test "[cli] - [interleave]: aliases reject before reading while distinct copies pass" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const copy = try cli.tempPath(allocator, &tmp.sub_path, "copy.fastq");
    const aliases = [_][]const u8{
        path,
        try cli.tempPath(allocator, &tmp.sub_path, "./input.fastq"),
        try cli.tempPath(allocator, &tmp.sub_path, "hard.fastq"),
        try cli.tempPath(allocator, &tmp.sub_path, "link.fastq"),
    };
    const payload = "@SRR1.1 1 length=4\nACGT\n+\nIIII\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });
    try tmp.dir.hardLink("input.fastq", tmp.dir, "hard.fastq", io, .{});
    try tmp.dir.symLink(io, "input.fastq", "link.fastq", .{});
    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, payload, .{});

    for ([_][]const u8{ payload, gzip.items, "", "invalid FASTQ" }) |bytes| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = bytes });
        try tmp.dir.writeFile(io, .{ .sub_path = "copy.fastq", .data = bytes });
        for (aliases) |alias| {
            const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{alias});
            try cli.expectResult(try cli.run(allocator, &.{ "interleave", path, alias }), 2, "", diagnostic);
        }
        for (0..2) |stdin_side| {
            const inputs: [2][]const u8 = if (stdin_side == 0) .{ "-", path } else .{ path, "-" };
            const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{inputs[1]});
            const file = try tmp.dir.openFile(io, "input.fastq", .{});
            defer file.close(io);
            try cli.expectResult(try cli.runWithStdinFile(allocator, &.{ "interleave", inputs[0], inputs[1] }, file), 2, "", diagnostic);
            if (bytes.len != 0) {
                var first: [1]u8 = undefined;
                try std.testing.expectEqual(@as(usize, 1), try file.readStreaming(io, &.{&first}));
                try std.testing.expectEqual(bytes[0], first[0]);
            }
        }
        if (std.mem.eql(u8, bytes, "invalid FASTQ")) continue;
        try cli.expectResult(try cli.run(allocator, &.{ "interleave", path, copy }), 0, if (bytes.len == 0) "" else payload ++ payload, "");
    }
}

test "[cli] - [interleave]: stdout aliases of either input reject before reading or writing" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const names = [_][]const u8{ "r1.fastq", "r2.fastq" };
    const paths = [_][]const u8{
        try cli.tempPath(allocator, &tmp.sub_path, names[0]),
        try cli.tempPath(allocator, &tmp.sub_path, names[1]),
    };
    const payload = "@same\nAC\n+\nII\n";
    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, payload, .{});
    for (0..2) |side| {
        try tmp.dir.writeFile(io, .{ .sub_path = names[side], .data = "" });
        try tmp.dir.hardLink(names[side], tmp.dir, "hard.fastq", io, .{});
        defer tmp.dir.deleteFile(io, "hard.fastq") catch {};
        try tmp.dir.symLink(io, names[side], "link.fastq", .{});
        defer tmp.dir.deleteFile(io, "link.fastq") catch {};
        const aliases = [_][]const u8{ names[side], "hard.fastq", "link.fastq" };
        for ([_][]const u8{ payload, gzip.items, "invalid FASTQ" }) |bytes| {
            for (aliases) |alias| {
                for ([_]bool{ false, true }) |append| {
                    for (names) |name| try tmp.dir.writeFile(io, .{ .sub_path = name, .data = bytes });
                    const output: std.Io.File = .{
                        .handle = try std.posix.openat(tmp.dir.handle, alias, .{
                            .ACCMODE = .WRONLY,
                            .TRUNC = !append,
                            .APPEND = append,
                            .CLOEXEC = true,
                        }, 0),
                        .flags = .{ .nonblocking = false },
                    };
                    defer output.close(io);
                    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: input file is output file\n", .{paths[side]});
                    try cli.expectResult(try cli.runWithStdoutFile(allocator, &.{ "interleave", paths[0], paths[1] }, output, null), 2, "", diagnostic);
                    for (names, 0..) |name, index| {
                        try std.testing.expectEqualStrings(if (!append and index == side) "" else bytes, try tmp.dir.readFileAlloc(io, name, allocator, .limited(1024)));
                    }
                    for (names) |stderr_name| {
                        const stderr_file: std.Io.File = .{
                            .handle = try std.posix.openat(tmp.dir.handle, stderr_name, .{
                                .ACCMODE = .WRONLY,
                                .APPEND = true,
                                .CLOEXEC = true,
                            }, 0),
                            .flags = .{ .nonblocking = false },
                        };
                        defer stderr_file.close(io);
                        try std.testing.expectEqual(@as(u8, 2), try cli.runWithOutputFiles(allocator, &.{ "interleave", paths[0], paths[1] }, output, stderr_file));
                        for (names, 0..) |name, index| {
                            try std.testing.expectEqualStrings(if (!append and index == side) "" else bytes, try tmp.dir.readFileAlloc(io, name, allocator, .limited(1024)));
                        }
                    }
                }
            }
        }
        for ([_]bool{ false, true }) |stdin_alias| {
            for (names) |name| try tmp.dir.writeFile(io, .{ .sub_path = name, .data = payload });
            const input = try tmp.dir.openFile(io, names[if (stdin_alias) side else 1 - side], .{});
            defer input.close(io);
            const output = try tmp.dir.openFile(io, names[side], .{ .mode = .write_only });
            defer output.close(io);
            var inputs = paths;
            inputs[if (stdin_alias) side else 1 - side] = "-";
            const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: input file is output file\n", .{inputs[side]});
            try cli.expectResult(try cli.runWithStdoutFile(allocator, &.{ "interleave", inputs[0], inputs[1] }, output, input), 2, "", diagnostic);
            var first: [1]u8 = undefined;
            try std.testing.expectEqual(@as(usize, 1), try input.readStreaming(io, &.{&first}));
            try std.testing.expectEqual(@as(u8, '@'), first[0]);
        }
    }
}

test "[cli] - [interleave]: shared input and output aliases preserve the paired-input error" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const payload = "@same\nAC\n+\nII\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });
    const output: std.Io.File = .{
        .handle = try std.posix.openat(tmp.dir.handle, "input.fastq", .{
            .ACCMODE = .WRONLY,
            .APPEND = true,
            .CLOEXEC = true,
        }, 0),
        .flags = .{ .nonblocking = false },
    };
    defer output.close(io);
    const args = [_][]const u8{ "interleave", path, path };
    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{path});
    try cli.expectResult(try cli.runWithStdoutFile(allocator, &args, output, null), 2, "", diagnostic);
    try std.testing.expectEqual(@as(u8, 2), try cli.runWithOutputFiles(allocator, &args, output, output));
    try std.testing.expectEqualStrings(payload, try tmp.dir.readFileAlloc(io, "input.fastq", allocator, .limited(1024)));
}

test "[cli] - [paired FIFO inputs]: a sequential writer completes without losing bytes" {
    const Producer = struct {
        io: std.Io,
        dir: std.Io.Dir,
        bytes: []const u8,
        err: ?anyerror = null,

        fn run(self: *@This()) void {
            self.write() catch |err| {
                self.err = err;
            };
        }

        fn write(self: *@This()) !void {
            for ([_][]const u8{ "r1", "r2" }) |name| {
                const file = try self.dir.openFile(self.io, name, .{ .mode = .write_only });
                defer file.close(self.io);
                // The first write must exceed the pipe capacity to exercise backpressure.
                try std.testing.expectEqual(.SUCCESS, std.os.linux.errno(
                    std.os.linux.fcntl(file.handle, std.os.linux.F.SETPIPE_SZ, 4096),
                ));
                try file.writeStreamingAll(self.io, self.bytes);
            }
        }
    };

    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const paths = [_][]const u8{
        try cli.tempPath(allocator, &tmp.sub_path, "r1"),
        try cli.tempPath(allocator, &tmp.sub_path, "r2"),
    };
    for ([_][:0]const u8{ "r1", "r2" }) |name| {
        try std.testing.expectEqual(.SUCCESS, std.os.linux.errno(std.os.linux.mknodat(
            tmp.dir.handle,
            name,
            std.os.linux.S.IFIFO | 0o600,
            0,
        )));
    }
    const payload = "@same\n" ++ "A" ** 3000 ++ "\n+\n" ++ "I" ** 3000 ++ "\n";
    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, payload, .{});
    for ([_][]const u8{ payload, gzip.items, "" }) |bytes| {
        var producer: Producer = .{ .io = io, .dir = tmp.dir, .bytes = bytes };
        var group: std.Io.Group = .init;
        defer group.cancel(io);
        try group.concurrent(io, Producer.run, .{&producer});
        try cli.expectResult(
            try cli.run(allocator, &.{ "interleave", paths[0], paths[1] }),
            0,
            if (bytes.len == 0) "" else payload ++ payload,
            "",
        );
        try group.await(io);
        if (producer.err) |err| return err;
    }
}

test "[cli] - [interleave]: terminal CR fields reject the pair before either mate is written" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const paths = [_][]const u8{
        try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq"),
        try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq"),
    };
    for (0..2) |bad_mate| {
        const r1 = if (bad_mate == 0) "@pair/1 x\r\r\nA\n+\n!\n" else "@pair/1\nA\n+\n!\n";
        const r2 = if (bad_mate == 1) "@pair/2\nT\n+\r\r\n#\n" else "@pair/2\nT\n+\n#\n";
        try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
        try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
        const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: record fields ending in CR cannot be written with LF endings\n", .{paths[bad_mate]});
        try cli.expectResult(try cli.run(allocator, &.{ "interleave", paths[0], paths[1] }), 1, "", diagnostic);
    }
}

test "[cli] - [interleave]: canonical and refill-spanning mates preserve exact output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const r1 = "@refill/1\nA\n+left\n!\n";

    var r2: std.ArrayList(u8) = .empty;
    try r2.appendSlice(allocator, "@refill/2 ");
    try r2.appendNTimes(allocator, 'x', 256 * 1024);
    try r2.appendSlice(allocator, "\nT\n+right\n#\n");
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2.items });

    const result = try cli.run(allocator, &.{ "interleave", r1_path, r2_path });
    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    try std.testing.expectEqual(r1.len + r2.items.len, result.stdout.len);
    try std.testing.expectEqualStrings(r1, result.stdout[0..r1.len]);
    try std.testing.expectEqualSlices(u8, r2.items, result.stdout[r1.len..]);
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
}

test "[cli] - [interleave]: stdin and mixed plain or gzip inputs preserve output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const r1 = "@gzip/1\nAC\n+one\n!!\n";
    const r2 = "@gzip/2\nGT\n+two\n##\n";
    const expected = r1 ++ r2;

    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
    try cli.expectResult(
        try cli.runWithStdin(allocator, &.{ "interleave", "-", r2_path }, r1, 1),
        0,
        expected,
        "",
    );
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try cli.expectResult(
        try cli.runWithStdin(allocator, &.{ "interleave", r1_path, "-" }, r2, 2),
        0,
        expected,
        "",
    );

    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, r1, .{});
    try cli.appendGzipMember(allocator, &gzip2, r2, .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        expected,
        "",
    );
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        expected,
        "",
    );

    gzip2.items[gzip2.items.len - 8] ^= 1;
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    const error_text = try std.fmt.allocPrint(
        allocator,
        "error: {s}: I/O error\n",
        .{r2_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        3,
        "",
        error_text,
    );
}

test "[cli] - [interleave]: validation precedence protects the failing pair" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@ok/1\nA\n+\n!\n@bad/1\n.\n+\n!\n",
    });
    try tmp.dir.writeFile(io, .{
        .sub_path = "r2.fastq",
        .data = "@ok/2\nT\n+\n#\n@other/2\n.\n+\n!\n",
    });
    const r1_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 1, line 2, offset 19)\n",
        .{r1_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        1,
        FIRST_PAIR,
        r1_error,
    );

    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@bad/1\n.\n+\n!\n",
    });
    try tmp.dir.writeFile(io, .{
        .sub_path = "r2.fastq",
        .data = "@bad/2\nA\nx\n!\n",
    });
    const r2_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S001: plus line must start with '+' " ++
            "(record 0, line 3, offset 9)\n",
        .{r2_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        1,
        "",
        r2_error,
    );
}

test "[cli] - [interleave]: later failure preserves complete pairs beyond the output buffer" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    var expected: std.ArrayList(u8) = .empty;
    for (0..4096) |_| {
        try r1.appendSlice(allocator, "@ok/1\nA\n+\n!\n");
        try r2.appendSlice(allocator, "@ok/2\nT\n+\n#\n");
        try expected.appendSlice(allocator, FIRST_PAIR);
    }
    try std.testing.expect(expected.items.len > 64 * 1024);
    const offset = r2.items.len + 7;
    try r1.appendSlice(allocator, "@bad/1\nA\n+\n!\n");
    try r2.appendSlice(allocator, "@bad/2\n.\n+\n#\n");
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2.items });
    const expected_stderr = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 4096, line 2, offset {d})\n",
        .{ r2_path, offset },
    );

    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        1,
        expected.items,
        expected_stderr,
    );
}

test "[cli] - [interleave]: pair mismatches and unequal counts are exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@left/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@right/2\nT\n+\n#\n" });
    const mismatch = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P001: paired identifiers or mate markers do not match " ++
            "(pair 0)\n" ++
            "  R1: input={s}, record=0, offset=0, first_token=left/1 " ++
            "[length=6, truncated=false], normalized_id=left " ++
            "[length=4, truncated=false], mate_markers=1\n" ++
            "  R2: input={s}, record=0, offset=0, first_token=right/2 " ++
            "[length=7, truncated=false], normalized_id=right " ++
            "[length=5, truncated=false], mate_markers=2\n",
        .{ r1_path, r2_path, r1_path, r2_path },
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        1,
        "",
        mismatch,
    );

    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@ok/1\nA\n+\n!\n@extra/1\nA\n+\n!\n",
    });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@ok/2\nT\n+\n#\n" });
    const unequal = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 1, last R2 record 0)\n",
        .{ r1_path, r2_path },
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        1,
        FIRST_PAIR,
        unequal,
    );
}

test "[cli] - [interleave]: arguments and limits fail before unsafe output" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try cli.expectResult(
        try cli.run(allocator, &.{"interleave"}),
        2,
        "",
        "error: interleave requires exactly two inputs\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "one.fastq" }),
        2,
        "",
        "error: interleave requires exactly two inputs\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "one", "two", "three" }),
        2,
        "",
        "error: interleave requires exactly two inputs\n",
    );
    try cli.expectResult(
        try cli.runWithStdin(allocator, &.{ "interleave", "-", "-" }, "unused", 1),
        2,
        "",
        "error: interleave inputs may contain standard input at most once\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--pair-names", "other", "a", "b" }),
        2,
        "",
        "error: --pair-names must be illumina or exact\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--pair-names" }),
        2,
        "",
        "error: --pair-names requires a value\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--alphabet", "dna", "a", "b" }),
        2,
        "",
        "error: --alphabet must be iupac or acgtn\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--json", "a", "b" }),
        2,
        "",
        "error: unknown interleave option: --json\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--paired", "a", "b" }),
        2,
        "",
        "error: unknown interleave option: --paired\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--", "--alphabet", "other" }),
        3,
        "",
        "error: --alphabet: file not found\n",
    );
}

test "[cli] - [interleave]: input, alphabet, and line limits identify their side" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const missing = try cli.tempPath(allocator, &tmp.sub_path, "missing.fastq");

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@bad/1\n.\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@bad/2\nA\n+\n!\n" });
    const alphabet_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 7)\n",
        .{r1_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{
            "interleave",
            "--alphabet",
            "acgtn",
            r1_path,
            r2_path,
        }),
        1,
        "",
        alphabet_error,
    );

    const line_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: line length limit exceeded\n",
        .{r1_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{
            "interleave",
            "--max-line-bytes",
            "4",
            r1_path,
            r2_path,
        }),
        4,
        "",
        line_error,
    );

    const missing_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: file not found\n",
        .{missing},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, missing }),
        3,
        "",
        missing_error,
    );
}

test "[cli] - [interleave]: write and flush failures return output I/O status" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@small/2\nA\n+\n!\n" });
    try cli.expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "interleave", "-", r2_path },
            "@small/1\nA\n+\n!\n",
        ),
        3,
        "",
        "error: standard output: I/O error\n",
    );

    const field_len = 70 * 1024;
    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    try r1.appendSlice(allocator, "@large/1\n");
    try r1.appendNTimes(allocator, 'A', field_len);
    try r1.appendSlice(allocator, "\n+\n");
    try r1.appendNTimes(allocator, '!', field_len);
    try r1.append(allocator, '\n');
    try r2.appendSlice(allocator, "@large/2\n");
    try r2.appendNTimes(allocator, 'T', field_len);
    try r2.appendSlice(allocator, "\n+\n");
    try r2.appendNTimes(allocator, '#', field_len);
    try r2.append(allocator, '\n');
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2.items });
    try cli.expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "interleave", "-", r2_path },
            r1.items,
        ),
        3,
        "",
        "error: standard output: I/O error\n",
    );
}
