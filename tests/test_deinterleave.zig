//! Installed CLI contracts for `z-fastq deinterleave`.

const std = @import("std");
const cli = @import("utilities.zig");

test "[cli] - [deinterleave]: fields, order, line endings, and round trip are exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try tempPath(allocator, &tmp.sub_path, "out2.fastq");

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
    const expected_interleaved =
        "@cluster 1:N:0:index-a\nAC\n+first annotation\n!~\n" ++
        "@cluster 2:Y:0:index-b\nGT\n+second annotation\n#$\n" ++
        "@legacy/1 opaque\nN\n+legacy\n#\n" ++
        "@legacy/2 other\nA\n+mate\n!\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = input });

    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            out1_path,
            "--out2",
            out2_path,
            input_path,
        }),
        0,
        "",
        "",
    );
    try expectFile(allocator, out1_path, expected1);
    try expectFile(allocator, out2_path, expected2);
    try expectResult(
        try cli.run(allocator, &.{ "check", "--paired", out1_path, out2_path }),
        0,
        "",
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", out1_path, out2_path }),
        0,
        expected_interleaved,
        "",
    );
}

test "[cli] - [deinterleave]: exact policy, empty input, stdin, and gzip are supported" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq.gz");
    const payload =
        "@same left\nA\n+one\n!\n@same right\nT\n+two\n#\n" ++
        "@empty one\n\n+empty-one\n\n@empty two\n\n+empty-two\n\n";
    const expected1 =
        "@same left\nA\n+one\n!\n" ++
        "@empty one\n\n+empty-one\n\n";
    const expected2 =
        "@same right\nT\n+two\n#\n" ++
        "@empty two\n\n+empty-two\n\n";

    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, payload, .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq.gz", .data = gzip.items });

    const gzip_out1 = try tempPath(allocator, &tmp.sub_path, "gzip-r1.fastq");
    const gzip_out2 = try tempPath(allocator, &tmp.sub_path, "gzip-r2.fastq");
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--pair-names",
            "exact",
            "--alphabet",
            "acgtn",
            "--out1",
            gzip_out1,
            "--out2",
            gzip_out2,
            input_path,
        }),
        0,
        "",
        "",
    );
    try expectFile(allocator, gzip_out1, expected1);
    try expectFile(allocator, gzip_out2, expected2);

    const stdin_out1 = try tempPath(allocator, &tmp.sub_path, "stdin-r1.fastq");
    const stdin_out2 = try tempPath(allocator, &tmp.sub_path, "stdin-r2.fastq");
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{
                "deinterleave",
                "--pair-names",
                "exact",
                "--alphabet",
                "acgtn",
                "--out1",
                stdin_out1,
                "--out2",
                stdin_out2,
                "-",
            },
            payload,
            1,
        ),
        0,
        "",
        "",
    );
    try expectFile(allocator, stdin_out1, expected1);
    try expectFile(allocator, stdin_out2, expected2);

    const empty_out1 = try tempPath(allocator, &tmp.sub_path, "empty-r1.fastq");
    const empty_out2 = try tempPath(allocator, &tmp.sub_path, "empty-r2.fastq");
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{
                "deinterleave",
                "--out1",
                empty_out1,
                "--out2",
                empty_out2,
                "-",
            },
            "",
            1,
        ),
        0,
        "",
        "",
    );
    try expectFile(allocator, empty_out1, "");
    try expectFile(allocator, empty_out2, "");
}

test "[cli] - [deinterleave]: structural, semantic, pair, and odd-count precedence is exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq");

    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@bad/1\n.\n+\n!\n@bad/2\nA\nx\n!\n",
    });
    var paths = try outputPaths(allocator, &tmp.sub_path, "structure");
    const structural_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S001: plus line must start with '+' " ++
            "(record 1, line 3, offset 22)\n",
        .{input_path},
    );
    try expectResult(
        try runDeinterleave(allocator, input_path, paths),
        1,
        "",
        structural_error,
    );
    try expectAbsent(paths[0]);
    try expectAbsent(paths[1]);

    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n",
    });
    paths = try outputPaths(allocator, &tmp.sub_path, "mismatch");
    const mismatch = try std.fmt.allocPrint(
        allocator,
        "error: {s}: P001: paired identifiers or mate markers do not match (pair 0)\n" ++
            "  R1: input={s}, record=0, offset=0, first_token=left/1 " ++
            "[length=6, truncated=false], normalized_id=left " ++
            "[length=4, truncated=false], mate_markers=1\n" ++
            "  R2: input={s}, record=1, offset=14, first_token=right/2 " ++
            "[length=7, truncated=false], normalized_id=right " ++
            "[length=5, truncated=false], mate_markers=2\n",
        .{ input_path, input_path, input_path },
    );
    try expectResult(
        try runDeinterleave(allocator, input_path, paths),
        1,
        "",
        mismatch,
    );
    try expectAbsent(paths[0]);
    try expectAbsent(paths[1]);

    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@odd/1\nA\n+\n!\n",
    });
    paths = try outputPaths(allocator, &tmp.sub_path, "odd");
    const odd = try std.fmt.allocPrint(
        allocator,
        "error: {s}: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 2, last R2 record 1)\n",
        .{input_path},
    );
    try expectResult(try runDeinterleave(allocator, input_path, paths), 1, "", odd);
    try expectAbsent(paths[0]);
    try expectAbsent(paths[1]);
}

test "[cli] - [deinterleave]: failed input removes both outputs" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq.gz");
    const payload = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, payload, .{});
    gzip.items[gzip.items.len - 8] ^= 1;
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq.gz", .data = gzip.items });
    const paths = try outputPaths(allocator, &tmp.sub_path, "corrupt");
    const expected = try std.fmt.allocPrint(allocator, "error: {s}: I/O error\n", .{input_path});

    try expectResult(
        try runDeinterleave(allocator, input_path, paths),
        3,
        "",
        expected,
    );
    try expectAbsent(paths[0]);
    try expectAbsent(paths[1]);
}

test "[cli] - [deinterleave]: validation failures on either mate remove both outputs" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq");

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
            .prefix = "final-structure",
            .input = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@odd/1\nA\n+\n",
            .exit_code = 1,
            .error_tail = "S004: unexpected end of file in quality line " ++
                "(record 2, line 4, offset 35)\n",
        },
    };
    for (cases) |case| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = case.input });
        const paths = try outputPaths(allocator, &tmp.sub_path, case.prefix);
        const expected = try std.fmt.allocPrint(
            allocator,
            "error: {s}: {s}",
            .{ input_path, case.error_tail },
        );
        try expectResult(
            try runDeinterleave(allocator, input_path, paths),
            case.exit_code,
            "",
            expected,
        );
        try expectAbsent(paths[0]);
        try expectAbsent(paths[1]);
    }

    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@ok/1\nA\n+\n!\n@ok/2\nAAAAAAAAA\n+\n!!!!!!!!!\n",
    });
    const line_paths = try outputPaths(allocator, &tmp.sub_path, "mate2-line");
    const line_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: line length limit exceeded\n",
        .{input_path},
    );
    try expectResult(
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
    try expectAbsent(line_paths[0]);
    try expectAbsent(line_paths[1]);
}

test "[cli] - [deinterleave]: output creation is exclusive and cleanup preserves prior paths" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq");
    const out1_path = try tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try tempPath(allocator, &tmp.sub_path, "out2.fastq");
    try tmp.dir.writeFile(io, .{
        .sub_path = "input.fastq",
        .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n",
    });

    try tmp.dir.writeFile(io, .{ .sub_path = "out1.fastq", .data = "keep-one" });
    const first_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: output path already exists\n",
        .{out1_path},
    );
    try expectResult(
        try runDeinterleave(allocator, input_path, .{ out1_path, out2_path }),
        3,
        "",
        first_error,
    );
    try expectFile(allocator, out1_path, "keep-one");
    try expectAbsent(out2_path);

    try tmp.dir.deleteFile(io, "out1.fastq");
    try tmp.dir.writeFile(io, .{ .sub_path = "out2.fastq", .data = "keep-two" });
    const second_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: output path already exists\n",
        .{out2_path},
    );
    try expectResult(
        try runDeinterleave(allocator, input_path, .{ out1_path, out2_path }),
        3,
        "",
        second_error,
    );
    try expectAbsent(out1_path);
    try expectFile(allocator, out2_path, "keep-two");

    try tmp.dir.deleteFile(io, "out2.fastq");
    const alias_path = try tempPath(allocator, &tmp.sub_path, "./out1.fastq");
    const alias_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: output path already exists\n",
        .{alias_path},
    );
    try expectResult(
        try runDeinterleave(allocator, input_path, .{ out1_path, alias_path }),
        3,
        "",
        alias_error,
    );
    try expectAbsent(out1_path);

    const input_alias_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: output path already exists\n",
        .{input_path},
    );
    try expectResult(
        try runDeinterleave(allocator, input_path, .{ input_path, out2_path }),
        3,
        "",
        input_alias_error,
    );
    try expectFile(
        allocator,
        input_path,
        "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n",
    );
    try expectAbsent(out2_path);

    try expectResult(
        try runDeinterleave(allocator, input_path, .{ out1_path, input_path }),
        3,
        "",
        input_alias_error,
    );
    try expectAbsent(out1_path);
    try expectFile(
        allocator,
        input_path,
        "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n",
    );

    try tmp.dir.writeFile(io, .{ .sub_path = "sentinel", .data = "keep-target" });
    try tmp.dir.symLink(io, "sentinel", "out1.fastq", .{});
    try expectResult(
        try runDeinterleave(allocator, input_path, .{ out1_path, out2_path }),
        3,
        "",
        first_error,
    );
    const link_stat = try tmp.dir.statFile(io, "out1.fastq", .{ .follow_symlinks = false });
    try std.testing.expectEqual(std.Io.File.Kind.sym_link, link_stat.kind);
    try expectFile(
        allocator,
        try tempPath(allocator, &tmp.sub_path, "sentinel"),
        "keep-target",
    );
    try expectAbsent(out2_path);
}

test "[cli] - [deinterleave]: replacement cleanup keeps diagnostic order and exit precedence" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const out1_path = try tempPath(allocator, &tmp.sub_path, "out1.fastq");
    const out2_path = try tempPath(allocator, &tmp.sub_path, "out2.fastq");

    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const child_io = threaded.io();
    var process = try std.process.spawn(child_io, .{
        .argv = &.{
            "zig-out/bin/z-fastq",
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

    var outputs_created = false;
    for (0..1000) |_| {
        const out1_exists = try pathExists(child_io, out1_path);
        const out2_exists = try pathExists(child_io, out2_path);
        if (out1_exists and out2_exists) {
            outputs_created = true;
            break;
        }
        try child_io.sleep(.fromMilliseconds(1), .awake);
    }
    try std.testing.expect(outputs_created);

    try std.Io.Dir.cwd().deleteFile(child_io, out1_path);
    try std.Io.Dir.cwd().writeFile(child_io, .{
        .sub_path = out1_path,
        .data = "replacement",
    });
    try std.Io.File.writeStreamingAll(
        process.stdin.?,
        child_io,
        "@bad/1\n.\n+\n!\n@bad/2\nA\n+\n!\n",
    );
    process.stdin.?.close(child_io);
    process.stdin = null;

    var stdout_buffer: [256]u8 = undefined;
    var stdout_reader = process.stdout.?.reader(child_io, &stdout_buffer);
    const stdout = try stdout_reader.interface.allocRemaining(allocator, .limited(4096));
    var stderr_buffer: [4096]u8 = undefined;
    var stderr_reader = process.stderr.?.reader(child_io, &stderr_buffer);
    const stderr = try stderr_reader.interface.allocRemaining(allocator, .limited(4096));
    const termination = try process.wait(child_io);
    const exit_code = switch (termination) {
        .exited => |code| code,
        else => return error.ChildProcessFailed,
    };
    const expected_stderr = try std.fmt.allocPrint(
        allocator,
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 7)\n" ++
            "error: {s}: cleanup failed: output path was replaced\n",
        .{out1_path},
    );

    try std.testing.expectEqual(@as(u8, 3), exit_code);
    try std.testing.expectEqual(@as(usize, 0), stdout.len);
    try std.testing.expectEqualStrings(expected_stderr, stderr);
    try expectFile(allocator, out1_path, "replacement");
    try expectAbsent(out2_path);
}

test "[cli] - [deinterleave]: argument and staging limits fail before output creation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectResult(
        try cli.run(allocator, &.{ "deinterleave", "input" }),
        2,
        "",
        "error: deinterleave requires --out1\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--out1", "one", "input" }),
        2,
        "",
        "error: deinterleave requires --out2\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            "one",
            "--out1",
            "two",
            "--out2",
            "three",
            "input",
        }),
        2,
        "",
        "error: --out1 may appear only once\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            "one",
            "--out2",
            "two",
            "--out2",
            "three",
            "input",
        }),
        2,
        "",
        "error: --out2 may appear only once\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--out1" }),
        2,
        "",
        "error: --out1 requires a value\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--out1", "one", "--out2" }),
        2,
        "",
        "error: --out2 requires a value\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            "-",
            "--out2",
            "two",
            "input",
        }),
        2,
        "",
        "error: deinterleave output paths cannot be standard output\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            "same",
            "--out2",
            "same",
            "input",
        }),
        2,
        "",
        "error: deinterleave output paths must differ\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            "one",
            "--out2",
            "two",
        }),
        2,
        "",
        "error: deinterleave requires exactly one input\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--out1",
            "one",
            "--out2",
            "two",
            "input",
            "extra",
        }),
        2,
        "",
        "error: deinterleave requires exactly one input\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "deinterleave", "--unknown" }),
        2,
        "",
        "error: unknown deinterleave option: --unknown\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "deinterleave",
            "--max-line-bytes",
            "18446744073709551615",
            "--out1",
            "one",
            "--out2",
            "two",
            "missing",
        }),
        4,
        "",
        "error: deinterleave record staging size exceeds supported limit\n",
    );

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const missing = try tempPath(allocator, &tmp.sub_path, "missing.fastq");
    const other = try tempPath(allocator, &tmp.sub_path, "other.fastq");
    const missing_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: file not found\n",
        .{missing},
    );
    try expectResult(
        try runDeinterleave(allocator, missing, .{ missing, other }),
        3,
        "",
        missing_error,
    );
    try expectAbsent(missing);
    try expectAbsent(other);
}

test "[cli] - [deinterleave]: refill-spanning records on both sides remain exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input_path = try tempPath(allocator, &tmp.sub_path, "input.fastq");
    const paths = try outputPaths(allocator, &tmp.sub_path, "large");
    const field_len = 512 * 1024;

    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    try appendRecord(allocator, &r1, "large-a/1", 'A', '!', field_len);
    try appendRecord(allocator, &r2, "large-a/2", 'T', '#', field_len);
    try appendRecord(allocator, &r1, "large-b/1", 'C', '$', 1);
    try appendRecord(allocator, &r2, "large-b/2", 'G', '%', 1);

    var interleaved: std.ArrayList(u8) = .empty;
    const first_record_len = field_len * 2 + "@large-a/1\n\n+\n\n".len;
    try interleaved.appendSlice(allocator, r1.items[0..first_record_len]);
    try interleaved.appendSlice(allocator, r2.items[0..first_record_len]);
    try interleaved.appendSlice(allocator, r1.items[first_record_len..]);
    try interleaved.appendSlice(allocator, r2.items[first_record_len..]);
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = interleaved.items });

    try expectResult(try runDeinterleave(allocator, input_path, paths), 0, "", "");
    try expectFile(allocator, paths[0], r1.items);
    try expectFile(allocator, paths[1], r2.items);
}

fn runDeinterleave(
    allocator: std.mem.Allocator,
    input_path: []const u8,
    paths: [2][]const u8,
) !cli.CommandResult {
    return cli.run(allocator, &.{
        "deinterleave",
        "--out1",
        paths[0],
        "--out2",
        paths[1],
        input_path,
    });
}

fn outputPaths(
    allocator: std.mem.Allocator,
    sub_path: []const u8,
    prefix: []const u8,
) ![2][]const u8 {
    return .{
        try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}-r1.fastq", .{
            sub_path,
            prefix,
        }),
        try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}-r2.fastq", .{
            sub_path,
            prefix,
        }),
    };
}

fn expectResult(
    result: cli.CommandResult,
    exit_code: u8,
    stdout: []const u8,
    stderr: []const u8,
) !void {
    try std.testing.expectEqual(exit_code, result.exit_code);
    try std.testing.expectEqualStrings(stdout, result.stdout);
    try std.testing.expectEqualStrings(stderr, result.stderr);
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
    try std.testing.expectEqualSlices(u8, expected, actual);
}

fn expectAbsent(path: []const u8) !void {
    try std.testing.expectError(
        error.FileNotFound,
        std.Io.Dir.cwd().access(std.testing.io, path, .{}),
    );
}

fn pathExists(io: std.Io, path: []const u8) !bool {
    std.Io.Dir.cwd().access(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return false,
        else => return err,
    };
    return true;
}

fn appendRecord(
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

fn tempPath(
    allocator: std.mem.Allocator,
    sub_path: []const u8,
    name: []const u8,
) ![]const u8 {
    return std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}", .{ sub_path, name });
}
