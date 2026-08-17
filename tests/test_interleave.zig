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
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");

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

    try expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        expected,
        "",
    );

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@same left\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@same right\nT\n+\n#\n" });
    try expectResult(
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
    try expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        "",
        "",
    );
}

test "[cli] - [interleave]: stdin and mixed plain or gzip inputs preserve output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const r1 = "@gzip/1\nAC\n+one\n!!\n";
    const r2 = "@gzip/2\nGT\n+two\n##\n";
    const expected = r1 ++ r2;

    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
    try expectResult(
        try cli.runWithStdin(allocator, &.{ "interleave", "-", r2_path }, r1, 1),
        0,
        expected,
        "",
    );
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try expectResult(
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
    try expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        0,
        expected,
        "",
    );
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    try expectResult(
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
    try expectResult(
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
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
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
    try expectResult(
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
    try expectResult(
        try cli.run(allocator, &.{ "interleave", r1_path, r2_path }),
        1,
        "",
        r2_error,
    );
}

test "[cli] - [interleave]: pair mismatches and unequal counts are exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");

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
    try expectResult(
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
    try expectResult(
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

    try expectResult(
        try cli.run(allocator, &.{"interleave"}),
        2,
        "",
        "error: interleave requires exactly two inputs\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "one.fastq" }),
        2,
        "",
        "error: interleave requires exactly two inputs\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "one", "two", "three" }),
        2,
        "",
        "error: interleave requires exactly two inputs\n",
    );
    try expectResult(
        try cli.runWithStdin(allocator, &.{ "interleave", "-", "-" }, "unused", 1),
        2,
        "",
        "error: interleave inputs may contain standard input at most once\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "--pair-names", "other", "a", "b" }),
        2,
        "",
        "error: --pair-names must be illumina or exact\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "--pair-names" }),
        2,
        "",
        "error: --pair-names requires a value\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "--alphabet", "dna", "a", "b" }),
        2,
        "",
        "error: --alphabet must be iupac or acgtn\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "--json", "a", "b" }),
        2,
        "",
        "error: unknown interleave option: --json\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "interleave", "--paired", "a", "b" }),
        2,
        "",
        "error: unknown interleave option: --paired\n",
    );
    try expectResult(
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
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const missing = try tempPath(allocator, &tmp.sub_path, "missing.fastq");

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@bad/1\n.\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@bad/2\nA\n+\n!\n" });
    const alphabet_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 7)\n",
        .{r1_path},
    );
    try expectResult(
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
    try expectResult(
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
    try expectResult(
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
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");

    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@small/2\nA\n+\n!\n" });
    try expectResult(
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
    try expectResult(
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

fn tempPath(
    allocator: std.mem.Allocator,
    sub_path: []const u8,
    name: []const u8,
) ![]const u8 {
    return std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}", .{ sub_path, name });
}
