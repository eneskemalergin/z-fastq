//! Installed CLI contracts for `z-fastq check`.

const std = @import("std");
const cli = @import("utilities.zig");

const FIXTURE_DIR = "tests/data/synthetic";

const FixtureExpect = struct {
    name: []const u8,
    code: ?[]const u8 = null,
    message: []const u8 = "",
    record_index: u64 = 0,
    line: u3 = 0,
    offset: u64 = 0,
};

const FIXTURES = [_]FixtureExpect{
    .{ .name = "acgtn_valid.fastq" },
    .{ .name = "basic_valid.fastq" },
    .{ .name = "crlf.fastq" },
    .{ .name = "missing_final_newline.fastq" },
    .{ .name = "empty_valid.fastq" },
    .{ .name = "iupac_valid.fastq" },
    .{
        .name = "bad_plus.fastq",
        .code = "S001",
        .message = "plus line must start with '+'",
        .line = 3,
        .offset = 15,
    },
    .{
        .name = "bad_alphabet.fastq",
        .code = "S002",
        .message = "sequence byte is outside the selected alphabet",
        .line = 2,
        .offset = 16,
    },
    .{
        .name = "bad_header.fastq",
        .code = "S003",
        .message = "header line must start with '@'",
        .line = 1,
        .offset = 0,
    },
    .{
        .name = "truncated_record.fastq",
        .code = "S004",
        .message = "unexpected end of file in quality line",
        .line = 4,
        .offset = 17,
    },
    .{
        .name = "bad_qual_length.fastq",
        .code = "S005",
        .message = "sequence and quality lengths differ",
        .line = 4,
        .offset = 19,
    },
    .{
        .name = "bad_quality_range.fastq",
        .code = "S006",
        .message = "quality byte must be ASCII 33 through 126",
        .line = 4,
        .offset = 22,
    },
};

test "[cli] - [check]: files and fragmented stdin produce exact fixture results" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (FIXTURES) |fixture| {
        const path = try std.fmt.allocPrint(allocator, "{s}/{s}", .{
            FIXTURE_DIR,
            fixture.name,
        });
        const data = try std.Io.Dir.cwd().readFileAlloc(
            std.testing.io,
            path,
            allocator,
            .limited(1024 * 1024),
        );
        const file_result = try cli.run(allocator, &.{ "check", path });
        const stdin_result = try cli.runWithStdin(allocator, &.{ "check", "-" }, data, 1);

        if (fixture.code) |code| {
            const file_error = try expectedError(
                allocator,
                path,
                code,
                fixture.message,
                fixture.record_index,
                fixture.line,
                fixture.offset,
            );
            const stdin_error = try expectedError(
                allocator,
                "-",
                code,
                fixture.message,
                fixture.record_index,
                fixture.line,
                fixture.offset,
            );
            try expectResult(file_result, 1, "", file_error);
            try expectResult(stdin_result, 1, "", stdin_error);
        } else {
            try expectResult(file_result, 0, "", "");
            try expectResult(stdin_result, 0, "", "");
        }
    }

    const empty_input = try cli.run(allocator, &.{ "check", "-" });
    try expectResult(empty_input, 0, "", "");
}

test "[cli] - [check]: alphabet policy and semantic precedence are exact" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const explicit_iupac = try cli.run(allocator, &.{
        "check",
        "--alphabet",
        "iupac",
        "tests/data/synthetic/iupac_valid.fastq",
    });
    try expectResult(explicit_iupac, 0, "", "");

    const narrow_valid = try cli.run(allocator, &.{
        "check",
        "--alphabet",
        "acgtn",
        "tests/data/synthetic/acgtn_valid.fastq",
    });
    try expectResult(narrow_valid, 0, "", "");

    const narrow = try cli.run(allocator, &.{
        "check",
        "--alphabet",
        "acgtn",
        "tests/data/synthetic/iupac_valid.fastq",
    });
    try expectResult(
        narrow,
        1,
        "",
        "error: tests/data/synthetic/iupac_valid.fastq: S002: " ++
            "sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 11)\n",
    );

    const structural = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\n.\nx\n \n",
        1,
    );
    try expectResult(
        structural,
        1,
        "",
        "error: -: S001: plus line must start with '+' " ++
            "(record 0, line 3, offset 5)\n",
    );

    const semantic = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\n.\n+\n\x7f\n",
        2,
    );
    try expectResult(
        semantic,
        1,
        "",
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 3)\n",
    );
}

test "[cli] - [check]: semantic locations survive every field position and line ending" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const original = "@r\nACGT\n+\n!!!!";

    for (0..4) |invalid_index| {
        const sequence_data = try allocator.dupe(u8, original);
        sequence_data[3 + invalid_index] = '.';
        const sequence_result = try cli.runWithStdin(
            allocator,
            &.{ "check", "-" },
            sequence_data,
            invalid_index + 1,
        );
        const sequence_error = try expectedError(
            allocator,
            "-",
            "S002",
            "sequence byte is outside the selected alphabet",
            0,
            2,
            3 + invalid_index,
        );
        try expectResult(sequence_result, 1, "", sequence_error);

        const quality_data = try allocator.dupe(u8, original);
        quality_data[10 + invalid_index] = if (invalid_index % 2 == 0) 32 else 127;
        const quality_result = try cli.runWithStdin(
            allocator,
            &.{ "check", "-" },
            quality_data,
            invalid_index + 1,
        );
        const quality_error = try expectedError(
            allocator,
            "-",
            "S006",
            "quality byte must be ASCII 33 through 126",
            0,
            4,
            10 + invalid_index,
        );
        try expectResult(quality_result, 1, "", quality_error);
    }

    const crlf = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\r\nAC.\r\n+\r\n!!!\r\n",
        3,
    );
    try expectResult(
        crlf,
        1,
        "",
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 6)\n",
    );

    const quality_boundaries = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\nAC\n+\n!~\n",
        1,
    );
    try expectResult(quality_boundaries, 0, "", "");

    const later_record = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@ok\nA\n+\n!\n@bad\n.\n+\n!\n",
        1,
    );
    try expectResult(
        later_record,
        1,
        "",
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 1, line 2, offset 15)\n",
    );
}

test "[cli] - [check]: gzip members validate and corrupt trailers exit as I/O" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var gzip: std.ArrayList(u8) = .empty;

    try cli.appendGzipMember(allocator, &gzip, "@one\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &gzip, "@two\nGT\n+\n!~\n", .{});
    const valid = try cli.runWithStdin(allocator, &.{ "check", "-" }, gzip.items, 1);
    try expectResult(valid, 0, "", "");

    var semantic_gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &semantic_gzip, "@one\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &semantic_gzip, "@two\nG.\n+\n!!\n", .{});
    const semantic = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        semantic_gzip.items,
        1,
    );
    try expectResult(
        semantic,
        1,
        "",
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 1, line 2, offset 19)\n",
    );

    const corrupt = try allocator.dupe(u8, gzip.items);
    corrupt[corrupt.len - 8] ^= 1;
    const damaged = try cli.runWithStdin(allocator, &.{ "check", "-" }, corrupt, 7);
    try expectResult(damaged, 3, "", "error: -: I/O error\n");
}

test "[cli] - [check]: argument failures occur before input is consumed" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct {
        args: []const []const u8,
        stderr: []const u8,
    }{
        .{
            .args = &.{"check"},
            .stderr = "error: check requires at least one input\n",
        },
        .{
            .args = &.{ "check", "--alphabet" },
            .stderr = "error: --alphabet requires a value\n",
        },
        .{
            .args = &.{ "check", "--alphabet", "dna", "-" },
            .stderr = "error: --alphabet must be iupac or acgtn\n",
        },
        .{
            .args = &.{ "check", "--json", "-" },
            .stderr = "error: unknown check option: --json\n",
        },
        .{
            .args = &.{ "check", "--bogus", "-" },
            .stderr = "error: unknown check option: --bogus\n",
        },
        .{
            .args = &.{ "check", "missing.fastq", "-", "-" },
            .stderr = "error: standard input may appear at most once\n",
        },
    };
    for (cases) |case| {
        const result = try cli.runWithStdin(allocator, case.args, "@r\nA\n+\n!\n", 1);
        try expectResult(result, 2, "", case.stderr);
    }

    const double_dash = try cli.run(allocator, &.{ "check", "--", "--alphabet" });
    try expectResult(
        double_dash,
        3,
        "",
        "error: --alphabet: file not found\n",
    );
}

test "[cli] - [check]: independent inputs continue and highest exit class wins" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "valid.fastq", .data = "@r\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "semantic.fastq", .data = "@r\n.\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "long.fastq", .data = "@r\nAAAAA\n+\n!!!!!\n" });

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const valid = try tempPath(allocator, &tmp.sub_path, "valid.fastq");
    const semantic = try tempPath(allocator, &tmp.sub_path, "semantic.fastq");
    const long = try tempPath(allocator, &tmp.sub_path, "long.fastq");
    const missing = try tempPath(allocator, &tmp.sub_path, "missing.fastq");

    const result = try cli.run(allocator, &.{
        "check",
        "--max-line-bytes",
        "4",
        valid,
        semantic,
        missing,
        long,
        valid,
        semantic,
    });
    const expected = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 3)\n" ++
            "error: {s}: file not found\n" ++
            "error: {s}: line length limit exceeded\n" ++
            "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 3)\n",
        .{ semantic, missing, long, semantic },
    );
    try expectResult(result, 4, "", expected);
}

fn expectedError(
    allocator: std.mem.Allocator,
    label: []const u8,
    code: []const u8,
    message: []const u8,
    record_index: u64,
    line: u3,
    offset: u64,
) ![]const u8 {
    return std.fmt.allocPrint(
        allocator,
        "error: {s}: {s}: {s} (record {d}, line {d}, offset {d})\n",
        .{ label, code, message, record_index, line, offset },
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
