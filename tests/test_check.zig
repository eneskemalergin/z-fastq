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
        .message = "header line must start with '@' and contain a nonempty identifier",
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
        const json_result = try cli.run(allocator, &.{ "check", "--json", path });

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
            try expectCheckJsonResult(allocator, json_result, path, fixture);
        } else {
            try expectResult(file_result, 0, "", "");
            try expectResult(stdin_result, 0, "", "");
            try expectCheckJsonResult(allocator, json_result, path, fixture);
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

    const empty_identifier = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@ description\nA\n+\n!\n",
        1,
    );
    try expectResult(
        empty_identifier,
        1,
        "",
        "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n",
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

    const truncated = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\n.\n+\n",
        1,
    );
    try expectResult(
        truncated,
        1,
        "",
        "error: -: S004: unexpected end of file in quality line " ++
            "(record 0, line 4, offset 7)\n",
    );

    const mismatch = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\n.\n+\n!!\n",
        1,
    );
    try expectResult(
        mismatch,
        1,
        "",
        "error: -: S005: sequence and quality lengths differ " ++
            "(record 0, line 4, offset 7)\n",
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

    const valid_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--json", "-" },
        gzip.items,
        1,
    );
    try std.testing.expectEqual(@as(u8, 0), valid_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), valid_json.stderr.len);
    var parsed_valid = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        valid_json.stdout,
        .{},
    );
    defer parsed_valid.deinit();
    const valid_results = try cli.expectJsonDocument(
        &parsed_valid.value,
        "z-fastq/check-v1",
    );
    try std.testing.expectEqual(@as(usize, 1), valid_results.len);
    try expectCheckStatus(valid_results[0], "-", "ok");

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

    const damaged_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--json", "-" },
        corrupt,
        7,
    );
    try std.testing.expectEqual(@as(u8, 3), damaged_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), damaged_json.stderr.len);
    var parsed_damaged = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        damaged_json.stdout,
        .{},
    );
    defer parsed_damaged.deinit();
    const damaged_results = try cli.expectJsonDocument(
        &parsed_damaged.value,
        "z-fastq/check-v1",
    );
    try std.testing.expectEqual(@as(usize, 1), damaged_results.len);
    try expectCheckFailureWithoutLocation(
        damaged_results[0],
        "-",
        "io_error",
        "I/O error",
    );
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
            .args = &.{ "check", "--bogus", "-" },
            .stderr = "error: unknown check option: --bogus\n",
        },
        .{
            .args = &.{ "check", "missing.fastq", "-", "-" },
            .stderr = "error: standard input may appear at most once\n",
        },
        .{
            .args = &.{ "check", "--json", "missing.fastq", "-", "-" },
            .stderr = "error: standard input may appear at most once\n",
        },
        .{
            .args = &.{ "check", "--paired", "missing.fastq" },
            .stderr = "error: check --paired requires exactly two inputs\n",
        },
        .{
            .args = &.{ "check", "--paired", "missing.fastq", "other.fastq", "extra.fastq" },
            .stderr = "error: check --paired requires exactly two inputs\n",
        },
        .{
            .args = &.{ "check", "--interleaved" },
            .stderr = "error: check --interleaved requires exactly one input\n",
        },
        .{
            .args = &.{ "check", "--interleaved", "missing.fastq", "extra.fastq" },
            .stderr = "error: check --interleaved requires exactly one input\n",
        },
        .{
            .args = &.{ "check", "--paired", "--interleaved", "-" },
            .stderr = "error: --paired and --interleaved are mutually exclusive\n",
        },
        .{
            .args = &.{ "check", "--pair-names" },
            .stderr = "error: --pair-names requires a value\n",
        },
        .{
            .args = &.{ "check", "--pair-names", "prefix", "-" },
            .stderr = "error: --pair-names must be illumina or exact\n",
        },
        .{
            .args = &.{ "check", "--pair-names", "exact", "-" },
            .stderr = "error: --pair-names requires --paired or --interleaved\n",
        },
        .{
            .args = &.{ "check", "--paired", "-", "-" },
            .stderr = "error: paired inputs may contain standard input at most once\n",
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

    const json_without_input = try cli.run(allocator, &.{ "check", "--json" });
    try expectResult(
        json_without_input,
        2,
        "",
        "error: check requires at least one input\n",
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

    const json_result = try cli.run(allocator, &.{
        "check",
        "--json",
        "--max-line-bytes",
        "4",
        valid,
        semantic,
        missing,
        long,
        valid,
    });
    try std.testing.expectEqual(@as(u8, 4), json_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), json_result.stderr.len);
    var parsed = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        json_result.stdout,
        .{},
    );
    defer parsed.deinit();
    const results = try cli.expectJsonDocument(&parsed.value, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 5), results.len);
    try expectCheckStatus(results[0], valid, "ok");
    try expectCheckFailure(results[1], semantic, "S002", 0, 3, 2);
    try expectCheckFailureWithoutLocation(results[2], missing, "io_error", "file not found");
    try expectCheckFailureWithoutLocation(
        results[3],
        long,
        "line_limit",
        "line length limit exceeded",
    );
    try expectCheckStatus(results[4], valid, "ok");
}

test "[cli] - [check-json]: options, stdin, escaped bytes, and output failure compose" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const stdin_result = try cli.runWithStdin(
        allocator,
        &.{ "check", "--alphabet", "acgtn", "--json", "--max-line-bytes", "4", "-" },
        "@r\nA\n+\n!\n",
        1,
    );
    try std.testing.expectEqual(@as(u8, 0), stdin_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), stdin_result.stderr.len);
    var stdin_json = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        stdin_result.stdout,
        .{},
    );
    defer stdin_json.deinit();
    const stdin_results = try cli.expectJsonDocument(&stdin_json.value, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), stdin_results.len);
    try expectCheckStatus(stdin_results[0], "-", "ok");

    const unsafe_path = "missing\n\\\"\xff.fastq";
    const escaped_result = try cli.run(allocator, &.{ "check", "--json", unsafe_path });
    try std.testing.expectEqual(@as(u8, 3), escaped_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), escaped_result.stderr.len);
    var escaped_json = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        escaped_result.stdout,
        .{},
    );
    defer escaped_json.deinit();
    const escaped_results = try cli.expectJsonDocument(&escaped_json.value, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), escaped_results.len);
    try expectCheckFailureWithoutLocation(
        escaped_results[0],
        "missing\\x0A\\\\\"\\xFF.fastq",
        "io_error",
        "file not found",
    );

    const double_dash = try cli.run(allocator, &.{ "check", "--json", "--", "--json" });
    try std.testing.expectEqual(@as(u8, 3), double_dash.exit_code);
    try std.testing.expectEqual(@as(usize, 0), double_dash.stderr.len);
    var double_dash_json = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        double_dash.stdout,
        .{},
    );
    defer double_dash_json.deinit();
    const double_dash_results = try cli.expectJsonDocument(
        &double_dash_json.value,
        "z-fastq/check-v1",
    );
    try std.testing.expectEqual(@as(usize, 1), double_dash_results.len);
    try expectCheckFailureWithoutLocation(
        double_dash_results[0],
        "--json",
        "io_error",
        "file not found",
    );

    const closed = try cli.runWithClosedStdout(
        allocator,
        &.{ "check", "--json", "-" },
        "@r\nA\n+\n!\n",
    );
    try std.testing.expectEqual(@as(u8, 3), closed.exit_code);
    try std.testing.expectEqual(@as(usize, 0), closed.stderr.len);
}

test "[cli] - [paired check]: documented name forms and input transports pass" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");

    const cases = [_]struct {
        header1: []const u8,
        header2: []const u8,
        exact: bool = false,
    }{
        .{ .header1 = "cluster 1:N:0:ACGT", .header2 = "cluster 2:N:0:ACGT" },
        .{ .header1 = "cluster/1 opaque", .header2 = "cluster/2 opaque" },
        .{ .header1 = "cluster 1/1 opaque", .header2 = "cluster 1/2 opaque" },
        .{ .header1 = "cluster instrument/1 opaque", .header2 = "cluster instrument/2 opaque" },
        .{
            .header1 = "instrument:run:flowcell:1:2:3:4:UMI 1:N:0:index-a",
            .header2 = "instrument:run:flowcell:1:2:3:4:UMI 2:Y:0:index-b",
        },
        .{ .header1 = "cluster/1 opaque", .header2 = "cluster 2:N:0:index" },
        .{ .header1 = "cluster\t1:N:0:index", .header2 = "cluster\t2:N:0:index" },
        .{ .header1 = "cluster opaque", .header2 = "cluster other" },
        .{
            .header1 = "cluster 1:N:0:index opaque/2",
            .header2 = "cluster 2:N:0:index opaque/1",
        },
        .{
            .header1 = "cluster/1 1:N:0:index/1 opaque",
            .header2 = "cluster/2 2:N:0:index/2 opaque",
        },
        .{
            .header1 = "cluster 1:N:0:index",
            .header2 = "cluster 2:N:0:index",
            .exact = true,
        },
    };
    for (cases) |case| {
        try writeFastq(&tmp, io, allocator, "r1.fastq", &.{case.header1});
        try writeFastq(&tmp, io, allocator, "r2.fastq", &.{case.header2});
        const result = if (case.exact)
            try cli.run(allocator, &.{
                "check",
                "--paired",
                "--pair-names",
                "exact",
                r1_path,
                r2_path,
            })
        else
            try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
        try expectResult(result, 0, "", "");
    }

    const double_dash = try cli.run(
        allocator,
        &.{ "check", "--paired", "--", r1_path, r2_path },
    );
    try expectResult(double_dash, 0, "", "");

    try writeFastqPayload(&tmp, io, "r1.fastq", "");
    try writeFastqPayload(&tmp, io, "r2.fastq", "");
    const empty_pair = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    try expectResult(empty_pair, 0, "", "");

    const empty_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        "",
        1,
    );
    try expectResult(empty_interleaved, 0, "", "");

    const interleaved =
        "@cluster 1:N:0:index\r\nAC\r\n+\r\n!!\r\n" ++
        "@cluster 2:N:0:index\r\nGT\r\n+\r\n!!\r\n";
    const stdin_result = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        interleaved,
        1,
    );
    try expectResult(stdin_result, 0, "", "");

    const interleaved_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        interleaved,
        3,
    );
    try std.testing.expectEqual(@as(u8, 0), interleaved_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), interleaved_json.stderr.len);
    var parsed_interleaved = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        interleaved_json.stdout,
        .{},
    );
    defer parsed_interleaved.deinit();
    const interleaved_results = try cli.expectJsonDocument(
        &parsed_interleaved.value,
        "z-fastq/check-v1",
    );
    try expectCheckStatus(interleaved_results[0], "-", "ok");

    const exact_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--pair-names", "exact", "-" },
        interleaved,
        5,
    );
    try expectResult(exact_interleaved, 0, "", "");

    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"stdin/2"});
    const stdin_r1 = try cli.runWithStdin(
        allocator,
        &.{ "check", "--paired", "-", r2_path },
        "@stdin/1\nA\n+\n!\n",
        1,
    );
    try expectResult(stdin_r1, 0, "", "");
    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{"stdin/1"});
    const stdin_r2 = try cli.runWithStdin(
        allocator,
        &.{ "check", "--paired", r1_path, "-" },
        "@stdin/2\nA\n+\n!\n",
        1,
    );
    try expectResult(stdin_r2, 0, "", "");

    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, "@gzip/1\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &gzip2, "@gzip/2\nGT\n+\n!!\n", .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    const gzip_result = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    try expectResult(gzip_result, 0, "", "");

    var interleaved_gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &interleaved_gzip, "@gzip/1\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &interleaved_gzip, "@gzip/2\nGT\n+\n!!\n", .{});
    const gzip_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        interleaved_gzip.items,
        3,
    );
    try expectResult(gzip_interleaved, 0, "", "");

    const corrupt_interleaved = try allocator.dupe(u8, interleaved_gzip.items);
    corrupt_interleaved[corrupt_interleaved.len - 8] ^= 1;
    const damaged_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        corrupt_interleaved,
        7,
    );
    try std.testing.expectEqual(@as(u8, 3), damaged_interleaved.exit_code);
    try std.testing.expectEqual(@as(usize, 0), damaged_interleaved.stderr.len);
    var parsed_damaged = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        damaged_interleaved.stdout,
        .{},
    );
    defer parsed_damaged.deinit();
    const damaged_results = try cli.expectJsonDocument(
        &parsed_damaged.value,
        "z-fastq/check-v1",
    );
    try expectCheckFailureWithoutLocation(
        damaged_results[0],
        "-",
        "io_error",
        "I/O error",
    );

    const json_result = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 0), json_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), json_result.stderr.len);
    var parsed = try std.json.parseFromSlice(std.json.Value, allocator, json_result.stdout, .{});
    defer parsed.deinit();
    const results = try cli.expectJsonDocument(&parsed.value, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), results.len);
    const result_object = results[0].object;
    try cli.expectJsonObjectKeys(result_object, &.{ "inputs", "status" });
    try expectJsonInputs(result_object.get("inputs"), r1_path, r2_path);
    try cli.expectJsonString(result_object.get("status"), "ok");

    const corrupt_r2 = try allocator.dupe(u8, gzip2.items);
    corrupt_r2[corrupt_r2.len - 8] ^= 1;
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = corrupt_r2 });
    const damaged_pair = try cli.run(
        allocator,
        &.{ "check", "--paired", r1_path, r2_path },
    );
    const damaged_pair_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: I/O error\n",
        .{r2_path},
    );
    try expectResult(damaged_pair, 3, "", damaged_pair_error);

    const damaged_pair_json = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 3), damaged_pair_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), damaged_pair_json.stderr.len);
    var parsed_damaged_pair = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        damaged_pair_json.stdout,
        .{},
    );
    defer parsed_damaged_pair.deinit();
    const damaged_pair_results = try cli.expectJsonDocument(
        &parsed_damaged_pair.value,
        "z-fastq/check-v1",
    );
    const damaged_pair_object = damaged_pair_results[0].object;
    try cli.expectJsonObjectKeys(
        damaged_pair_object,
        &.{ "inputs", "status", "failed_input", "error" },
    );
    try expectJsonInputs(damaged_pair_object.get("inputs"), r1_path, r2_path);
    try cli.expectJsonString(damaged_pair_object.get("status"), "error");
    try cli.expectJsonString(damaged_pair_object.get("failed_input"), r2_path);
    const damaged_pair_details = damaged_pair_object.get("error").?.object;
    try cli.expectJsonString(damaged_pair_details.get("code"), "io_error");
    try cli.expectJsonString(damaged_pair_details.get("message"), "I/O error");
    try cli.expectJsonNullErrorLocation(damaged_pair_details);
}

test "[cli] - [paired check]: P001 diagnostics retain exact bounded identity fields" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{"cluster/1"});
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"other/2"});

    const human = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    const expected = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P001: paired identifiers or mate markers do not match (pair 0)\n" ++
            "  R1: input={s}, record=0, offset=0, first_token=cluster/1 " ++
            "[length=9, truncated=false], normalized_id=cluster " ++
            "[length=7, truncated=false], mate_markers=1\n" ++
            "  R2: input={s}, record=0, offset=0, first_token=other/2 " ++
            "[length=7, truncated=false], normalized_id=other " ++
            "[length=5, truncated=false], mate_markers=2\n",
        .{ r1_path, r2_path, r1_path, r2_path },
    );
    try expectResult(human, 1, "", expected);

    const json_result = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 1), json_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), json_result.stderr.len);
    var parsed = try std.json.parseFromSlice(std.json.Value, allocator, json_result.stdout, .{});
    defer parsed.deinit();
    const results = try cli.expectJsonDocument(&parsed.value, "z-fastq/check-v1");
    const result_object = results[0].object;
    try cli.expectJsonObjectKeys(result_object, &.{ "inputs", "status", "error" });
    try expectJsonInputs(result_object.get("inputs"), r1_path, r2_path);
    try cli.expectJsonString(result_object.get("status"), "error");
    try expectPairNameJson(result_object.get("error"), "cluster/1", "cluster", "other/2", "other");

    const mismatch_cases = [_]struct {
        header1: []const u8,
        header2: []const u8,
        exact: bool = false,
    }{
        .{ .header1 = "cluster/2", .header2 = "cluster/1" },
        .{ .header1 = "cluster/1", .header2 = "cluster/1" },
        .{ .header1 = "cluster/1", .header2 = "cluster" },
        .{ .header1 = "cluster/1 2:N:0:index", .header2 = "cluster/2" },
        .{ .header1 = "cluster 1/2", .header2 = "cluster 1/1" },
        .{ .header1 = "/1", .header2 = "/2" },
        .{
            .header1 = "instrument:run:flowcell:1:2:3:4:UMI-A 1:N:0:index",
            .header2 = "instrument:run:flowcell:1:2:3:4:UMI-B 2:N:0:index",
        },
        .{ .header1 = "cluster/1", .header2 = "cluster/2", .exact = true },
    };
    for (mismatch_cases) |case| {
        try writeFastq(&tmp, io, allocator, "r1.fastq", &.{case.header1});
        try writeFastq(&tmp, io, allocator, "r2.fastq", &.{case.header2});
        const result = if (case.exact)
            try cli.run(allocator, &.{
                "check",
                "--paired",
                "--pair-names",
                "exact",
                "--json",
                r1_path,
                r2_path,
            })
        else
            try cli.run(allocator, &.{
                "check",
                "--paired",
                "--json",
                r1_path,
                r2_path,
            });
        try std.testing.expectEqual(@as(u8, 1), result.exit_code);
        try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
        var mismatch_json = try std.json.parseFromSlice(
            std.json.Value,
            allocator,
            result.stdout,
            .{},
        );
        defer mismatch_json.deinit();
        const mismatch_results = try cli.expectJsonDocument(
            &mismatch_json.value,
            "z-fastq/check-v1",
        );
        const mismatch_object = mismatch_results[0].object;
        try cli.expectJsonString(mismatch_object.get("status"), "error");
        try cli.expectJsonString(mismatch_object.get("error").?.object.get("code"), "P001");
    }

    try writeFastqPayload(&tmp, io, "r1.fastq", "@clu\x00ster/1\nA\n+\n!\n");
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"other/2"});
    const escaped = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    var parsed_escaped = try std.json.parseFromSlice(std.json.Value, allocator, escaped.stdout, .{});
    defer parsed_escaped.deinit();
    const escaped_results = try cli.expectJsonDocument(
        &parsed_escaped.value,
        "z-fastq/check-v1",
    );
    const escaped_tokens = escaped_results[0].object.get("error").?.object
        .get("first_tokens").?.array.items;
    try cli.expectJsonString(escaped_tokens[0].object.get("prefix"), "clu\\x00ster/1");

    const interleaved_payload =
        "@cluster/1\nA\n+\n!\n" ++
        "@other/2\nA\n+\n!\n";
    const interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        interleaved_payload,
        1,
    );
    try expectResult(
        interleaved,
        1,
        "",
        "error: -: P001: paired identifiers or mate markers do not match (pair 0)\n" ++
            "  R1: input=-, record=0, offset=0, first_token=cluster/1 " ++
            "[length=9, truncated=false], normalized_id=cluster " ++
            "[length=7, truncated=false], mate_markers=1\n" ++
            "  R2: input=-, record=1, offset=17, first_token=other/2 " ++
            "[length=7, truncated=false], normalized_id=other " ++
            "[length=5, truncated=false], mate_markers=2\n",
    );
    const interleaved_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        interleaved_payload,
        3,
    );
    var parsed_interleaved = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        interleaved_json.stdout,
        .{},
    );
    defer parsed_interleaved.deinit();
    const interleaved_results = try cli.expectJsonDocument(
        &parsed_interleaved.value,
        "z-fastq/check-v1",
    );
    const interleaved_error = interleaved_results[0].object.get("error").?.object;
    try expectJsonIntegerPair(interleaved_error.get("record_indexes"), 0, 1);
    try expectJsonIntegerPair(interleaved_error.get("byte_offsets"), 0, 17);
    const interleaved_tokens = interleaved_error.get("first_tokens").?.array.items;
    try expectBoundedJson(interleaved_tokens[0], "cluster/1", 9, false);

    const later_interleaved_payload =
        "@long-cluster-name/1\nA\n+\n!\n" ++
        "@long-cluster-name/2\nA\n+\n!\n" ++
        "@x/1\nA\n+\n!\n" ++
        "@y/2\nA\n+\n!\n";
    const later_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        later_interleaved_payload,
        7,
    );
    try std.testing.expectEqual(@as(u8, 1), later_interleaved.exit_code);
    try std.testing.expectEqual(@as(usize, 0), later_interleaved.stderr.len);
    var parsed_later = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        later_interleaved.stdout,
        .{},
    );
    defer parsed_later.deinit();
    const later_results = try cli.expectJsonDocument(
        &parsed_later.value,
        "z-fastq/check-v1",
    );
    const later_error = later_results[0].object.get("error").?.object;
    try cli.expectJsonInteger(later_error.get("pair_index"), 1);
    try expectJsonIntegerPair(later_error.get("record_indexes"), 2, 3);

    const long1 = try allocator.alloc(u8, 160);
    const long2 = try allocator.alloc(u8, 160);
    @memset(long1, 'A');
    @memset(long2, 'A');
    long2[long2.len - 1] = 'B';
    const header1 = try std.fmt.allocPrint(allocator, "{s}/1", .{long1});
    const header2 = try std.fmt.allocPrint(allocator, "{s}/2", .{long2});
    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{header1});
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{header2});
    const bounded = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 1), bounded.exit_code);
    try std.testing.expect(bounded.stdout.len < 2048);
    var parsed_bounded = try std.json.parseFromSlice(std.json.Value, allocator, bounded.stdout, .{});
    defer parsed_bounded.deinit();
    const bounded_results = try cli.expectJsonDocument(&parsed_bounded.value, "z-fastq/check-v1");
    const bounded_error = bounded_results[0].object.get("error").?.object;
    const normalized = bounded_error.get("normalized_ids").?.array.items;
    try expectBoundedJson(normalized[0], long1[0..128], 160, true);
    try expectBoundedJson(normalized[1], long2[0..128], 160, true);

    const long_interleaved_payload = try std.fmt.allocPrint(
        allocator,
        "@{s}/1\nA\n+\n!\n@{s}/2\nA\n+\n!\n",
        .{ long1, long2 },
    );
    const long_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        long_interleaved_payload,
        17,
    );
    try std.testing.expectEqual(@as(u8, 1), long_interleaved.exit_code);
    try std.testing.expectEqual(@as(usize, 0), long_interleaved.stderr.len);
    try std.testing.expect(long_interleaved.stdout.len < 2048);
    var parsed_long_interleaved = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        long_interleaved.stdout,
        .{},
    );
    defer parsed_long_interleaved.deinit();
    const long_interleaved_results = try cli.expectJsonDocument(
        &parsed_long_interleaved.value,
        "z-fastq/check-v1",
    );
    const long_interleaved_error = long_interleaved_results[0].object
        .get("error").?.object;
    const long_first_tokens = long_interleaved_error.get("first_tokens").?.array.items;
    const long_normalized = long_interleaved_error.get("normalized_ids").?.array.items;
    try expectBoundedJson(long_first_tokens[0], long1[0..128], 162, true);
    try expectBoundedJson(long_first_tokens[1], long2[0..128], 162, true);
    try expectBoundedJson(long_normalized[0], long1[0..128], 160, true);
    try expectBoundedJson(long_normalized[1], long2[0..128], 160, true);
}

test "[cli] - [paired check]: P002 and semantic precedence are exact" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{ "first/1", "second/1" });
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"first/2"});

    const count = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    const count_error = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 1, last R2 record 0)\n",
        .{ r1_path, r2_path },
    );
    try expectResult(count, 1, "", count_error);

    const count_json = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 1), count_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), count_json.stderr.len);
    var parsed_count = try std.json.parseFromSlice(std.json.Value, allocator, count_json.stdout, .{});
    defer parsed_count.deinit();
    const count_results = try cli.expectJsonDocument(&parsed_count.value, "z-fastq/check-v1");
    const count_object = count_results[0].object;
    try cli.expectJsonObjectKeys(count_object, &.{ "inputs", "status", "error" });
    try expectJsonInputs(count_object.get("inputs"), r1_path, r2_path);
    try cli.expectJsonString(count_object.get("status"), "error");
    const pair_count_error = count_object.get("error").?.object;
    try cli.expectJsonObjectKeys(pair_count_error, &.{
        "code",
        "message",
        "pair_index",
        "remaining_side",
        "record_indexes",
    });
    try cli.expectJsonString(pair_count_error.get("code"), "P002");
    try cli.expectJsonString(pair_count_error.get("message"), "paired input is missing a mate");
    try cli.expectJsonInteger(pair_count_error.get("pair_index"), 1);
    try cli.expectJsonString(pair_count_error.get("remaining_side"), "R1");
    try expectJsonIntegerPair(pair_count_error.get("record_indexes"), 1, 0);

    try writeFastqPayload(&tmp, io, "r1.fastq", "");
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"only/2"});
    const r2_remaining = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    const r2_remaining_error = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P002: paired input is missing a mate " ++
            "(pair 0, remaining R2, last R1 record none, last R2 record 0)\n",
        .{ r1_path, r2_path },
    );
    try expectResult(r2_remaining, 1, "", r2_remaining_error);

    const odd = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        "@odd/1\n.\n+\n!\n",
        1,
    );
    try expectResult(
        odd,
        1,
        "",
        "error: -: P002: paired input is missing a mate " ++
            "(pair 0, remaining R1, last R1 record 0, last R2 record none)\n",
    );

    const odd_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        "@odd/1\nA\n+\n!\n",
        2,
    );
    try std.testing.expectEqual(@as(u8, 1), odd_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), odd_json.stderr.len);
    var parsed_odd = try std.json.parseFromSlice(std.json.Value, allocator, odd_json.stdout, .{});
    defer parsed_odd.deinit();
    const odd_results = try cli.expectJsonDocument(&parsed_odd.value, "z-fastq/check-v1");
    const odd_object = odd_results[0].object;
    try cli.expectJsonObjectKeys(odd_object, &.{ "input", "status", "error" });
    try cli.expectJsonString(odd_object.get("input"), "-");
    try cli.expectJsonString(odd_object.get("status"), "error");
    const odd_error = odd_object.get("error").?.object;
    try cli.expectJsonString(odd_error.get("code"), "P002");
    const odd_indexes = odd_error.get("record_indexes").?.array.items;
    try cli.expectJsonInteger(odd_indexes[0], 0);
    try std.testing.expect(odd_indexes[1] == .null);

    const later_odd = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        "@a/1\nA\n+\n!\n@a/2\nA\n+\n!\n@b/1\nA\n+\n!\n",
        3,
    );
    try expectResult(
        later_odd,
        1,
        "",
        "error: -: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 2, last R2 record 1)\n",
    );

    try writeFastqPayload(&tmp, io, "r1.fastq", "@a/1\n.\n+\n!\n");
    try writeFastqPayload(&tmp, io, "r2.fastq", "@b/2\nA\n+\n!\n");
    const r1_semantic = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    const r1_error = try expectedError(
        allocator,
        r1_path,
        "S002",
        "sequence byte is outside the selected alphabet",
        0,
        2,
        5,
    );
    try expectResult(r1_semantic, 1, "", r1_error);

    try writeFastqPayload(&tmp, io, "r2.fastq", "@b/2\n.\n+\n!\n");
    const both_semantic = try cli.run(
        allocator,
        &.{ "check", "--paired", r1_path, r2_path },
    );
    try expectResult(both_semantic, 1, "", r1_error);

    try writeFastqPayload(&tmp, io, "r1.fastq", "@a/1\nR\n+\n!\n");
    try writeFastqPayload(&tmp, io, "r2.fastq", "@a/2\nA\n+\n!\n");
    const alphabet = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--alphabet",
        "acgtn",
        r1_path,
        r2_path,
    });
    const alphabet_error = try expectedError(
        allocator,
        r1_path,
        "S002",
        "sequence byte is outside the selected alphabet",
        0,
        2,
        5,
    );
    try expectResult(alphabet, 1, "", alphabet_error);

    const line_limit = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--max-line-bytes",
        "3",
        r1_path,
        r2_path,
    });
    const line_error = try std.fmt.allocPrint(
        allocator,
        "error: {s}: line length limit exceeded\n",
        .{r1_path},
    );
    try expectResult(line_limit, 4, "", line_error);

    try writeFastqPayload(&tmp, io, "r1.fastq", "@a/1\nA\n+\n!\n");
    try writeFastqPayload(&tmp, io, "r2.fastq", "@b/2\n.\n+\n!\n");
    const r2_semantic = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    const r2_error = try expectedError(
        allocator,
        r2_path,
        "S002",
        "sequence byte is outside the selected alphabet",
        0,
        2,
        5,
    );
    try expectResult(r2_semantic, 1, "", r2_error);

    const r2_semantic_json = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 1), r2_semantic_json.exit_code);
    var parsed_semantic = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        r2_semantic_json.stdout,
        .{},
    );
    defer parsed_semantic.deinit();
    const semantic_results = try cli.expectJsonDocument(
        &parsed_semantic.value,
        "z-fastq/check-v1",
    );
    const semantic_object = semantic_results[0].object;
    try cli.expectJsonObjectKeys(
        semantic_object,
        &.{ "inputs", "status", "failed_input", "error" },
    );
    try cli.expectJsonString(semantic_object.get("failed_input"), r2_path);
    const semantic_error = semantic_object.get("error").?.object;
    try cli.expectJsonString(semantic_error.get("code"), "S002");
    try cli.expectJsonErrorLocation(
        semantic_error,
        0,
        5,
        2,
    );

    try writeFastqPayload(&tmp, io, "r1.fastq", "@a/1\n.\n+\n!\n");
    try writeFastqPayload(&tmp, io, "r2.fastq", "@a/2\nA\n+\n");
    const r2_truncated = try cli.run(
        allocator,
        &.{ "check", "--paired", r1_path, r2_path },
    );
    const r2_truncated_error = try expectedError(
        allocator,
        r2_path,
        "S004",
        "unexpected end of file in quality line",
        0,
        4,
        9,
    );
    try expectResult(r2_truncated, 1, "", r2_truncated_error);

    const structural = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        "@a/1\n.\n+\n!\n@a/2\nA\nx\n!\n",
        1,
    );
    try expectResult(
        structural,
        1,
        "",
        "error: -: S001: plus line must start with '+' " ++
            "(record 1, line 3, offset 18)\n",
    );

    const empty_identifier = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--pair-names", "exact", "-" },
        "@ description\nA\n+\n!\n@ unrelated\nT\n+\n#\n",
        1,
    );
    try expectResult(
        empty_identifier,
        1,
        "",
        "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n",
    );
}

fn writeFastq(
    tmp: *std.testing.TmpDir,
    io: std.Io,
    allocator: std.mem.Allocator,
    name: []const u8,
    headers: []const []const u8,
) !void {
    var data: std.ArrayList(u8) = .empty;
    defer data.deinit(allocator);
    for (headers) |header| {
        try data.append(allocator, '@');
        try data.appendSlice(allocator, header);
        try data.appendSlice(allocator, "\nA\n+\n!\n");
    }
    try tmp.dir.writeFile(io, .{ .sub_path = name, .data = data.items });
}

fn writeFastqPayload(
    tmp: *std.testing.TmpDir,
    io: std.Io,
    name: []const u8,
    payload: []const u8,
) !void {
    try tmp.dir.writeFile(io, .{ .sub_path = name, .data = payload });
}

fn expectJsonInputs(value: ?std.json.Value, input1: []const u8, input2: []const u8) !void {
    const inputs = switch (value orelse return error.UnexpectedJsonShape) {
        .array => |array| array.items,
        else => return error.UnexpectedJsonShape,
    };
    try std.testing.expectEqual(@as(usize, 2), inputs.len);
    try cli.expectJsonString(inputs[0], input1);
    try cli.expectJsonString(inputs[1], input2);
}

fn expectPairNameJson(
    value: ?std.json.Value,
    first1: []const u8,
    normalized1: []const u8,
    first2: []const u8,
    normalized2: []const u8,
) !void {
    const object = switch (value orelse return error.UnexpectedJsonShape) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try cli.expectJsonObjectKeys(object, &.{
        "code",
        "message",
        "pair_index",
        "record_indexes",
        "byte_offsets",
        "first_tokens",
        "normalized_ids",
        "mate_markers",
    });
    try cli.expectJsonString(object.get("code"), "P001");
    try cli.expectJsonString(
        object.get("message"),
        "paired identifiers or mate markers do not match",
    );
    try cli.expectJsonInteger(object.get("pair_index"), 0);
    try expectJsonIntegerPair(object.get("record_indexes"), 0, 0);
    try expectJsonIntegerPair(object.get("byte_offsets"), 0, 0);

    const first_tokens = object.get("first_tokens").?.array.items;
    const normalized_ids = object.get("normalized_ids").?.array.items;
    try expectBoundedJson(first_tokens[0], first1, first1.len, false);
    try expectBoundedJson(first_tokens[1], first2, first2.len, false);
    try expectBoundedJson(normalized_ids[0], normalized1, normalized1.len, false);
    try expectBoundedJson(normalized_ids[1], normalized2, normalized2.len, false);

    const markers = object.get("mate_markers").?.array.items;
    try std.testing.expectEqual(@as(usize, 2), markers.len);
    try std.testing.expectEqual(@as(usize, 1), markers[0].array.items.len);
    try std.testing.expectEqual(@as(usize, 1), markers[1].array.items.len);
    try cli.expectJsonInteger(markers[0].array.items[0], 1);
    try cli.expectJsonInteger(markers[1].array.items[0], 2);
}

fn expectJsonIntegerPair(value: ?std.json.Value, first: u64, second: u64) !void {
    const values = switch (value orelse return error.UnexpectedJsonShape) {
        .array => |array| array.items,
        else => return error.UnexpectedJsonShape,
    };
    try std.testing.expectEqual(@as(usize, 2), values.len);
    try cli.expectJsonInteger(values[0], first);
    try cli.expectJsonInteger(values[1], second);
}

fn expectBoundedJson(
    value: std.json.Value,
    prefix: []const u8,
    length: usize,
    truncated: bool,
) !void {
    const object = switch (value) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try cli.expectJsonObjectKeys(object, &.{ "prefix", "length", "truncated" });
    try cli.expectJsonString(object.get("prefix"), prefix);
    try cli.expectJsonInteger(object.get("length"), length);
    const actual_truncated = switch (object.get("truncated") orelse
        return error.UnexpectedJsonShape) {
        .bool => |actual| actual,
        else => return error.UnexpectedJsonShape,
    };
    try std.testing.expectEqual(truncated, actual_truncated);
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

fn expectCheckJsonResult(
    allocator: std.mem.Allocator,
    command_result: cli.CommandResult,
    input: []const u8,
    fixture: FixtureExpect,
) !void {
    try std.testing.expectEqual(
        @as(u8, if (fixture.code == null) 0 else 1),
        command_result.exit_code,
    );
    try std.testing.expectEqual(@as(usize, 0), command_result.stderr.len);
    var parsed = try std.json.parseFromSlice(
        std.json.Value,
        allocator,
        command_result.stdout,
        .{},
    );
    defer parsed.deinit();
    const results = try cli.expectJsonDocument(&parsed.value, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), results.len);
    if (fixture.code) |code| {
        try expectCheckFailure(
            results[0],
            input,
            code,
            fixture.record_index,
            fixture.offset,
            fixture.line,
        );
        const result = results[0].object;
        const error_object = result.get("error").?.object;
        try cli.expectJsonString(error_object.get("message"), fixture.message);
    } else {
        try expectCheckStatus(results[0], input, "ok");
    }
}

fn expectCheckStatus(value: std.json.Value, input: []const u8, status: []const u8) !void {
    const object = switch (value) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try cli.expectJsonObjectKeys(object, &.{ "input", "status" });
    try cli.expectJsonString(object.get("input"), input);
    try cli.expectJsonString(object.get("status"), status);
}

fn expectCheckFailure(
    value: std.json.Value,
    input: []const u8,
    code: []const u8,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u3,
) !void {
    const error_object = try cli.expectJsonError(value, input, code);
    try cli.expectJsonErrorLocation(
        error_object,
        record_index,
        byte_offset,
        line_in_record,
    );
}

fn expectCheckFailureWithoutLocation(
    value: std.json.Value,
    input: []const u8,
    code: []const u8,
    message: []const u8,
) !void {
    const error_object = try cli.expectJsonError(value, input, code);
    try cli.expectJsonString(error_object.get("message"), message);
    try cli.expectJsonNullErrorLocation(error_object);
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
