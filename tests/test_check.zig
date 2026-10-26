//! Installed CLI contracts for `z-fastq check`.

const std = @import("std");
const cli = @import("utilities.zig");

const FIXTURE_DIR = cli.FIXTURE_DIR;
/// The identifier prefix kept in P001 diagnostics (`PAIR_DIAGNOSTIC_PREFIX_BYTES` in src/main.zig).
const DIAGNOSTIC_PREFIX_BYTES = 128;
/// An identifier longer than the diagnostic prefix, so the prefix is the only part shown.
const LONG_IDENTIFIER_BYTES = DIAGNOSTIC_PREFIX_BYTES + 32;
/// Upper bound on the size of a bounded P001 report.
const BOUNDED_REPORT_BYTES = 2048;
/// Length of each directory component in the long-path test; a path longer than the JSON buffer.
const LONG_COMPONENT_BYTES = 200;

const FixtureExpect = struct {
    name: []const u8,
    code: ?[]const u8 = null,
    message: []const u8 = "",
    record_index: u64 = 0,
    line: u3 = 0,
    offset: u64 = 0,
};

const FIXTURES = [_]FixtureExpect{
    .{ .name = "crlf.fastq" },
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
        const path = try cli.fixturePath(allocator, fixture.name);
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
            try cli.expectResult(file_result, 1, "", file_error);
            try cli.expectResult(stdin_result, 1, "", stdin_error);
            try expectCheckJsonResult(allocator, json_result, path, fixture);
        } else {
            try cli.expectResult(file_result, 0, "", "");
            try cli.expectResult(stdin_result, 0, "", "");
            try expectCheckJsonResult(allocator, json_result, path, fixture);
        }
    }
}

test "[integration] - [check]: counts agree with count and stats under explicit validation rules" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "records");
    const cases = [_]struct {
        input: []const u8,
        reads: u64,
        alphabet: []const u8 = "acgtn",
        check_error: []const u8 = "",
        stats_ok: bool = true,
    }{
        .{ .input = "", .reads = 0 },
        .{ .input = "@r\r\nAcgTn\r\n+note\r\n!5?I~\r\n", .reads = 1 },
        .{
            .input = "@r\n?\n+\n!\n",
            .reads = 1,
            .check_error = "S002: sequence byte is outside the selected alphabet (record 0, line 2, offset 3)\n",
        },
        .{
            .input = "@r\nA\n+\n \n",
            .reads = 1,
            .check_error = "S006: quality byte must be ASCII 33 through 126 (record 0, line 4, offset 7)\n",
            .stats_ok = false,
        },
    };
    for (cases, 0..) |case, index| {
        var gzip: std.ArrayList(u8) = .empty;
        try cli.appendGzipMember(allocator, &gzip, case.input, .{});
        for ([_][]const u8{ case.input, gzip.items }, 0..) |bytes, encoding| {
            errdefer std.debug.print("count comparison case={d} encoding={d}\n", .{ index, encoding });
            try tmp.dir.writeFile(io, .{ .sub_path = "records", .data = bytes });
            const diagnostic = if (case.check_error.len == 0) "" else try std.fmt.allocPrint(allocator, "error: {s}: {s}", .{ path, case.check_error });
            try cli.expectResult(
                try cli.run(allocator, &.{ "check", "--alphabet", case.alphabet, "--max-line-bytes", "64", path }),
                if (case.check_error.len == 0) 0 else 1,
                "",
                diagnostic,
            );
            try cli.expectResult(
                try cli.run(allocator, &.{ "count", "--max-line-bytes", "64", path }),
                0,
                try std.fmt.allocPrint(allocator, "{d}\n", .{case.reads}),
                "",
            );
            if (!case.stats_ok) {
                try cli.expectResult(
                    try cli.run(allocator, &.{ "stats", "--max-line-bytes", "64", path }),
                    1,
                    "",
                    diagnostic,
                );
                continue;
            }
            const stats = try cli.run(allocator, &.{ "stats", "--json", "--max-line-bytes", "64", path });
            try std.testing.expectEqual(@as(u8, 0), stats.exit_code);
            try std.testing.expectEqualStrings("", stats.stderr);
            const results = try cli.parseJsonResults(allocator, stats.stdout, "z-fastq/stats-v1");
            try std.testing.expectEqual(@as(usize, 1), results.len);
            try cli.expectJsonString(results[0].object.get("input"), path);
            try cli.expectJsonString(results[0].object.get("status"), "ok");
            try cli.expectJsonInteger(results[0].object.get("reads"), case.reads);
        }
    }
}

test "[cli] - [check]: the alphabet option accepts only its selected sequence symbols" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const explicit_iupac = try cli.run(allocator, &.{
        "check",
        "--alphabet",
        "iupac",
        "tests/data/synthetic/iupac_valid.fastq",
    });
    try cli.expectResult(explicit_iupac, 0, "", "");

    const narrow_valid = try cli.run(allocator, &.{
        "check",
        "--alphabet",
        "acgtn",
        "tests/data/synthetic/acgtn_valid.fastq",
    });
    try cli.expectResult(narrow_valid, 0, "", "");

    const narrow = try cli.run(allocator, &.{
        "check",
        "--alphabet",
        "acgtn",
        "tests/data/synthetic/iupac_valid.fastq",
    });
    try cli.expectResult(
        narrow,
        1,
        "",
        "error: tests/data/synthetic/iupac_valid.fastq: S002: " ++
            "sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 11)\n",
    );
}

test "[cli] - [check]: a record reports the first rule it breaks" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    // Each input breaks several rules; the reported error is the one earliest in the record.
    const cases = [_]struct {
        args: []const []const u8 = &.{ "check", "-" },
        input: []const u8,
        chunk_len: usize = 1,
        stderr: []const u8,
    }{
        .{
            .input = "@r\n.\nx\n \n",
            .stderr = "error: -: S001: plus line must start with '+' (record 0, line 3, offset 5)\n",
        },
        .{
            .args = &.{ "check", "--interleaved", "-" },
            .input = "@a/1\nR\n+\n \n@a/2\nAAA\n+\n!\n!\n",
            .chunk_len = 26,
            .stderr = "error: -: S005: sequence and quality lengths differ (record 1, line 4, offset 22)\n",
        },
        .{
            .input = "@ description\nA\n+\n!\n",
            .stderr = "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
                "(record 0, line 1, offset 0)\n",
        },
        .{
            .input = "@r\n.\n+\n\x7f\n",
            .chunk_len = 2,
            .stderr = "error: -: S002: sequence byte is outside the selected alphabet (record 0, line 2, offset 3)\n",
        },
        .{
            .input = "@r\n.\n+\n",
            .stderr = "error: -: S004: unexpected end of file in quality line (record 0, line 4, offset 7)\n",
        },
        .{
            .input = "@r\n.\n+\n!!\n",
            .stderr = "error: -: S005: sequence and quality lengths differ (record 0, line 4, offset 7)\n",
        },
    };
    for (cases) |case| {
        const result = try cli.runWithStdin(allocator, case.args, case.input, case.chunk_len);
        try cli.expectResult(result, 1, "", case.stderr);
    }
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
        try cli.expectResult(sequence_result, 1, "", sequence_error);
    }

    const crlf = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@r\r\nAC.\r\n+\r\n!!!\r\n",
        3,
    );
    try cli.expectResult(
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
    try cli.expectResult(quality_boundaries, 0, "", "");

    const later_record = try cli.runWithStdin(
        allocator,
        &.{ "check", "-" },
        "@ok\nA\n+\n!\n@bad\n.\n+\n!\n",
        1,
    );
    try cli.expectResult(
        later_record,
        1,
        "",
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 1, line 2, offset 15)\n",
    );
}

test "[cli] - [check]: files and stdin preserve prefix detection and gzip failures" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input");
    const plain = "@one\nAC\n+\n!!\n@two\nGT\n+\n!~\n";
    var gzip: std.ArrayList(u8) = .empty;

    try cli.appendGzipMember(allocator, &gzip, "@one\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &gzip, "@two\nGT\n+\n!~\n", .{});
    var semantic_gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &semantic_gzip, "@one\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &semantic_gzip, "@two\nG.\n+\n!!\n", .{});
    const corrupt = try allocator.dupe(u8, gzip.items);
    corrupt[cli.gzipCrcOffset(corrupt.len)] ^= 1;

    const cases = [_]struct {
        bytes: []const u8,
        exit_code: u8 = 0,
        message: []const u8 = "",
    }{
        .{ .bytes = plain },
        .{ .bytes = gzip.items },
        .{ .bytes = gzip.items[0..9], .exit_code = 3, .message = "I/O error\n" },
        .{ .bytes = corrupt, .exit_code = 3, .message = "I/O error\n" },
        .{
            .bytes = semantic_gzip.items,
            .exit_code = 1,
            .message = "S002: sequence byte is outside the selected alphabet " ++
                "(record 1, line 2, offset 19)\n",
        },
    };
    for (cases, 0..) |case, index| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input", .data = case.bytes });
        for ([_][]const u8{ path, "-" }, 0..) |label, transport| {
            errdefer std.debug.print("check input case {d}, input {s}\n", .{ index, label });
            const result = if (transport == 0)
                try cli.run(allocator, &.{ "check", path })
            else
                try cli.runWithStdin(allocator, &.{ "check", "-" }, case.bytes, 1);
            const stderr = if (case.message.len == 0) "" else try std.fmt.allocPrint(
                allocator,
                "error: {s}: {s}",
                .{ label, case.message },
            );
            try cli.expectResult(result, case.exit_code, "", stderr);
        }
    }

    const valid_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--json", "-" },
        gzip.items,
        1,
    );
    try std.testing.expectEqual(@as(u8, 0), valid_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), valid_json.stderr.len);
    const valid_results = try cli.parseJsonResults(allocator, valid_json.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), valid_results.len);
    try expectCheckStatus(valid_results[0], "-", "ok");

    const damaged_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--json", "-" },
        corrupt,
        7,
    );
    try std.testing.expectEqual(@as(u8, 3), damaged_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), damaged_json.stderr.len);
    const damaged_results = try cli.parseJsonResults(allocator, damaged_json.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), damaged_results.len);
    try cli.expectJsonFailureWithoutLocation(
        damaged_results[0],
        "-",
        "io_error",
        "I/O error",
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
    const valid = try cli.tempPath(allocator, &tmp.sub_path, "valid.fastq");
    const semantic = try cli.tempPath(allocator, &tmp.sub_path, "semantic.fastq");
    const long = try cli.tempPath(allocator, &tmp.sub_path, "long.fastq");
    const missing = try cli.tempPath(allocator, &tmp.sub_path, "missing.fastq");

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
    try cli.expectResult(result, 4, "", expected);

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
    const results = try cli.parseJsonResults(allocator, json_result.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 5), results.len);
    try expectCheckStatus(results[0], valid, "ok");
    try cli.expectJsonFailure(results[1], semantic, "S002", 0, 3, 2);
    try cli.expectJsonFailureWithoutLocation(results[2], missing, "io_error", "file not found");
    try cli.expectJsonFailureWithoutLocation(
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
    const stdin_results = try cli.parseJsonResults(allocator, stdin_result.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), stdin_results.len);
    try expectCheckStatus(stdin_results[0], "-", "ok");

    const unsafe_path = "missing\n\\\"\xff.fastq";
    const escaped_result = try cli.run(allocator, &.{ "check", "--json", unsafe_path });
    try std.testing.expectEqual(@as(u8, 3), escaped_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), escaped_result.stderr.len);
    const escaped_results = try cli.parseJsonResults(allocator, escaped_result.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), escaped_results.len);
    try cli.expectJsonFailureWithoutLocation(
        escaped_results[0],
        "missing\\x0A\\\\\"\\xFF.fastq",
        "io_error",
        "file not found",
    );
}

test "[cli] - [check-json]: a double dash makes a leading-hyphen name a file path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const double_dash = try cli.run(allocator, &.{ "check", "--json", "--", "--json" });
    try std.testing.expectEqual(@as(u8, 3), double_dash.exit_code);
    try std.testing.expectEqual(@as(usize, 0), double_dash.stderr.len);
    const double_dash_results = try cli.parseJsonResults(allocator, double_dash.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), double_dash_results.len);
    try cli.expectJsonFailureWithoutLocation(
        double_dash_results[0],
        "--json",
        "io_error",
        "file not found",
    );
}

test "[cli] - [check-json]: output failure preserves an observed limit status in every mode" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const long_directory = ("\x01" ** LONG_COMPONENT_BYTES ++ "/") ** 17;
    try tmp.dir.createDirPath(io, long_directory);
    const empty = FIXTURE_DIR ++ "/empty_valid.fastq";

    for ([_][]const u8{ "limit.fastq", long_directory ++ "limit.fastq" }) |name| {
        try tmp.dir.writeFile(io, .{ .sub_path = name, .data = "@abcde\nA\n+\n!\n" });
        const path = try cli.tempPath(allocator, &tmp.sub_path, name);
        for ([_][]const []const u8{
            &.{ "check", "--json", "--max-line-bytes", "4", path, empty },
            &.{ "check", "--paired", "--json", "--max-line-bytes", "4", path, empty },
            &.{ "check", "--interleaved", "--json", "--max-line-bytes", "4", path },
        }) |args| {
            const captured = try cli.run(allocator, args);
            try std.testing.expectEqual(@as(u8, 4), captured.exit_code);
            try std.testing.expectEqualStrings("", captured.stderr);
            try std.testing.expectEqual(name.len > LONG_COMPONENT_BYTES, captured.stdout.len > cli.JSON_OUTPUT_BUFFER_BYTES);
            const results = try cli.parseJsonResults(allocator, captured.stdout, "z-fastq/check-v1");
            const limit_error = results[0].object.get("error").?.object;
            try cli.expectJsonString(limit_error.get("code"), "line_limit");
            try cli.expectJsonString(limit_error.get("message"), "line length limit exceeded");
            const failed = try cli.runWithClosedStdout(allocator, args, "");
            try std.testing.expectEqual(@as(u8, 4), failed.exit_code);
            try std.testing.expectEqualStrings("", failed.stdout);
            try std.testing.expectEqualStrings("error: standard output: I/O error\n", failed.stderr);
        }
    }
}

test "[cli] - [paired check]: aliases reject before reading or writing JSON" {
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

    for ([_][]const u8{payload}) |bytes| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = bytes });
        try tmp.dir.writeFile(io, .{ .sub_path = "copy.fastq", .data = bytes });
        for ([_][]const []const u8{ &.{ "check", "--paired" }, &.{ "check", "--paired", "--json" } }) |command| {
            for (aliases) |alias| {
                const args = try std.mem.concat(allocator, []const u8, &.{ command, &.{ path, alias } });
                const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{alias});
                try cli.expectResult(try cli.run(allocator, args), 2, "", diagnostic);
            }
            for (0..2) |stdin_side| {
                const inputs: [2][]const u8 = if (stdin_side == 0) .{ "-", path } else .{ path, "-" };
                const args = try std.mem.concat(allocator, []const u8, &.{ command, &inputs });
                const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{inputs[1]});
                const file = try tmp.dir.openFile(io, "input.fastq", .{});
                defer file.close(io);
                try cli.expectResult(try cli.runWithStdinFile(allocator, args, file), 2, "", diagnostic);
                if (bytes.len != 0) {
                    var first: [1]u8 = undefined;
                    try std.testing.expectEqual(@as(usize, 1), try file.readStreaming(io, &.{&first}));
                    try std.testing.expectEqual(bytes[0], first[0]);
                }
            }
            const args = try std.mem.concat(allocator, []const u8, &.{ command, &.{ path, copy } });
            const result = try cli.run(allocator, args);
            try std.testing.expectEqual(@as(u8, 0), result.exit_code);
            try std.testing.expectEqualStrings("", result.stderr);
            if (command.len == 2) {
                try std.testing.expectEqualStrings("", result.stdout);
            } else {
                const results = try cli.parseJsonResults(allocator, result.stdout, "z-fastq/check-v1");
                try std.testing.expectEqual(@as(usize, 1), results.len);
                try cli.expectJsonObjectKeys(results[0].object, &.{ "inputs", "status" });
                try expectJsonInputs(results[0].object.get("inputs"), path, copy);
                try cli.expectJsonString(results[0].object.get("status"), "ok");
            }
        }
    }
}

test "[cli] - [paired check]: identifier pairs pass under their policy or report P001" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    const Expected = enum { pass, p001 };
    const cases = [_]struct {
        header1: []const u8,
        header2: []const u8,
        policy: []const u8 = "illumina",
        expected: Expected = .pass,
    }{
        .{ .header1 = "cluster 1:N:0:ACGT", .header2 = "cluster 2:N:0:ACGT" },
        .{ .header1 = "cluster/1 opaque", .header2 = "cluster/2 opaque" },
        .{
            .header1 = "instrument:run:flowcell:1:2:3:4:UMI 1:N:0:index-a",
            .header2 = "instrument:run:flowcell:1:2:3:4:UMI 2:Y:0:index-b",
        },
        .{ .header1 = "cluster 1:N:0:index", .header2 = "cluster 2:N:0:index", .policy = "exact" },
        .{ .header1 = "cluster/2", .header2 = "cluster/1", .expected = .p001 },
        .{ .header1 = "cluster/1", .header2 = "cluster", .expected = .p001 },
        .{ .header1 = "cluster/1", .header2 = "cluster/2", .policy = "exact", .expected = .p001 },
    };
    for (cases) |case| {
        try writeFastq(&tmp, io, allocator, "r1.fastq", &.{case.header1});
        try writeFastq(&tmp, io, allocator, "r2.fastq", &.{case.header2});
        switch (case.expected) {
            .pass => try cli.expectResult(
                try cli.run(allocator, &.{ "check", "--paired", "--pair-names", case.policy, r1_path, r2_path }),
                0,
                "",
                "",
            ),
            .p001 => {
                const result = try cli.run(allocator, &.{ "check", "--paired", "--pair-names", case.policy, "--json", r1_path, r2_path });
                try std.testing.expectEqual(@as(u8, 1), result.exit_code);
                try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
                const results = try cli.parseJsonResults(allocator, result.stdout, "z-fastq/check-v1");
                try cli.expectJsonString(results[0].object.get("status"), "error");
                try cli.expectJsonString(results[0].object.get("error").?.object.get("code"), "P001");
            },
        }
    }
}

test "[cli] - [paired check]: double dash, empty inputs, interleaved input, and standard input pairs pass" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{"cluster 1:N:0:ACGT"});
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"cluster 2:N:0:ACGT"});
    const double_dash = try cli.run(
        allocator,
        &.{ "check", "--paired", "--", r1_path, r2_path },
    );
    try cli.expectResult(double_dash, 0, "", "");

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "" });
    const empty_pair = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    try cli.expectResult(empty_pair, 0, "", "");

    const empty_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        "",
        1,
    );
    try cli.expectResult(empty_interleaved, 0, "", "");

    const interleaved =
        "@cluster 1:N:0:index\r\nAC\r\n+\r\n!!\r\n" ++
        "@cluster 2:N:0:index\r\nGT\r\n+\r\n!!\r\n";
    const stdin_result = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        interleaved,
        1,
    );
    try cli.expectResult(stdin_result, 0, "", "");

    const interleaved_json = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        interleaved,
        3,
    );
    try std.testing.expectEqual(@as(u8, 0), interleaved_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), interleaved_json.stderr.len);
    const interleaved_results = try cli.parseJsonResults(allocator, interleaved_json.stdout, "z-fastq/check-v1");
    try expectCheckStatus(interleaved_results[0], "-", "ok");

    const exact_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--pair-names", "exact", "-" },
        interleaved,
        5,
    );
    try cli.expectResult(exact_interleaved, 0, "", "");

    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"stdin/2"});
    const stdin_r1 = try cli.runWithStdin(
        allocator,
        &.{ "check", "--paired", "-", r2_path },
        "@stdin/1\nA\n+\n!\n",
        1,
    );
    try cli.expectResult(stdin_r1, 0, "", "");
    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{"stdin/1"});
    const stdin_r2 = try cli.runWithStdin(
        allocator,
        &.{ "check", "--paired", r1_path, "-" },
        "@stdin/2\nA\n+\n!\n",
        1,
    );
    try cli.expectResult(stdin_r2, 0, "", "");
}

test "[cli] - [paired check]: gzip pairs and interleaved members pass with exact JSON" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, "@gzip/1\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &gzip2, "@gzip/2\nGT\n+\n!!\n", .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    const gzip_result = try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path });
    try cli.expectResult(gzip_result, 0, "", "");

    var interleaved_gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &interleaved_gzip, "@gzip/1\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &interleaved_gzip, "@gzip/2\nGT\n+\n!!\n", .{});
    const gzip_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        interleaved_gzip.items,
        3,
    );
    try cli.expectResult(gzip_interleaved, 0, "", "");
}

test "[cli] - [paired check]: a damaged gzip member is an I/O failure in text and JSON" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, "@gzip/1\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &gzip2, "@gzip/2\nGT\n+\n!!\n", .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    var interleaved_gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &interleaved_gzip, "@gzip/1\nAC\n+\n!!\n", .{});
    try cli.appendGzipMember(allocator, &interleaved_gzip, "@gzip/2\nGT\n+\n!!\n", .{});
    const corrupt_interleaved = try allocator.dupe(u8, interleaved_gzip.items);
    corrupt_interleaved[cli.gzipCrcOffset(corrupt_interleaved.len)] ^= 1;
    const damaged_interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "--json", "-" },
        corrupt_interleaved,
        7,
    );
    try std.testing.expectEqual(@as(u8, 3), damaged_interleaved.exit_code);
    try std.testing.expectEqual(@as(usize, 0), damaged_interleaved.stderr.len);
    const damaged_results = try cli.parseJsonResults(allocator, damaged_interleaved.stdout, "z-fastq/check-v1");
    try cli.expectJsonFailureWithoutLocation(
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
    const results = try cli.parseJsonResults(allocator, json_result.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), results.len);
    const result_object = results[0].object;
    try cli.expectJsonObjectKeys(result_object, &.{ "inputs", "status" });
    try expectJsonInputs(result_object.get("inputs"), r1_path, r2_path);
    try cli.expectJsonString(result_object.get("status"), "ok");

    const corrupt_r2 = try allocator.dupe(u8, gzip2.items);
    corrupt_r2[cli.gzipCrcOffset(corrupt_r2.len)] ^= 1;
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
    try cli.expectResult(damaged_pair, 3, "", damaged_pair_error);

    const damaged_pair_json = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 3), damaged_pair_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), damaged_pair_json.stderr.len);
    const damaged_pair_results = try cli.parseJsonResults(allocator, damaged_pair_json.stdout, "z-fastq/check-v1");
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

test "[cli] - [paired check]: P001 text and JSON report both identifiers with their normalized forms" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
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
    try cli.expectResult(human, 1, "", expected);

    const json_result = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    try std.testing.expectEqual(@as(u8, 1), json_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), json_result.stderr.len);
    const results = try cli.parseJsonResults(allocator, json_result.stdout, "z-fastq/check-v1");
    const result_object = results[0].object;
    try cli.expectJsonObjectKeys(result_object, &.{ "inputs", "status", "error" });
    try expectJsonInputs(result_object.get("inputs"), r1_path, r2_path);
    try cli.expectJsonString(result_object.get("status"), "error");
    try expectPairNameJson(result_object.get("error"), "cluster/1", "cluster", "other/2", "other");
}

test "[cli] - [paired check]: a NUL byte in an identifier is escaped in the P001 prefix" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@clu\x00ster/1\nA\n+\n!\n" });
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"other/2"});
    const escaped = try cli.run(allocator, &.{
        "check",
        "--paired",
        "--json",
        r1_path,
        r2_path,
    });
    const escaped_results = try cli.parseJsonResults(allocator, escaped.stdout, "z-fastq/check-v1");
    const escaped_tokens = escaped_results[0].object.get("error").?.object
        .get("first_tokens").?.array.items;
    try cli.expectJsonString(escaped_tokens[0].object.get("prefix"), "clu\\x00ster/1");
}

test "[cli] - [check]: interleaved identifier mismatch reports the pair record indexes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const interleaved_payload =
        "@cluster/1\nA\n+\n!\n" ++
        "@other/2\nA\n+\n!\n";
    const interleaved = try cli.runWithStdin(
        allocator,
        &.{ "check", "--interleaved", "-" },
        interleaved_payload,
        1,
    );
    try cli.expectResult(
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
    const interleaved_results = try cli.parseJsonResults(allocator, interleaved_json.stdout, "z-fastq/check-v1");
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
    const later_results = try cli.parseJsonResults(allocator, later_interleaved.stdout, "z-fastq/check-v1");
    const later_error = later_results[0].object.get("error").?.object;
    try cli.expectJsonInteger(later_error.get("pair_index"), 1);
    try expectJsonIntegerPair(later_error.get("record_indexes"), 2, 3);
}

test "[cli] - [paired check]: long identifiers are bounded in P001 reports" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const long1 = try allocator.alloc(u8, LONG_IDENTIFIER_BYTES);
    const long2 = try allocator.alloc(u8, LONG_IDENTIFIER_BYTES);
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
    try std.testing.expect(bounded.stdout.len < BOUNDED_REPORT_BYTES);
    const bounded_results = try cli.parseJsonResults(allocator, bounded.stdout, "z-fastq/check-v1");
    const bounded_error = bounded_results[0].object.get("error").?.object;
    const normalized = bounded_error.get("normalized_ids").?.array.items;
    try expectBoundedJson(normalized[0], long1[0..DIAGNOSTIC_PREFIX_BYTES], LONG_IDENTIFIER_BYTES, true);
    try expectBoundedJson(normalized[1], long2[0..DIAGNOSTIC_PREFIX_BYTES], LONG_IDENTIFIER_BYTES, true);

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
    try std.testing.expect(long_interleaved.stdout.len < BOUNDED_REPORT_BYTES);
    const long_interleaved_results = try cli.parseJsonResults(allocator, long_interleaved.stdout, "z-fastq/check-v1");
    const long_interleaved_error = long_interleaved_results[0].object
        .get("error").?.object;
    const long_first_tokens = long_interleaved_error.get("first_tokens").?.array.items;
    const long_normalized = long_interleaved_error.get("normalized_ids").?.array.items;
    try expectBoundedJson(long_first_tokens[0], long1[0..DIAGNOSTIC_PREFIX_BYTES], LONG_IDENTIFIER_BYTES + 2, true);
    try expectBoundedJson(long_first_tokens[1], long2[0..DIAGNOSTIC_PREFIX_BYTES], LONG_IDENTIFIER_BYTES + 2, true);
    try expectBoundedJson(long_normalized[0], long1[0..DIAGNOSTIC_PREFIX_BYTES], LONG_IDENTIFIER_BYTES, true);
    try expectBoundedJson(long_normalized[1], long2[0..DIAGNOSTIC_PREFIX_BYTES], LONG_IDENTIFIER_BYTES, true);
}

test "[cli] - [paired check]: missing mates and semantic failures report the exact first error" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    // A structural failure in R2 is reported before a semantic failure in R1.
    const Reported = enum { pair, r1, r2 };
    const cases = [_]struct {
        r1: []const u8,
        r2: []const u8,
        reported: Reported,
        message: []const u8,
    }{
        .{
            .r1 = "@first/1\nA\n+\n!\n@second/1\nA\n+\n!\n",
            .r2 = "@first/2\nA\n+\n!\n",
            .reported = .pair,
            .message = "P002: paired input is missing a mate " ++
                "(pair 1, remaining R1, last R1 record 1, last R2 record 0)",
        },
        .{
            .r1 = "",
            .r2 = "@only/2\nA\n+\n!\n",
            .reported = .pair,
            .message = "P002: paired input is missing a mate " ++
                "(pair 0, remaining R2, last R1 record none, last R2 record 0)",
        },
        .{
            .r1 = "@a/1\n.\n+\n!\n",
            .r2 = "@b/2\nA\n+\n!\n",
            .reported = .r1,
            .message = "S002: sequence byte is outside the selected alphabet (record 0, line 2, offset 5)",
        },
        .{
            .r1 = "@a/1\n.\n+\n!\n",
            .r2 = "@b/2\n.\n+\n!\n",
            .reported = .r1,
            .message = "S002: sequence byte is outside the selected alphabet (record 0, line 2, offset 5)",
        },
        .{
            .r1 = "@a/1\n.\n+\n!\n",
            .r2 = "@a/2\nA\n+\n",
            .reported = .r2,
            .message = "S004: unexpected end of file in quality line (record 0, line 4, offset 9)",
        },
    };
    for (cases) |case| {
        try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = case.r1 });
        try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = case.r2 });
        const expected = switch (case.reported) {
            .pair => try std.fmt.allocPrint(allocator, "error: {s} + {s}: {s}\n", .{ r1_path, r2_path, case.message }),
            .r1 => try std.fmt.allocPrint(allocator, "error: {s}: {s}\n", .{ r1_path, case.message }),
            .r2 => try std.fmt.allocPrint(allocator, "error: {s}: {s}\n", .{ r2_path, case.message }),
        };
        try cli.expectResult(try cli.run(allocator, &.{ "check", "--paired", r1_path, r2_path }), 1, "", expected);
    }
}

test "[cli] - [paired check]: a missing mate in JSON names the pair, side, and last records" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    try writeFastq(&tmp, io, allocator, "r1.fastq", &.{ "first/1", "second/1" });
    try writeFastq(&tmp, io, allocator, "r2.fastq", &.{"first/2"});

    const count_json = try cli.run(allocator, &.{ "check", "--paired", "--json", r1_path, r2_path });
    try std.testing.expectEqual(@as(u8, 1), count_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), count_json.stderr.len);
    const count_results = try cli.parseJsonResults(allocator, count_json.stdout, "z-fastq/check-v1");
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
    const results = try cli.parseJsonResults(allocator, command_result.stdout, "z-fastq/check-v1");
    try std.testing.expectEqual(@as(usize, 1), results.len);
    if (fixture.code) |code| {
        try cli.expectJsonFailure(
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
