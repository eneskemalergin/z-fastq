//! Installed CLI contracts for `z-fastq count`.

const std = @import("std");
const builtin = @import("builtin");

const ZFASTQ_BIN = if (builtin.os.tag == .windows)
    "zig-out\\bin\\z-fastq.exe"
else
    "zig-out/bin/z-fastq";
const FIXTURE_DIR = "tests/data/synthetic";
const EXPECTED_USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count    Count records in plain FASTQ inputs
    \\
    \\General options:
    \\  -h, --help           Show this help message
    \\  -V, --version        Print version
    \\
    \\Count options:
    \\  --max-line-bytes N   Override default line length limit
    \\
    \\Count usage:
    \\  z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
;

const FixtureExpect = struct {
    path: []const u8,
    exit_code: u8,
    stdout: []const u8 = "",
    stderr: []const u8 = "",
    stdin_stderr: []const u8 = "",
};

const FIXTURES = [_]FixtureExpect{
    .{ .path = "basic_valid.fastq", .exit_code = 0, .stdout = "5\n" },
    .{ .path = "crlf.fastq", .exit_code = 0, .stdout = "2\n" },
    .{ .path = "missing_final_newline.fastq", .exit_code = 0, .stdout = "1\n" },
    .{
        .path = "bad_plus.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/bad_plus.fastq: S001: " ++
            "plus line must start with '+' (record 0, line 3, offset 15)\n",
        .stdin_stderr = "error: -: S001: " ++
            "plus line must start with '+' (record 0, line 3, offset 15)\n",
    },
    .{
        .path = "bad_qual_length.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/bad_qual_length.fastq: S005: " ++
            "sequence and quality lengths differ (record 0, line 4, offset 19)\n",
        .stdin_stderr = "error: -: S005: " ++
            "sequence and quality lengths differ (record 0, line 4, offset 19)\n",
    },
    .{
        .path = "bad_header.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/bad_header.fastq: S003: " ++
            "header line must start with '@' (record 0, line 1, offset 0)\n",
        .stdin_stderr = "error: -: S003: " ++
            "header line must start with '@' (record 0, line 1, offset 0)\n",
    },
    .{
        .path = "truncated_record.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/truncated_record.fastq: S004: " ++
            "unexpected end of file in quality line (record 0, line 4, offset 17)\n",
        .stdin_stderr = "error: -: S004: " ++
            "unexpected end of file in quality line (record 0, line 4, offset 17)\n",
    },
};

fn fixturePath(allocator: std.mem.Allocator, name: []const u8) ![]const u8 {
    return std.fmt.allocPrint(allocator, "{s}/{s}", .{ FIXTURE_DIR, name });
}

fn runCount(
    allocator: std.mem.Allocator,
    path: []const u8,
) !CommandResult {
    return runCli(allocator, &.{ "count", path });
}

const CommandResult = struct {
    exit_code: u8,
    stdout: []u8,
    stderr: []u8,
};

fn runCli(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
) !CommandResult {
    return runCliWithStdin(allocator, command_args, "", 1);
}

fn runCliWithStdin(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdin_data: []const u8,
    chunk_len: usize,
) !CommandResult {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var argv: std.ArrayList([]const u8) = .empty;
    defer argv.deinit(allocator);
    try argv.append(allocator, ZFASTQ_BIN);
    try argv.appendSlice(allocator, command_args);

    var proc = try std.process.spawn(io, .{
        .argv = argv.items,
        .stdin = .pipe,
        .stdout = .pipe,
        .stderr = .pipe,
    });
    defer proc.kill(io);

    try writeAndCloseStdin(io, &proc, stdin_data, chunk_len);
    return collectCommandResult(allocator, io, &proc);
}

fn collectCommandResult(
    allocator: std.mem.Allocator,
    io: std.Io,
    proc: *std.process.Child,
) !CommandResult {
    var stdout_buf: [4096]u8 = undefined;
    var stdout_reader = proc.stdout.?.reader(io, &stdout_buf);
    const stdout = try stdout_reader.interface.allocRemaining(allocator, .limited(1024 * 1024));
    errdefer allocator.free(stdout);
    var stderr_buf: [4096]u8 = undefined;
    var stderr_reader = proc.stderr.?.reader(io, &stderr_buf);
    const stderr = try stderr_reader.interface.allocRemaining(allocator, .limited(1024 * 1024));
    errdefer allocator.free(stderr);

    const wait = try proc.wait(io);
    const exit_code: u8 = switch (wait) {
        .exited => |code| @intCast(code),
        else => return error.ChildProcessFailed,
    };

    return .{ .exit_code = exit_code, .stdout = stdout, .stderr = stderr };
}

fn runCliWithClosedStdout(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdin_data: []const u8,
) !CommandResult {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var argv: std.ArrayList([]const u8) = .empty;
    defer argv.deinit(allocator);
    try argv.append(allocator, ZFASTQ_BIN);
    try argv.appendSlice(allocator, command_args);

    var proc = try std.process.spawn(io, .{
        .argv = argv.items,
        .stdin = .pipe,
        .stdout = .close,
        .stderr = .pipe,
    });
    defer proc.kill(io);

    try writeAndCloseStdin(io, &proc, stdin_data, @max(stdin_data.len, 1));
    var stderr_buf: [4096]u8 = undefined;
    var stderr_reader = proc.stderr.?.reader(io, &stderr_buf);
    const stderr = try stderr_reader.interface.allocRemaining(allocator, .limited(1024 * 1024));

    const wait = try proc.wait(io);
    const exit_code: u8 = switch (wait) {
        .exited => |code| @intCast(code),
        else => return error.ChildProcessFailed,
    };

    return .{ .exit_code = exit_code, .stdout = try allocator.alloc(u8, 0), .stderr = stderr };
}

fn runCliWithClosedStdin(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
) !CommandResult {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var argv: std.ArrayList([]const u8) = .empty;
    defer argv.deinit(allocator);
    try argv.append(allocator, ZFASTQ_BIN);
    try argv.appendSlice(allocator, command_args);

    var proc = try std.process.spawn(io, .{
        .argv = argv.items,
        .stdin = .close,
        .stdout = .pipe,
        .stderr = .pipe,
    });
    defer proc.kill(io);

    return collectCommandResult(allocator, io, &proc);
}

fn writeAndCloseStdin(
    io: std.Io,
    proc: *std.process.Child,
    data: []const u8,
    chunk_len: usize,
) !void {
    std.debug.assert(chunk_len > 0);
    var pos: usize = 0;
    while (pos < data.len) {
        const end = pos + @min(chunk_len, data.len - pos);
        try std.Io.File.writeStreamingAll(proc.stdin.?, io, data[pos..end]);
        pos = end;
    }
    proc.stdin.?.close(io);
    proc.stdin = null;
}

test "[cli] - [count]: files and fragmented stdin produce exact fixture results" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (FIXTURES) |fixture| {
        const path = try fixturePath(allocator, fixture.path);
        const file_result = try runCount(allocator, path);
        const data = try std.Io.Dir.cwd().readFileAlloc(
            std.testing.io,
            path,
            allocator,
            .limited(1024 * 1024),
        );
        const stdin_result = try runCliWithStdin(allocator, &.{ "count", "-" }, data, 1);

        try std.testing.expectEqual(fixture.exit_code, file_result.exit_code);
        try std.testing.expectEqualStrings(fixture.stdout, file_result.stdout);
        try std.testing.expectEqualStrings(fixture.stderr, file_result.stderr);
        try std.testing.expectEqual(fixture.exit_code, stdin_result.exit_code);
        try std.testing.expectEqualStrings(fixture.stdout, stdin_result.stdout);
        try std.testing.expectEqualStrings(fixture.stdin_stderr, stdin_result.stderr);
    }
}

test "[cli] - [count]: empty and over-limit stdin preserve exit classes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const empty = try runCliWithStdin(allocator, &.{ "count", "-" }, "", 1);
    try std.testing.expectEqual(@as(u8, 0), empty.exit_code);
    try std.testing.expectEqualStrings("0\n", empty.stdout);
    try std.testing.expectEqual(@as(usize, 0), empty.stderr.len);

    const limited = try runCliWithStdin(
        allocator,
        &.{ "count", "--max-line-bytes", "4", "-" },
        "@r\nAAAAA\n+\n!!!!!\n",
        2,
    );
    try std.testing.expectEqual(@as(u8, 4), limited.exit_code);
    try std.testing.expectEqual(@as(usize, 0), limited.stdout.len);
    try std.testing.expectEqualStrings(
        "error: -: line length limit exceeded\n",
        limited.stderr,
    );
}

test "[cli] - [count]: stdin is explicit and duplicate use fails before input processing" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const implicit = try runCliWithStdin(allocator, &.{"count"}, "", 1);
    try std.testing.expectEqual(@as(u8, 2), implicit.exit_code);
    try std.testing.expectEqual(@as(usize, 0), implicit.stdout.len);
    try std.testing.expectEqualStrings(
        "error: count requires at least one input\n",
        implicit.stderr,
    );

    const duplicate = try runCliWithStdin(
        allocator,
        &.{
            "count",
            "tests/data/synthetic/does_not_exist.fastq",
            "-",
            "-",
        },
        "",
        1,
    );
    try std.testing.expectEqual(@as(u8, 2), duplicate.exit_code);
    try std.testing.expectEqual(@as(usize, 0), duplicate.stdout.len);
    try std.testing.expectEqualStrings(
        "error: standard input may appear at most once\n",
        duplicate.stderr,
    );

    const after_double_dash = try runCliWithStdin(
        allocator,
        &.{ "count", "--", "-" },
        "@r\nA\n+\n!\n",
        1,
    );
    try std.testing.expectEqual(@as(u8, 0), after_double_dash.exit_code);
    try std.testing.expectEqualStrings("1\n", after_double_dash.stdout);
    try std.testing.expectEqual(@as(usize, 0), after_double_dash.stderr.len);
}

test "[cli] - [count]: a missing input exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const result = try runCount(allocator, "tests/data/synthetic/does_not_exist.fastq");
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(
        "error: tests/data/synthetic/does_not_exist.fastq: file not found\n",
        result.stderr,
    );
}

test "[cli] - [count]: a long missing path is reported without truncation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try allocator.alloc(u8, 600);
    @memset(path, 'x');

    const result = try runCount(allocator, path);
    const expected = try std.fmt.allocPrint(
        allocator,
        "error: {s}: failed to open file\n",
        .{path},
    );

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(expected, result.stderr);
}

test "[cli] - [count]: all paths run and the highest exit class wins" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const result = try runCli(allocator, &.{
        "count",
        "tests/data/synthetic/missing_final_newline.fastq",
        "tests/data/synthetic/bad_header.fastq",
        "tests/data/synthetic/does_not_exist.fastq",
    });
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqualStrings("1\n", result.stdout);
    try std.testing.expectEqualStrings(
        "error: tests/data/synthetic/bad_header.fastq: S003: " ++
            "header line must start with '@' (record 0, line 1, offset 0)\n" ++
            "error: tests/data/synthetic/does_not_exist.fastq: file not found\n",
        result.stderr,
    );
}

test "[cli] - [count]: files around stdin retain output order and exit precedence" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const result = try runCliWithStdin(
        arena.allocator(),
        &.{
            "count",
            "tests/data/synthetic/missing_final_newline.fastq",
            "-",
            "tests/data/synthetic/does_not_exist.fastq",
            "tests/data/synthetic/basic_valid.fastq",
        },
        "not-a-header\nA\n+\n!\n",
        3,
    );

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqualStrings("1\n5\n", result.stdout);
    try std.testing.expectEqualStrings(
        "error: -: S003: header line must start with '@' " ++
            "(record 0, line 1, offset 0)\n" ++
            "error: tests/data/synthetic/does_not_exist.fastq: file not found\n",
        result.stderr,
    );
}

test "[cli] - [root]: help, version, and usage failures are exact" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const help = try runCli(allocator, &.{"--help"});
    try std.testing.expectEqual(@as(u8, 0), help.exit_code);
    try std.testing.expectEqualStrings(EXPECTED_USAGE, help.stdout);
    try std.testing.expectEqual(@as(usize, 0), help.stderr.len);

    const short_help = try runCli(allocator, &.{"-h"});
    try std.testing.expectEqual(@as(u8, 0), short_help.exit_code);
    try std.testing.expectEqualStrings(help.stdout, short_help.stdout);
    try std.testing.expectEqual(@as(usize, 0), short_help.stderr.len);

    const version = try runCli(allocator, &.{"--version"});
    try std.testing.expectEqual(@as(u8, 0), version.exit_code);
    try std.testing.expectEqualStrings("z-fastq 0.0.3\n", version.stdout);
    try std.testing.expectEqual(@as(usize, 0), version.stderr.len);

    const short_version = try runCli(allocator, &.{"-V"});
    try std.testing.expectEqual(@as(u8, 0), short_version.exit_code);
    try std.testing.expectEqualStrings(version.stdout, short_version.stdout);
    try std.testing.expectEqual(@as(usize, 0), short_version.stderr.len);

    const invalid = try runCli(allocator, &.{ "count", "--max-line-bytes", "nope" });
    try std.testing.expectEqual(@as(u8, 2), invalid.exit_code);
    try std.testing.expectEqual(@as(usize, 0), invalid.stdout.len);
    try std.testing.expectEqualStrings(
        "error: invalid --max-line-bytes value\n",
        invalid.stderr,
    );

    const missing_command = try runCli(allocator, &.{});
    try std.testing.expectEqual(@as(u8, 2), missing_command.exit_code);
    try std.testing.expectEqual(@as(usize, 0), missing_command.stdout.len);
    try std.testing.expectEqualStrings(EXPECTED_USAGE, missing_command.stderr);

    const unknown_command = try runCli(allocator, &.{"unknown"});
    try std.testing.expectEqual(@as(u8, 2), unknown_command.exit_code);
    try std.testing.expectEqual(@as(usize, 0), unknown_command.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unknown command: unknown\n" ++ EXPECTED_USAGE,
        unknown_command.stderr,
    );

    const missing_path = try runCli(allocator, &.{"count"});
    try std.testing.expectEqual(@as(u8, 2), missing_path.exit_code);
    try std.testing.expectEqual(@as(usize, 0), missing_path.stdout.len);
    try std.testing.expectEqualStrings(
        "error: count requires at least one input\n",
        missing_path.stderr,
    );

    const missing_limit = try runCli(allocator, &.{ "count", "--max-line-bytes" });
    try std.testing.expectEqual(@as(u8, 2), missing_limit.exit_code);
    try std.testing.expectEqual(@as(usize, 0), missing_limit.stdout.len);
    try std.testing.expectEqualStrings(
        "error: --max-line-bytes requires a value\n",
        missing_limit.stderr,
    );

    const overflow = try runCli(allocator, &.{
        "count",
        "--max-line-bytes",
        "340282366920938463463374607431768211456",
    });
    try std.testing.expectEqual(@as(u8, 4), overflow.exit_code);
    try std.testing.expectEqual(@as(usize, 0), overflow.stdout.len);
    try std.testing.expectEqualStrings(
        "error: --max-line-bytes exceeds supported limit\n",
        overflow.stderr,
    );
}

test "[cli] - [output]: a closed stdout exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for ([_]struct {
        args: []const []const u8,
        stdin: []const u8 = "",
        exit_code: u8 = 3,
        stderr: []const u8 = "",
    }{
        .{ .args = &.{"--help"} },
        .{ .args = &.{"--version"} },
        .{ .args = &.{ "count", "tests/data/synthetic/basic_valid.fastq" } },
        .{ .args = &.{ "count", "-" }, .stdin = "@r\nA\n+\n!\n" },
        .{
            .args = &.{
                "count",
                "--max-line-bytes",
                "16",
                "-",
                "tests/data/synthetic/missing_final_newline.fastq",
            },
            .stdin = "@r\nAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!\n",
            .exit_code = 4,
            .stderr = "error: -: line length limit exceeded\n",
        },
    }) |case| {
        const result = try runCliWithClosedStdout(allocator, case.args, case.stdin);
        try std.testing.expectEqual(case.exit_code, result.exit_code);
        try std.testing.expectEqualStrings(case.stderr, result.stderr);
    }
}

test "[cli] - [count]: a closed stdin exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const result = try runCliWithClosedStdin(arena.allocator(), &.{ "count", "-" });

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings("error: -: I/O error\n", result.stderr);
}

test "[cli] - [diagnostics]: untrusted command, option, and path bytes use escaped ASCII" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const command = try runCli(allocator, &.{"bad\n\\\x1b"});
    try std.testing.expectEqual(@as(u8, 2), command.exit_code);
    try std.testing.expectEqual(@as(usize, 0), command.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unknown command: bad\\x0A\\\\\\x1B\n" ++ EXPECTED_USAGE,
        command.stderr,
    );

    const option = try runCli(allocator, &.{ "count", "-bad\t\\\x1b" });
    try std.testing.expectEqual(@as(u8, 2), option.exit_code);
    try std.testing.expectEqual(@as(usize, 0), option.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unknown count option: -bad\\x09\\\\\\x1B\n",
        option.stderr,
    );

    const path = try runCount(allocator, "unsafe\n\x1b.fastq");
    try std.testing.expectEqual(@as(u8, 3), path.exit_code);
    try std.testing.expectEqual(@as(usize, 0), path.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unsafe\\x0A\\x1B.fastq: file not found\n",
        path.stderr,
    );
}

test "[cli] - [count]: an unknown option exits with usage status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const result = try runCli(arena.allocator(), &.{ "count", "--bogus" });

    try std.testing.expectEqual(@as(u8, 2), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings("error: unknown count option: --bogus\n", result.stderr);
}

test "[cli] - [count]: double dash treats a leading-hyphen argument as a path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const path = "-z-fastq-v002-missing/record.fastq";
    const result = try runCli(arena.allocator(), &.{
        "count",
        "--",
        path,
    });

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(
        "error: " ++ path ++ ": file not found\n",
        result.stderr,
    );
}

test "[cli] - [count]: a line may exceed the read buffer within its limit" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const line_len: usize = 262145;

    {
        const file = try tmp.dir.createFile(io, "long.fastq", .{});
        defer file.close(io);
        try std.Io.File.writeStreamingAll(file, io, "@long\n");
        var chunk: [4096]u8 = undefined;
        @memset(&chunk, 'A');
        var remaining = line_len;
        while (remaining > 0) {
            const n = @min(remaining, chunk.len);
            try std.Io.File.writeStreamingAll(file, io, chunk[0..n]);
            remaining -= n;
        }
        try std.Io.File.writeStreamingAll(file, io, "\n+\n");
        @memset(&chunk, '!');
        remaining = line_len;
        while (remaining > 0) {
            const n = @min(remaining, chunk.len);
            try std.Io.File.writeStreamingAll(file, io, chunk[0..n]);
            remaining -= n;
        }
        try std.Io.File.writeStreamingAll(file, io, "\n");
    }

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const path = try std.fmt.allocPrint(
        arena.allocator(),
        ".zig-cache/tmp/{s}/long.fastq",
        .{tmp.sub_path},
    );
    const result = try runCount(arena.allocator(), path);
    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    try std.testing.expectEqualStrings("1\n", result.stdout);
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);

    const limited = try runCli(arena.allocator(), &.{
        "count",
        "--max-line-bytes",
        "262144",
        path,
    });
    const limited_expected = try std.fmt.allocPrint(
        arena.allocator(),
        "error: {s}: line length limit exceeded\n",
        .{path},
    );
    try std.testing.expectEqual(@as(u8, 4), limited.exit_code);
    try std.testing.expectEqual(@as(usize, 0), limited.stdout.len);
    try std.testing.expectEqualStrings(limited_expected, limited.stderr);
}

test "[cli] - [count]: a lone CR at EOF remains quality content" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const file = try tmp.dir.createFile(io, "lone-cr.fastq", .{});
    defer file.close(io);
    try file.writePositionalAll(io, "@r\nA\n+\n!\r", 0);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const path = try std.fmt.allocPrint(
        arena.allocator(),
        ".zig-cache/tmp/{s}/lone-cr.fastq",
        .{tmp.sub_path},
    );

    const result = try runCount(arena.allocator(), path);
    const expected = try std.fmt.allocPrint(
        arena.allocator(),
        "error: {s}: S005: sequence and quality lengths differ " ++
            "(record 0, line 4, offset 7)\n",
        .{path},
    );

    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(expected, result.stderr);
}
