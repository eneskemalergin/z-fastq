//! Installed CLI contracts for `z-fastq count`.

const std = @import("std");
const builtin = @import("builtin");

const ZFASTQ_BIN = if (builtin.os.tag == .windows)
    "zig-out\\bin\\z-fastq.exe"
else
    "zig-out/bin/z-fastq";
const FIXTURE_DIR = "tests/data/synthetic";

const FixtureExpect = struct {
    path: []const u8,
    exit_code: u8,
    stdout: ?[]const u8 = null,
    stderr_code: ?[]const u8 = null,
    stderr: ?[]const u8 = null,
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
    },
    .{ .path = "bad_qual_length.fastq", .exit_code = 1, .stderr_code = "S005" },
    .{ .path = "bad_header.fastq", .exit_code = 1, .stderr_code = "S003" },
    .{ .path = "truncated_record.fastq", .exit_code = 1, .stderr_code = "S004" },
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
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var argv: std.ArrayList([]const u8) = .empty;
    defer argv.deinit(allocator);
    try argv.append(allocator, ZFASTQ_BIN);
    try argv.appendSlice(allocator, command_args);

    var proc = try std.process.spawn(io, .{
        .argv = argv.items,
        .stdout = .pipe,
        .stderr = .pipe,
    });
    defer proc.kill(io);

    var stdout_buf: [4096]u8 = undefined;
    var stdout_reader = proc.stdout.?.reader(io, &stdout_buf);
    const stdout = try stdout_reader.interface.allocRemaining(allocator, .limited(1024 * 1024));
    var stderr_buf: [4096]u8 = undefined;
    var stderr_reader = proc.stderr.?.reader(io, &stderr_buf);
    const stderr = try stderr_reader.interface.allocRemaining(allocator, .limited(1024 * 1024));

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
        .stdout = .close,
        .stderr = .pipe,
    });
    defer proc.kill(io);

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

test "[cli] - [count]: valid fixture files produce expected counts" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (FIXTURES) |fixture| {
        const path = try fixturePath(allocator, fixture.path);
        const result = try runCount(allocator, path);

        try std.testing.expectEqual(fixture.exit_code, result.exit_code);
        if (fixture.stdout) |want| {
            try std.testing.expectEqualStrings(want, result.stdout);
        }
        if (fixture.stderr_code) |code| {
            try std.testing.expect(std.mem.indexOf(u8, result.stderr, code) != null);
        }
        if (fixture.stderr) |want| {
            try std.testing.expectEqualStrings(want, result.stderr);
        }
    }
}
test "[cli] - [count]: a missing input exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const result = try runCount(allocator, "tests/data/synthetic/does_not_exist.fastq");
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
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

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expect(std.mem.indexOf(u8, result.stderr, path) != null);
    try std.testing.expect(std.mem.endsWith(u8, result.stderr, ": failed to open file\n"));
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
    try std.testing.expect(std.mem.indexOf(u8, result.stderr, "S003") != null);
    try std.testing.expect(std.mem.indexOf(u8, result.stderr, "file not found") != null);
}

test "[cli] - [root]: help, version, and usage failures are exact" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const help = try runCli(allocator, &.{"--help"});
    try std.testing.expectEqual(@as(u8, 0), help.exit_code);
    try std.testing.expect(std.mem.indexOf(u8, help.stdout, "usage: z-fastq") != null);
    try std.testing.expectEqual(@as(usize, 0), help.stderr.len);

    const short_help = try runCli(allocator, &.{"-h"});
    try std.testing.expectEqual(@as(u8, 0), short_help.exit_code);
    try std.testing.expectEqualStrings(help.stdout, short_help.stdout);

    const version = try runCli(allocator, &.{"--version"});
    try std.testing.expectEqual(@as(u8, 0), version.exit_code);
    try std.testing.expectEqualStrings("z-fastq 0.0.2\n", version.stdout);
    try std.testing.expectEqual(@as(usize, 0), version.stderr.len);

    const short_version = try runCli(allocator, &.{"-V"});
    try std.testing.expectEqual(@as(u8, 0), short_version.exit_code);
    try std.testing.expectEqualStrings(version.stdout, short_version.stdout);

    const invalid = try runCli(allocator, &.{ "count", "--max-line-bytes", "nope" });
    try std.testing.expectEqual(@as(u8, 2), invalid.exit_code);
    try std.testing.expectEqualStrings(
        "error: invalid --max-line-bytes value\n",
        invalid.stderr,
    );

    const missing_command = try runCli(allocator, &.{});
    try std.testing.expectEqual(@as(u8, 2), missing_command.exit_code);
    try std.testing.expect(std.mem.startsWith(u8, missing_command.stderr, "usage: z-fastq"));

    const unknown_command = try runCli(allocator, &.{"unknown"});
    try std.testing.expectEqual(@as(u8, 2), unknown_command.exit_code);
    try std.testing.expect(std.mem.startsWith(
        u8,
        unknown_command.stderr,
        "error: unknown command: unknown\nusage: z-fastq",
    ));

    const missing_path = try runCli(allocator, &.{"count"});
    try std.testing.expectEqual(@as(u8, 2), missing_path.exit_code);
    try std.testing.expectEqualStrings(
        "error: count requires at least one file path\n",
        missing_path.stderr,
    );

    const missing_limit = try runCli(allocator, &.{ "count", "--max-line-bytes" });
    try std.testing.expectEqual(@as(u8, 2), missing_limit.exit_code);
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
    try std.testing.expectEqualStrings(
        "error: --max-line-bytes exceeds supported limit\n",
        overflow.stderr,
    );
}

test "[cli] - [output]: a closed stdout exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for ([_][]const []const u8{
        &.{"--help"},
        &.{"--version"},
        &.{ "count", "tests/data/synthetic/basic_valid.fastq" },
    }) |args| {
        const result = try runCliWithClosedStdout(allocator, args);
        try std.testing.expectEqual(@as(u8, 3), result.exit_code);
        try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    }
}

test "[cli] - [diagnostics]: untrusted command, option, and path bytes use escaped ASCII" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const command = try runCli(allocator, &.{"bad\n\\\x1b"});
    try std.testing.expectEqual(@as(u8, 2), command.exit_code);
    try std.testing.expect(std.mem.startsWith(
        u8,
        command.stderr,
        "error: unknown command: bad\\x0A\\\\\\x1B\nusage: z-fastq",
    ));

    const option = try runCli(allocator, &.{ "count", "-bad\t\\\x1b" });
    try std.testing.expectEqual(@as(u8, 2), option.exit_code);
    try std.testing.expectEqualStrings(
        "error: unknown count option: -bad\\x09\\\\\\x1B\n",
        option.stderr,
    );

    const path = try runCount(allocator, "unsafe\n\x1b.fastq");
    try std.testing.expectEqual(@as(u8, 3), path.exit_code);
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
    try std.testing.expectEqualStrings("error: unknown count option: --bogus\n", result.stderr);
}

test "[cli] - [count]: double dash treats a leading-hyphen argument as a path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const result = try runCli(arena.allocator(), &.{
        "count",
        "--",
        "tests/data/synthetic/basic_valid.fastq",
    });

    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    try std.testing.expectEqualStrings("5\n", result.stdout);
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

    const limited = try runCli(arena.allocator(), &.{
        "count",
        "--max-line-bytes",
        "262144",
        path,
    });
    try std.testing.expectEqual(@as(u8, 4), limited.exit_code);
    try std.testing.expect(
        std.mem.indexOf(u8, limited.stderr, "line length limit exceeded") != null,
    );
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
