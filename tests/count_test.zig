//! CLI integration tests for `z-fastq count`.

const std = @import("std");
const builtin = @import("builtin");

const ZFASTQ_BIN = if (builtin.os.tag == .windows) "zig-out\\bin\\z-fastq.exe" else "zig-out/bin/z-fastq";
const FIXTURE_DIR = "tests/data/synthetic";

const FixtureExpect = struct {
    path: []const u8,
    exit_code: u8,
    stdout: ?[]const u8 = null,
};

const fixtures = [_]FixtureExpect{
    .{ .path = "basic_valid.fastq", .exit_code = 0, .stdout = "5\n" },
    .{ .path = "crlf.fastq", .exit_code = 0, .stdout = "2\n" },
    .{ .path = "missing_final_newline.fastq", .exit_code = 0, .stdout = "1\n" },
    .{ .path = "bad_qual_length.fastq", .exit_code = 1 },
    .{ .path = "bad_header.fastq", .exit_code = 1 },
    .{ .path = "truncated_record.fastq", .exit_code = 1 },
};

fn fixturePath(allocator: std.mem.Allocator, name: []const u8) ![]const u8 {
    return std.fmt.allocPrint(allocator, "{s}/{s}", .{ FIXTURE_DIR, name });
}

fn runCount(
    allocator: std.mem.Allocator,
    path: []const u8,
) !struct { exit_code: u8, stdout: []u8 } {
    var threaded = std.Io.Threaded.init(allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var argv: std.ArrayList([]const u8) = .empty;
    defer argv.deinit(allocator);
    try argv.append(allocator, ZFASTQ_BIN);
    try argv.append(allocator, "count");
    try argv.append(allocator, path);

    var proc = try std.process.spawn(io, .{
        .argv = argv.items,
        .stdout = .pipe,
        .stderr = .pipe,
    });
    defer proc.kill(io);

    var read_buf: [4096]u8 = undefined;
    var stdout_reader = proc.stdout.?.reader(io, &read_buf);
    const stdout = try stdout_reader.interface.allocRemaining(allocator, .limited(1024 * 1024));

    const wait = try proc.wait(io);
    const exit_code: u8 = switch (wait) {
        .exited => |code| @intCast(code),
        else => return error.ChildProcessFailed,
    };

    return .{ .exit_code = exit_code, .stdout = stdout };
}

test "count: file expectations" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (fixtures) |fixture| {
        const path = try fixturePath(allocator, fixture.path);
        const result = try runCount(allocator, path);

        try std.testing.expectEqual(fixture.exit_code, result.exit_code);
        if (fixture.stdout) |want| {
            try std.testing.expectEqualStrings(want, result.stdout);
        }
    }
}
test "count: exits 3 when input file is missing" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const result = try runCount(allocator, "tests/data/synthetic/does_not_exist.fastq");
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
}
