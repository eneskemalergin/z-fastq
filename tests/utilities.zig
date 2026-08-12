//! Shared utilities for installed CLI tests.

const std = @import("std");
const builtin = @import("builtin");

const ZFASTQ_BIN = if (builtin.os.tag == .windows)
    "zig-out\\bin\\z-fastq.exe"
else
    "zig-out/bin/z-fastq";

pub const CommandResult = struct {
    exit_code: u8,
    stdout: []u8,
    stderr: []u8,
};

pub const GzipOptions = struct {
    extra: []const u8 = "",
    name: ?[]const u8 = null,
    comment: ?[]const u8 = null,
    header_crc: bool = false,
};

pub fn appendGzipMember(
    allocator: std.mem.Allocator,
    output: *std.ArrayList(u8),
    payload: []const u8,
    options: GzipOptions,
) !void {
    if (payload.len > std.math.maxInt(u16) or options.extra.len > std.math.maxInt(u16)) {
        return error.TestFixtureTooLarge;
    }

    var flags: u8 = 0;
    if (options.header_crc) flags |= 0x02;
    if (options.extra.len != 0) flags |= 0x04;
    if (options.name != null) flags |= 0x08;
    if (options.comment != null) flags |= 0x10;

    const header_start = output.items.len;
    try output.appendSlice(allocator, &.{
        0x1f, 0x8b, 0x08, flags, 0, 0, 0, 0, 0, 0xff,
    });
    if (options.extra.len != 0) {
        try appendInt(u16, allocator, output, @intCast(options.extra.len));
        try output.appendSlice(allocator, options.extra);
    }
    if (options.name) |name| {
        try output.appendSlice(allocator, name);
        try output.append(allocator, 0);
    }
    if (options.comment) |comment| {
        try output.appendSlice(allocator, comment);
        try output.append(allocator, 0);
    }
    if (options.header_crc) {
        var crc: std.hash.Crc32 = .init();
        crc.update(output.items[header_start..]);
        try appendInt(u16, allocator, output, @truncate(crc.final()));
    }

    try output.append(allocator, 0x01);
    try appendInt(u16, allocator, output, @intCast(payload.len));
    try appendInt(u16, allocator, output, ~@as(u16, @intCast(payload.len)));
    try output.appendSlice(allocator, payload);

    var payload_crc: std.hash.Crc32 = .init();
    payload_crc.update(payload);
    try appendInt(u32, allocator, output, payload_crc.final());
    try appendInt(u32, allocator, output, @truncate(payload.len));
}

pub fn run(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
) !CommandResult {
    return runWithStdin(allocator, command_args, "", 1);
}

pub fn runWithStdin(
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

pub fn runWithClosedStdout(
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
    const stderr = try stderr_reader.interface.allocRemaining(
        allocator,
        .limited(1024 * 1024),
    );
    const wait = try proc.wait(io);
    const exit_code: u8 = switch (wait) {
        .exited => |code| @intCast(code),
        else => return error.ChildProcessFailed,
    };
    return .{
        .exit_code = exit_code,
        .stdout = try allocator.alloc(u8, 0),
        .stderr = stderr,
    };
}

pub fn runWithClosedStdin(
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

fn collectCommandResult(
    allocator: std.mem.Allocator,
    io: std.Io,
    proc: *std.process.Child,
) !CommandResult {
    var stdout_buf: [4096]u8 = undefined;
    var stdout_reader = proc.stdout.?.reader(io, &stdout_buf);
    const stdout = try stdout_reader.interface.allocRemaining(
        allocator,
        .limited(1024 * 1024),
    );
    errdefer allocator.free(stdout);
    var stderr_buf: [4096]u8 = undefined;
    var stderr_reader = proc.stderr.?.reader(io, &stderr_buf);
    const stderr = try stderr_reader.interface.allocRemaining(
        allocator,
        .limited(1024 * 1024),
    );
    errdefer allocator.free(stderr);

    const wait = try proc.wait(io);
    const exit_code: u8 = switch (wait) {
        .exited => |code| @intCast(code),
        else => return error.ChildProcessFailed,
    };
    return .{ .exit_code = exit_code, .stdout = stdout, .stderr = stderr };
}

fn appendInt(
    comptime T: type,
    allocator: std.mem.Allocator,
    output: *std.ArrayList(u8),
    value: T,
) !void {
    var bytes: [@sizeOf(T)]u8 = undefined;
    std.mem.writeInt(T, &bytes, value, .little);
    try output.appendSlice(allocator, &bytes);
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
