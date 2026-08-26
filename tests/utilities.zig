//! Shared utilities for installed CLI tests.

const std = @import("std");

const ZFASTQ_BIN = "zig-out/bin/z-fastq";
const PROCESS_OUTPUT_LIMIT = 1024 * 1024;

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

pub fn expectJsonObjectKeys(
    object: std.json.ObjectMap,
    expected: []const []const u8,
) !void {
    try std.testing.expectEqual(expected.len, object.count());
    var iterator = object.iterator();
    var index: usize = 0;
    while (iterator.next()) |entry| : (index += 1) {
        try std.testing.expectEqualStrings(expected[index], entry.key_ptr.*);
    }
}

pub fn expectJsonDocument(
    value: *const std.json.Value,
    schema: []const u8,
) ![]const std.json.Value {
    const object = switch (value.*) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try expectJsonObjectKeys(object, &.{ "schema", "tool", "byte_strings", "results" });
    try expectJsonString(object.get("schema"), schema);
    try expectJsonString(object.get("byte_strings"), "escaped-bytes-v1");

    const tool_value = object.get("tool") orelse return error.UnexpectedJsonShape;
    const tool = switch (tool_value) {
        .object => |tool| tool,
        else => return error.UnexpectedJsonShape,
    };
    try expectJsonObjectKeys(tool, &.{ "name", "version" });
    try expectJsonString(tool.get("name"), "z-fastq");
    const version = switch (tool.get("version") orelse return error.UnexpectedJsonShape) {
        .string => |string| string,
        else => return error.UnexpectedJsonShape,
    };
    _ = std.SemanticVersion.parse(version) catch return error.UnexpectedJsonShape;

    const results_value = object.get("results") orelse return error.UnexpectedJsonShape;
    return switch (results_value) {
        .array => |results| results.items,
        else => error.UnexpectedJsonShape,
    };
}

pub fn expectJsonString(value: ?std.json.Value, expected: []const u8) !void {
    const actual = switch (value orelse return error.UnexpectedJsonShape) {
        .string => |string| string,
        else => return error.UnexpectedJsonShape,
    };
    try std.testing.expectEqualStrings(expected, actual);
}

pub fn expectJsonInteger(value: ?std.json.Value, expected: u64) !void {
    const integer = switch (value orelse return error.UnexpectedJsonShape) {
        .integer => |integer| integer,
        else => return error.UnexpectedJsonShape,
    };
    try std.testing.expectEqual(@as(i64, @intCast(expected)), integer);
}

pub fn expectJsonError(
    value: std.json.Value,
    input: []const u8,
    code: []const u8,
) !std.json.ObjectMap {
    const object = switch (value) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try expectJsonObjectKeys(object, &.{ "input", "status", "error" });
    try expectJsonString(object.get("input"), input);
    try expectJsonString(object.get("status"), "error");
    const error_value = object.get("error") orelse return error.UnexpectedJsonShape;
    const error_object = switch (error_value) {
        .object => |error_object| error_object,
        else => return error.UnexpectedJsonShape,
    };
    try expectJsonObjectKeys(
        error_object,
        &.{ "code", "message", "record_index", "byte_offset", "line_in_record" },
    );
    try expectJsonString(error_object.get("code"), code);
    return error_object;
}

pub fn expectJsonErrorLocation(
    object: std.json.ObjectMap,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u64,
) !void {
    try expectJsonInteger(object.get("record_index"), record_index);
    try expectJsonInteger(object.get("byte_offset"), byte_offset);
    try expectJsonInteger(object.get("line_in_record"), line_in_record);
}

pub fn expectJsonNullErrorLocation(object: std.json.ObjectMap) !void {
    inline for (&.{ "record_index", "byte_offset", "line_in_record" }) |field| {
        const value = object.get(field) orelse return error.UnexpectedJsonShape;
        try std.testing.expect(value == .null);
    }
}

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

    const child_stdin = proc.stdin.?;
    proc.stdin = null;
    var stdin_writer: StdinWriter = .{
        .io = io,
        .file = child_stdin,
        .data = stdin_data,
        .chunk_len = chunk_len,
    };
    var group: std.Io.Group = .init;
    defer group.cancel(io);
    group.async(io, writeChildStdin, .{&stdin_writer});

    const output = try collectCommandOutput(allocator, io, &proc);
    errdefer {
        allocator.free(output.stdout);
        allocator.free(output.stderr);
    }
    try group.await(io);
    const exit_code = try waitForExit(io, &proc);
    if (stdin_writer.err) |err| return err;
    return .{
        .exit_code = exit_code,
        .stdout = output.stdout,
        .stderr = output.stderr,
    };
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

    const child_stdin = proc.stdin.?;
    proc.stdin = null;
    var stdin_writer: StdinWriter = .{
        .io = io,
        .file = child_stdin,
        .data = stdin_data,
        .chunk_len = @max(stdin_data.len, 1),
    };
    var group: std.Io.Group = .init;
    defer group.cancel(io);
    group.async(io, writeChildStdin, .{&stdin_writer});

    var stderr_buf: [4096]u8 = undefined;
    var stderr_reader = proc.stderr.?.reader(io, &stderr_buf);
    const stderr = try stderr_reader.interface.allocRemaining(
        allocator,
        .limited(PROCESS_OUTPUT_LIMIT),
    );
    errdefer allocator.free(stderr);
    try group.await(io);
    const exit_code = try waitForExit(io, &proc);
    if (stdin_writer.err) |err| return err;
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

    const output = try collectCommandOutput(allocator, io, &proc);
    errdefer {
        allocator.free(output.stdout);
        allocator.free(output.stderr);
    }
    return .{
        .exit_code = try waitForExit(io, &proc),
        .stdout = output.stdout,
        .stderr = output.stderr,
    };
}

const CommandOutput = struct {
    stdout: []u8,
    stderr: []u8,
};

fn collectCommandOutput(
    allocator: std.mem.Allocator,
    io: std.Io,
    proc: *std.process.Child,
) !CommandOutput {
    var buffer: std.Io.File.MultiReader.Buffer(2) = undefined;
    var multi_reader: std.Io.File.MultiReader = undefined;
    multi_reader.init(
        allocator,
        io,
        buffer.toStreams(),
        &.{ proc.stdout.?, proc.stderr.? },
    );
    defer multi_reader.deinit();

    const stdout_reader = multi_reader.reader(0);
    const stderr_reader = multi_reader.reader(1);
    while (multi_reader.fill(4096, .none)) |_| {
        if (stdout_reader.buffered().len > PROCESS_OUTPUT_LIMIT or
            stderr_reader.buffered().len > PROCESS_OUTPUT_LIMIT)
        {
            return error.StreamTooLong;
        }
    } else |err| switch (err) {
        error.EndOfStream => {},
        else => |other| return other,
    }
    try multi_reader.checkAnyError();

    const stdout = try multi_reader.toOwnedSlice(0);
    errdefer allocator.free(stdout);
    const stderr = try multi_reader.toOwnedSlice(1);
    errdefer allocator.free(stderr);

    return .{ .stdout = stdout, .stderr = stderr };
}

fn waitForExit(io: std.Io, proc: *std.process.Child) !u8 {
    return switch (try proc.wait(io)) {
        .exited => |code| @intCast(code),
        else => return error.ChildProcessFailed,
    };
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

const StdinWriter = struct {
    io: std.Io,
    file: std.Io.File,
    data: []const u8,
    chunk_len: usize,
    err: ?std.Io.File.Writer.Error = null,
};

fn writeChildStdin(writer: *StdinWriter) void {
    writeAndCloseStdin(
        writer.io,
        writer.file,
        writer.data,
        writer.chunk_len,
    ) catch |err| {
        writer.err = err;
    };
}

fn writeAndCloseStdin(
    io: std.Io,
    file: std.Io.File,
    data: []const u8,
    chunk_len: usize,
) std.Io.File.Writer.Error!void {
    std.debug.assert(chunk_len > 0);
    defer file.close(io);
    var pos: usize = 0;
    write_stdin: while (pos < data.len) {
        const end = pos + @min(chunk_len, data.len - pos);
        std.Io.File.writeStreamingAll(file, io, data[pos..end]) catch |err| switch (err) {
            // Usage failures may close stdin before the test finishes feeding it.
            error.BrokenPipe => break :write_stdin,
            else => return err,
        };
        pos = end;
    }
}
