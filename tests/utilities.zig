//! Shared utilities for installed CLI tests.

const std = @import("std");

const ZFASTQ_BIN = "zig-out/bin/z-fastq";
const PROCESS_OUTPUT_LIMIT = 1024 * 1024;
pub const PROCESS_TIMEOUT = std.Io.Duration.fromSeconds(30);

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

pub const ProcessDeadline = struct {
    io: std.Io,
    process: *std.process.Child,
    deadline: std.Io.Clock.Timestamp,
    pid_file: std.Io.File,
    finished: bool = false,

    pub fn init(
        io: std.Io,
        process: *std.process.Child,
        duration: std.Io.Duration,
    ) !ProcessDeadline {
        const raw_pid_fd = std.os.linux.pidfd_open(process.id.?, 0);
        const pid_fd = switch (std.os.linux.errno(raw_pid_fd)) {
            .SUCCESS => std.math.cast(std.posix.fd_t, raw_pid_fd) orelse
                return error.ChildProcessControlFailed,
            else => return error.ChildProcessControlFailed,
        };
        return .{
            .io = io,
            .process = process,
            .deadline = std.Io.Clock.Timestamp.now(io, .awake).addDuration(.{
                .raw = duration,
                .clock = .awake,
            }),
            .pid_file = .{ .handle = pid_fd, .flags = .{ .nonblocking = false } },
        };
    }

    pub fn deinit(self: *ProcessDeadline) void {
        if (!self.finished) self.process.kill(self.io);
        self.pid_file.close(self.io);
        self.* = undefined;
    }

    pub fn timeout(self: *const ProcessDeadline) std.Io.Timeout {
        return .{ .deadline = self.deadline };
    }

    pub fn expired(self: *const ProcessDeadline) bool {
        return std.Io.Clock.Timestamp.now(self.io, .awake).raw.nanoseconds >=
            self.deadline.raw.nanoseconds;
    }

    pub fn waitForExit(self: *ProcessDeadline) !u8 {
        var poll_fd = [1]std.posix.pollfd{.{
            .fd = self.pid_file.handle,
            .events = std.posix.POLL.IN,
            .revents = 0,
        }};
        const timeout_ms = self.remainingMilliseconds() orelse return error.Timeout;
        if (try std.posix.poll(&poll_fd, timeout_ms) == 0) return error.Timeout;
        if (poll_fd[0].revents & std.posix.POLL.NVAL != 0) {
            return error.ChildProcessControlFailed;
        }

        const termination = try self.process.wait(self.io);
        self.finished = true;
        return switch (termination) {
            .exited => |code| code,
            else => error.ChildProcessFailed,
        };
    }

    pub fn failTimeout(
        self: *ProcessDeadline,
    ) error{ ChildProcessCleanupFailed, ChildProcessTimedOut } {
        self.process.kill(self.io);
        self.finished = true;
        if (self.process.id != null or
            self.process.stdin != null or
            self.process.stdout != null or
            self.process.stderr != null)
        {
            return error.ChildProcessCleanupFailed;
        }
        return error.ChildProcessTimedOut;
    }

    fn remainingMilliseconds(self: *const ProcessDeadline) ?i32 {
        const now = std.Io.Clock.Timestamp.now(self.io, .awake);
        const remaining_ns = self.deadline.raw.nanoseconds - now.raw.nanoseconds;
        if (remaining_ns <= 0) return null;
        const milliseconds = @divFloor(
            remaining_ns + std.time.ns_per_ms - 1,
            std.time.ns_per_ms,
        );
        return @intCast(@min(milliseconds, std.math.maxInt(i32)));
    }
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
    return runInstalled(
        allocator,
        command_args,
        .{ .data = .{ .bytes = stdin_data, .chunk_len = chunk_len } },
        .capture,
        PROCESS_TIMEOUT,
    );
}

pub fn runWithClosedStdout(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdin_data: []const u8,
) !CommandResult {
    return runInstalled(
        allocator,
        command_args,
        .{ .data = .{
            .bytes = stdin_data,
            .chunk_len = @max(stdin_data.len, 1),
        } },
        .closed,
        PROCESS_TIMEOUT,
    );
}

pub fn runWithClosedStdin(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
) !CommandResult {
    return runInstalled(
        allocator,
        command_args,
        .closed,
        .capture,
        PROCESS_TIMEOUT,
    );
}

pub fn runStalledStdinForTest(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    timeout: std.Io.Duration,
) !void {
    const result = try runInstalled(
        allocator,
        command_args,
        .stalled,
        .capture,
        timeout,
    );
    defer allocator.free(result.stdout);
    defer allocator.free(result.stderr);
    return error.ExpectedChildProcessTimeout;
}

pub fn finishSpawned(
    allocator: std.mem.Allocator,
    io: std.Io,
    process: *std.process.Child,
    process_deadline: *ProcessDeadline,
) !CommandResult {
    const output = collectCommandOutput(
        allocator,
        io,
        process,
        process_deadline.timeout(),
    ) catch |err| switch (err) {
        error.Timeout => return process_deadline.failTimeout(),
        else => |other| return other,
    };
    errdefer {
        allocator.free(output.stdout);
        allocator.free(output.stderr);
    }
    return .{
        .exit_code = process_deadline.waitForExit() catch |err| switch (err) {
            error.Timeout => return process_deadline.failTimeout(),
            else => |other| return other,
        },
        .stdout = output.stdout,
        .stderr = output.stderr,
    };
}

const StdinMode = union(enum) {
    closed,
    stalled,
    data: struct {
        bytes: []const u8,
        chunk_len: usize,
    },
};

const StdoutMode = enum {
    capture,
    closed,
};

fn runInstalled(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdin_mode: StdinMode,
    stdout_mode: StdoutMode,
    timeout: std.Io.Duration,
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
        .stdin = switch (stdin_mode) {
            .closed => .close,
            .stalled, .data => .pipe,
        },
        .stdout = switch (stdout_mode) {
            .capture => .pipe,
            .closed => .close,
        },
        .stderr = .pipe,
    });
    defer proc.kill(io);

    var process_deadline = try ProcessDeadline.init(io, &proc, timeout);
    defer process_deadline.deinit();

    var stdin_writer: StdinWriter = undefined;
    var stdin_group: std.Io.Group = .init;
    defer stdin_group.cancel(io);
    const has_stdin_writer = switch (stdin_mode) {
        .closed, .stalled => false,
        .data => |data| writer: {
            const child_stdin = proc.stdin.?;
            proc.stdin = null;
            stdin_writer = .{
                .io = io,
                .file = child_stdin,
                .data = data.bytes,
                .chunk_len = data.chunk_len,
            };
            stdin_group.async(io, writeChildStdin, .{&stdin_writer});
            break :writer true;
        },
    };

    const output = switch (stdout_mode) {
        .capture => collectCommandOutput(
            allocator,
            io,
            &proc,
            process_deadline.timeout(),
        ),
        .closed => collectCommandStderr(
            allocator,
            io,
            &proc,
            process_deadline.timeout(),
        ),
    } catch |err| switch (err) {
        error.Timeout => {
            stdin_group.cancel(io);
            return process_deadline.failTimeout();
        },
        else => |other| return other,
    };
    errdefer {
        allocator.free(output.stdout);
        allocator.free(output.stderr);
    }
    const exit_code = process_deadline.waitForExit() catch |err| switch (err) {
        error.Timeout => {
            stdin_group.cancel(io);
            return process_deadline.failTimeout();
        },
        else => |other| return other,
    };
    stdin_group.cancel(io);
    if (has_stdin_writer) {
        if (stdin_writer.err) |err| return err;
    }
    return .{
        .exit_code = exit_code,
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
    timeout: std.Io.Timeout,
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
    while (multi_reader.fill(4096, timeout)) |_| {
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

fn collectCommandStderr(
    allocator: std.mem.Allocator,
    io: std.Io,
    proc: *std.process.Child,
    timeout: std.Io.Timeout,
) !CommandOutput {
    var buffer: std.Io.File.MultiReader.Buffer(1) = undefined;
    var multi_reader: std.Io.File.MultiReader = undefined;
    multi_reader.init(
        allocator,
        io,
        buffer.toStreams(),
        &.{proc.stderr.?},
    );
    defer multi_reader.deinit();

    const stderr_reader = multi_reader.reader(0);
    while (multi_reader.fill(4096, timeout)) |_| {
        if (stderr_reader.buffered().len > PROCESS_OUTPUT_LIMIT) {
            return error.StreamTooLong;
        }
    } else |err| switch (err) {
        error.EndOfStream => {},
        else => |other| return other,
    }
    try multi_reader.checkAnyError();

    const stdout = try allocator.alloc(u8, 0);
    errdefer allocator.free(stdout);
    const stderr = try multi_reader.toOwnedSlice(0);
    errdefer allocator.free(stderr);
    return .{ .stdout = stdout, .stderr = stderr };
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

fn writeChildStdin(writer: *StdinWriter) std.Io.Cancelable!void {
    writeAndCloseStdin(
        writer.io,
        writer.file,
        writer.data,
        writer.chunk_len,
    ) catch |err| {
        if (err == error.Canceled) return error.Canceled;
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
