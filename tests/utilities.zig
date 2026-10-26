//! Shared utilities for installed CLI tests.

const std = @import("std");

/// The CLI's JSON output buffer: a longer document is flushed before it ends.
pub const JSON_OUTPUT_BUFFER_BYTES = 16 * 1024;

/// Tracked synthetic fixtures, relative to the repository root.
pub const FIXTURE_DIR = "tests/data/synthetic";

pub const EXPECTED_VERSION = @import("test_options").package_version;
pub const ZFASTQ_BIN = "zig-out/bin/z-fastq";
/// Captured stream limit of the harness: a longer stream is a test failure, not a truncated result.
pub const PROCESS_OUTPUT_LIMIT = 1024 * 1024;
pub const PROCESS_TIMEOUT = std.Io.Duration.fromSeconds(30);

pub const CommandResult = struct {
    exit_code: u8,
    stdout: []u8,
    stderr: []u8,
};

pub fn expectResult(
    result: CommandResult,
    exit_code: u8,
    stdout: []const u8,
    stderr: []const u8,
) !void {
    try std.testing.expectEqual(exit_code, result.exit_code);
    try std.testing.expectEqualStrings(stdout, result.stdout);
    try std.testing.expectEqualStrings(stderr, result.stderr);
}

pub fn expectFailedPrefix(
    result: CommandResult,
    exit_code: u8,
    valid_stdout: []const u8,
    stderr: []const u8,
) !void {
    try std.testing.expectEqual(exit_code, result.exit_code);
    try expectPrefix(valid_stdout, result.stdout);
    try std.testing.expectEqualStrings(stderr, result.stderr);
}

pub fn expectPrefix(valid: []const u8, actual: []const u8) !void {
    try std.testing.expect(actual.len <= valid.len);
    try std.testing.expectEqualSlices(u8, valid[0..actual.len], actual);
}

pub fn tempPath(
    allocator: std.mem.Allocator,
    sub_path: []const u8,
    name: []const u8,
) ![]const u8 {
    return std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}", .{ sub_path, name });
}

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
    try expectJsonString(tool.get("version"), EXPECTED_VERSION);

    const results_value = object.get("results") orelse return error.UnexpectedJsonShape;
    return switch (results_value) {
        .array => |results| results.items,
        else => error.UnexpectedJsonShape,
    };
}

/// Parses a versioned JSON document and returns its results. The parsed tree is
/// allocated from `allocator` and not freed separately, so pass an arena.
pub fn parseJsonResults(
    allocator: std.mem.Allocator,
    bytes: []const u8,
    schema: []const u8,
) ![]const std.json.Value {
    const value = try std.json.parseFromSliceLeaky(std.json.Value, allocator, bytes, .{});
    return expectJsonDocument(&value, schema);
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

fn expectJsonError(
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

/// A failed result with the given code and record location.
pub fn expectJsonFailure(
    value: std.json.Value,
    input: []const u8,
    code: []const u8,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u64,
) !void {
    const error_object = try expectJsonError(value, input, code);
    try expectJsonErrorLocation(error_object, record_index, byte_offset, line_in_record);
}

/// A failed result with the given code and message and no record location.
pub fn expectJsonFailureWithoutLocation(
    value: std.json.Value,
    input: []const u8,
    code: []const u8,
    message: []const u8,
) !void {
    const error_object = try expectJsonError(value, input, code);
    try expectJsonString(error_object.get("message"), message);
    try expectJsonNullErrorLocation(error_object);
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

/// Offset of the CRC-32 field in the eight-byte gzip trailer that ends a member
/// of `member_len` bytes; the input size follows it.
pub fn gzipCrcOffset(member_len: usize) usize {
    return member_len - GZIP_TRAILER_BYTES;
}

/// Offset of the little-endian input size (ISIZE) field in the gzip trailer.
pub fn gzipSizeOffset(member_len: usize) usize {
    return member_len - GZIP_TRAILER_BYTES + 4;
}

const GZIP_TRAILER_BYTES = 8;

pub fn fixturePath(allocator: std.mem.Allocator, name: []const u8) ![]const u8 {
    return std.fmt.allocPrint(allocator, "{s}/{s}", .{ FIXTURE_DIR, name });
}

/// Opens an existing file for writing in append mode, the way an output
/// redirection does, so a write lands after any bytes already present.
pub fn openAppendOutput(dir: std.Io.Dir, sub_path: []const u8) !std.Io.File {
    return .{
        .handle = try std.posix.openat(dir.handle, sub_path, .{
            .ACCMODE = .WRONLY,
            .APPEND = true,
            .CLOEXEC = true,
        }, 0),
        .flags = .{ .nonblocking = false },
    };
}

pub fn run(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
) !CommandResult {
    return runWith(allocator, command_args, .{});
}

pub fn runWithStdin(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdin_data: []const u8,
    chunk_len: usize,
) !CommandResult {
    return runWith(allocator, command_args, .{
        .stdin = .{ .data = .{ .bytes = stdin_data, .chunk_len = chunk_len } },
    });
}

pub fn runWithClosedStdout(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdin_data: []const u8,
) !CommandResult {
    return runWith(allocator, command_args, .{
        .stdin = .{ .data = .{ .bytes = stdin_data, .chunk_len = @max(stdin_data.len, 1) } },
        .stdout = .closed,
    });
}

pub fn runWithStdinFile(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    file: std.Io.File,
) !CommandResult {
    return runWith(allocator, command_args, .{ .stdin = .{ .file = file } });
}

pub fn runWithStdoutFile(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdout_file: std.Io.File,
    stdin_file: ?std.Io.File,
) !CommandResult {
    return runWith(allocator, command_args, .{
        .stdin = if (stdin_file) |file| .{ .file = file } else .{ .data = .{ .bytes = "", .chunk_len = 1 } },
        .stdout = .{ .file = stdout_file },
    });
}

/// Runs the command with standard output and error redirected to files and
/// standard input closed, and returns only the exit status.
pub fn runWithOutputFiles(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    stdout_file: std.Io.File,
    stderr_file: std.Io.File,
) !u8 {
    const result = try runWith(allocator, command_args, .{
        .stdin = .closed,
        .stdout = .{ .file = stdout_file },
        .stderr = .{ .file = stderr_file },
    });
    return result.exit_code;
}

pub fn runWithClosedStdin(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
) !CommandResult {
    return runWith(allocator, command_args, .{ .stdin = .closed });
}

/// Completes a child that the caller spawned with piped stdout and stderr.
pub fn finishSpawned(
    allocator: std.mem.Allocator,
    io: std.Io,
    process: *std.process.Child,
    process_deadline: *ProcessDeadline,
) !CommandResult {
    const outputs = collectPipes(
        2,
        allocator,
        io,
        .{ process.stdout.?, process.stderr.? },
        process_deadline.timeout(),
        PROCESS_OUTPUT_LIMIT,
    ) catch |err| switch (err) {
        error.Timeout => return process_deadline.failTimeout(),
        else => |other| return other,
    };
    errdefer {
        allocator.free(outputs[0]);
        allocator.free(outputs[1]);
    }
    return .{
        .exit_code = process_deadline.waitForExit() catch |err| switch (err) {
            error.Timeout => return process_deadline.failTimeout(),
            else => |other| return other,
        },
        .stdout = outputs[0],
        .stderr = outputs[1],
    };
}

pub const StdinMode = union(enum) {
    closed,
    stalled,
    file: std.Io.File,
    data: struct {
        bytes: []const u8,
        chunk_len: usize,
    },
};

pub const StdoutMode = union(enum) {
    capture,
    closed,
    file: std.Io.File,
};

pub const StderrMode = union(enum) {
    capture,
    file: std.Io.File,
};

pub const RunOptions = struct {
    stdin: StdinMode = .{ .data = .{ .bytes = "", .chunk_len = 1 } },
    stdout: StdoutMode = .capture,
    stderr: StderrMode = .capture,
    /// Time allowed for the child to exit; a longer run is killed and reported as a timeout.
    timeout: std.Io.Duration = PROCESS_TIMEOUT,
    /// A captured stream longer than this many bytes aborts the run with `error.StreamTooLong`.
    output_limit: usize = PROCESS_OUTPUT_LIMIT,
};

/// The single spawn path for installed CLI runs. Captured streams are read
/// until end of stream, so a child is never blocked on a full pipe.
pub fn runWith(
    allocator: std.mem.Allocator,
    command_args: []const []const u8,
    options: RunOptions,
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
        .stdin = switch (options.stdin) {
            .closed => .close,
            .stalled, .data => .pipe,
            .file => |file| .{ .file = file },
        },
        .stdout = switch (options.stdout) {
            .capture => .pipe,
            .closed => .close,
            .file => |file| .{ .file = file },
        },
        .stderr = switch (options.stderr) {
            .capture => .pipe,
            .file => |file| .{ .file = file },
        },
    });
    defer proc.kill(io);

    var process_deadline = try ProcessDeadline.init(io, &proc, options.timeout);
    defer process_deadline.deinit();

    var stdin_writer: StdinWriter = undefined;
    var stdin_group: std.Io.Group = .init;
    defer stdin_group.cancel(io);
    const has_stdin_writer = switch (options.stdin) {
        .closed, .stalled, .file => false,
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

    const output = collectOutput(allocator, io, &proc, process_deadline.timeout(), options) catch |err| switch (err) {
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

/// Reads the piped streams of the child; a stream that is not piped yields an empty slice.
fn collectOutput(
    allocator: std.mem.Allocator,
    io: std.Io,
    proc: *std.process.Child,
    timeout: std.Io.Timeout,
    options: RunOptions,
) !CommandOutput {
    const stdout_piped = options.stdout == .capture;
    const stderr_piped = options.stderr == .capture;
    if (stdout_piped and stderr_piped) {
        const outputs = try collectPipes(
            2,
            allocator,
            io,
            .{ proc.stdout.?, proc.stderr.? },
            timeout,
            options.output_limit,
        );
        return .{ .stdout = outputs[0], .stderr = outputs[1] };
    }
    if (stdout_piped) {
        const outputs = try collectPipes(1, allocator, io, .{proc.stdout.?}, timeout, options.output_limit);
        return .{ .stdout = outputs[0], .stderr = try allocator.alloc(u8, 0) };
    }
    if (stderr_piped) {
        const outputs = try collectPipes(1, allocator, io, .{proc.stderr.?}, timeout, options.output_limit);
        return .{ .stdout = try allocator.alloc(u8, 0), .stderr = outputs[0] };
    }
    return .{ .stdout = try allocator.alloc(u8, 0), .stderr = try allocator.alloc(u8, 0) };
}

fn collectPipes(
    comptime n: usize,
    allocator: std.mem.Allocator,
    io: std.Io,
    files: [n]std.Io.File,
    timeout: std.Io.Timeout,
    output_limit: usize,
) ![n][]u8 {
    var buffer: std.Io.File.MultiReader.Buffer(n) = undefined;
    var multi_reader: std.Io.File.MultiReader = undefined;
    multi_reader.init(allocator, io, buffer.toStreams(), &files);
    defer multi_reader.deinit();

    while (multi_reader.fill(4096, timeout)) |_| {
        for (0..n) |index| {
            if (multi_reader.reader(index).buffered().len > output_limit) {
                return error.StreamTooLong;
            }
        }
    } else |err| switch (err) {
        error.EndOfStream => {},
        else => |other| return other,
    }
    try multi_reader.checkAnyError();

    var outputs: [n][]u8 = undefined;
    var owned: usize = 0;
    errdefer for (outputs[0..owned]) |output| allocator.free(output);
    while (owned < n) : (owned += 1) {
        outputs[owned] = try multi_reader.toOwnedSlice(owned);
    }
    return outputs;
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
