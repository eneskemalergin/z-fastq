//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("root.zig");
const fastq = @import("fastq.zig");
const io_layer = @import("io.zig");

const INVALID_QUALITY_MESSAGE = "quality byte must be ASCII 33 through 126";

const USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count    Count records in plain or gzip FASTQ inputs
    \\  stats    Report aggregate FASTQ statistics
    \\  check    Validate FASTQ structure, sequence alphabet, and quality range
    \\
    \\General options:
    \\  -h, --help           Show this help message
    \\  -V, --version        Print version
    \\
    \\Input options:
    \\  --max-line-bytes N   Override default line length limit
    \\
    \\Machine output:
    \\  --json               Emit versioned JSON (stats and check only)
    \\
    \\Check options:
    \\  --alphabet POLICY    Select iupac (default) or acgtn sequence symbols
    \\
    \\Count usage:
    \\  z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
    \\Stats usage:
    \\  z-fastq stats [--json] [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
    \\Check usage:
    \\  z-fastq check [--json] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
;

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;

    var args = std.process.Args.Iterator.initAllocator(init.minimal.args, gpa) catch {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
        std.process.exit(3);
    };
    defer args.deinit();
    _ = args.skip();

    const cmd = args.next() orelse {
        printUsageAndExit(io);
    };

    if (std.mem.eql(u8, cmd, "--help") or std.mem.eql(u8, cmd, "-h")) {
        printHelpAndExit(io);
    }
    if (std.mem.eql(u8, cmd, "--version") or std.mem.eql(u8, cmd, "-V")) {
        printVersionAndExit(io);
    }

    const command: Command = if (std.mem.eql(u8, cmd, "count"))
        .count
    else if (std.mem.eql(u8, cmd, "stats"))
        .stats
    else if (std.mem.eql(u8, cmd, "check"))
        .check
    else {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown command: ") catch {};
        writeEscaped(.stderr(), io, cmd);
        std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
        printUsageAndExit(io);
    };

    var max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES;
    var alphabet: zfastq.Alphabet = .iupac;
    var json_output = false;
    var positional = std.ArrayList([]const u8).empty;
    defer positional.deinit(gpa);

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--")) {
            while (args.next()) |path| {
                positional.append(gpa, path) catch {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: out of memory\n",
                    ) catch {};
                    std.process.exit(3);
                };
            }
            break;
        }
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        }
        if (std.mem.eql(u8, arg, "--max-line-bytes")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --max-line-bytes requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            max_line_bytes = std.fmt.parseInt(usize, value, 10) catch |err| switch (err) {
                error.Overflow => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: --max-line-bytes exceeds supported limit\n",
                    ) catch {};
                    std.process.exit(4);
                },
                error.InvalidCharacter => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: invalid --max-line-bytes value\n",
                    ) catch {};
                    std.process.exit(2);
                },
            };
            continue;
        }
        if (command != .count and std.mem.eql(u8, arg, "--json")) {
            json_output = true;
            continue;
        }
        if (command == .check and std.mem.eql(u8, arg, "--alphabet")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --alphabet requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            alphabet = if (std.mem.eql(u8, value, "iupac"))
                .iupac
            else if (std.mem.eql(u8, value, "acgtn"))
                .acgtn
            else {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --alphabet must be iupac or acgtn\n",
                ) catch {};
                std.process.exit(2);
            };
            continue;
        }
        if (arg.len > 1 and std.mem.startsWith(u8, arg, "-")) {
            std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown ") catch {};
            std.Io.File.writeStreamingAll(.stderr(), io, @tagName(command)) catch {};
            std.Io.File.writeStreamingAll(.stderr(), io, " option: ") catch {};
            writeEscaped(.stderr(), io, arg);
            std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
            std.process.exit(2);
        }
        positional.append(gpa, arg) catch {
            std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
            std.process.exit(3);
        };
    }

    const options = InputOptions{ .max_line_bytes = max_line_bytes };
    const code = switch (command) {
        .count => runCount(io, positional.items, options),
        .stats => runStats(io, gpa, positional.items, options, json_output),
        .check => runCheck(io, positional.items, .{
            .max_line_bytes = max_line_bytes,
            .alphabet = alphabet,
        }, json_output),
    };
    std.process.exit(code);
}

const Command = enum { count, stats, check };

fn printUsageAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, USAGE) catch {};
    std.process.exit(2);
}

fn printHelpAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, USAGE) catch std.process.exit(3);
    std.process.exit(0);
}

fn printVersionAndExit(io: std.Io) noreturn {
    var buf: [64]u8 = undefined;
    const line = std.fmt.bufPrint(&buf, "z-fastq {s}\n", .{zfastq.VERSION}) catch "z-fastq\n";
    std.Io.File.writeStreamingAll(.stdout(), io, line) catch std.process.exit(3);
    std.process.exit(0);
}

// --- Count command ---

const InputOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
};

fn runCount(
    io: std.Io,
    inputs: []const []const u8,
    options: InputOptions,
) u8 {
    if (inputs.len == 0) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: count requires at least one input\n",
        ) catch {};
        return 2;
    }

    var has_stdin = false;
    for (inputs) |input| {
        if (!std.mem.eql(u8, input, "-")) continue;
        if (has_stdin) {
            std.Io.File.writeStreamingAll(
                .stderr(),
                io,
                "error: standard input may appear at most once\n",
            ) catch {};
            return 2;
        }
        has_stdin = true;
    }

    var exit_code: u8 = 0;
    for (inputs) |input| {
        const count_result = if (std.mem.eql(u8, input, "-"))
            countStdin(io, options)
        else
            countFile(io, input, options);
        const count = count_result catch |err| switch (err) {
            error.Io => {
                exit_code = @max(exit_code, 3);
                continue;
            },
            error.Format => {
                exit_code = @max(exit_code, 1);
                continue;
            },
            error.Limit => {
                printPathError(io, input, "line length limit exceeded");
                exit_code = @max(exit_code, 4);
                continue;
            },
            error.OutOfMemory => {
                std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                return 3;
            },
        };

        printCount(io, count) catch return @max(exit_code, 3);
    }

    return exit_code;
}

const CountError = error{
    Io,
    Format,
    Limit,
    OutOfMemory,
};

fn countFile(
    io: std.Io,
    path: []const u8,
    options: InputOptions,
) CountError!u64 {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => {
            printPathError(io, path, "file not found");
            return error.Io;
        },
        else => {
            printPathError(io, path, "failed to open file");
            return error.Io;
        },
    };
    defer file.close(io);

    return countInput(io, path, file, options);
}

fn countStdin(io: std.Io, options: InputOptions) CountError!u64 {
    return countInput(io, "-", std.Io.File.stdin(), options);
}

fn countInput(
    io: std.Io,
    label: []const u8,
    file: std.Io.File,
    options: InputOptions,
) CountError!u64 {
    var sniff_buffer: [2]u8 = undefined;
    var sniff_reader = file.readerStreaming(io, &sniff_buffer);
    const prefix = sniff_reader.interface.peek(2) catch |err| switch (err) {
        error.EndOfStream => null,
        error.ReadFailed => {
            printPathError(io, label, "I/O error");
            return error.Io;
        },
    };

    if (prefix) |bytes| {
        if (std.mem.eql(u8, bytes, &.{ 0x1f, 0x8b })) {
            var gzip_buffer: [64 * 1024]u8 = undefined;
            // The gzip reader must replay the magic consumed by the sniff reader.
            @memcpy(gzip_buffer[0..2], bytes);
            var gzip_reader = file.readerStreaming(io, &gzip_buffer);
            gzip_reader.pos = 2;
            gzip_reader.interface.end = 2;
            var gzip_source = io_layer.GzipSource.init(&gzip_reader.interface);
            return countGzipSource(io, label, &gzip_source, options);
        }
    }

    var plain_source = zfastq.io.plain.ReaderSource.init(&sniff_reader.interface);
    return countSource(io, label, plain_source.byteSource(), options);
}

fn countGzipSource(
    io: std.Io,
    label: []const u8,
    source: *io_layer.GzipSource,
    options: InputOptions,
) CountError!u64 {
    const scan_options = zfastq.count_scan.Options{
        .max_line_bytes = options.max_line_bytes,
    };
    var scanner = zfastq.count_scan.Scanner.init(scan_options);
    var decompressor_buffer: [io_layer.COUNT_DECOMPRESS_BUFFER_BYTES]u8 = undefined;
    while (io_layer.readGzipChunk(source, &decompressor_buffer) catch {
        printPathError(io, label, "I/O error");
        return error.Io;
    }) |decoded| {
        _ = scanner.feed(decoded) catch |err| {
            return mapScanError(io, label, &scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapScanError(io, label, &scanner, err);
    };
    return scanner.record_index;
}

fn countSource(
    io: std.Io,
    label: []const u8,
    source: zfastq.io.ByteSource,
    options: InputOptions,
) CountError!u64 {
    const scan_options = zfastq.count_scan.Options{
        .max_line_bytes = options.max_line_bytes,
    };
    var scanner = zfastq.count_scan.Scanner.init(scan_options);
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const n = source.read(&buf) catch {
            printPathError(io, label, "I/O error");
            return error.Io;
        };
        if (n == 0) break;
        _ = scanner.feed(buf[0..n]) catch |err| {
            return mapScanError(io, label, &scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapScanError(io, label, &scanner, err);
    };
    return scanner.record_index;
}

fn mapScanError(
    io: std.Io,
    path: []const u8,
    scanner: *zfastq.count_scan.Scanner,
    err: zfastq.ReaderError,
) CountError {
    return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => {
            if (scanner.takeLastError()) |details| {
                printParseError(io, path, details);
            }
            return error.Format;
        },
        error.LineTooLong => error.Limit,
        error.OutOfMemory => error.OutOfMemory,
        error.Io => error.Io,
    };
}

fn printCount(io: std.Io, n: u64) !void {
    var buf: [32]u8 = undefined;
    const text = try std.fmt.bufPrint(&buf, "{d}\n", .{n});
    try std.Io.File.writeStreamingAll(.stdout(), io, text);
}

const CommandFailure = struct {
    code: []const u8,
    message: []const u8,
    exit_code: u8,
    record_index: ?u64 = null,
    byte_offset: ?u64 = null,
    line_in_record: ?u3 = null,

    fn lint(details: zfastq.ParseError) CommandFailure {
        return .{
            .code = zfastq.codeTag(details.code),
            .message = details.message,
            .exit_code = 1,
            .record_index = details.record_index,
            .byte_offset = details.byte_offset,
            .line_in_record = details.line_in_record,
        };
    }

    fn plain(code: []const u8, message: []const u8, exit_code: u8) CommandFailure {
        return .{ .code = code, .message = message, .exit_code = exit_code };
    }
};

const StatsOutcome = union(enum) {
    success: zfastq.Stats,
    failure: CommandFailure,
};

fn validateRecordCommandInputs(
    io: std.Io,
    command: Command,
    inputs: []const []const u8,
) ?u8 {
    if (inputs.len == 0) {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
        std.Io.File.writeStreamingAll(.stderr(), io, @tagName(command)) catch {};
        std.Io.File.writeStreamingAll(.stderr(), io, " requires at least one input\n") catch {};
        return 2;
    }

    var has_stdin = false;
    for (inputs) |input| {
        if (!std.mem.eql(u8, input, "-")) continue;
        if (has_stdin) {
            std.Io.File.writeStreamingAll(
                .stderr(),
                io,
                "error: standard input may appear at most once\n",
            ) catch {};
            return 2;
        }
        has_stdin = true;
    }
    return null;
}

// --- Check command ---

const CheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
};

// Keep check code generation out of the measured count dispatch path.
noinline fn runCheck(
    io: std.Io,
    inputs: []const []const u8,
    options: CheckOptions,
    json_output: bool,
) u8 {
    if (validateRecordCommandInputs(io, .check, inputs)) |exit_code| return exit_code;
    if (json_output) return runCheckJson(io, inputs, options);

    var exit_code: u8 = 0;
    for (inputs) |input| {
        const failure = if (std.mem.eql(u8, input, "-"))
            checkStdin(io, options)
        else
            checkFile(io, input, options);
        if (failure) |details| {
            printCommandFailure(io, input, details);
            exit_code = @max(exit_code, details.exit_code);
        }
    }
    return exit_code;
}

fn runCheckJson(io: std.Io, inputs: []const []const u8, options: CheckOptions) u8 {
    var stdout_buffer: [16 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/check-v1") catch return 3;
    var exit_code: u8 = 0;
    for (inputs) |input| {
        const failure = if (std.mem.eql(u8, input, "-"))
            checkStdin(io, options)
        else
            checkFile(io, input, options);
        writeCheckJsonResult(&json, input, failure) catch return 3;
        if (failure) |details| exit_code = @max(exit_code, details.exit_code);
    }
    finishJsonDocument(&json) catch return 3;
    stdout_writer.interface.flush() catch return 3;
    return exit_code;
}

fn checkFile(
    io: std.Io,
    path: []const u8,
    options: CheckOptions,
) ?CommandFailure {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return CommandFailure.plain("io_error", "file not found", 3),
        else => return CommandFailure.plain("io_error", "failed to open file", 3),
    };

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);
    return checkSource(input.byteSource(), options);
}

fn checkStdin(
    io: std.Io,
    options: CheckOptions,
) ?CommandFailure {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);
    return checkSource(input.byteSource(), options);
}

fn checkSource(
    source: zfastq.io.ByteSource,
    options: CheckOptions,
) ?CommandFailure {
    var scanner = fastq.CheckScanner.init(
        .{ .max_line_bytes = options.max_line_bytes },
        .{ .alphabet = options.alphabet },
    );
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const n = source.read(&buf) catch {
            return CommandFailure.plain("io_error", "I/O error", 3);
        };
        if (n == 0) break;
        _ = scanner.feed(buf[0..n]) catch |err| {
            return mapCheckScannerFailure(&scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapCheckScannerFailure(&scanner, err);
    };
    return null;
}

fn mapCheckScannerFailure(
    scanner: *fastq.CheckScanner,
    err: fastq.CheckScannerError,
) CommandFailure {
    return switch (err) {
        error.Format => if (scanner.takeLastError()) |details|
            CommandFailure.lint(details)
        else
            CommandFailure.plain("io_error", "validation failed without details", 3),
        error.LineTooLong => CommandFailure.plain(
            "line_limit",
            "line length limit exceeded",
            4,
        ),
        error.ArithmeticLimit => CommandFailure.plain(
            "arithmetic_limit",
            "input location exceeds supported limit",
            4,
        ),
    };
}

// --- Stats command ---

// Keep stats code generation out of the measured count dispatch path.
noinline fn runStats(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: InputOptions,
    json_output: bool,
) u8 {
    if (validateRecordCommandInputs(io, .stats, inputs)) |exit_code| return exit_code;
    if (json_output) return runStatsJson(io, allocator, inputs, options);

    var exit_code: u8 = 0;
    var printed_block = false;
    for (inputs) |input| {
        const outcome = if (std.mem.eql(u8, input, "-"))
            statsStdin(io, allocator, options)
        else
            statsFile(io, allocator, input, options);
        const stats = switch (outcome) {
            .success => |stats| stats,
            .failure => |failure| {
                if (std.mem.eql(u8, failure.code, "out_of_memory")) {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: out of memory\n",
                    ) catch {};
                    return 3;
                }
                printCommandFailure(io, input, failure);
                exit_code = @max(exit_code, failure.exit_code);
                continue;
            },
        };

        printStats(io, input, stats.result(), printed_block) catch
            return @max(exit_code, 3);
        printed_block = true;
    }
    return exit_code;
}

fn runStatsJson(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: InputOptions,
) u8 {
    var stdout_buffer: [16 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/stats-v1") catch return 3;
    var exit_code: u8 = 0;
    for (inputs) |input| {
        const outcome = if (std.mem.eql(u8, input, "-"))
            statsStdin(io, allocator, options)
        else
            statsFile(io, allocator, input, options);
        writeStatsJsonResult(&json, input, outcome) catch return 3;
        switch (outcome) {
            .success => {},
            .failure => |failure| exit_code = @max(exit_code, failure.exit_code),
        }
    }
    finishJsonDocument(&json) catch return 3;
    stdout_writer.interface.flush() catch return 3;
    return exit_code;
}

fn statsFile(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    options: InputOptions,
) StatsOutcome {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .{ .failure = CommandFailure.plain(
            "io_error",
            "file not found",
            3,
        ) },
        else => return .{ .failure = CommandFailure.plain(
            "io_error",
            "failed to open file",
            3,
        ) },
    };

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return .{ .failure = CommandFailure.plain("io_error", "I/O error", 3) };
    };
    defer input.deinit(io);
    return collectStats(allocator, input.byteSource(), options);
}

fn statsStdin(
    io: std.Io,
    allocator: std.mem.Allocator,
    options: InputOptions,
) StatsOutcome {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        return .{ .failure = CommandFailure.plain("io_error", "I/O error", 3) };
    };
    defer input.deinit(io);
    return collectStats(allocator, input.byteSource(), options);
}

const RecordInput = struct {
    file: std.Io.File,
    owns_file: bool,
    read_buffer: [64 * 1024]u8,
    file_reader: std.Io.File.Reader,
    source: union(enum) {
        plain: io_layer.ReaderSource,
        gzip: io_layer.GzipSource,
    },

    fn init(
        self: *RecordInput,
        io: std.Io,
        file: std.Io.File,
        owns_file: bool,
    ) error{Io}!void {
        self.file = file;
        self.owns_file = owns_file;
        self.file_reader = file.readerStreaming(io, &self.read_buffer);
        const prefix = self.file_reader.interface.peek(2) catch |err| switch (err) {
            error.EndOfStream => null,
            error.ReadFailed => return error.Io,
        };
        if (prefix) |bytes| {
            if (std.mem.eql(u8, bytes, &.{ 0x1f, 0x8b })) {
                self.source = .{ .gzip = io_layer.GzipSource.init(&self.file_reader.interface) };
                return;
            }
        }
        self.source = .{ .plain = io_layer.ReaderSource.init(&self.file_reader.interface) };
    }

    fn deinit(self: *RecordInput, io: std.Io) void {
        if (self.owns_file) self.file.close(io);
        self.* = undefined;
    }

    fn byteSource(self: *RecordInput) zfastq.io.ByteSource {
        return switch (self.source) {
            .plain => |*source| source.byteSource(),
            .gzip => |*source| source.byteSource(),
        };
    }
};

fn collectStats(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    options: InputOptions,
) StatsOutcome {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return .{ .failure = CommandFailure.plain(
        "out_of_memory",
        "out of memory",
        3,
    ) };
    defer reader.deinit();

    var stats: zfastq.Stats = .{};
    while (reader.next() catch |err| {
        return .{ .failure = mapReaderFailure(&reader, err) };
    }) |record| {
        stats.addRecord(record) catch |err| switch (err) {
            error.S006InvalidQuality => {
                const quality_error = stats.takeLastQualityError() orelse {
                    return .{ .failure = CommandFailure.plain(
                        "io_error",
                        "quality validation failed without details",
                        3,
                    ) };
                };
                const offsets = reader.currentRecordOffsets() orelse {
                    return .{ .failure = CommandFailure.plain(
                        "io_error",
                        "record location is unavailable",
                        3,
                    ) };
                };
                const relative_offset = std.math.cast(u64, quality_error.byte_index) orelse
                    return .{ .failure = CommandFailure.plain(
                        "arithmetic_limit",
                        "statistics arithmetic limit exceeded",
                        4,
                    ) };
                const byte_offset = std.math.add(
                    u64,
                    offsets.quality,
                    relative_offset,
                ) catch return .{ .failure = CommandFailure.plain(
                    "arithmetic_limit",
                    "statistics arithmetic limit exceeded",
                    4,
                ) };
                return .{ .failure = CommandFailure.lint(.{
                    .code = .s006_invalid_quality_range,
                    .message = INVALID_QUALITY_MESSAGE,
                    .record_index = reader.recordIndex() - 1,
                    .byte_offset = byte_offset,
                    .line_in_record = 4,
                }) };
            },
            error.S005LengthMismatch => {
                const offsets = reader.currentRecordOffsets() orelse {
                    return .{ .failure = CommandFailure.plain(
                        "io_error",
                        "record location is unavailable",
                        3,
                    ) };
                };
                return .{ .failure = CommandFailure.lint(.{
                    .code = .s005_length_mismatch,
                    .message = "sequence and quality lengths differ",
                    .record_index = reader.recordIndex() - 1,
                    .byte_offset = offsets.quality,
                    .line_in_record = 4,
                }) };
            },
            error.Overflow => return .{ .failure = CommandFailure.plain(
                "arithmetic_limit",
                "statistics arithmetic limit exceeded",
                4,
            ) },
        };
    }
    return .{ .success = stats };
}

fn mapReaderFailure(
    reader: *zfastq.Reader,
    err: zfastq.ReaderError,
) CommandFailure {
    return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => if (reader.takeLastError()) |details|
            CommandFailure.lint(details)
        else
            CommandFailure.plain("io_error", "validation failed without details", 3),
        error.LineTooLong => CommandFailure.plain(
            "line_limit",
            "line length limit exceeded",
            4,
        ),
        error.OutOfMemory => CommandFailure.plain("out_of_memory", "out of memory", 3),
        error.Io => CommandFailure.plain("io_error", "I/O error", 3),
    };
}

// --- Machine output ---

fn beginJsonDocument(json: *std.json.Stringify, schema: []const u8) !void {
    try json.beginObject();
    try json.objectField("schema");
    try json.write(schema);
    try json.objectField("tool");
    try json.beginObject();
    try json.objectField("name");
    try json.write("z-fastq");
    try json.objectField("version");
    try json.write(zfastq.VERSION);
    try json.endObject();
    try json.objectField("byte_strings");
    try json.write("escaped-bytes-v1");
    try json.objectField("results");
    try json.beginArray();
}

fn finishJsonDocument(json: *std.json.Stringify) !void {
    try json.endArray();
    try json.endObject();
    try json.writer.writeByte('\n');
}

fn writeCheckJsonResult(
    json: *std.json.Stringify,
    input: []const u8,
    failure: ?CommandFailure,
) !void {
    try json.beginObject();
    try json.objectField("input");
    try writeEscapedJsonString(json, input);
    try json.objectField("status");
    try json.write(if (failure == null) "ok" else "error");
    if (failure) |details| try writeJsonFailure(json, details);
    try json.endObject();
}

fn writeStatsJsonResult(
    json: *std.json.Stringify,
    input: []const u8,
    outcome: StatsOutcome,
) !void {
    try json.beginObject();
    try json.objectField("input");
    try writeEscapedJsonString(json, input);
    switch (outcome) {
        .failure => |failure| {
            try json.objectField("status");
            try json.write("error");
            try writeJsonFailure(json, failure);
        },
        .success => |stats| {
            const result = stats.result();
            try json.objectField("status");
            try json.write("ok");
            try json.objectField("reads");
            try json.write(result.reads);
            try json.objectField("bases");
            try json.write(result.bases);
            try json.objectField("min_length");
            try json.write(result.min_length);
            try json.objectField("max_length");
            try json.write(result.max_length);
            try json.objectField("mean_length");
            try json.write(result.mean_length);
            try json.objectField("a");
            try json.write(result.a);
            try json.objectField("c");
            try json.write(result.c);
            try json.objectField("g");
            try json.write(result.g);
            try json.objectField("t");
            try json.write(result.t);
            try json.objectField("n");
            try json.write(result.n);
            try json.objectField("other_bases");
            try json.write(result.other_bases);
            try json.objectField("gc_fraction");
            try json.write(result.gc_fraction);
            try json.objectField("quality_sum");
            try json.write(result.quality_sum);
            try json.objectField("mean_quality");
            try json.write(result.mean_quality);
            try json.objectField("q20_bases");
            try json.write(result.q20_bases);
            try json.objectField("q20_fraction");
            try json.write(result.q20_fraction);
            try json.objectField("q30_bases");
            try json.write(result.q30_bases);
            try json.objectField("q30_fraction");
            try json.write(result.q30_fraction);
        },
    }
    try json.endObject();
}

fn writeJsonFailure(json: *std.json.Stringify, failure: CommandFailure) !void {
    try json.objectField("error");
    try json.beginObject();
    try json.objectField("code");
    try json.write(failure.code);
    try json.objectField("message");
    try json.write(failure.message);
    try json.objectField("record_index");
    try json.write(failure.record_index);
    try json.objectField("byte_offset");
    try json.write(failure.byte_offset);
    try json.objectField("line_in_record");
    try json.write(failure.line_in_record);
    try json.endObject();
}

fn writeEscapedJsonString(json: *std.json.Stringify, bytes: []const u8) !void {
    try json.beginWriteRaw();
    try json.writer.writeByte('"');
    var run_start: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte >= 0x20 and byte <= 0x7e and byte != '\\') continue;

        try std.json.Stringify.encodeJsonStringChars(
            bytes[run_start..index],
            .{},
            json.writer,
        );
        if (byte == '\\') {
            try std.json.Stringify.encodeJsonStringChars("\\\\", .{}, json.writer);
        } else {
            const escaped = [4]u8{ '\\', 'x', HEX[byte >> 4], HEX[byte & 0x0f] };
            try std.json.Stringify.encodeJsonStringChars(&escaped, .{}, json.writer);
        }
        run_start = index + 1;
    }
    try std.json.Stringify.encodeJsonStringChars(bytes[run_start..], .{}, json.writer);
    try json.writer.writeByte('"');
    json.endWriteRaw();
}

fn printStats(
    io: std.Io,
    label: []const u8,
    result: zfastq.StatsResult,
    separator: bool,
) !void {
    if (separator) try std.Io.File.writeStreamingAll(.stdout(), io, "\n");
    try std.Io.File.writeStreamingAll(.stdout(), io, "input: ");
    try writeEscapedAll(.stdout(), io, label);
    try std.Io.File.writeStreamingAll(.stdout(), io, "\n");
    try writeUnsignedField(io, "reads", result.reads);
    try writeUnsignedField(io, "bases", result.bases);
    try writeOptionalUnsignedField(io, "min_length", result.min_length);
    try writeOptionalUnsignedField(io, "max_length", result.max_length);
    try writeRatioField(io, "mean_length", result.bases, result.reads);
    try writeUnsignedField(io, "a", result.a);
    try writeUnsignedField(io, "c", result.c);
    try writeUnsignedField(io, "g", result.g);
    try writeUnsignedField(io, "t", result.t);
    try writeUnsignedField(io, "n", result.n);
    try writeUnsignedField(io, "other_bases", result.other_bases);
    try writeRatioField(
        io,
        "gc_fraction",
        @as(u128, result.g) + result.c,
        @as(u128, result.a) + result.c + result.g + result.t,
    );
    try writeUnsignedField(io, "quality_sum", result.quality_sum);
    try writeRatioField(io, "mean_quality", result.quality_sum, result.bases);
    try writeUnsignedField(io, "q20_bases", result.q20_bases);
    try writeRatioField(io, "q20_fraction", result.q20_bases, result.bases);
    try writeUnsignedField(io, "q30_bases", result.q30_bases);
    try writeRatioField(io, "q30_fraction", result.q30_bases, result.bases);
}

fn writeUnsignedField(io: std.Io, name: []const u8, value: u64) !void {
    var buf: [96]u8 = undefined;
    const line = try std.fmt.bufPrint(&buf, "{s}: {d}\n", .{ name, value });
    try std.Io.File.writeStreamingAll(.stdout(), io, line);
}

fn writeOptionalUnsignedField(io: std.Io, name: []const u8, value: ?u64) !void {
    if (value) |number| return writeUnsignedField(io, name, number);
    var buf: [64]u8 = undefined;
    const line = try std.fmt.bufPrint(&buf, "{s}: -\n", .{name});
    try std.Io.File.writeStreamingAll(.stdout(), io, line);
}

fn writeRatioField(
    io: std.Io,
    name: []const u8,
    numerator: u128,
    denominator: u128,
) !void {
    var buf: [96]u8 = undefined;
    if (denominator == 0) {
        const line = try std.fmt.bufPrint(&buf, "{s}: -\n", .{name});
        return std.Io.File.writeStreamingAll(.stdout(), io, line);
    }
    const scale = 1_000_000;
    const rounded = (numerator * scale + denominator / 2) / denominator;
    const line = try std.fmt.bufPrint(
        &buf,
        "{s}: {d}.{d:0>6}\n",
        .{ name, rounded / scale, rounded % scale },
    );
    try std.Io.File.writeStreamingAll(.stdout(), io, line);
}

fn printPathError(io: std.Io, path: []const u8, message: []const u8) void {
    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
}

fn printCommandFailure(io: std.Io, path: []const u8, failure: CommandFailure) void {
    if (failure.record_index == null or
        failure.byte_offset == null or
        failure.line_in_record == null)
    {
        return printPathError(io, path, failure.message);
    }

    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, failure.code) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, failure.message) catch {};

    var buf: [96]u8 = undefined;
    const suffix = std.fmt.bufPrint(
        &buf,
        " (record {d}, line {d}, offset {d})\n",
        .{
            failure.record_index.?,
            failure.line_in_record.?,
            failure.byte_offset.?,
        },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, suffix) catch {};
}

fn printParseError(io: std.Io, path: []const u8, details: zfastq.ParseError) void {
    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, zfastq.codeTag(details.code)) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, details.message) catch {};

    var buf: [96]u8 = undefined;
    const suffix = std.fmt.bufPrint(
        &buf,
        " (record {d}, line {d}, offset {d})\n",
        .{
            details.record_index,
            details.line_in_record,
            details.byte_offset,
        },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, suffix) catch {};
}

const HEX = "0123456789ABCDEF";

fn writeEscaped(file: std.Io.File, io: std.Io, bytes: []const u8) void {
    writeEscapedAll(file, io, bytes) catch {};
}

fn writeEscapedAll(file: std.Io.File, io: std.Io, bytes: []const u8) !void {
    var run_start: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte >= 0x20 and byte <= 0x7e and byte != '\\') continue;

        try std.Io.File.writeStreamingAll(file, io, bytes[run_start..index]);
        if (byte == '\\') {
            try std.Io.File.writeStreamingAll(file, io, "\\\\");
        } else {
            const escaped = [4]u8{ '\\', 'x', HEX[byte >> 4], HEX[byte & 0x0f] };
            try std.Io.File.writeStreamingAll(file, io, &escaped);
        }
        run_start = index + 1;
    }
    try std.Io.File.writeStreamingAll(file, io, bytes[run_start..]);
}

test "[unit] - [machine output]: handled non-lint errors share one exact shape" {
    const cases = [_]struct {
        code: []const u8,
        message: []const u8,
        exit_code: u8,
    }{
        .{ .code = "io_error", .message = "I/O error", .exit_code = 3 },
        .{ .code = "line_limit", .message = "line length limit exceeded", .exit_code = 4 },
        .{ .code = "arithmetic_limit", .message = "arithmetic limit exceeded", .exit_code = 4 },
        .{ .code = "out_of_memory", .message = "out of memory", .exit_code = 3 },
    };

    for (cases) |case| {
        var storage: [256]u8 = undefined;
        var writer = std.Io.Writer.fixed(&storage);
        var json: std.json.Stringify = .{ .writer = &writer };
        try json.beginObject();
        try writeJsonFailure(
            &json,
            CommandFailure.plain(case.code, case.message, case.exit_code),
        );
        try json.endObject();

        var parsed = try std.json.parseFromSlice(
            std.json.Value,
            std.testing.allocator,
            writer.buffered(),
            .{},
        );
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqual(@as(usize, 1), object.count());
        const failure = object.get("error").?.object;
        try std.testing.expectEqual(@as(usize, 5), failure.count());
        try std.testing.expectEqualStrings(case.code, failure.get("code").?.string);
        try std.testing.expectEqualStrings(case.message, failure.get("message").?.string);
        try std.testing.expect(failure.get("record_index").? == .null);
        try std.testing.expect(failure.get("byte_offset").? == .null);
        try std.testing.expect(failure.get("line_in_record").? == .null);
    }
}

test "[edge] - [stats-json]: preserves the maximum u64 counter" {
    var storage: [1024]u8 = undefined;
    var writer = std.Io.Writer.fixed(&storage);
    var json: std.json.Stringify = .{ .writer = &writer };
    try writeStatsJsonResult(&json, "-", .{ .success = .{
        .reads = std.math.maxInt(u64),
    } });

    var parsed = try std.json.parseFromSlice(
        std.json.Value,
        std.testing.allocator,
        writer.buffered(),
        .{ .parse_numbers = false },
    );
    defer parsed.deinit();
    const result = parsed.value.object.get("reads").?;
    try std.testing.expectEqualStrings("18446744073709551615", result.number_string);
}

test "[failure] - [stats command]: reader allocation failure becomes a handled result" {
    var source = zfastq.io.plain.SliceSource.init("");
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 0,
    });

    const outcome = collectStats(failing.allocator(), source.byteSource(), .{});
    const failure = switch (outcome) {
        .success => return error.ExpectedFailure,
        .failure => |failure| failure,
    };
    try std.testing.expect(failing.has_induced_failure);
    try std.testing.expectEqualStrings("out_of_memory", failure.code);
    try std.testing.expectEqualStrings("out of memory", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expect(failure.record_index == null);
    try std.testing.expect(failure.byte_offset == null);
    try std.testing.expect(failure.line_in_record == null);
}
