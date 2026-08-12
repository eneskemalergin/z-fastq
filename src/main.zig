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
    \\Check options:
    \\  --alphabet POLICY    Select iupac (default) or acgtn sequence symbols
    \\
    \\Count usage:
    \\  z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
    \\Stats usage:
    \\  z-fastq stats [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
    \\Check usage:
    \\  z-fastq check [--alphabet iupac|acgtn] [--max-line-bytes N] <path|-> [<path|-> ...]
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
        .stats => runStats(io, gpa, positional.items, options),
        .check => runCheck(io, positional.items, .{
            .max_line_bytes = max_line_bytes,
            .alphabet = alphabet,
        }),
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

const RecordCommandError = error{
    Io,
    Format,
    LineLimit,
    ArithmeticLimit,
    OutOfMemory,
};

const CheckCommandError = error{
    Io,
    Format,
    LineLimit,
    ArithmeticLimit,
};

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
) u8 {
    if (inputs.len == 0) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: check requires at least one input\n",
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
        const result = if (std.mem.eql(u8, input, "-"))
            checkStdin(io, options)
        else
            checkFile(io, input, options);
        result catch |err| switch (err) {
            error.Io => exit_code = @max(exit_code, 3),
            error.Format => exit_code = @max(exit_code, 1),
            error.LineLimit => {
                printPathError(io, input, "line length limit exceeded");
                exit_code = @max(exit_code, 4);
            },
            error.ArithmeticLimit => {
                printPathError(io, input, "input location exceeds supported limit");
                exit_code = @max(exit_code, 4);
            },
        };
    }
    return exit_code;
}

fn checkFile(
    io: std.Io,
    path: []const u8,
    options: CheckOptions,
) CheckCommandError!void {
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

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        printPathError(io, path, "I/O error");
        return error.Io;
    };
    defer input.deinit(io);
    return checkSource(io, path, input.byteSource(), options);
}

fn checkStdin(
    io: std.Io,
    options: CheckOptions,
) CheckCommandError!void {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        printPathError(io, "-", "I/O error");
        return error.Io;
    };
    defer input.deinit(io);
    return checkSource(io, "-", input.byteSource(), options);
}

fn checkSource(
    io: std.Io,
    label: []const u8,
    source: zfastq.io.ByteSource,
    options: CheckOptions,
) CheckCommandError!void {
    var scanner = fastq.CheckScanner.init(
        .{ .max_line_bytes = options.max_line_bytes },
        .{ .alphabet = options.alphabet },
    );
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const n = source.read(&buf) catch {
            printPathError(io, label, "I/O error");
            return error.Io;
        };
        if (n == 0) break;
        _ = scanner.feed(buf[0..n]) catch |err| {
            return mapCheckScannerError(io, label, &scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapCheckScannerError(io, label, &scanner, err);
    };
}

fn mapCheckScannerError(
    io: std.Io,
    path: []const u8,
    scanner: *fastq.CheckScanner,
    err: fastq.CheckScannerError,
) CheckCommandError {
    return switch (err) {
        error.Format => {
            if (scanner.takeLastError()) |details| printParseError(io, path, details);
            return error.Format;
        },
        error.LineTooLong => error.LineLimit,
        error.ArithmeticLimit => error.ArithmeticLimit,
    };
}

// --- Stats command ---

// Keep stats code generation out of the measured count dispatch path.
noinline fn runStats(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: InputOptions,
) u8 {
    if (inputs.len == 0) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: stats requires at least one input\n",
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
    var printed_block = false;
    for (inputs) |input| {
        const stats_result = if (std.mem.eql(u8, input, "-"))
            statsStdin(io, allocator, options)
        else
            statsFile(io, allocator, input, options);
        const stats = stats_result catch |err| switch (err) {
            error.Io => {
                exit_code = @max(exit_code, 3);
                continue;
            },
            error.Format => {
                exit_code = @max(exit_code, 1);
                continue;
            },
            error.LineLimit => {
                printPathError(io, input, "line length limit exceeded");
                exit_code = @max(exit_code, 4);
                continue;
            },
            error.ArithmeticLimit => {
                printPathError(io, input, "statistics arithmetic limit exceeded");
                exit_code = @max(exit_code, 4);
                continue;
            },
            error.OutOfMemory => {
                std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                return 3;
            },
        };

        printStats(io, input, stats.result(), printed_block) catch
            return @max(exit_code, 3);
        printed_block = true;
    }
    return exit_code;
}

fn statsFile(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    options: InputOptions,
) RecordCommandError!zfastq.Stats {
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

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        printPathError(io, path, "I/O error");
        return error.Io;
    };
    defer input.deinit(io);
    return collectStats(io, allocator, path, input.byteSource(), options);
}

fn statsStdin(
    io: std.Io,
    allocator: std.mem.Allocator,
    options: InputOptions,
) RecordCommandError!zfastq.Stats {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        printPathError(io, "-", "I/O error");
        return error.Io;
    };
    defer input.deinit(io);
    return collectStats(io, allocator, "-", input.byteSource(), options);
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
    io: std.Io,
    allocator: std.mem.Allocator,
    label: []const u8,
    source: zfastq.io.ByteSource,
    options: InputOptions,
) RecordCommandError!zfastq.Stats {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return error.OutOfMemory;
    defer reader.deinit();

    var stats: zfastq.Stats = .{};
    while (reader.next() catch |err| {
        return mapReaderError(io, label, &reader, err);
    }) |record| {
        stats.addRecord(record) catch |err| switch (err) {
            error.S006InvalidQuality => {
                const quality_error = stats.takeLastQualityError() orelse {
                    printPathError(io, label, "quality validation failed without details");
                    return error.Io;
                };
                const offsets = reader.currentRecordOffsets() orelse {
                    printPathError(io, label, "record location is unavailable");
                    return error.Io;
                };
                const relative_offset = std.math.cast(u64, quality_error.byte_index) orelse
                    return error.ArithmeticLimit;
                const byte_offset = std.math.add(
                    u64,
                    offsets.quality,
                    relative_offset,
                ) catch return error.ArithmeticLimit;
                printParseError(io, label, .{
                    .code = .s006_invalid_quality_range,
                    .message = INVALID_QUALITY_MESSAGE,
                    .record_index = reader.recordIndex() - 1,
                    .byte_offset = byte_offset,
                    .line_in_record = 4,
                });
                return error.Format;
            },
            error.S005LengthMismatch => {
                const offsets = reader.currentRecordOffsets() orelse {
                    printPathError(io, label, "record location is unavailable");
                    return error.Io;
                };
                printParseError(io, label, .{
                    .code = .s005_length_mismatch,
                    .message = "sequence and quality lengths differ",
                    .record_index = reader.recordIndex() - 1,
                    .byte_offset = offsets.quality,
                    .line_in_record = 4,
                });
                return error.Format;
            },
            error.Overflow => return error.ArithmeticLimit,
        };
    }
    return stats;
}

fn mapReaderError(
    io: std.Io,
    path: []const u8,
    reader: *zfastq.Reader,
    err: zfastq.ReaderError,
) RecordCommandError {
    return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => {
            if (reader.takeLastError()) |details| {
                printParseError(io, path, details);
            }
            return error.Format;
        },
        error.LineTooLong => error.LineLimit,
        error.OutOfMemory => error.OutOfMemory,
        error.Io => {
            printPathError(io, path, "I/O error");
            return error.Io;
        },
    };
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
