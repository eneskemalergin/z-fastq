//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("z-fastq");

const USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count    Count records in plain or gzip FASTQ inputs
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

    var max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES;
    var positional = std.ArrayList([]const u8).empty;
    defer positional.deinit(gpa);

    if (std.mem.eql(u8, cmd, "count")) {
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
            if (arg.len > 1 and std.mem.startsWith(u8, arg, "-")) {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: unknown count option: ",
                ) catch {};
                writeEscaped(.stderr(), io, arg);
                std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
                std.process.exit(2);
            }
            positional.append(gpa, arg) catch {
                std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                std.process.exit(3);
            };
        }

        const code = runCount(io, positional.items, .{
            .max_line_bytes = max_line_bytes,
        });
        std.process.exit(code);
    }

    std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown command: ") catch {};
    writeEscaped(.stderr(), io, cmd);
    std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
    printUsageAndExit(io);
}

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

const CountOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
};

fn runCount(
    io: std.Io,
    inputs: []const []const u8,
    options: CountOptions,
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
    options: CountOptions,
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

fn countStdin(io: std.Io, options: CountOptions) CountError!u64 {
    return countInput(io, "-", std.Io.File.stdin(), options);
}

fn countInput(
    io: std.Io,
    label: []const u8,
    file: std.Io.File,
    options: CountOptions,
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
            var gzip_source = zfastq.io.gzip.ReaderSource.init(&gzip_reader.interface);
            return countSource(io, label, gzip_source.byteSource(), options);
        }
    }

    var plain_source = zfastq.io.plain.ReaderSource.init(&sniff_reader.interface);
    return countSource(io, label, plain_source.byteSource(), options);
}

fn countSource(
    io: std.Io,
    label: []const u8,
    source: zfastq.io.ByteSource,
    options: CountOptions,
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
    var run_start: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte >= 0x20 and byte <= 0x7e and byte != '\\') continue;

        std.Io.File.writeStreamingAll(file, io, bytes[run_start..index]) catch return;
        if (byte == '\\') {
            std.Io.File.writeStreamingAll(file, io, "\\\\") catch return;
        } else {
            const escaped = [4]u8{ '\\', 'x', HEX[byte >> 4], HEX[byte & 0x0f] };
            std.Io.File.writeStreamingAll(file, io, &escaped) catch return;
        }
        run_start = index + 1;
    }
    std.Io.File.writeStreamingAll(file, io, bytes[run_start..]) catch {};
}
