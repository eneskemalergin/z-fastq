//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("z-fastq");

const USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count    Count records in plain FASTQ files
    \\
    \\General options:
    \\  -h, --help           Show this help message
    \\  -V, --version        Print version
    \\
    \\Count options:
    \\  --max-line-bytes N   Override default line length limit
    \\
    \\Count usage:
    \\  z-fastq count <file.fastq> [file2.fastq ...]
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
            if (std.mem.startsWith(u8, arg, "-")) {
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
    paths: []const []const u8,
    options: CountOptions,
) u8 {
    if (paths.len == 0) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: count requires at least one file path\n",
        ) catch {};
        return 2;
    }

    var exit_code: u8 = 0;
    for (paths) |path| {
        const count_result = countFile(io, path, options) catch |err| switch (err) {
            error.Io => {
                exit_code = @max(exit_code, 3);
                continue;
            },
            error.Format => {
                exit_code = @max(exit_code, 1);
                continue;
            },
            error.Limit => {
                printPathError(io, path, "line length limit exceeded");
                exit_code = @max(exit_code, 4);
                continue;
            },
            error.OutOfMemory => {
                std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                return 3;
            },
        };

        printCount(io, count_result) catch return 3;
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

    const scan_options = zfastq.count_scan.Options{
        .max_line_bytes = options.max_line_bytes,
    };
    var scanner = zfastq.count_scan.Scanner.init(scan_options);
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const n = file.readStreaming(io, &.{&buf}) catch |err| switch (err) {
            error.EndOfStream => 0,
            else => {
                printPathError(io, path, "I/O error");
                return error.Io;
            },
        };
        if (n == 0) break;
        _ = scanner.feed(buf[0..n]) catch |err| {
            return mapScanError(io, path, &scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapScanError(io, path, &scanner, err);
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
