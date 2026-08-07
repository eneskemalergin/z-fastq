//! `z-fastq count`: count records in plain FASTQ files.

const std = @import("std");
const zfastq = @import("z-fastq");

pub const Options = struct {
    max_line_bytes: usize = zfastq.limits.default_max_line_bytes,
};

pub fn run(
    io: std.Io,
    paths: []const []const u8,
    options: Options,
) u8 {
    if (paths.len == 0) {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: count requires at least one file path\n") catch {};
        return 2;
    }

    var exit_code: u8 = 0;
    for (paths) |path| {
        const count_result = countFile(io, path, options) catch |err| switch (err) {
            error.Io => {
                printPathError(io, path, "I/O error");
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

// --- Private helpers ---

const CountError = error{
    Io,
    Format,
    Limit,
    OutOfMemory,
};

fn countFile(
    io: std.Io,
    path: []const u8,
    options: Options,
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

    var scanner: zfastq.count_scan.Scanner = undefined;
    const scan_options = zfastq.count_scan.Options{
        .max_line_bytes = options.max_line_bytes,
    };
    const count = zfastq.count_scan.countPlainFile(io, file, scan_options, &scanner) catch |err| {
        return mapScanError(io, path, &scanner, err);
    };
    return count;
}

fn mapScanError(
    io: std.Io,
    path: []const u8,
    scanner: *zfastq.count_scan.Scanner,
    err: zfastq.ReaderError,
) CountError {
    return switch (err) {
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
    var buf: [512]u8 = undefined;
    const line = std.fmt.bufPrint(&buf, "error: {s}: {s}\n", .{ path, message }) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, line) catch {};
}

fn printParseError(io: std.Io, path: []const u8, details: zfastq.ParseError) void {
    var buf: [512]u8 = undefined;
    const line = std.fmt.bufPrint(
        &buf,
        "error: {s}: {s}: {s} (record {d}, line {d}, offset {d})\n",
        .{
            path,
            zfastq.codeTag(details.code),
            details.message,
            details.record_index,
            details.line_in_record,
            details.byte_offset,
        },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, line) catch {};
}
