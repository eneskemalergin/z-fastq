//! Public aggregate-statistics contracts.

const std = @import("std");
const builtin = @import("builtin");
const zfastq = @import("z-fastq");
const cli = @import("utilities.zig");
const CommandResult = cli.CommandResult;
const runCli = cli.runWithStdin;
const runCliWithClosedStdout = cli.runWithClosedStdout;
const runCliWithClosedStdin = cli.runWithClosedStdin;

const SAMPLE_FASTQ = "@r\nACGTN?\n+\n!5?I~!\n";
const SAMPLE_GZIP = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0xff,
    0x01, 0x13, 0x00, 0xec, 0xff, 0x40, 0x72, 0x0a, 0x41, 0x43,
    0x47, 0x54, 0x4e, 0x3f, 0x0a, 0x2b, 0x0a, 0x21, 0x35, 0x3f,
    0x49, 0x7e, 0x21, 0x0a, 0x56, 0x1d, 0xf4, 0xf1, 0x13, 0x00,
    0x00, 0x00,
};
const INVALID_QUALITY_GZIP = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0xff,
    0x01, 0x11, 0x00, 0xee, 0xff, 0x40, 0x62, 0x61, 0x64, 0x0d,
    0x0a, 0x41, 0x43, 0x0d, 0x0a, 0x2b, 0x0d, 0x0a, 0x21, 0x7f,
    0x0d, 0x0a, 0xe4, 0xbd, 0xef, 0xc2, 0x11, 0x00, 0x00, 0x00,
};
const SAMPLE_FIELDS =
    \\reads: 1
    \\bases: 6
    \\min_length: 6
    \\max_length: 6
    \\mean_length: 6.000000
    \\a: 1
    \\c: 1
    \\g: 1
    \\t: 1
    \\n: 1
    \\other_bases: 1
    \\gc_fraction: 0.500000
    \\quality_sum: 183
    \\mean_quality: 30.500000
    \\q20_bases: 4
    \\q20_fraction: 0.666667
    \\q30_bases: 3
    \\q30_fraction: 0.500000
    \\
;

test "[unit] - [statistics]: records accumulate exact public results" {
    var stats: zfastq.Stats = .{};
    try stats.addRecord(record("AaCcGgTtNnR?", "!5?I~!!!!!!!"));
    try stats.addRecord(record("", ""));

    const result = stats.result();
    try std.testing.expectEqual(@as(u64, 2), result.reads);
    try std.testing.expectEqual(@as(u64, 12), result.bases);
    try std.testing.expectEqual(@as(u64, 0), result.min_length.?);
    try std.testing.expectEqual(@as(u64, 12), result.max_length.?);
    try expectOptionalApprox(6.0, result.mean_length);
    try std.testing.expectEqual(@as(u64, 2), result.a);
    try std.testing.expectEqual(@as(u64, 2), result.c);
    try std.testing.expectEqual(@as(u64, 2), result.g);
    try std.testing.expectEqual(@as(u64, 2), result.t);
    try std.testing.expectEqual(@as(u64, 2), result.n);
    try std.testing.expectEqual(@as(u64, 2), result.other_bases);
    try expectOptionalApprox(0.5, result.gc_fraction);
    try std.testing.expectEqual(@as(u64, 183), result.quality_sum);
    try expectOptionalApprox(15.25, result.mean_quality);
    try std.testing.expectEqual(@as(u64, 4), result.q20_bases);
    try expectOptionalApprox(1.0 / 3.0, result.q20_fraction);
    try std.testing.expectEqual(@as(u64, 3), result.q30_bases);
    try expectOptionalApprox(0.25, result.q30_fraction);
}

test "[unit] - [statistics]: empty and ambiguous inputs expose undefined denominators" {
    var empty: zfastq.Stats = .{};
    const empty_result = empty.result();
    try std.testing.expectEqual(@as(u64, 0), empty_result.reads);
    try std.testing.expectEqual(@as(u64, 0), empty_result.bases);
    try std.testing.expect(empty_result.min_length == null);
    try std.testing.expect(empty_result.max_length == null);
    try std.testing.expect(empty_result.mean_length == null);
    try std.testing.expect(empty_result.gc_fraction == null);
    try std.testing.expect(empty_result.mean_quality == null);
    try std.testing.expect(empty_result.q20_fraction == null);
    try std.testing.expect(empty_result.q30_fraction == null);

    var ambiguous: zfastq.Stats = .{};
    try ambiguous.addRecord(record("NnRY", "!!!!"));
    const ambiguous_result = ambiguous.result();
    try std.testing.expect(ambiguous_result.gc_fraction == null);
    try std.testing.expectEqual(@as(u64, 2), ambiguous_result.n);
    try std.testing.expectEqual(@as(u64, 2), ambiguous_result.other_bases);
}

test "[failure] - [statistics]: invalid quality rejects the entire record" {
    var stats: zfastq.Stats = .{};
    try std.testing.expectError(
        error.S006InvalidQuality,
        stats.addRecord(record("ACG", "! \x7f")),
    );

    const details = stats.takeLastQualityError().?;
    try std.testing.expectEqual(@as(usize, 1), details.byte_index);
    try std.testing.expectEqual(@as(u8, 32), details.byte);
    try std.testing.expect(stats.takeLastQualityError() == null);
    try std.testing.expectEqualDeep(zfastq.StatsResult{
        .reads = 0,
        .bases = 0,
        .min_length = null,
        .max_length = null,
        .mean_length = null,
        .a = 0,
        .c = 0,
        .g = 0,
        .t = 0,
        .n = 0,
        .other_bases = 0,
        .gc_fraction = null,
        .quality_sum = 0,
        .mean_quality = null,
        .q20_bases = 0,
        .q20_fraction = null,
        .q30_bases = 0,
        .q30_fraction = null,
    }, stats.result());
}

test "[failure] - [statistics]: length mismatch precedes quality validation" {
    var stats: zfastq.Stats = .{};
    try stats.addRecord(record("A", "!"));
    const before = stats.result();

    try std.testing.expectError(
        error.S005LengthMismatch,
        stats.addRecord(record("AC", " ")),
    );
    try std.testing.expectEqualDeep(before, stats.result());
}

test "[failure] - [statistics]: every public counter addition is checked" {
    const maximum = std.math.maxInt(u64);
    const cases = [_]zfastq.Stats{
        .{ .reads = maximum },
        .{ .bases = maximum },
        .{ .a = maximum },
        .{ .c = maximum },
        .{ .g = maximum },
        .{ .t = maximum },
        .{ .n = maximum },
        .{ .other_bases = maximum },
        .{ .quality_sum = maximum },
        .{ .q20_bases = maximum },
        .{ .q30_bases = maximum },
    };
    const overflow_record = record("ACGTNX", "~~~~~~");

    for (cases) |initial| {
        var stats = initial;
        try std.testing.expectError(error.Overflow, stats.addRecord(overflow_record));
        try std.testing.expectEqualDeep(initial.result(), stats.result());
    }
}

test "[unit] - [statistics]: maximum materialized counters do not wrap" {
    const maximum = std.math.maxInt(u64);
    const stats = zfastq.Stats{
        .reads = maximum,
        .bases = maximum,
        .min_length = maximum,
        .max_length = maximum,
        .a = maximum,
        .c = maximum,
        .g = maximum,
        .t = maximum,
        .n = maximum,
        .other_bases = maximum,
        .quality_sum = maximum,
        .q20_bases = maximum,
        .q30_bases = maximum,
    };

    const result = stats.result();
    try expectOptionalApprox(1.0, result.mean_length);
    try expectOptionalApprox(0.5, result.gc_fraction);
    try expectOptionalApprox(1.0, result.mean_quality);
    try expectOptionalApprox(1.0, result.q20_fraction);
    try expectOptionalApprox(1.0, result.q30_fraction);
}

test "[unit] - [quality]: Phred+33 accepts both boundaries only" {
    try std.testing.expectEqual(@as(u8, 0), try zfastq.decodePhred33('!'));
    try std.testing.expectEqual(@as(u8, 93), try zfastq.decodePhred33('~'));
    try std.testing.expectError(error.InvalidQuality, zfastq.decodePhred33(' '));
    try std.testing.expectError(error.InvalidQuality, zfastq.decodePhred33(127));
}

test "[cli] - [stats]: plain, gzip, and fragmented stdin print identical fields" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    try writeTempFile(io, &tmp, "reads.fastq.gz", SAMPLE_FASTQ);
    try writeTempFile(io, &tmp, "reads.bin", &SAMPLE_GZIP);
    const plain_path = try tempPath(allocator, &tmp, "reads.fastq.gz");
    const gzip_path = try tempPath(allocator, &tmp, "reads.bin");

    const plain = try runCli(allocator, &.{ "stats", plain_path }, "", 1);
    const gzip = try runCli(allocator, &.{ "stats", gzip_path }, "", 1);
    const plain_stdin = try runCli(allocator, &.{ "stats", "-" }, SAMPLE_FASTQ, 1);
    const gzip_stdin = try runCli(allocator, &.{ "stats", "-" }, &SAMPLE_GZIP, 1);

    for ([_]struct {
        result: CommandResult,
        label: []const u8,
    }{
        .{ .result = plain, .label = plain_path },
        .{ .result = gzip, .label = gzip_path },
        .{ .result = plain_stdin, .label = "-" },
        .{ .result = gzip_stdin, .label = "-" },
    }) |case| {
        try std.testing.expectEqual(@as(u8, 0), case.result.exit_code);
        try std.testing.expectEqualStrings(
            try expectedStats(allocator, case.label),
            case.result.stdout,
        );
        try std.testing.expectEqual(@as(usize, 0), case.result.stderr.len);
    }
}

test "[cli] - [stats]: empty input and a zero-length record remain distinct" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const empty = try runCli(allocator, &.{ "stats", "-" }, "", 1);
    try std.testing.expectEqual(@as(u8, 0), empty.exit_code);
    try std.testing.expectEqualStrings(
        \\input: -
        \\reads: 0
        \\bases: 0
        \\min_length: -
        \\max_length: -
        \\mean_length: -
        \\a: 0
        \\c: 0
        \\g: 0
        \\t: 0
        \\n: 0
        \\other_bases: 0
        \\gc_fraction: -
        \\quality_sum: 0
        \\mean_quality: -
        \\q20_bases: 0
        \\q20_fraction: -
        \\q30_bases: 0
        \\q30_fraction: -
        \\
    , empty.stdout);
    try std.testing.expectEqual(@as(usize, 0), empty.stderr.len);

    const zero = try runCli(allocator, &.{ "stats", "-" }, "@z\n\n+\n\n", 1);
    try std.testing.expectEqual(@as(u8, 0), zero.exit_code);
    try std.testing.expectEqualStrings(
        \\input: -
        \\reads: 1
        \\bases: 0
        \\min_length: 0
        \\max_length: 0
        \\mean_length: 0.000000
        \\a: 0
        \\c: 0
        \\g: 0
        \\t: 0
        \\n: 0
        \\other_bases: 0
        \\gc_fraction: -
        \\quality_sum: 0
        \\mean_quality: -
        \\q20_bases: 0
        \\q20_fraction: -
        \\q30_bases: 0
        \\q30_fraction: -
        \\
    , zero.stdout);
    try std.testing.expectEqual(@as(usize, 0), zero.stderr.len);
}

test "[cli] - [stats]: S006 reports the exact decompressed quality-byte offset" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const invalid = "@bad\r\nAC\r\n+\r\n!\x7f\r\n";

    for ([_][]const u8{ invalid, &INVALID_QUALITY_GZIP }) |input| {
        const result = try runCli(allocator, &.{ "stats", "-" }, input, 1);

        try std.testing.expectEqual(@as(u8, 1), result.exit_code);
        try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
        try std.testing.expectEqualStrings(
            "error: -: S006: quality byte must be ASCII 33 through 126 " ++
                "(record 0, line 4, offset 14)\n",
            result.stderr,
        );
    }

    const later = try runCli(
        allocator,
        &.{ "stats", "-" },
        "@ok\nA\n+\n!\n" ++ invalid,
        1,
    );
    try std.testing.expectEqual(@as(u8, 1), later.exit_code);
    try std.testing.expectEqual(@as(usize, 0), later.stdout.len);
    try std.testing.expectEqualStrings(
        "error: -: S006: quality byte must be ASCII 33 through 126 " ++
            "(record 1, line 4, offset 24)\n",
        later.stderr,
    );
}

test "[cli] - [stats]: successful blocks survive independent higher-class failures" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    try writeTempFile(io, &tmp, "valid.fastq", SAMPLE_FASTQ);
    try writeTempFile(io, &tmp, "invalid.fastq", "@bad\nA\n+\n \n");
    const valid_path = try tempPath(allocator, &tmp, "valid.fastq");
    const invalid_path = try tempPath(allocator, &tmp, "invalid.fastq");
    const missing_path = try tempPath(allocator, &tmp, "missing.fastq");

    const result = try runCli(
        allocator,
        &.{ "stats", valid_path, invalid_path, missing_path, valid_path },
        "",
        1,
    );
    const one_block = try expectedStats(allocator, valid_path);
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqualStrings(
        try std.fmt.allocPrint(allocator, "{s}\n{s}", .{ one_block, one_block }),
        result.stdout,
    );
    try std.testing.expectEqualStrings(
        try std.fmt.allocPrint(
            allocator,
            "error: {s}: S006: quality byte must be ASCII 33 through 126 " ++
                "(record 0, line 4, offset 9)\n" ++
                "error: {s}: file not found\n",
            .{ invalid_path, missing_path },
        ),
        result.stderr,
    );
}

test "[cli] - [stats]: input labels use escaped ASCII" {
    if (builtin.os.tag == .windows) return error.SkipZigTest;

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const name = "unsafe\n\\.fastq";
    try writeTempFile(io, &tmp, name, SAMPLE_FASTQ);
    const path = try tempPath(allocator, &tmp, name);
    const result = try runCli(allocator, &.{ "stats", path }, "", 1);
    const escaped_path = try std.mem.replaceOwned(
        u8,
        allocator,
        path,
        "\n\\",
        "\\x0A\\\\",
    );

    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    try std.testing.expectEqualStrings(
        try expectedStats(allocator, escaped_path),
        result.stdout,
    );
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
}

test "[cli] - [stats]: arguments, line limits, damaged gzip, and output I/O are explicit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const missing = try runCli(allocator, &.{"stats"}, "", 1);
    try expectCommand(missing, 2, "", "error: stats requires at least one input\n");

    const duplicate_stdin = try runCli(allocator, &.{ "stats", "-", "-" }, "", 1);
    try expectCommand(
        duplicate_stdin,
        2,
        "",
        "error: standard input may appear at most once\n",
    );

    const unknown = try runCli(allocator, &.{ "stats", "--json" }, "", 1);
    try expectCommand(unknown, 2, "", "error: unknown stats option: --json\n");

    const line_limit = try runCli(
        allocator,
        &.{ "stats", "--max-line-bytes", "1", "-" },
        "@r\nAA\n+\n!!\n",
        1,
    );
    try expectCommand(line_limit, 4, "", "error: -: line length limit exceeded\n");

    const structural = try runCli(allocator, &.{ "stats", "-" }, "@r\nAA\n+\n!\n", 1);
    try expectCommand(
        structural,
        1,
        "",
        "error: -: S005: sequence and quality lengths differ " ++
            "(record 0, line 4, offset 8)\n",
    );

    const one_byte = try runCli(allocator, &.{ "stats", "-" }, "x", 1);
    try expectCommand(
        one_byte,
        1,
        "",
        "error: -: S003: header line must start with '@' " ++
            "(record 0, line 1, offset 0)\n",
    );

    var damaged = SAMPLE_GZIP;
    damaged[damaged.len - 8] ^= 1;
    const corrupt = try runCli(allocator, &.{ "stats", "-" }, &damaged, 1);
    try expectCommand(corrupt, 3, "", "error: -: I/O error\n");

    const closed = try runCliWithClosedStdout(
        allocator,
        &.{ "stats", "-" },
        SAMPLE_FASTQ,
    );
    try std.testing.expectEqual(@as(u8, 3), closed.exit_code);
    try std.testing.expectEqual(@as(usize, 0), closed.stderr.len);

    const closed_stdin = try runCliWithClosedStdin(allocator, &.{ "stats", "-" });
    try expectCommand(closed_stdin, 3, "", "error: -: I/O error\n");
}

fn record(sequence: []const u8, quality: []const u8) zfastq.Record {
    return .{
        .header = "r",
        .id = "r",
        .sequence = sequence,
        .plus = "",
        .quality = quality,
    };
}

fn expectOptionalApprox(expected: f64, actual: ?f64) !void {
    try std.testing.expectApproxEqAbs(expected, actual.?, 1e-12);
}

fn writeTempFile(
    io: std.Io,
    tmp: *std.testing.TmpDir,
    name: []const u8,
    bytes: []const u8,
) !void {
    const file = try tmp.dir.createFile(io, name, .{});
    defer file.close(io);
    try std.Io.File.writeStreamingAll(file, io, bytes);
}

fn tempPath(
    allocator: std.mem.Allocator,
    tmp: *const std.testing.TmpDir,
    name: []const u8,
) ![]const u8 {
    return std.fmt.allocPrint(
        allocator,
        ".zig-cache/tmp/{s}/{s}",
        .{ tmp.sub_path, name },
    );
}

fn expectedStats(allocator: std.mem.Allocator, label: []const u8) ![]const u8 {
    return std.fmt.allocPrint(allocator, "input: {s}\n{s}", .{ label, SAMPLE_FIELDS });
}

fn expectCommand(
    result: CommandResult,
    exit_code: u8,
    stdout: []const u8,
    stderr: []const u8,
) !void {
    try std.testing.expectEqual(exit_code, result.exit_code);
    try std.testing.expectEqualStrings(stdout, result.stdout);
    try std.testing.expectEqualStrings(stderr, result.stderr);
}
