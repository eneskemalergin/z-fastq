//! Public aggregate-statistics contracts.

const std = @import("std");
const zfastq = @import("z-fastq");
const cli = @import("utilities.zig");
const CommandResult = cli.CommandResult;
const runCliWithStdin = cli.runWithStdin;
const runCliWithClosedStdout = cli.runWithClosedStdout;

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

const EMPTY_FIELDS =
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

test "[edge] - [statistics]: record updates cross 32 bits without losing counters or ratios" {
    const crossing: u64 = 1 << 32;
    // Each seeded history consists of ACGTNX records with six Phred-93 qualities.
    for ([_]u64{ crossing / 558, crossing / 6, crossing - 1 }) |initial_reads| {
        errdefer std.debug.print("stats 32-bit crossing: initial reads {d}\n", .{initial_reads});
        var stats = zfastq.Stats{
            .reads = initial_reads,
            .bases = initial_reads * 6,
            .min_length = 6,
            .max_length = 6,
            .a = initial_reads,
            .c = initial_reads,
            .g = initial_reads,
            .t = initial_reads,
            .n = initial_reads,
            .other_bases = initial_reads,
            .quality_sum = initial_reads * 558,
            .q20_bases = initial_reads * 6,
            .q30_bases = initial_reads * 6,
        };
        for (1..4) |added| {
            try stats.addRecord(record("ACGTNX", "~~~~~~"));
            const reads = initial_reads + added;
            try std.testing.expectEqualDeep(zfastq.StatsResult{
                .reads = reads,
                .bases = reads * 6,
                .min_length = 6,
                .max_length = 6,
                .mean_length = 6,
                .a = reads,
                .c = reads,
                .g = reads,
                .t = reads,
                .n = reads,
                .other_bases = reads,
                .gc_fraction = 0.5,
                .quality_sum = reads * 558,
                .mean_quality = 93,
                .q20_bases = reads * 6,
                .q20_fraction = 1,
                .q30_bases = reads * 6,
                .q30_fraction = 1,
            }, stats.result());
        }

        const before = stats.result();
        try std.testing.expectError(error.S006InvalidQuality, stats.addRecord(record("ACGTNX", "~~~ ~~")));
        try std.testing.expectEqualDeep(before, stats.result());
        const details = stats.takeLastQualityError().?;
        try std.testing.expectEqual(@as(usize, 3), details.byte_index);
        try std.testing.expectEqual(@as(u8, 32), details.byte);
        try std.testing.expect(stats.takeLastQualityError() == null);
    }
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

test "[cli] - [stats]: files and stdin preserve prefixes, gzip members, and exact fields" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    var members: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &members, "", .{});
    const empty_member_end = members.items.len;
    try cli.appendGzipMember(allocator, &members, SAMPLE_FASTQ[0..1], .{});
    try cli.appendGzipMember(allocator, &members, SAMPLE_FASTQ[1..], .{});
    const later_header = SAMPLE_GZIP ++ [_]u8{ 0x1f, 0x8b };

    const cases = [_]struct {
        name: []const u8 = "input",
        bytes: []const u8,
        fields: []const u8 = SAMPLE_FIELDS,
        exit_code: u8 = 0,
        message: []const u8 = "",
    }{
        .{ .name = "reads.fastq.gz", .bytes = SAMPLE_FASTQ },
        .{ .name = "reads.bin", .bytes = &SAMPLE_GZIP },
        .{ .bytes = members.items },
        .{ .bytes = "", .fields = EMPTY_FIELDS },
        .{ .bytes = members.items[0..empty_member_end], .fields = EMPTY_FIELDS },
        .{
            .bytes = "@r",
            .exit_code = 1,
            .message = "S004: unexpected end of file in sequence line (record 0, line 2, offset 2)\n",
        },
        .{ .bytes = SAMPLE_GZIP[0..9], .exit_code = 3, .message = "I/O error\n" },
        .{ .bytes = &later_header, .exit_code = 3, .message = "I/O error\n" },
    };
    for (cases, 0..) |case, index| {
        try tmp.dir.writeFile(io, .{ .sub_path = case.name, .data = case.bytes });
        const path = try cli.tempPath(allocator, &tmp.sub_path, case.name);
        for ([_][]const u8{ path, "-" }, 0..) |label, transport| {
            errdefer std.debug.print("stats input case {d}, input {s}\n", .{ index, label });
            const result = if (transport == 0)
                try cli.run(allocator, &.{ "stats", path })
            else
                try runCliWithStdin(allocator, &.{ "stats", "-" }, case.bytes, 1);
            const stdout = if (case.exit_code == 0) try std.fmt.allocPrint(
                allocator,
                "input: {s}\n{s}",
                .{ label, case.fields },
            ) else "";
            const stderr = if (case.message.len == 0) "" else try std.fmt.allocPrint(
                allocator,
                "error: {s}: {s}",
                .{ label, case.message },
            );
            try cli.expectResult(result, case.exit_code, stdout, stderr);
        }
    }
}

test "[cli] - [stats-json]: JSON shapes match the sample, empty, and zero-length inputs" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct {
        args: []const []const u8,
        bytes: []const u8,
        shape: enum { sample, empty, zero_length },
    }{
        .{ .args = &.{ "stats", "--json", "--max-line-bytes", "64", "-" }, .bytes = SAMPLE_FASTQ, .shape = .sample },
        .{ .args = &.{ "stats", "--json", "-" }, .bytes = &SAMPLE_GZIP, .shape = .sample },
        .{ .args = &.{ "stats", "--json", "-" }, .bytes = "", .shape = .empty },
        .{ .args = &.{ "stats", "--json", "-" }, .bytes = "@z\n\n+\n\n", .shape = .zero_length },
    };
    for (cases) |case| {
        const result = try runCliWithStdin(allocator, case.args, case.bytes, 1);
        try std.testing.expectEqual(@as(u8, 0), result.exit_code);
        try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
        const results = try cli.parseJsonResults(allocator, result.stdout, "z-fastq/stats-v1");
        try std.testing.expectEqual(@as(usize, 1), results.len);
        switch (case.shape) {
            .sample => try expectSampleJsonStats(results[0], "-"),
            .empty => try expectEmptyJsonStats(results[0], "-"),
            .zero_length => try expectZeroLengthJsonStats(results[0]),
        }
    }
}

test "[cli] - [stats]: empty input and a zero-length record remain distinct" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const empty = try runCliWithStdin(allocator, &.{ "stats", "-" }, "", 1);
    try std.testing.expectEqual(@as(u8, 0), empty.exit_code);
    try std.testing.expectEqualStrings("input: -\n" ++ EMPTY_FIELDS, empty.stdout);
    try std.testing.expectEqual(@as(usize, 0), empty.stderr.len);

    const zero = try runCliWithStdin(allocator, &.{ "stats", "-" }, "@z\n\n+\n\n", 1);
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
test "[cli] - [stats]: source counts and exact decimal ties agree with human and JSON fields" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const cases = [_]struct {
        reads: usize,
        bases: usize,
        mean_length: []const u8,
        fraction: []const u8,
        mean_quality: []const u8,
    }{
        .{ .reads = 1, .bases = 127, .mean_length = "127.000000", .fraction = "0.007874", .mean_quality = "0.244094" },
        .{ .reads = 1, .bases = 128, .mean_length = "128.000000", .fraction = "0.007813", .mean_quality = "0.242188" },
        .{ .reads = 1, .bases = 129, .mean_length = "129.000000", .fraction = "0.007752", .mean_quality = "0.240310" },
        .{ .reads = 128, .bases = 1, .mean_length = "0.007813", .fraction = "1.000000", .mean_quality = "31.000000" },
    };
    for (cases) |case| {
        errdefer std.debug.print("stats reads={d} bases={d}\n", .{ case.reads, case.bases });
        var input: std.ArrayList(u8) = .empty;
        try input.ensureTotalCapacity(allocator, case.bases * 2 + case.reads * 10);
        // One C and one Phred-31 byte make the ratios exact decimal ties at 128 bases.
        try input.appendSlice(allocator, "@r\nC");
        try input.appendNTimes(allocator, 'A', case.bases - 1);
        try input.appendSlice(allocator, "\n+\n@");
        try input.appendNTimes(allocator, '!', case.bases - 1);
        try input.append(allocator, '\n');
        for (1..case.reads) |_| try input.appendSlice(allocator, "@e\n\n+\n\n");
        const min_length = if (case.reads == 1) case.bases else 0;
        const expected = try std.fmt.allocPrint(
            allocator,
            "input: -\nreads: {d}\nbases: {d}\nmin_length: {d}\nmax_length: {d}\nmean_length: {s}\n" ++
                "a: {d}\nc: 1\ng: 0\nt: 0\nn: 0\nother_bases: 0\ngc_fraction: {s}\n" ++
                "quality_sum: 31\nmean_quality: {s}\nq20_bases: 1\nq20_fraction: {s}\n" ++
                "q30_bases: 1\nq30_fraction: {s}\n",
            .{ case.reads, case.bases, min_length, case.bases, case.mean_length, case.bases - 1, case.fraction, case.mean_quality, case.fraction, case.fraction },
        );
        var gzip: std.ArrayList(u8) = .empty;
        try cli.appendGzipMember(allocator, &gzip, input.items, .{});
        for ([_][]const u8{ input.items, gzip.items }) |bytes| {
            try cli.expectResult(try runCliWithStdin(allocator, &.{ "stats", "-" }, bytes, bytes.len), 0, expected, "");
            const json = try runCliWithStdin(allocator, &.{ "stats", "--json", "-" }, bytes, bytes.len);
            try std.testing.expectEqual(@as(u8, 0), json.exit_code);
            try std.testing.expectEqualStrings("", json.stderr);
            const results = try cli.parseJsonResults(allocator, json.stdout, "z-fastq/stats-v1");
            try std.testing.expectEqual(@as(usize, 1), results.len);
            const object = results[0].object;
            try cli.expectJsonObjectKeys(object, &STATS_RESULT_KEYS);
            try cli.expectJsonString(object.get("input"), "-");
            try cli.expectJsonString(object.get("status"), "ok");
            inline for (
                .{ "reads", "bases", "min_length", "max_length", "a", "c", "g", "t", "n", "other_bases", "quality_sum", "q20_bases", "q30_bases" },
                .{ case.reads, case.bases, min_length, case.bases, case.bases - 1, 1, 0, 0, 0, 0, 31, 1, 1 },
            ) |field, count| {
                try cli.expectJsonInteger(object.get(field), count);
            }
            const bases: f64 = @floatFromInt(case.bases);
            const reads: f64 = @floatFromInt(case.reads);
            try expectJsonNumber(object.get("mean_length"), bases / reads);
            try expectJsonNumber(object.get("gc_fraction"), 1.0 / bases);
            try expectJsonNumber(object.get("mean_quality"), 31.0 / bases);
            try expectJsonNumber(object.get("q20_fraction"), 1.0 / bases);
            try expectJsonNumber(object.get("q30_fraction"), 1.0 / bases);
        }
    }
}

test "[cli] - [stats]: S006 reports the exact decompressed quality-byte offset" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const invalid = "@bad\r\nAC\r\n+\r\n!\x7f\r\n";
    const invalid_lf = try runCliWithStdin(
        allocator,
        &.{ "stats", "-" },
        "@bad\nAC\n+\n!\x7f\n",
        1,
    );
    try std.testing.expectEqual(@as(u8, 1), invalid_lf.exit_code);
    try std.testing.expectEqual(@as(usize, 0), invalid_lf.stdout.len);
    try std.testing.expectEqualStrings(
        "error: -: S006: quality byte must be ASCII 33 through 126 " ++
            "(record 0, line 4, offset 11)\n",
        invalid_lf.stderr,
    );

    for ([_][]const u8{ invalid, &INVALID_QUALITY_GZIP }) |input| {
        const result = try runCliWithStdin(allocator, &.{ "stats", "-" }, input, 1);

        try std.testing.expectEqual(@as(u8, 1), result.exit_code);
        try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
        try std.testing.expectEqualStrings(
            "error: -: S006: quality byte must be ASCII 33 through 126 " ++
                "(record 0, line 4, offset 14)\n",
            result.stderr,
        );
    }

    const later = try runCliWithStdin(
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

test "[cli] - [stats]: embedded quality endings preserve structural error precedence" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const cases = [_]struct { input: []const u8, code: []const u8, offset: usize }{
        .{ .input = "@r\nACG\n+\n!\n!\n", .code = "S005", .offset = 9 },
        .{ .input = "@r\nACG\n+\n \n!\n", .code = "S005", .offset = 9 },
        .{ .input = "@r\nACG\n+\n!\n!\r\n", .code = "S005", .offset = 9 },
        .{ .input = "@r\nACG\n+\n!\n\r\n", .code = "S005", .offset = 9 },
        .{ .input = "@r\nAC\n+\n!\r\r\n", .code = "S006", .offset = 9 },
        .{ .input = "@r\r\nACG\r\n+\r\n!\x7f!\r\n", .code = "S006", .offset = 13 },
    };
    for (cases) |case| {
        for ([_][]const u8{ "", "@ok\nA\n+\n!\n" }) |prefix| {
            const input = try std.mem.concat(allocator, u8, &.{ prefix, case.input });
            const record_index: u64 = @intFromBool(prefix.len != 0);
            const offset = prefix.len + case.offset;
            const message = if (std.mem.eql(u8, case.code, "S005"))
                "sequence and quality lengths differ"
            else
                "quality byte must be ASCII 33 through 126";
            const expected = try std.fmt.allocPrint(
                allocator,
                "error: -: {s}: {s} (record {d}, line 4, offset {d})\n",
                .{ case.code, message, record_index, offset },
            );
            const human = try runCliWithStdin(allocator, &.{ "stats", "-" }, input, input.len);
            try cli.expectResult(human, 1, "", expected);
            const json = try runCliWithStdin(allocator, &.{ "stats", "--json", "-" }, input, input.len);
            try std.testing.expectEqual(@as(u8, 1), json.exit_code);
            try std.testing.expectEqualStrings("", json.stderr);
            const results = try cli.parseJsonResults(allocator, json.stdout, "z-fastq/stats-v1");
            try std.testing.expectEqual(@as(usize, 1), results.len);
            try cli.expectJsonFailure(results[0], "-", case.code, record_index, offset, 4);
        }
    }
}

test "[cli] - [stats]: every input runs and the highest exit class wins" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    // The limit file has an 8-byte sequence line; the sample's lines fit the 6-byte limit used below.
    try tmp.dir.writeFile(io, .{ .sub_path = "valid.fastq", .data = SAMPLE_FASTQ });
    try tmp.dir.writeFile(io, .{ .sub_path = "invalid.fastq", .data = "@bad\nA\n+\n \n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "limit.fastq", .data = "@r\nACGTACGT\n+\n!!!!!!!!\n" });
    const valid_path = try cli.tempPath(allocator, &tmp.sub_path, "valid.fastq");
    const invalid_path = try cli.tempPath(allocator, &tmp.sub_path, "invalid.fastq");
    const limit_path = try cli.tempPath(allocator, &tmp.sub_path, "limit.fastq");
    const missing_path = try cli.tempPath(allocator, &tmp.sub_path, "missing.fastq");
    const one_block = try expectedStats(allocator, valid_path);
    const two_blocks = try std.fmt.allocPrint(allocator, "{s}\n{s}", .{ one_block, one_block });

    const cases = [_]struct {
        args: []const []const u8,
        exit_code: u8,
        stdout: []const u8,
        stderr: []const u8,
    }{
        .{
            .args = &.{ "stats", valid_path, invalid_path, missing_path, valid_path },
            .exit_code = 3,
            .stdout = two_blocks,
            .stderr = try std.fmt.allocPrint(
                allocator,
                "error: {s}: S006: quality byte must be ASCII 33 through 126 " ++
                    "(record 0, line 4, offset 9)\n" ++
                    "error: {s}: file not found\n",
                .{ invalid_path, missing_path },
            ),
        },
        .{
            .args = &.{ "stats", valid_path, missing_path, valid_path },
            .exit_code = 3,
            .stdout = two_blocks,
            .stderr = try std.fmt.allocPrint(allocator, "error: {s}: file not found\n", .{missing_path}),
        },
        .{
            .args = &.{ "stats", "--max-line-bytes", "6", limit_path, missing_path, valid_path },
            .exit_code = 4,
            .stdout = one_block,
            .stderr = try std.fmt.allocPrint(
                allocator,
                "error: {s}: line length limit exceeded\n" ++
                    "error: {s}: file not found\n",
                .{ limit_path, missing_path },
            ),
        },
    };
    for (cases) |case| {
        const result = try cli.run(allocator, case.args);
        try cli.expectResult(result, case.exit_code, case.stdout, case.stderr);
    }
}

test "[cli] - [stats]: input labels use escaped ASCII" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const name = "unsafe\n\\.fastq";
    try tmp.dir.writeFile(io, .{ .sub_path = name, .data = SAMPLE_FASTQ });
    const path = try cli.tempPath(allocator, &tmp.sub_path, name);
    const result = try runCliWithStdin(allocator, &.{ "stats", path }, "", 1);
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

    const json_result = try runCliWithStdin(allocator, &.{ "stats", "--json", path }, "", 1);
    try std.testing.expectEqual(@as(u8, 0), json_result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), json_result.stderr.len);
    const results = try cli.parseJsonResults(allocator, json_result.stdout, "z-fastq/stats-v1");
    try std.testing.expectEqual(@as(usize, 1), results.len);
    const object = results[0].object;
    try cli.expectJsonString(object.get("input"), escaped_path);
}

test "[cli] - [stats]: line limits and structural errors exit with their documented status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const line_limit = try runCliWithStdin(
        allocator,
        &.{ "stats", "--max-line-bytes", "1", "-" },
        "@r\nAA\n+\n!!\n",
        1,
    );
    try cli.expectResult(line_limit, 4, "", "error: -: line length limit exceeded\n");

    const structural = try runCliWithStdin(allocator, &.{ "stats", "-" }, "@r\nAA\n+\n!\n", 1);
    try cli.expectResult(
        structural,
        1,
        "",
        "error: -: S005: sequence and quality lengths differ " ++
            "(record 0, line 4, offset 8)\n",
    );

    const mismatch_and_invalid = try runCliWithStdin(
        allocator,
        &.{ "stats", "-" },
        "@r\nAA\n+\n \n",
        1,
    );
    try cli.expectResult(
        mismatch_and_invalid,
        1,
        "",
        "error: -: S005: sequence and quality lengths differ " ++
            "(record 0, line 4, offset 8)\n",
    );

    const one_byte = try runCliWithStdin(allocator, &.{ "stats", "-" }, "x", 1);
    try cli.expectResult(
        one_byte,
        1,
        "",
        "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n",
    );

    const empty_identifier = try runCliWithStdin(
        allocator,
        &.{ "stats", "-" },
        "@ description\nA\n+\n!\n",
        1,
    );
    try cli.expectResult(
        empty_identifier,
        1,
        "",
        "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n",
    );
}

test "[cli] - [stats]: a damaged gzip member is an I/O failure in human and JSON output" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var damaged = SAMPLE_GZIP;
    damaged[cli.gzipCrcOffset(damaged.len)] ^= 1;
    const corrupt = try runCliWithStdin(allocator, &.{ "stats", "-" }, &damaged, 1);
    try cli.expectResult(corrupt, 3, "", "error: -: I/O error\n");

    const corrupt_json = try runCliWithStdin(allocator, &.{ "stats", "--json", "-" }, &damaged, 1);
    try std.testing.expectEqual(@as(u8, 3), corrupt_json.exit_code);
    try std.testing.expectEqual(@as(usize, 0), corrupt_json.stderr.len);
    const corrupt_results = try cli.parseJsonResults(allocator, corrupt_json.stdout, "z-fastq/stats-v1");
    try std.testing.expectEqual(@as(usize, 1), corrupt_results.len);
    try cli.expectJsonFailureWithoutLocation(
        corrupt_results[0],
        "-",
        "io_error",
        "I/O error",
    );
}

test "[cli] - [stats]: output failure preserves an observed limit status" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const long_directory = ("\x01" ** 200 ++ "/") ** 17;
    try tmp.dir.createDirPath(io, long_directory);
    try tmp.dir.writeFile(io, .{ .sub_path = "valid.fastq", .data = "@a\nA\n+\n!\n" });
    const valid = try cli.tempPath(allocator, &tmp.sub_path, "valid.fastq");

    for ([_][]const u8{ "limit.fastq", long_directory ++ "limit.fastq" }) |name| {
        try tmp.dir.writeFile(io, .{ .sub_path = name, .data = "@abcde\nA\n+\n!\n" });
        const path = try cli.tempPath(allocator, &tmp.sub_path, name);
        const args: []const []const u8 = &.{ "stats", "--json", "--max-line-bytes", "4", path, valid };
        const captured = try cli.run(allocator, args);
        try std.testing.expectEqual(@as(u8, 4), captured.exit_code);
        try std.testing.expectEqualStrings("", captured.stderr);
        try std.testing.expectEqual(name.len > 200, captured.stdout.len > cli.JSON_OUTPUT_BUFFER_BYTES);
        const results = try cli.parseJsonResults(allocator, captured.stdout, "z-fastq/stats-v1");
        try std.testing.expectEqual(@as(usize, 2), results.len);
        try cli.expectJsonString(results[0].object.get("error").?.object.get("code"), "line_limit");
        try cli.expectJsonString(results[1].object.get("status"), "ok");
        try cli.expectResult(try runCliWithClosedStdout(allocator, args, ""), 4, "", "error: standard output: I/O error\n");
    }

    const path = try cli.tempPath(allocator, &tmp.sub_path, "limit.fastq");
    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: line length limit exceeded\nerror: standard output: I/O error\n", .{path});
    try cli.expectResult(
        try runCliWithClosedStdout(allocator, &.{ "stats", "--max-line-bytes", "4", path, valid }, ""),
        4,
        "",
        diagnostic,
    );
}

test "[cli] - [stats-json]: mixed results preserve order, variants, and stream separation" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "valid.fastq", .data = SAMPLE_FASTQ });
    try tmp.dir.writeFile(io, .{ .sub_path = "invalid.fastq", .data = "@bad\nA\n+\n \n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "long.fastq", .data = "@r\nAAAAAAA\n+\n!!!!!!!\n" });

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const valid = try cli.tempPath(allocator, &tmp.sub_path, "valid.fastq");
    const invalid = try cli.tempPath(allocator, &tmp.sub_path, "invalid.fastq");
    const long = try cli.tempPath(allocator, &tmp.sub_path, "long.fastq");
    const missing = try cli.tempPath(allocator, &tmp.sub_path, "missing.fastq");

    const result = try runCliWithStdin(
        allocator,
        &.{
            "stats",
            "--json",
            "--max-line-bytes",
            "6",
            valid,
            invalid,
            missing,
            long,
            valid,
        },
        "",
        1,
    );
    try std.testing.expectEqual(@as(u8, 4), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    const results = try cli.parseJsonResults(allocator, result.stdout, "z-fastq/stats-v1");
    try std.testing.expectEqual(@as(usize, 5), results.len);
    try expectSampleJsonStats(results[0], valid);
    try cli.expectJsonFailure(results[1], invalid, "S006", 0, 9, 4);
    try cli.expectJsonFailureWithoutLocation(results[2], missing, "io_error", "file not found");
    try cli.expectJsonFailureWithoutLocation(
        results[3],
        long,
        "line_limit",
        "line length limit exceeded",
    );
    try expectSampleJsonStats(results[4], valid);
}

fn expectZeroLengthJsonStats(value: std.json.Value) !void {
    const object = value.object;
    try cli.expectJsonObjectKeys(object, &STATS_RESULT_KEYS);
    try cli.expectJsonString(object.get("input"), "-");
    try cli.expectJsonString(object.get("status"), "ok");
    try cli.expectJsonInteger(object.get("reads"), 1);
    inline for (.{ "bases", "min_length", "max_length", "a", "c", "g", "t", "n", "other_bases", "quality_sum", "q20_bases", "q30_bases" }) |field| {
        try cli.expectJsonInteger(object.get(field), 0);
    }
    try expectJsonNumber(object.get("mean_length"), 0.0);
    inline for (.{ "gc_fraction", "mean_quality", "q20_fraction", "q30_fraction" }) |field| {
        try std.testing.expect(object.get(field).? == .null);
    }
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

fn expectedStats(allocator: std.mem.Allocator, label: []const u8) ![]const u8 {
    return std.fmt.allocPrint(allocator, "input: {s}\n{s}", .{ label, SAMPLE_FIELDS });
}

const STATS_RESULT_KEYS = [_][]const u8{
    "input",
    "status",
    "reads",
    "bases",
    "min_length",
    "max_length",
    "mean_length",
    "a",
    "c",
    "g",
    "t",
    "n",
    "other_bases",
    "gc_fraction",
    "quality_sum",
    "mean_quality",
    "q20_bases",
    "q20_fraction",
    "q30_bases",
    "q30_fraction",
};

fn expectSampleJsonStats(value: std.json.Value, input: []const u8) !void {
    const object = switch (value) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try cli.expectJsonObjectKeys(object, &STATS_RESULT_KEYS);
    try cli.expectJsonString(object.get("input"), input);
    try cli.expectJsonString(object.get("status"), "ok");
    try cli.expectJsonInteger(object.get("reads"), 1);
    try cli.expectJsonInteger(object.get("bases"), 6);
    try cli.expectJsonInteger(object.get("min_length"), 6);
    try cli.expectJsonInteger(object.get("max_length"), 6);
    try expectJsonNumber(object.get("mean_length"), 6.0);
    inline for (&.{ "a", "c", "g", "t", "n", "other_bases" }) |field| {
        try cli.expectJsonInteger(object.get(field), 1);
    }
    try expectJsonNumber(object.get("gc_fraction"), 0.5);
    try cli.expectJsonInteger(object.get("quality_sum"), 183);
    try expectJsonNumber(object.get("mean_quality"), 30.5);
    try cli.expectJsonInteger(object.get("q20_bases"), 4);
    try expectJsonNumber(object.get("q20_fraction"), 2.0 / 3.0);
    try cli.expectJsonInteger(object.get("q30_bases"), 3);
    try expectJsonNumber(object.get("q30_fraction"), 0.5);
}

fn expectEmptyJsonStats(value: std.json.Value, input: []const u8) !void {
    const object = switch (value) {
        .object => |object| object,
        else => return error.UnexpectedJsonShape,
    };
    try cli.expectJsonObjectKeys(object, &STATS_RESULT_KEYS);
    try cli.expectJsonString(object.get("input"), input);
    try cli.expectJsonString(object.get("status"), "ok");
    inline for (&.{
        "reads",
        "bases",
        "a",
        "c",
        "g",
        "t",
        "n",
        "other_bases",
        "quality_sum",
        "q20_bases",
        "q30_bases",
    }) |field| {
        try cli.expectJsonInteger(object.get(field), 0);
    }
    inline for (&.{
        "min_length",
        "max_length",
        "mean_length",
        "gc_fraction",
        "mean_quality",
        "q20_fraction",
        "q30_fraction",
    }) |field| {
        const field_value = object.get(field) orelse return error.UnexpectedJsonShape;
        try std.testing.expect(field_value == .null);
    }
}

fn expectJsonNumber(value: ?std.json.Value, expected: f64) !void {
    const number = switch (value orelse return error.UnexpectedJsonShape) {
        .float => |number| number,
        .integer => |number| @as(f64, @floatFromInt(number)),
        else => return error.UnexpectedJsonShape,
    };
    try std.testing.expectApproxEqAbs(expected, number, 1e-12);
}
