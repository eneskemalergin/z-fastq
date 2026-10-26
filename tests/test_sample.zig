//! Private `sample` RNG vectors and installed `z-fastq sample` behavior.

const std = @import("std");
const sampling = @import("sample_internal");
const cli = @import("utilities.zig");

const BASIC_PATH = "tests/data/synthetic/basic_valid.fastq";
const EMPTY_PATH = "tests/data/synthetic/empty_valid.fastq";
const EMPTY_RECORD = "@empty\n\n+\n\n";
const OUTPUT_BUFFER_BYTES = 64 * 1024;
const READ_BUFFER_BYTES = 256 * 1024;

const BASIC_SELECTED_HALF =
    \\@read1 run=1 lane=1
    \\ACGTACGTACGT
    \\+
    \\IIIIIIIIIIII
    \\@read3_with_plus_header
    \\AAAA
    \\+read3_with_plus_header
    \\!!!!
    \\@read5_phred_zero
    \\A
    \\+
    \\!
    \\
;
const BASIC_EXACT_ONE =
    \\@read5_phred_zero
    \\A
    \\+
    \\!
    \\
;
const BASIC_EXACT_TWO =
    \\@read3_with_plus_header
    \\AAAA
    \\+read3_with_plus_header
    \\!!!!
    \\@read5_phred_zero
    \\A
    \\+
    \\!
    \\
;

fn appendPair(
    allocator: std.mem.Allocator,
    r1: ?*std.ArrayList(u8),
    r2: ?*std.ArrayList(u8),
    interleaved: *std.ArrayList(u8),
    pair_index: usize,
) !void {
    var record1_buffer: [128]u8 = undefined;
    const record1 = try std.fmt.bufPrint(
        &record1_buffer,
        "@pair{d}/1 lane=left\nAC\n+left-{d}\n!!\n",
        .{ pair_index, pair_index },
    );
    var record2_buffer: [128]u8 = undefined;
    const record2 = try std.fmt.bufPrint(
        &record2_buffer,
        "@pair{d}/2 lane=right\nGT\n+right-{d}\n##\n",
        .{ pair_index, pair_index },
    );
    if (r1) |output| try output.appendSlice(allocator, record1);
    if (r2) |output| try output.appendSlice(allocator, record2);
    try interleaved.appendSlice(allocator, record1);
    try interleaved.appendSlice(allocator, record2);
}

fn appendSelectedPairs(
    allocator: std.mem.Allocator,
    output: *std.ArrayList(u8),
    indexes: []const u8,
) !void {
    for (indexes) |pair_index| {
        try appendPair(allocator, null, null, output, pair_index);
    }
}

test "[cli] - [sample]: stdout aliases of the input reject before reading or writing" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const input = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const payload = "@same\nAC\n+\nII\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = payload });
    try tmp.dir.hardLink("r1.fastq", tmp.dir, "r1-hard.fastq", io, .{});
    try tmp.dir.symLink(io, "r1.fastq", "r1-link.fastq", .{});
    const args = [_][]const u8{ "sample", "--count", "1", input };
    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: input file is output file\n", .{input});

    for ([_][]const u8{ "r1-hard.fastq", "r1-link.fastq" }) |alias| {
        for ([_]bool{ false, true }) |append| {
            try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = payload });
            const output: std.Io.File = .{
                .handle = try std.posix.openat(tmp.dir.handle, alias, .{
                    .ACCMODE = .WRONLY,
                    .TRUNC = !append,
                    .APPEND = append,
                    .CLOEXEC = true,
                }, 0),
                .flags = .{ .nonblocking = false },
            };
            defer output.close(io);

            try cli.expectResult(try cli.runWithStdoutFile(allocator, &args, output, null), 2, "", diagnostic);
            try std.testing.expectEqualStrings(if (append) payload else "", try tmp.dir.readFileAlloc(io, "r1.fastq", allocator, .limited(1024)));

            const stderr_file = try cli.openAppendOutput(tmp.dir, "r1.fastq");
            defer stderr_file.close(io);
            const result = try cli.runWith(allocator, &args, .{ .stdout = .{ .file = output }, .stderr = .{ .file = stderr_file } });
            try std.testing.expectEqual(@as(u8, 2), result.exit_code);
            try std.testing.expectEqualStrings(if (append) payload else "", try tmp.dir.readFileAlloc(io, "r1.fastq", allocator, .limited(1024)));
        }
    }
}

test "[cli] - [sample]: shared input and output aliases preserve the paired-input error" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const payload = "@same\nAC\n+\nII\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });
    const output = try cli.openAppendOutput(tmp.dir, "input.fastq");
    defer output.close(io);
    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{path});
    for ([_][2][]const u8{ .{ "--fraction", "0" }, .{ "--fraction", "1" }, .{ "--count", "0" }, .{ "--count", "1" } }) |selection| {
        const args = [_][]const u8{ "sample", "--paired", selection[0], selection[1], path, path };
        try cli.expectResult(try cli.runWithStdoutFile(allocator, &args, output, null), 2, "", diagnostic);
        const result = try cli.runWith(allocator, &args, .{ .stdout = .{ .file = output }, .stderr = .{ .file = output } });
        try std.testing.expectEqual(@as(u8, 2), result.exit_code);
        try std.testing.expectEqualStrings(payload, try tmp.dir.readFileAlloc(io, "input.fastq", allocator, .limited(1024)));
    }
}

test "[cli] - [sample]: redirected stdin aliases reject without consuming input" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const other = try cli.tempPath(allocator, &tmp.sub_path, "other.fastq");
    const payload = "@same\nAC\n+\nII\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });
    try tmp.dir.writeFile(io, .{ .sub_path = "other.fastq", .data = payload });
    const output = try tmp.dir.openFile(io, "input.fastq", .{ .mode = .write_only });
    defer output.close(io);
    for ([_][]const []const u8{
        &.{ "sample", "--fraction", "1", "-" },
        &.{ "sample", "--interleaved", "--fraction", "0", "-" },
        &.{ "sample", "--paired", "--fraction", "1", "-", other },
        &.{ "sample", "--paired", "--fraction", "0", other, "-" },
    }) |args| {
        const input = try tmp.dir.openFile(io, "input.fastq", .{});
        defer input.close(io);
        try cli.expectResult(try cli.runWithStdoutFile(allocator, args, output, input), 2, "", "error: -: input file is output file\n");
        var first: [1]u8 = undefined;
        try std.testing.expectEqual(@as(usize, 1), try input.readStreaming(io, &.{&first}));
        try std.testing.expectEqual(@as(u8, '@'), first[0]);
        try std.testing.expectEqualStrings(payload, try tmp.dir.readFileAlloc(io, "input.fastq", allocator, .limited(1024)));
    }
}

test "[unit] - [MT19937-64]: scalar seed 5489 matches the published vector" {
    const expected = [_]u64{
        14514284786278117030,
        4620546740167642908,
        13109570281517897720,
        17462938647148434322,
        355488278567739596,
        7469126240319926998,
        4635995468481642529,
        418970542659199878,
        9604170989252516556,
        6358044926049913402,
        5058016125798318033,
        10349215569089701407,
    };
    var generator = sampling.Mt19937_64.init(5489);
    for (expected) |value| try std.testing.expectEqual(value, generator.nextU64());
}

test "[edge] - [MT19937-64]: maximum seed matches the frozen sequence" {
    const expected = [_]u64{
        478026398904862820,
        13243134898385798468,
        709236020254955927,
        9482188692832154854,
        17279096482229114326,
        9673544723405839539,
        5170222943873112136,
        7179667202894142212,
    };
    var generator = sampling.Mt19937_64.init(std.math.maxInt(u64));
    for (expected) |value| try std.testing.expectEqual(value, generator.nextU64());
}

test "[unit] - [fraction selector]: seed 11 matches frozen seqtk indexes" {
    const selected = [_]usize{ 0, 2, 4, 5, 10, 11, 13, 14, 15 };
    var selector = sampling.Selector.init(
        try sampling.Fraction.parse("0.5"),
        11,
    );
    var selected_pos: usize = 0;
    for (0..16) |record_index| {
        const expected = selected_pos < selected.len and selected[selected_pos] == record_index;
        try std.testing.expectEqual(expected, selector.selectRecord());
        if (expected) selected_pos += 1;
    }
    try std.testing.expectEqual(selected.len, selected_pos);
}

test "[edge] - [fraction selector]: a draw equal to the fraction is not selected" {
    var generator = sampling.Mt19937_64.init(11);
    const draw = generator.nextUnitFloat();
    var equal = sampling.Selector.init(.{ .probability = draw }, 11);
    var above = sampling.Selector.init(.{
        .probability = std.math.nextAfter(f64, draw, 1.0),
    }, 11);

    try std.testing.expect(!equal.selectRecord());
    try std.testing.expect(above.selectRecord());
}

test "[edge] - [fraction selector]: parsed zero and one avoid random state" {
    var tiny = [_]u8{'0'} ** 403;
    tiny[1] = '.';
    tiny[tiny.len - 1] = '1';
    var none = sampling.Selector.init(try sampling.Fraction.parse("0.000"), 11);
    var all = sampling.Selector.init(try sampling.Fraction.parse("1.000"), 11);
    const rounded_none = try sampling.Fraction.parse(&tiny);
    const rounded_all = try sampling.Fraction.parse("0.99999999999999999");

    try std.testing.expectEqual(.none, std.meta.activeTag(none));
    try std.testing.expectEqual(.all, std.meta.activeTag(all));
    try std.testing.expectEqual(.none, std.meta.activeTag(rounded_none));
    try std.testing.expectEqual(.all, std.meta.activeTag(rounded_all));
    for (0..4) |_| {
        try std.testing.expect(!none.selectRecord());
        try std.testing.expect(all.selectRecord());
    }
}

test "[unit] - [sample numbers]: fraction and seed grammars are exact" {
    for ([_][]const u8{ "0", "1", "0.0", "0.000001", "0.5", "0.999999", "1.0", "1.000" }) |text| {
        _ = try sampling.Fraction.parse(text);
    }
    for ([_][]const u8{
        "",
        "00",
        "01",
        "0.",
        ".5",
        "1.",
        "1.1",
        "2",
        "+0.5",
        "-0.5",
        "0e0",
        "NaN",
        "inf",
        " 0.5",
        "0.5 ",
    }) |text| {
        try std.testing.expectError(error.InvalidFraction, sampling.Fraction.parse(text));
    }

    try std.testing.expectEqual(@as(u64, 0), try sampling.parseSeed("0"));
    try std.testing.expectEqual(@as(u64, 11), try sampling.parseSeed("00011"));
    try std.testing.expectEqual(std.math.maxInt(u64), try sampling.parseSeed("18446744073709551615"));
    for ([_][]const u8{ "", "+1", "-1", "1.0", "1e2", " 1", "1 " }) |text| {
        try std.testing.expectError(error.InvalidSeed, sampling.parseSeed(text));
        try std.testing.expectError(error.InvalidCount, sampling.parseCount(text));
    }
    try std.testing.expectError(error.Overflow, sampling.parseSeed("18446744073709551616"));

    try std.testing.expectEqual(@as(u64, 0), try sampling.parseCount("0"));
    try std.testing.expectEqual(@as(u64, 11), try sampling.parseCount("00011"));
    try std.testing.expectEqual(std.math.maxInt(u64), try sampling.parseCount("18446744073709551615"));
    try std.testing.expectError(error.Overflow, sampling.parseCount("18446744073709551616"));
}

test "[unit] - [exact selector]: seed 11 matches frozen seqtk indexes" {
    var zero = sampling.ExactSelector.init(0, 11);
    defer zero.deinit(std.testing.allocator);
    var untouched = sampling.Mt19937_64.init(11);
    try zero.considerRecord(std.testing.allocator, 1);
    try std.testing.expectEqual(untouched.nextU64(), zero.generator.nextU64());

    const Expected = union(enum) {
        none,
        all,
        indexes: []const u64,
    };
    const cases = [_]struct {
        target: u64,
        expected: Expected,
    }{
        .{ .target = 0, .expected = .none },
        .{ .target = 1, .expected = .{ .indexes = &.{5} } },
        .{ .target = 2, .expected = .{ .indexes = &.{ 3, 5 } } },
        .{ .target = 3, .expected = .{ .indexes = &.{ 2, 4, 5 } } },
        .{ .target = 5, .expected = .all },
        .{ .target = 8, .expected = .all },
    };

    for (cases) |case| {
        var selector = sampling.ExactSelector.init(case.target, 11);
        defer selector.deinit(std.testing.allocator);
        var batched = sampling.ExactSelector.init(case.target, 11);
        defer batched.deinit(std.testing.allocator);
        for (1..6) |record_number| {
            try selector.considerRecord(std.testing.allocator, record_number);
        }
        try batched.considerRecordsThrough(std.testing.allocator, 2);
        try batched.considerRecordsThrough(std.testing.allocator, 5);
        const actual = selector.finish();
        const batched_actual = batched.finish();
        try std.testing.expectEqual(
            std.meta.activeTag(actual),
            std.meta.activeTag(batched_actual),
        );
        switch (case.expected) {
            .indexes => |expected| {
                const indexes = switch (actual) {
                    .indexes => |indexes| indexes,
                    else => return error.TestExpectedEqual,
                };
                const batched_indexes = batched_actual.indexes;
                try std.testing.expectEqual(expected.len, indexes.len());
                for (expected, 0..) |expected_index, index| {
                    try std.testing.expectEqual(expected_index, indexes.at(index));
                    try std.testing.expectEqual(expected_index, batched_indexes.at(index));
                }
                try std.testing.expectEqual(
                    expected.len,
                    selector.indexes.items.len,
                );
                try std.testing.expectEqual(
                    @as(usize, 0),
                    selector.index_high_words.items.len,
                );
            },
            .none => {
                try std.testing.expect(actual == .none);
                try std.testing.expectEqual(@as(usize, 0), selector.indexes.items.len);
            },
            .all => {
                try std.testing.expect(actual == .all);
                try std.testing.expectEqual(@as(usize, 0), selector.indexes.items.len);
            },
        }
    }
}

fn exerciseExactSelectorAllocations(allocator: std.mem.Allocator) !void {
    var selector = sampling.ExactSelector.init(37, 11);
    defer selector.deinit(allocator);
    for (1..101) |record_number| {
        try selector.considerRecord(allocator, record_number);
    }
    const selection = selector.finish();
    try std.testing.expect(selection == .indexes);
    try std.testing.expectEqual(@as(usize, 37), selection.indexes.len());
}

test "[failure] - [exact selector]: first index allocation failure is reported and releases state" {
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    try std.testing.expectError(error.OutOfMemory, exerciseExactSelectorAllocations(failing.allocator()));
}

test "[cli] - [sample]: fractions 0 and 1 keep fields and canonicalize LF" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const basic = try std.Io.Dir.cwd().readFileAlloc(
        std.testing.io,
        BASIC_PATH,
        allocator,
        .limited(1024 * 1024),
    );
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "0", BASIC_PATH }), 0, "", "");
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "1", BASIC_PATH }), 0, basic, "");
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "0.5", BASIC_PATH }), 0, BASIC_SELECTED_HALF, "");
    try cli.expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", "tests/data/synthetic/crlf.fastq" }),
        0,
        "@crlf_read1\nACGT\n+\nIIII\n@crlf_read2\nTGCA\n+\nJJJJ\n",
        "",
    );
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "1", "-" }), 0, "", "");
}

test "[cli] - [sample]: selecting all records canonicalizes LF and CRLF line endings" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "records");
    const cases = [_]struct {
        input: []const u8,
        canonical: []const u8,
        count: []const u8,
        gzip: bool,
    }{
        .{
            .input = "@a\nRyuN\n+first\n!5?~\n",
            .canonical = "@a\nRyuN\n+first\n!5?~\n",
            .count = "1",
            .gzip = false,
        },
        .{
            .input = "@a\r\nRyuN\r\n+first\r\n!5?~\r\n",
            .canonical = "@a\nRyuN\n+first\n!5?~\n",
            .count = "1",
            .gzip = true,
        },
    };
    for (cases) |case| {
        var gzip: std.ArrayList(u8) = .empty;
        try cli.appendGzipMember(allocator, &gzip, case.input, .{});
        try tmp.dir.writeFile(io, .{ .sub_path = "records", .data = if (case.gzip) gzip.items else case.input });
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "1", path }), 0, case.canonical, "");
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--count", case.count, path }), 0, case.canonical, "");
    }
}

test "[cli] - [sample]: a selected record with a terminal CR field fails before output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: record fields ending in CR cannot be written with LF endings\n", .{path});
    const input = "@r\r\r\nA\n+\n!\n";
    var gzip: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip, input, .{});
    for ([_][]const u8{ input, gzip.items }) |bytes| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = bytes });
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "1", path }), 1, "", diagnostic);
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--fraction", "0", path }), 0, "", "");
    }
}

test "[cli] - [sample]: a terminal CR in a pair fails after the complete earlier pair" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const pair_path = try cli.tempPath(allocator, &tmp.sub_path, "pair.fastq");
    const first1 = "@ok/1\nA\n+\n!\n";
    const first2 = "@ok/2\nT\n+\n#\n";
    const padding = try allocator.alloc(u8, OUTPUT_BUFFER_BYTES * 4);
    @memset(padding, 'x');

    const shapes = [_]usize{ 0, padding.len };
    for (shapes) |pad| {
        for (0..2) |bad_mate| {
            const mate1 = try std.fmt.allocPrint(allocator, "@pair/1 {s}\nA\n+{s}\n!\n", .{
                padding[0..pad], if (bad_mate == 0) "\r\r" else "",
            });
            const mate2 = try std.fmt.allocPrint(allocator, "@pair/2 {s}\nT\n+{s}\n#\n", .{
                padding[0..pad], if (bad_mate == 1) "\r\r" else "",
            });
            const r1 = try std.mem.concat(allocator, u8, &.{ first1, mate1 });
            const r2 = try std.mem.concat(allocator, u8, &.{ first2, mate2 });
            const pairs = try std.mem.concat(allocator, u8, &.{ first1, first2, mate1, mate2 });
            try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
            try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
            try tmp.dir.writeFile(io, .{ .sub_path = "pair.fastq", .data = pairs });
            const bad_path = if (bad_mate == 0) r1_path else r2_path;
            const paired_error = try std.fmt.allocPrint(allocator, "error: {s}: record fields ending in CR cannot be written with LF endings\n", .{bad_path});
            const interleaved_error = try std.fmt.allocPrint(allocator, "error: {s}: record fields ending in CR cannot be written with LF endings\n", .{pair_path});
            for ([_]struct { flag: []const u8, amount: []const u8 }{ .{ .flag = "--fraction", .amount = "1" }, .{ .flag = "--count", .amount = "2" } }) |mode| {
                try cli.expectResult(try cli.run(allocator, &.{ "sample", "--paired", mode.flag, mode.amount, r1_path, r2_path }), 1, first1 ++ first2, paired_error);
                try cli.expectResult(try cli.run(allocator, &.{ "sample", "--interleaved", mode.flag, mode.amount, pair_path }), 1, first1 ++ first2, interleaved_error);
            }
        }
    }
}

test "[cli] - [sample]: input failures precede unwritable fields in split gzip members" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const first = "@p/1 x\r\r\nA\n+\n!\n";
    const cases = [_]struct { second: []const u8, diagnostic: []const u8 }{
        .{ .second = "bad\nT\n+\n#\n", .diagnostic = "error: -: S003: header line must start with '@' and contain a nonempty identifier (record 1, line 1, offset 15)\n" },
        .{ .second = "@p/2\nT\n+\n \n", .diagnostic = "error: -: S006: quality byte must be ASCII 33 through 126 (record 1, line 4, offset 24)\n" },
        .{ .second = "", .diagnostic = "error: -: P002: paired input is missing a mate (pair 0, remaining R1, last R1 record 0, last R2 record none)\n" },
    };
    for (cases) |case| {
        const input = try std.mem.concat(allocator, u8, &.{ first, case.second });
        for ([_]usize{ 0, first.len, @min(first.len + 1, input.len), input.len }) |split| {
            var gzip: std.ArrayList(u8) = .empty;
            try cli.appendGzipMember(allocator, &gzip, input[0..split], .{});
            try cli.appendGzipMember(allocator, &gzip, input[split..], .{});
            try cli.expectResult(try cli.runWithStdin(allocator, &.{ "sample", "--interleaved", "--fraction", "1", "-" }, gzip.items, gzip.items.len), 1, "", case.diagnostic);
        }
    }
}

test "[cli] - [sample]: a FIFO is rejected by paired exact-count runs without waiting for a writer" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try std.testing.expectEqual(.SUCCESS, std.os.linux.errno(std.os.linux.mknodat(
        tmp.dir.handle,
        "input.fifo",
        std.os.linux.S.IFIFO | 0o600,
        0,
    )));
    try tmp.dir.symLink(io, "input.fifo", "link.fastq", .{});

    for ([_][]const u8{ "input.fifo", "link.fastq" }) |name| {
        const path = try cli.tempPath(allocator, &tmp.sub_path, name);
        const diagnostic = try std.fmt.allocPrint(
            allocator,
            "error: {s}: exact-count sampling requires a regular file\n",
            .{path},
        );
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--paired", "--count", "1", path, BASIC_PATH }), 2, "", diagnostic);
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--paired", "--count", "1", BASIC_PATH, path }), 2, "", diagnostic);
    }
}

test "[cli] - [sample]: follows regular-file symlinks in both passes" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const mate1 = "@p/1\nA\n+\n!\n";
    const mate2 = "@p/2\nT\n+\n#\n";
    const files = [_]struct { name: []const u8, link: []const u8, bytes: []const u8 }{
        .{ .name = "r1.fastq", .link = "r1.link", .bytes = mate1 },
        .{ .name = "r2.fastq", .link = "r2.link", .bytes = mate2 },
        .{ .name = "pairs.fastq", .link = "pairs.link", .bytes = mate1 ++ mate2 },
    };
    var paths: [files.len][]const u8 = undefined;
    for (files, &paths) |file, *path| {
        try tmp.dir.symLink(io, file.name, file.link, .{});
        path.* = try cli.tempPath(allocator, &tmp.sub_path, file.link);
    }
    for ([_]bool{ false, true }) |gzip| {
        for (files) |file| {
            var compressed: std.ArrayList(u8) = .empty;
            if (gzip) try cli.appendGzipMember(allocator, &compressed, file.bytes, .{});
            try tmp.dir.writeFile(io, .{
                .sub_path = file.name,
                .data = if (gzip) compressed.items else file.bytes,
            });
        }
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--count", "1", paths[0] }), 0, mate1, "");
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--paired", "--count", "1", paths[0], paths[1] }), 0, mate1 ++ mate2, "");
        try cli.expectResult(try cli.run(allocator, &.{ "sample", "--interleaved", "--count", "1", paths[2] }), 0, mate1 ++ mate2, "");
    }
}

test "[cli] - [sample]: count boundaries preserve records in input order" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const basic = try std.Io.Dir.cwd().readFileAlloc(
        std.testing.io,
        BASIC_PATH,
        allocator,
        .limited(1024 * 1024),
    );
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const empty_file = try tmp.dir.createFile(std.testing.io, "empty.fastq", .{});
    empty_file.close(std.testing.io);
    const empty_path = try cli.tempPath(allocator, &tmp.sub_path, "empty.fastq");

    const cases = [_]struct {
        count: []const u8,
        path: []const u8,
        expected: []const u8,
    }{
        .{ .count = "0", .path = BASIC_PATH, .expected = "" },
        .{ .count = "1", .path = BASIC_PATH, .expected = BASIC_EXACT_ONE },
        .{ .count = "2", .path = BASIC_PATH, .expected = BASIC_EXACT_TWO },
        .{ .count = "5", .path = BASIC_PATH, .expected = basic },
        .{ .count = "8", .path = BASIC_PATH, .expected = basic },
        .{
            .count = "2",
            .path = "tests/data/synthetic/crlf.fastq",
            .expected = "@crlf_read1\nACGT\n+\nIIII\n@crlf_read2\nTGCA\n+\nJJJJ\n",
        },
        .{ .count = "18446744073709551615", .path = EMPTY_PATH, .expected = EMPTY_RECORD },
        .{ .count = "18446744073709551615", .path = empty_path, .expected = "" },
    };
    for (cases) |case| {
        try cli.expectResult(
            try cli.run(allocator, &.{ "sample", "--count", case.count, case.path }),
            0,
            case.expected,
            "",
        );
    }
}

test "[cli] - [sample]: corrupt final gzip member keeps pairs from checked members" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const record = "@r\nAAA\n+\n!!!\n";
    const input = try allocator.alloc(u8, 40329 * record.len);
    for (0..input.len / record.len) |index| {
        @memcpy(input[index * record.len ..][0..record.len], record);
    }
    var gzip: std.ArrayList(u8) = .empty;
    var offset: usize = 0;
    while (offset < input.len) {
        const end = @min(input.len, offset + std.math.maxInt(u16));
        try cli.appendGzipMember(allocator, &gzip, input[offset..end], .{});
        offset = end;
    }
    gzip.items[cli.gzipCrcOffset(gzip.items.len)] ^= 1;
    const final_member_start = (input.len - 1) / std.math.maxInt(u16) * std.math.maxInt(u16);
    const checked_pairs_len = final_member_start / (2 * record.len) * (2 * record.len);

    const result = try cli.runWithStdin(
        allocator,
        &.{ "sample", "--interleaved", "--fraction", "1", "-" },
        gzip.items,
        gzip.items.len,
    );
    try cli.expectFailedPrefix(result, 3, input, "error: -: I/O error\n");
    try std.testing.expect(result.stdout.len >= checked_pairs_len);
}

test "[integration] - [sample]: plain and gzip file and stdin select identical records" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const basic = try std.Io.Dir.cwd().readFileAlloc(
        std.testing.io,
        BASIC_PATH,
        allocator,
        .limited(1024 * 1024),
    );
    var gzip: std.ArrayList(u8) = .empty;
    defer gzip.deinit(allocator);
    try cli.appendGzipMember(allocator, &gzip, basic, .{});

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const gzip_name = "sample.fastq.gz";
    {
        const file = try tmp.dir.createFile(std.testing.io, gzip_name, .{});
        defer file.close(std.testing.io);
        try std.Io.File.writeStreamingAll(file, std.testing.io, gzip.items);
    }
    const gzip_path = try cli.tempPath(allocator, &tmp.sub_path, gzip_name);

    const file_plain = try cli.run(allocator, &.{ "sample", "--fraction", "0.5", "--seed", "11", BASIC_PATH });
    const file_gzip = try cli.run(allocator, &.{ "sample", "--fraction", "0.5", "--seed", "11", gzip_path });
    const stdin_plain = try cli.runWithStdin(allocator, &.{ "sample", "--fraction", "0.5", "--seed", "11", "-" }, basic, 1);
    const stdin_gzip = try cli.runWithStdin(allocator, &.{ "sample", "--fraction", "0.5", "--seed", "11", "-" }, gzip.items, 1);
    for ([_]cli.CommandResult{ file_plain, file_gzip, stdin_plain, stdin_gzip }) |result| {
        try cli.expectResult(result, 0, BASIC_SELECTED_HALF, "");
    }

    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--count", "2", "--seed", "11", BASIC_PATH }), 0, BASIC_EXACT_TWO, "");
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--count", "2", "--seed", "11", gzip_path }), 0, BASIC_EXACT_TWO, "");
}

test "[cli] - [sample]: pair indexes match frozen fraction and exact-count decisions" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    var interleaved: std.ArrayList(u8) = .empty;
    for (0..16) |pair_index| {
        try appendPair(allocator, &r1, &r2, &interleaved, pair_index);
    }
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2.items });
    try tmp.dir.writeFile(io, .{
        .sub_path = "interleaved.fastq",
        .data = interleaved.items,
    });
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const interleaved_path = try cli.tempPath(
        allocator,
        &tmp.sub_path,
        "interleaved.fastq",
    );

    const all = &[_]u8{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const cases = [_]struct {
        flag: []const u8,
        amount: []const u8,
        seed: []const u8,
        indexes: []const u8,
    }{
        .{ .flag = "--fraction", .amount = "0", .seed = "0", .indexes = &.{} },
        .{ .flag = "--fraction", .amount = "0.1", .seed = "11", .indexes = &.{ 4, 11 } },
        .{ .flag = "--fraction", .amount = "0.5", .seed = "0", .indexes = &.{ 0, 2, 5, 7, 10, 13, 14, 15 } },
        .{ .flag = "--fraction", .amount = "0.999999", .seed = "0", .indexes = all },
        .{ .flag = "--fraction", .amount = "1", .seed = "0", .indexes = all },
        .{ .flag = "--count", .amount = "0", .seed = "11", .indexes = &.{} },
        .{ .flag = "--count", .amount = "3", .seed = "11", .indexes = &.{ 4, 10, 11 } },
        .{ .flag = "--count", .amount = "2", .seed = "18446744073709551615", .indexes = &.{ 8, 15 } },
        .{ .flag = "--count", .amount = "18446744073709551615", .seed = "11", .indexes = all },
    };

    var expected: std.ArrayList(u8) = .empty;
    for (cases) |case| {
        expected.clearRetainingCapacity();
        try appendSelectedPairs(allocator, &expected, case.indexes);
        try cli.expectResult(try cli.run(allocator, &.{
            "sample", "--paired", case.flag, case.amount, "--seed", case.seed, r1_path, r2_path,
        }), 0, expected.items, "");
        try cli.expectResult(try cli.run(allocator, &.{
            "sample", "--interleaved", case.flag, case.amount, "--seed", case.seed, interleaved_path,
        }), 0, expected.items, "");
    }
}

test "[integration] - [sample]: gzip and exact-name fields match plain input" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const r1 = "@same left opaque\r\nRY\r\n+left annotation\r\n!~\r\n";
    const r2 = "@same right opaque\r\nKM\r\n+right annotation\r\n#$\r\n";
    const interleaved = r1 ++ r2;
    const expected =
        "@same left opaque\nRY\n+left annotation\n!~\n" ++
        "@same right opaque\nKM\n+right annotation\n#$\n";
    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    var gzip_interleaved: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, r1, .{});
    try cli.appendGzipMember(allocator, &gzip2, r2, .{});
    try cli.appendGzipMember(allocator, &gzip_interleaved, interleaved, .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq.gz", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq.gz", .data = gzip2.items });
    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq.gz",
        .data = gzip_interleaved.items,
    });
    try tmp.dir.writeFile(io, .{ .sub_path = "empty1.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "empty2.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "empty-pairs.fastq", .data = "" });

    const gzip1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq.gz");
    const gzip2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq.gz");
    const gzip_pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "pairs.fastq.gz");
    const empty1_path = try cli.tempPath(allocator, &tmp.sub_path, "empty1.fastq");
    const empty2_path = try cli.tempPath(allocator, &tmp.sub_path, "empty2.fastq");
    const empty_pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "empty-pairs.fastq");

    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--paired", "--count", "1", "--pair-names", "exact", gzip1_path, gzip2_path,
    }), 0, expected, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--interleaved", "--count", "1", "--pair-names", "exact", gzip_pairs_path,
    }), 0, expected, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--paired", "--count", "18446744073709551615", empty1_path, empty2_path,
    }), 0, "", "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--interleaved", "--count", "18446744073709551615", empty_pairs_path,
    }), 0, "", "");
}

test "[cli] - [sample]: select-all modes agree with interleave and canonical fields" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    const r1 = "@same left opaque\r\nAC\r\n+left annotation\r\n!~\r\n";
    const r2 = "@same right opaque\r\nGT\r\n+right annotation\r\n#$\r\n";
    const interleaved = r1 ++ r2;
    const expected =
        "@same left opaque\nAC\n+left annotation\n!~\n" ++
        "@same right opaque\nGT\n+right annotation\n#$\n";
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "pairs.fastq");

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = interleaved });
    try cli.expectResult(
        try cli.run(allocator, &.{ "interleave", "--pair-names", "exact", r1_path, r2_path }),
        0,
        expected,
        "",
    );
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--paired", "--count", "1", "--pair-names", "exact", r1_path, r2_path,
    }), 0, expected, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--interleaved", "--count", "1", "--pair-names", "exact", pairs_path,
    }), 0, expected, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--paired", "--fraction", "1", "--pair-names", "exact", r1_path, r2_path,
    }), 0, expected, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--interleaved", "--fraction", "1", "--pair-names", "exact", pairs_path,
    }), 0, expected, "");

    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    var gzip_pairs: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, r1, .{});
    try cli.appendGzipMember(allocator, &gzip2, r2, .{});
    try cli.appendGzipMember(allocator, &gzip_pairs, interleaved, .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = gzip2.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = gzip_pairs.items });
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--paired", "--count", "2", "--pair-names", "exact", r1_path, r2_path,
    }), 0, expected, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--interleaved", "--count", "2", "--pair-names", "exact", pairs_path,
    }), 0, expected, "");

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = "" });
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path }), 0, "", "");
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--interleaved", "--fraction", "1", pairs_path }), 0, "", "");
}

test "[cli] - [sample]: gzip and stdin layouts preserve selected pairs" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    var interleaved: std.ArrayList(u8) = .empty;
    for (0..3) |pair_index| {
        try appendPair(allocator, &r1, &r2, &interleaved, pair_index);
    }
    var expected: std.ArrayList(u8) = .empty;
    try appendSelectedPairs(allocator, &expected, &.{ 0, 2 });

    var gzip1: std.ArrayList(u8) = .empty;
    var gzip2: std.ArrayList(u8) = .empty;
    var gzip_interleaved: std.ArrayList(u8) = .empty;
    try cli.appendGzipMember(allocator, &gzip1, r1.items, .{});
    try cli.appendGzipMember(allocator, &gzip2, r2.items, .{});
    try cli.appendGzipMember(allocator, &gzip_interleaved, interleaved.items, .{});
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq.gz", .data = gzip1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq.gz", .data = gzip2.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq.gz", .data = gzip_interleaved.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2.items });
    const gzip1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq.gz");
    const gzip2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq.gz");
    const gzip_interleaved_path = try cli.tempPath(
        allocator,
        &tmp.sub_path,
        "pairs.fastq.gz",
    );
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--paired", "--fraction", "0.5", gzip1_path, gzip2_path,
    }), 0, expected.items, "");
    try cli.expectResult(try cli.run(allocator, &.{
        "sample", "--interleaved", "--fraction", "0.5", gzip_interleaved_path,
    }), 0, expected.items, "");
    try cli.expectResult(try cli.runWithStdin(
        allocator,
        &.{ "sample", "--paired", "--fraction", "0.5", "-", r2_path },
        r1.items,
        1,
    ), 0, expected.items, "");
    try cli.expectResult(try cli.runWithStdin(
        allocator,
        &.{ "sample", "--paired", "--fraction", "0.5", r1_path, "-" },
        r2.items,
        2,
    ), 0, expected.items, "");
    try cli.expectResult(try cli.runWithStdin(
        allocator,
        &.{ "sample", "--interleaved", "--fraction", "0.5", "-" },
        gzip_interleaved.items,
        1,
    ), 0, expected.items, "");
    for ([_]usize{ 1, 7, interleaved.items.len }) |chunk_len| {
        try cli.expectResult(
            try cli.runWithStdin(
                allocator,
                &.{ "sample", "--interleaved", "--fraction", "0.5", "-" },
                interleaved.items,
                chunk_len,
            ),
            0,
            expected.items,
            "",
        );
    }

    var first_member: std.ArrayList(u8) = .empty;
    var second_member: std.ArrayList(u8) = .empty;
    var damaged_member: std.ArrayList(u8) = .empty;
    for (0..3) |pair_index| {
        try appendPair(allocator, null, null, &first_member, pair_index);
        try appendPair(allocator, null, null, &second_member, pair_index);
    }
    try appendPair(allocator, null, null, &damaged_member, 3);
    var corrupt_chain: std.ArrayList(u8) = .empty;
    for (0..4) |_| {
        try cli.appendGzipMember(allocator, &corrupt_chain, first_member.items, .{});
        try cli.appendGzipMember(allocator, &corrupt_chain, second_member.items, .{});
    }
    try cli.appendGzipMember(allocator, &corrupt_chain, damaged_member.items, .{});
    corrupt_chain.items[cli.gzipCrcOffset(corrupt_chain.items.len)] ^= 1;
    const damaged = try cli.runWithStdin(
        allocator,
        &.{ "sample", "--interleaved", "--fraction", "1", "-" },
        corrupt_chain.items,
        4096,
    );
    try std.testing.expectEqual(@as(u8, 3), damaged.exit_code);
    try std.testing.expectEqualStrings("error: -: I/O error\n", damaged.stderr);
    var intact: std.ArrayList(u8) = .empty;
    for (0..8) |_| try intact.appendSlice(allocator, first_member.items);
    try std.testing.expect(damaged.stdout.len >= intact.items.len);
    try std.testing.expectEqualStrings(intact.items, damaged.stdout[0..intact.items.len]);
    try cli.expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "0", "-" },
            damaged.stdout,
            4096,
        ),
        0,
        "",
        "",
    );
}

test "[cli] - [sample]: a sequence and quality length mismatch fails before any pair is written" {
    const input = "@a/1\nR\n+\n \n@a/2\nAAA\n+\n!\n!\n";
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    try cli.expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
            input,
            input.len,
        ),
        1,
        "",
        "error: -: S005: sequence and quality lengths differ " ++
            "(record 1, line 4, offset 22)\n",
    );
}

test "[cli] - [sample]: an identifier mismatch after one complete pair keeps that pair" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const first_pair = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@ok/1\nA\n+\n!\n@left/1\nC\n+\n$\n",
    });
    try tmp.dir.writeFile(io, .{
        .sub_path = "r2.fastq",
        .data = "@ok/2\nT\n+\n#\n@right/2\nG\n+\n%\n",
    });
    const mismatch = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P001: paired identifiers or mate markers do not match (pair 1)\n",
        .{ r1_path, r2_path },
    );
    const result = try cli.run(allocator, &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path });
    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqualStrings(first_pair, result.stdout);
    try cli.expectFirstLine(mismatch, result.stderr);
}

test "[cli] - [sample]: a missing mate after one complete pair keeps that pair" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const first_pair = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@ok/1\nA\n+\n!\n@extra/1\nC\n+\n$\n",
    });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@ok/2\nT\n+\n#\n" });
    const unequal = try std.fmt.allocPrint(
        allocator,
        "error: {s} + {s}: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 1, last R2 record 0)\n",
        .{ r1_path, r2_path },
    );
    const result = try cli.run(allocator, &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path });
    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqualStrings(first_pair, result.stdout);
    try cli.expectFirstLine(unequal, result.stderr);
}

test "[cli] - [sample]: a sequence error in mate two keeps the complete pair before it" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const first_pair = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@ok/1\nA\n+\n!\n@bad/1\nC\n+\n$\n",
    });
    try tmp.dir.writeFile(io, .{
        .sub_path = "r2.fastq",
        .data = "@ok/2\nT\n+\n#\n@bad/2\n.\n+\n%\n",
    });
    const mate2_semantic = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 1, line 2, offset 19)\n",
        .{r2_path},
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path }),
        1,
        first_pair,
        mate2_semantic,
    );
}

test "[cli] - [sample]: a plus-line error after complete pairs keeps them" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const late_structural = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@bad/1\nC\n+\n$\n@bad/2\nG\nx\n%\n";
    try cli.expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
            late_structural,
            1,
        ),
        1,
        "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n",
        "error: -: S001: plus line must start with '+' " ++
            "(record 3, line 3, offset 46)\n",
    );
}

test "[cli] - [sample]: a closed standard output after a complete pair reports an I/O error" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "pairs.fastq");
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@ok/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@ok/2\nT\n+\n#\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n" });
    const cases = [_][]const []const u8{
        &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path },
        &.{ "sample", "--paired", "--count", "1", r1_path, r2_path },
        &.{ "sample", "--interleaved", "--count", "1", pairs_path },
    };
    for (cases) |args| {
        try cli.expectResult(try cli.runWithClosedStdout(allocator, args, ""), 3, "", "error: standard output: I/O error\n");
    }
}

test "[edge] - [sample]: interleaved records spanning a read-buffer refill keep their fields" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "pairs.fastq");

    const mate_len = READ_BUFFER_BYTES + 1;
    var large: std.ArrayList(u8) = .empty;
    try large.appendSlice(allocator, "@large/1 ");
    try large.appendNTimes(allocator, 'x', mate_len);
    try large.appendSlice(allocator, "\nA\n+first\n!\n@large/2 ");
    try large.appendNTimes(allocator, 'y', mate_len);
    try large.appendSlice(allocator, "\nT\n+second\n#\n");
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = large.items });
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--interleaved", "--fraction", "1", pairs_path }), 0, large.items, "");
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--interleaved", "--count", "1", pairs_path }), 0, large.items, "");
    try cli.expectResult(try cli.run(allocator, &.{ "sample", "--interleaved", "--fraction", "0", pairs_path }), 0, "", "");
}

test "[cli] - [sample]: an identifier mismatch across a read-buffer refill fails the first pair" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "pairs.fastq");

    var input: std.ArrayList(u8) = .empty;
    try input.appendSlice(allocator, "@large/1 ");
    try input.appendNTimes(allocator, 'x', 2 * READ_BUFFER_BYTES);
    try input.appendSlice(allocator, "\nA\n+first\n!\n@other/2\nT\n+second\n#\n");
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = input.items });

    const mismatch = try std.fmt.allocPrint(
        allocator,
        "error: {s}: P001: paired identifiers or mate markers do not match (pair 0)\n",
        .{pairs_path},
    );
    const result = try cli.run(allocator, &.{ "sample", "--interleaved", "--fraction", "0", pairs_path });
    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqualStrings("", result.stdout);
    try cli.expectFirstLine(mismatch, result.stderr);
}

test "[cli] - [sample]: an alias of a paired input rejects zero-selection runs before reading" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const link = try cli.tempPath(allocator, &tmp.sub_path, "link.fastq");
    const payload = "@SRR1.1 1 length=4\nACGT\n+\nIIII\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });
    try tmp.dir.symLink(io, "input.fastq", "link.fastq", .{});

    const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{link});
    for ([_][]const u8{ "--fraction", "--count" }) |mode| {
        for ([_][]const u8{ "0", "1" }) |amount| {
            try cli.expectResult(try cli.run(allocator, &.{ "sample", "--paired", mode, amount, path, link }), 2, "", diagnostic);
        }
    }

    const stdin_diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: paired inputs refer to the same file\n", .{path});
    const file = try tmp.dir.openFile(io, "input.fastq", .{});
    defer file.close(io);
    try cli.expectResult(try cli.runWithStdinFile(allocator, &.{ "sample", "--paired", "--fraction", "0", "-", path }, file), 2, "", stdin_diagnostic);
    var first: [1]u8 = undefined;
    try std.testing.expectEqual(@as(usize, 1), try file.readStreaming(io, &.{&first}));
    try std.testing.expectEqual(payload[0], first[0]);
}

test "[cli] - [sample]: paired and interleaved argument shapes are rejected with exit 2" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct {
        args: []const []const u8,
        stderr: []const u8,
    }{
        .{
            .args = &.{ "sample", "--fraction", "1", "--pair-names", "exact", EMPTY_PATH },
            .stderr = "error: --pair-names requires --paired or --interleaved\n",
        },
        .{
            .args = &.{ "sample", "--paired", "--interleaved", "--fraction", "1", EMPTY_PATH },
            .stderr = "error: --paired and --interleaved are mutually exclusive\n",
        },
        .{
            .args = &.{ "sample", "--paired", "--fraction", "1", EMPTY_PATH },
            .stderr = "error: sample --paired requires exactly two inputs\n",
        },
        .{
            .args = &.{ "sample", "--interleaved", "--fraction", "1", EMPTY_PATH, BASIC_PATH },
            .stderr = "error: sample --interleaved requires exactly one input\n",
        },
        .{
            .args = &.{ "sample", "--paired", "--fraction", "1", "-", "-" },
            .stderr = "error: paired sample inputs may contain standard input at most once\n",
        },
        .{
            .args = &.{ "sample", "--paired", "--count", "1", "-", EMPTY_PATH },
            .stderr = "error: paired exact-count sampling requires file paths\n",
        },
        .{
            .args = &.{ "sample", "--interleaved", "--count", "1", "-" },
            .stderr = "error: paired exact-count sampling requires file paths\n",
        },
    };
    for (cases) |case| {
        try cli.expectResult(try cli.run(allocator, case.args), 2, "", case.stderr);
    }
}

test "[cli] - [sample]: a closed standard output or input reports its I/O error with exit 3" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const pair = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    try cli.expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
            pair,
        ),
        3,
        "",
        "error: standard output: I/O error\n",
    );
    try cli.expectResult(
        try cli.runWithClosedStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
        ),
        3,
        "",
        "error: -: I/O error\n",
    );
}

test "[cli] - [sample]: every record is validated before selection in fraction and count modes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct {
        args: []const []const u8,
        stderr: []const u8,
    }{
        .{
            .args = &.{ "sample", "--fraction", "0", "tests/data/synthetic/bad_alphabet.fastq" },
            .stderr = "error: tests/data/synthetic/bad_alphabet.fastq: S002: " ++
                "sequence byte is outside the selected alphabet (record 0, line 2, offset 16)\n",
        },
        .{
            .args = &.{ "sample", "--count", "0", "tests/data/synthetic/bad_plus.fastq" },
            .stderr = "error: tests/data/synthetic/bad_plus.fastq: S001: " ++
                "plus line must start with '+' (record 0, line 3, offset 15)\n",
        },
        .{
            .args = &.{ "sample", "--fraction", "1", "--alphabet", "acgtn", "tests/data/synthetic/iupac_valid.fastq" },
            .stderr = "error: tests/data/synthetic/iupac_valid.fastq: S002: " ++
                "sequence byte is outside the selected alphabet (record 0, line 2, offset 11)\n",
        },
    };
    for (cases) |case| {
        try cli.expectResult(try cli.run(allocator, case.args), 1, "", case.stderr);
    }
}

test "[cli] - [sample]: a later malformed record preserves output below and above 64 KiB" {
    for ([_]usize{ 1, 8192 }) |records| {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        defer arena.deinit();
        const allocator = arena.allocator();
        var input: std.ArrayList(u8) = .empty;
        for (0..records) |_| try input.appendSlice(allocator, "@ok\nA\n+\n!\n");
        const prefix_len = input.items.len;
        try input.appendSlice(allocator, "@bad\nA\nx\n!\n");
        const expected_stderr = try std.fmt.allocPrint(
            allocator,
            "error: -: S001: plus line must start with '+' " ++
                "(record {d}, line 3, offset {d})\n",
            .{ records, prefix_len + 7 },
        );

        try cli.expectResult(
            try cli.runWithStdin(
                allocator,
                &.{ "sample", "--fraction", "1", "-" },
                input.items,
                if (records == 1) 2 else 4096,
            ),
            1,
            input.items[0..prefix_len],
            expected_stderr,
        );
    }
}

test "[cli] - [sample]: a later structural failure keeps complete pairs beyond the output buffer" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1_path = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const pairs_path = try cli.tempPath(allocator, &tmp.sub_path, "pairs.fastq");
    var r1: std.ArrayList(u8) = .empty;
    var r2: std.ArrayList(u8) = .empty;
    var pairs: std.ArrayList(u8) = .empty;
    var pair_count: usize = 0;
    while (pairs.items.len <= OUTPUT_BUFFER_BYTES) : (pair_count += 1) {
        try appendPair(allocator, &r1, &r2, &pairs, pair_count);
    }
    const r2_prefix_len = r2.items.len;
    const prefix_len = pairs.items.len;
    const tail1 = "@bad/1\nA\n+\n!\n";
    const tail2 = "@bad/2\n.\n+\n#\n";
    try r1.appendSlice(allocator, tail1);
    try r2.appendSlice(allocator, tail2);
    try pairs.appendSlice(allocator, tail1 ++ tail2);
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2.items });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = pairs.items });

    for ([_]bool{ false, true }) |interleaved| {
        const record_index = if (interleaved) 2 * pair_count + 1 else pair_count;
        const offset = (if (interleaved) prefix_len + tail1.len else r2_prefix_len) + 7;
        const expected_stderr = try std.fmt.allocPrint(
            allocator,
            "error: {s}: S002: sequence byte is outside the selected alphabet " ++
                "(record {d}, line 2, offset {d})\n",
            .{ if (interleaved) pairs_path else r2_path, record_index, offset },
        );
        for ([_]bool{ false, true }) |exact| {
            const selection = if (exact) "--count" else "--fraction";
            const amount = if (exact) try std.fmt.allocPrint(allocator, "{d}", .{pair_count + 1}) else "1";
            const args: []const []const u8 = if (interleaved)
                &.{ "sample", "--interleaved", selection, amount, pairs_path }
            else
                &.{ "sample", "--paired", selection, amount, r1_path, r2_path };
            try cli.expectResult(
                try cli.run(allocator, args),
                1,
                if (exact) "" else pairs.items[0..prefix_len],
                expected_stderr,
            );
        }
    }
}

test "[cli] - [sample]: piped input and output make progress together" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var input: std.ArrayList(u8) = .empty;
    for (0..8192) |_| try input.appendSlice(allocator, "@r\nACGT\n+\n!!!!\n");
    try cli.expectResult(try cli.runWithStdin(allocator, &.{ "sample", "--fraction", "1", "-" }, input.items, 4096), 0, input.items, "");
}

test "[cli] - [sample]: an invalid fraction and an overflowing count exit with their classes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try cli.expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1.1", EMPTY_PATH }),
        2,
        "",
        "error: invalid --fraction value\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "sample", "--count", "18446744073709551616", EMPTY_PATH }),
        4,
        "",
        "error: --count exceeds supported limit\n",
    );
}

test "[cli] - [sample]: invocation shape is validated before input" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const cases = [_]struct {
        args: []const []const u8,
        stderr: []const u8,
    }{
        .{ .args = &.{ "sample", EMPTY_PATH }, .stderr = "error: sample requires --fraction P or --count K\n" },
        .{ .args = &.{ "sample", "--count", "1" }, .stderr = "error: sample requires exactly one input\n" },
        .{ .args = &.{ "sample", "--fraction", "0", EMPTY_PATH, BASIC_PATH }, .stderr = "error: sample requires exactly one input\n" },
        .{
            .args = &.{ "sample", "--fraction", "0.5", "--count", "1", BASIC_PATH },
            .stderr = "error: --fraction and --count are mutually exclusive\n",
        },
        .{
            .args = &.{ "sample", "--count", "1", "-" },
            .stderr = "error: exact-count sampling requires a file path\n",
        },
    };
    for (cases) |case| {
        try cli.expectResult(try cli.run(allocator, case.args), 2, "", case.stderr);
    }

    try cli.expectResult(
        try cli.runWith(allocator, &.{ "sample", "--count", "1", "-" }, .{ .stdin = .stalled }),
        2,
        "",
        "error: exact-count sampling requires a file path\n",
    );
}

test "[cli] - [sample]: the largest line limits accept valid input in both modes" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input.fastq");
    const valid = "@p/1\nA\n+\n!\n@p/2\nT\n+\n!\n";
    const invalid = "@p/1\nA\n+\n!\n@p/2\nX\n+\n!\n";
    const lint_error = try std.fmt.allocPrint(allocator, "error: {s}: S002: sequence byte is outside the selected alphabet (record 1, line 2, offset 16)\n", .{path});

    for ([_][]const u8{ valid, invalid }) |payload| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = payload });
        const rejects_input = std.mem.eql(u8, payload, invalid);
        for ([_][]const u8{ "--count", "--fraction" }) |mode| {
            for ([_][]const u8{ "4611686018427387903", "18446744073709551615" }) |limit| {
                try cli.expectResult(
                    try cli.run(allocator, &.{ "sample", "--interleaved", mode, "1", "--max-line-bytes", limit, path }),
                    if (rejects_input) 1 else 0,
                    if (rejects_input) "" else payload,
                    if (rejects_input) lint_error else "",
                );
                try std.testing.expectEqualStrings(payload, try tmp.dir.readFileAlloc(io, "input.fastq", allocator, .limited(1024)));
            }
        }
    }
}

test "[cli] - [sample]: a line over the limit exits 4 and a missing input exits 3" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try cli.expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", "--max-line-bytes", "3", BASIC_PATH }),
        4,
        "",
        "error: tests/data/synthetic/basic_valid.fastq: line length limit exceeded\n",
    );
    try cli.expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", "missing.fastq" }),
        3,
        "",
        "error: missing.fastq: file not found\n",
    );
}

test "[cli] - [sample]: closed standard output reports its I/O error with exit 3" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try cli.expectResult(
        try cli.runWithClosedStdout(allocator, &.{ "sample", "--fraction", "1", BASIC_PATH }, ""),
        3,
        "",
        "error: standard output: I/O error\n",
    );
    try cli.expectResult(
        try cli.runWithClosedStdout(allocator, &.{ "sample", "--fraction", "1", "-" }, "@ok\nA\n+\n!\n@bad\nA\nx\n!\n"),
        3,
        "",
        "error: -: S001: plus line must start with '+' " ++
            "(record 1, line 3, offset 17)\n" ++
            "error: standard output: I/O error\n",
    );
    try cli.expectResult(
        try cli.runWithClosedStdout(allocator, &.{ "sample", "--count", "1", BASIC_PATH }, ""),
        3,
        "",
        "error: standard output: I/O error\n",
    );
}
