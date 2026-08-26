//! Private RNG and installed CLI contracts for `z-fastq sample`.

const std = @import("std");
const sampling = @import("sample_internal");
const cli = @import("utilities.zig");

const BASIC_PATH = "tests/data/synthetic/basic_valid.fastq";
const EMPTY_PATH = "tests/data/synthetic/empty_valid.fastq";
const EMPTY_RECORD = "@empty\n\n+\n\n";
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
    r1: *std.ArrayList(u8),
    r2: *std.ArrayList(u8),
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
    try r1.appendSlice(allocator, record1);
    try r2.appendSlice(allocator, record2);
    try interleaved.appendSlice(allocator, record1);
    try interleaved.appendSlice(allocator, record2);
}

fn appendSelectedPairs(
    allocator: std.mem.Allocator,
    output: *std.ArrayList(u8),
    indexes: []const u8,
) !void {
    var unused_r1: std.ArrayList(u8) = .empty;
    defer unused_r1.deinit(allocator);
    var unused_r2: std.ArrayList(u8) = .empty;
    defer unused_r2.deinit(allocator);
    for (indexes) |pair_index| {
        try appendPair(allocator, &unused_r1, &unused_r2, output, pair_index);
    }
}

fn appendExactSelectedPairs(
    allocator: std.mem.Allocator,
    output: *std.ArrayList(u8),
    pair_count: u64,
    target: u64,
    seed: u64,
) !void {
    var selector = sampling.ExactSelector.init(target, seed);
    defer selector.deinit(allocator);
    var pair_number: u64 = 1;
    while (pair_number <= pair_count) : (pair_number += 1) {
        try selector.considerRecord(allocator, pair_number);
    }
    switch (selector.finish()) {
        .none => {},
        .all => for (0..pair_count) |pair_index| {
            try appendSelectedPairs(allocator, output, &.{@intCast(pair_index)});
        },
        .indexes => |indexes| for (indexes) |selected| {
            try appendSelectedPairs(allocator, output, &.{@intCast(selected - 1)});
        },
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

test "[edge] - [fraction selector]: equality is excluded by the strict boundary" {
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
    try std.testing.expectEqual(std.math.maxInt(u64), try sampling.parseSeed(
        "18446744073709551615",
    ));
    for ([_][]const u8{ "", "+1", "-1", "1.0", "1e2", " 1", "1 " }) |text| {
        try std.testing.expectError(error.InvalidSeed, sampling.parseSeed(text));
    }
    try std.testing.expectError(
        error.Overflow,
        sampling.parseSeed("18446744073709551616"),
    );

    try std.testing.expectEqual(@as(u64, 0), try sampling.parseCount("0"));
    try std.testing.expectEqual(@as(u64, 11), try sampling.parseCount("00011"));
    try std.testing.expectEqual(std.math.maxInt(u64), try sampling.parseCount(
        "18446744073709551615",
    ));
    for ([_][]const u8{ "", "+1", "-1", "1.0", "1e2", " 1", "1 " }) |text| {
        try std.testing.expectError(error.InvalidCount, sampling.parseCount(text));
    }
    try std.testing.expectError(
        error.Overflow,
        sampling.parseCount("18446744073709551616"),
    );
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
        for (1..6) |record_number| {
            try selector.considerRecord(std.testing.allocator, record_number);
        }
        const actual = selector.finish();
        switch (case.expected) {
            .indexes => |expected| {
                const indexes = switch (actual) {
                    .indexes => |indexes| indexes,
                    else => return error.TestExpectedEqual,
                };
                try std.testing.expectEqualSlices(u64, expected, indexes);
                try std.testing.expectEqual(
                    selector.indexes.items.len,
                    selector.indexes.capacity,
                );
            },
            .none => {
                try std.testing.expect(actual == .none);
                try std.testing.expectEqual(@as(usize, 0), selector.indexes.capacity);
            },
            .all => {
                try std.testing.expect(actual == .all);
                try std.testing.expectEqual(@as(usize, 0), selector.indexes.capacity);
            },
        }
    }
}

test "[unit] - [exact selector]: indexes materialize only after the target is exceeded" {
    var selector = sampling.ExactSelector.init(5, 11);
    defer selector.deinit(std.testing.allocator);

    for (1..6) |record_number| {
        try selector.considerRecord(std.testing.allocator, record_number);
        try std.testing.expectEqual(@as(usize, 0), selector.indexes.capacity);
    }
    try selector.considerRecord(std.testing.allocator, 6);
    try std.testing.expectEqual(@as(usize, 5), selector.indexes.items.len);
    try std.testing.expectEqual(@as(usize, 5), selector.indexes.capacity);
    try std.testing.expect(selector.finish() == .indexes);
}

fn exerciseExactSelectorAllocations(allocator: std.mem.Allocator) !void {
    var selector = sampling.ExactSelector.init(37, 11);
    defer selector.deinit(allocator);
    for (1..101) |record_number| {
        try selector.considerRecord(allocator, record_number);
    }
    const selection = selector.finish();
    try std.testing.expect(selection == .indexes);
    try std.testing.expectEqual(@as(usize, 37), selection.indexes.len);
}

test "[failure] - [exact selector]: materialized index allocation is released" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseExactSelectorAllocations,
        .{},
    );
}

test "[cli] - [sample]: boundary fractions preserve fields and canonicalize LF" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const basic = try std.Io.Dir.cwd().readFileAlloc(
        std.testing.io,
        BASIC_PATH,
        allocator,
        .limited(1024 * 1024),
    );
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "0", BASIC_PATH }),
        0,
        "",
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", BASIC_PATH }),
        0,
        basic,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "0.5", BASIC_PATH }),
        0,
        BASIC_SELECTED_HALF,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", "tests/data/synthetic/crlf.fastq" }),
        0,
        "@crlf_read1\nACGT\n+\nIIII\n@crlf_read2\nTGCA\n+\nJJJJ\n",
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", "-" }),
        0,
        "",
        "",
    );
}

test "[cli] - [exact sample]: count boundaries preserve records in input order" {
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
    const empty_path = try std.fmt.allocPrint(
        allocator,
        ".zig-cache/tmp/{s}/empty.fastq",
        .{tmp.sub_path},
    );

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
        .{ .count = "18446744073709551615", .path = EMPTY_PATH, .expected = EMPTY_RECORD },
        .{ .count = "18446744073709551615", .path = empty_path, .expected = "" },
    };
    for (cases) |case| {
        try expectResult(
            try cli.run(allocator, &.{ "sample", "--count", case.count, case.path }),
            0,
            case.expected,
            "",
        );
    }
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
    const gzip_path = try std.fmt.allocPrint(
        allocator,
        ".zig-cache/tmp/{s}/{s}",
        .{ tmp.sub_path, gzip_name },
    );

    const file_plain = try cli.run(
        allocator,
        &.{ "sample", "--fraction", "0.5", "--seed", "11", BASIC_PATH },
    );
    const file_gzip = try cli.run(
        allocator,
        &.{ "sample", "--fraction", "0.5", "--seed", "11", gzip_path },
    );
    const stdin_plain = try cli.runWithStdin(
        allocator,
        &.{ "sample", "--fraction", "0.5", "--seed", "11", "-" },
        basic,
        1,
    );
    const stdin_gzip = try cli.runWithStdin(
        allocator,
        &.{ "sample", "--fraction", "0.5", "--seed", "11", "-" },
        gzip.items,
        1,
    );
    for ([_]cli.CommandResult{ file_plain, file_gzip, stdin_plain, stdin_gzip }) |result| {
        try expectResult(result, 0, BASIC_SELECTED_HALF, "");
    }

    try expectResult(
        try cli.run(
            allocator,
            &.{ "sample", "--count", "2", "--seed", "11", BASIC_PATH },
        ),
        0,
        BASIC_EXACT_TWO,
        "",
    );
    try expectResult(
        try cli.run(
            allocator,
            &.{ "sample", "--count", "2", "--seed", "11", gzip_path },
        ),
        0,
        BASIC_EXACT_TWO,
        "",
    );
}

test "[cli] - [paired fraction sample]: frozen pair indexes match logical record decisions" {
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
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const interleaved_path = try tempPath(
        allocator,
        &tmp.sub_path,
        "interleaved.fastq",
    );

    const all = &[_]u8{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const Case = struct {
        seed: []const u8,
        fraction: []const u8,
        indexes: []const u8,
    };
    const cases = [_]Case{
        .{ .seed = "0", .fraction = "0", .indexes = &.{} },
        .{ .seed = "0", .fraction = "0.000001", .indexes = &.{} },
        .{ .seed = "0", .fraction = "0.1", .indexes = &.{ 2, 5 } },
        .{ .seed = "0", .fraction = "0.5", .indexes = &.{ 0, 2, 5, 7, 10, 13, 14, 15 } },
        .{ .seed = "0", .fraction = "0.999999", .indexes = all },
        .{ .seed = "0", .fraction = "1", .indexes = all },
        .{ .seed = "1", .fraction = "0", .indexes = &.{} },
        .{ .seed = "1", .fraction = "0.000001", .indexes = &.{} },
        .{ .seed = "1", .fraction = "0.1", .indexes = &.{ 3, 7, 10 } },
        .{ .seed = "1", .fraction = "0.5", .indexes = &.{ 0, 1, 2, 3, 4, 6, 7, 10, 13, 14, 15 } },
        .{ .seed = "1", .fraction = "0.999999", .indexes = all },
        .{ .seed = "1", .fraction = "1", .indexes = all },
        .{ .seed = "11", .fraction = "0", .indexes = &.{} },
        .{ .seed = "11", .fraction = "0.000001", .indexes = &.{} },
        .{ .seed = "11", .fraction = "0.1", .indexes = &.{ 4, 11 } },
        .{ .seed = "11", .fraction = "0.5", .indexes = &.{ 0, 2, 4, 5, 10, 11, 13, 14, 15 } },
        .{ .seed = "11", .fraction = "0.999999", .indexes = all },
        .{ .seed = "11", .fraction = "1", .indexes = all },
        .{ .seed = "4294967295", .fraction = "0", .indexes = &.{} },
        .{ .seed = "4294967295", .fraction = "0.000001", .indexes = &.{} },
        .{ .seed = "4294967295", .fraction = "0.1", .indexes = &.{ 7, 15 } },
        .{ .seed = "4294967295", .fraction = "0.5", .indexes = &.{ 0, 4, 5, 6, 7, 8, 9, 11, 15 } },
        .{ .seed = "4294967295", .fraction = "0.999999", .indexes = all },
        .{ .seed = "4294967295", .fraction = "1", .indexes = all },
        .{ .seed = "18446744073709551615", .fraction = "0", .indexes = &.{} },
        .{ .seed = "18446744073709551615", .fraction = "0.000001", .indexes = &.{} },
        .{ .seed = "18446744073709551615", .fraction = "0.1", .indexes = &.{ 0, 2, 8, 15 } },
        .{ .seed = "18446744073709551615", .fraction = "0.5", .indexes = &.{ 0, 2, 6, 7, 8, 9, 12, 14, 15 } },
        .{ .seed = "18446744073709551615", .fraction = "0.999999", .indexes = all },
        .{ .seed = "18446744073709551615", .fraction = "1", .indexes = all },
    };

    var expected: std.ArrayList(u8) = .empty;
    for (cases) |case| {
        expected.clearRetainingCapacity();
        try appendSelectedPairs(allocator, &expected, case.indexes);
        try expectResult(
            try cli.run(allocator, &.{
                "sample",
                "--paired",
                "--fraction",
                case.fraction,
                "--seed",
                case.seed,
                r1_path,
                r2_path,
            }),
            0,
            expected.items,
            "",
        );
        try expectResult(
            try cli.run(allocator, &.{
                "sample",
                "--interleaved",
                "--fraction",
                case.fraction,
                "--seed",
                case.seed,
                interleaved_path,
            }),
            0,
            expected.items,
            "",
        );
    }
}

test "[cli] - [paired exact sample]: pair indexes match the exact selector" {
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
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const interleaved_path = try tempPath(
        allocator,
        &tmp.sub_path,
        "interleaved.fastq",
    );

    const cases = [_]struct {
        count: []const u8,
        seed: []const u8,
        target: u64,
        seed_value: u64,
    }{
        .{ .count = "0", .seed = "11", .target = 0, .seed_value = 11 },
        .{ .count = "1", .seed = "11", .target = 1, .seed_value = 11 },
        .{ .count = "3", .seed = "11", .target = 3, .seed_value = 11 },
        .{ .count = "16", .seed = "11", .target = 16, .seed_value = 11 },
        .{ .count = "20", .seed = "11", .target = 20, .seed_value = 11 },
        .{
            .count = "2",
            .seed = "18446744073709551615",
            .target = 2,
            .seed_value = std.math.maxInt(u64),
        },
        .{
            .count = "18446744073709551615",
            .seed = "11",
            .target = std.math.maxInt(u64),
            .seed_value = 11,
        },
    };

    var expected: std.ArrayList(u8) = .empty;
    for (cases) |case| {
        expected.clearRetainingCapacity();
        try appendExactSelectedPairs(
            allocator,
            &expected,
            16,
            case.target,
            case.seed_value,
        );
        try expectResult(
            try cli.run(allocator, &.{
                "sample",
                "--paired",
                "--count",
                case.count,
                "--seed",
                case.seed,
                r1_path,
                r2_path,
            }),
            0,
            expected.items,
            "",
        );
        try expectResult(
            try cli.run(allocator, &.{
                "sample",
                "--interleaved",
                "--count",
                case.count,
                "--seed",
                case.seed,
                interleaved_path,
            }),
            0,
            expected.items,
            "",
        );
    }
}

test "[integration] - [paired exact sample]: gzip and exact-name fields match plain input" {
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

    const gzip1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq.gz");
    const gzip2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq.gz");
    const gzip_pairs_path = try tempPath(allocator, &tmp.sub_path, "pairs.fastq.gz");
    const empty1_path = try tempPath(allocator, &tmp.sub_path, "empty1.fastq");
    const empty2_path = try tempPath(allocator, &tmp.sub_path, "empty2.fastq");
    const empty_pairs_path = try tempPath(allocator, &tmp.sub_path, "empty-pairs.fastq");

    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--count",
            "1",
            "--pair-names",
            "exact",
            gzip1_path,
            gzip2_path,
        }),
        0,
        expected,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "1",
            "--pair-names",
            "exact",
            gzip_pairs_path,
        }),
        0,
        expected,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--count",
            "18446744073709551615",
            empty1_path,
            empty2_path,
        }),
        0,
        "",
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "18446744073709551615",
            empty_pairs_path,
        }),
        0,
        "",
        "",
    );
}

test "[cli] - [paired fraction sample]: fields, exact names, empty input, and LF output are preserved" {
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
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = r2 });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = interleaved });
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const pairs_path = try tempPath(allocator, &tmp.sub_path, "pairs.fastq");

    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--fraction",
            "1",
            "--pair-names",
            "exact",
            r1_path,
            r2_path,
        }),
        0,
        expected,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "1",
            "--pair-names",
            "exact",
            pairs_path,
        }),
        0,
        expected,
        "",
    );

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "" });
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = "" });
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path }),
        0,
        "",
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "1",
            pairs_path,
        }),
        0,
        "",
        "",
    );
}

test "[integration] - [paired fraction sample]: gzip, stdin layouts, and source chunks agree" {
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
    const gzip1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq.gz");
    const gzip2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq.gz");
    const gzip_interleaved_path = try tempPath(
        allocator,
        &tmp.sub_path,
        "pairs.fastq.gz",
    );
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");

    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--fraction",
            "0.5",
            gzip1_path,
            gzip2_path,
        }),
        0,
        expected.items,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "0.5",
            gzip_interleaved_path,
        }),
        0,
        expected.items,
        "",
    );
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--paired", "--fraction", "0.5", "-", r2_path },
            r1.items,
            1,
        ),
        0,
        expected.items,
        "",
    );
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--paired", "--fraction", "0.5", r1_path, "-" },
            r2.items,
            2,
        ),
        0,
        expected.items,
        "",
    );
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "0.5", "-" },
            gzip_interleaved.items,
            1,
        ),
        0,
        expected.items,
        "",
    );
    for (1..interleaved.items.len + 1) |chunk_len| {
        try expectResult(
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
    var discard1: std.ArrayList(u8) = .empty;
    var discard2: std.ArrayList(u8) = .empty;
    for (0..500) |pair_index| {
        try appendPair(allocator, &discard1, &discard2, &first_member, pair_index);
        try appendPair(allocator, &discard1, &discard2, &second_member, pair_index);
    }
    try appendPair(allocator, &discard1, &discard2, &damaged_member, 500);
    var corrupt_chain: std.ArrayList(u8) = .empty;
    for (0..4) |_| {
        try cli.appendGzipMember(allocator, &corrupt_chain, first_member.items, .{});
        try cli.appendGzipMember(allocator, &corrupt_chain, second_member.items, .{});
    }
    try cli.appendGzipMember(allocator, &corrupt_chain, damaged_member.items, .{});
    corrupt_chain.items[corrupt_chain.items.len - 8] ^= 1;
    const damaged = try cli.runWithStdin(
        allocator,
        &.{ "sample", "--interleaved", "--fraction", "1", "-" },
        corrupt_chain.items,
        4096,
    );
    try std.testing.expectEqual(@as(u8, 3), damaged.exit_code);
    try std.testing.expectEqualStrings("error: -: I/O error\n", damaged.stderr);
    try std.testing.expect(damaged.stdout.len != 0);
    try std.testing.expect(std.mem.startsWith(u8, damaged.stdout, "@pair0/1 lane=left\n"));
    try expectResult(
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

test "[cli] - [paired sample]: failing pairs leave only earlier complete pairs" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const r1_path = try tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2_path = try tempPath(allocator, &tmp.sub_path, "r2.fastq");
    const pairs_path = try tempPath(allocator, &tmp.sub_path, "pairs.fastq");
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
        "error: {s} + {s}: P001: paired identifiers or mate markers do not match " ++
            "(pair 1)\n" ++
            "  R1: input={s}, record=1, offset=12, first_token=left/1 " ++
            "[length=6, truncated=false], normalized_id=left " ++
            "[length=4, truncated=false], mate_markers=1\n" ++
            "  R2: input={s}, record=1, offset=12, first_token=right/2 " ++
            "[length=7, truncated=false], normalized_id=right " ++
            "[length=5, truncated=false], mate_markers=2\n",
        .{ r1_path, r2_path, r1_path, r2_path },
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--fraction",
            "1",
            r1_path,
            r2_path,
        }),
        1,
        first_pair,
        mismatch,
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--count",
            "1",
            r1_path,
            r2_path,
        }),
        1,
        "",
        mismatch,
    );
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "0", "-" },
            "@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n",
            1,
        ),
        1,
        "",
        "error: -: P001: paired identifiers or mate markers do not match (pair 0)\n" ++
            "  R1: input=-, record=0, offset=0, first_token=left/1 " ++
            "[length=6, truncated=false], normalized_id=left " ++
            "[length=4, truncated=false], mate_markers=1\n" ++
            "  R2: input=-, record=1, offset=14, first_token=right/2 " ++
            "[length=7, truncated=false], normalized_id=right " ++
            "[length=5, truncated=false], mate_markers=2\n",
    );

    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = "@left/1\nA\n+\n!\n@right/2\nT\n+\n#\n",
    });
    const buffered_mismatch = try std.fmt.allocPrint(
        allocator,
        "error: {s}: P001: paired identifiers or mate markers do not match (pair 0)\n" ++
            "  R1: input={s}, record=0, offset=0, first_token=left/1 " ++
            "[length=6, truncated=false], normalized_id=left " ++
            "[length=4, truncated=false], mate_markers=1\n" ++
            "  R2: input={s}, record=1, offset=14, first_token=right/2 " ++
            "[length=7, truncated=false], normalized_id=right " ++
            "[length=5, truncated=false], mate_markers=2\n",
        .{ pairs_path, pairs_path, pairs_path },
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "0",
            pairs_path,
        }),
        1,
        "",
        buffered_mismatch,
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "1",
            pairs_path,
        }),
        1,
        "",
        buffered_mismatch,
    );

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
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--fraction",
            "1",
            r1_path,
            r2_path,
        }),
        1,
        first_pair,
        mate2_semantic,
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--count",
            "1",
            r1_path,
            r2_path,
        }),
        1,
        "",
        mate2_semantic,
    );
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "0", "-" },
            "@bad/1\nA\n+\n!\n@bad/2\n.\n+\n!\n",
            1,
        ),
        1,
        "",
        "error: -: S002: sequence byte is outside the selected alphabet " ++
            "(record 1, line 2, offset 20)\n",
    );

    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = "@bad/1\n.\n+\n!\n@bad/2\nA\nx\n!\n",
    });
    const buffered_precedence = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S001: plus line must start with '+' " ++
            "(record 1, line 3, offset 22)\n",
        .{pairs_path},
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "0",
            pairs_path,
        }),
        1,
        "",
        buffered_precedence,
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "1",
            pairs_path,
        }),
        1,
        "",
        buffered_precedence,
    );

    try tmp.dir.writeFile(io, .{
        .sub_path = "r1.fastq",
        .data = "@bad/1\n.\n+\n!\n",
    });
    try tmp.dir.writeFile(io, .{
        .sub_path = "r2.fastq",
        .data = "@bad/2\nA\n+\n!\n",
    });
    const mate1_semantic = try std.fmt.allocPrint(
        allocator,
        "error: {s}: S002: sequence byte is outside the selected alphabet " ++
            "(record 0, line 2, offset 7)\n",
        .{r1_path},
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--fraction",
            "0",
            r1_path,
            r2_path,
        }),
        1,
        "",
        mate1_semantic,
    );

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
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--fraction",
            "1",
            r1_path,
            r2_path,
        }),
        1,
        first_pair,
        unequal,
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--paired",
            "--count",
            "1",
            r1_path,
            r2_path,
        }),
        1,
        "",
        unequal,
    );

    const late_structural = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@bad/1\nC\n+\n$\n@bad/2\nG\nx\n%\n";
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
            late_structural,
            1,
        ),
        1,
        first_pair,
        "error: -: S001: plus line must start with '+' " ++
            "(record 3, line 3, offset 46)\n",
    );
    try expectResult(
        try cli.runWithStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "0", "-" },
            "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@odd/1\nC\n+\n$\n",
            1,
        ),
        1,
        "",
        "error: -: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 2, last R2 record 1)\n",
    );
    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n@odd/1\nC\n+\n$\n",
    });
    const odd = try std.fmt.allocPrint(
        allocator,
        "error: {s}: P002: paired input is missing a mate " ++
            "(pair 1, remaining R1, last R1 record 2, last R2 record 1)\n",
        .{pairs_path},
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "1",
            pairs_path,
        }),
        1,
        "",
        odd,
    );

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@ok/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@ok/2\nT\n+\n#\n" });
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--paired", "--fraction", "1", r1_path, r2_path },
            "",
        ),
        3,
        "",
        "",
    );
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--paired", "--count", "1", r1_path, r2_path },
            "",
        ),
        3,
        "",
        "",
    );
    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n",
    });
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--interleaved", "--count", "1", pairs_path },
            "",
        ),
        3,
        "",
        "",
    );
}

test "[edge] - [paired sample]: refill-spanning 512 KiB mates preserve output" {
    const io = std.testing.io;
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const pairs_path = try tempPath(allocator, &tmp.sub_path, "pairs.fastq");

    var first_large: std.ArrayList(u8) = .empty;
    try first_large.appendSlice(allocator, "@large/1 ");
    try first_large.appendNTimes(allocator, 'x', 512 * 1024);
    try first_large.appendSlice(allocator, "\nA\n+first\n!\n@large/2 ");
    try first_large.appendNTimes(allocator, 'y', 70 * 1024);
    try first_large.appendSlice(allocator, "\nT\n+second\n#\n");
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = first_large.items });
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "1",
            pairs_path,
        }),
        0,
        first_large.items,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "1",
            pairs_path,
        }),
        0,
        first_large.items,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "0",
            pairs_path,
        }),
        0,
        "",
        "",
    );

    var second_large: std.ArrayList(u8) = .empty;
    try second_large.appendSlice(allocator, "@large/1\nA\n+first\n!\n@large/2 ");
    try second_large.appendNTimes(allocator, 'z', 512 * 1024);
    try second_large.appendSlice(allocator, "\nT\n+second\n#\n");
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = second_large.items });
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--fraction",
            "1",
            pairs_path,
        }),
        0,
        second_large.items,
        "",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--interleaved",
            "--count",
            "1",
            pairs_path,
        }),
        0,
        second_large.items,
        "",
    );
}

test "[cli] - [paired sample]: arguments and output failures retain their classes" {
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
        try expectResult(try cli.run(allocator, case.args), 2, "", case.stderr);
    }

    const pair = "@ok/1\nA\n+\n!\n@ok/2\nT\n+\n#\n";
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
            pair,
        ),
        3,
        "",
        "",
    );
    try expectResult(
        try cli.runWithClosedStdin(
            allocator,
            &.{ "sample", "--interleaved", "--fraction", "1", "-" },
        ),
        3,
        "",
        "error: -: I/O error\n",
    );
}

test "[cli] - [sample]: every record is validated before selection" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "0",
            "tests/data/synthetic/bad_alphabet.fastq",
        }),
        1,
        "",
        "error: tests/data/synthetic/bad_alphabet.fastq: S002: " ++
            "sequence byte is outside the selected alphabet (record 0, line 2, offset 16)\n",
    );
    const exact_cases = [_]struct {
        path: []const u8,
        stderr: []const u8,
    }{
        .{
            .path = "tests/data/synthetic/bad_plus.fastq",
            .stderr = "error: tests/data/synthetic/bad_plus.fastq: S001: " ++
                "plus line must start with '+' (record 0, line 3, offset 15)\n",
        },
        .{
            .path = "tests/data/synthetic/bad_alphabet.fastq",
            .stderr = "error: tests/data/synthetic/bad_alphabet.fastq: S002: " ++
                "sequence byte is outside the selected alphabet (record 0, line 2, offset 16)\n",
        },
        .{
            .path = "tests/data/synthetic/bad_header.fastq",
            .stderr = "error: tests/data/synthetic/bad_header.fastq: S003: " ++
                "header line must start with '@' and contain a nonempty identifier " ++
                "(record 0, line 1, offset 0)\n",
        },
        .{
            .path = "tests/data/synthetic/truncated_record.fastq",
            .stderr = "error: tests/data/synthetic/truncated_record.fastq: S004: " ++
                "unexpected end of file in quality line (record 0, line 4, offset 17)\n",
        },
        .{
            .path = "tests/data/synthetic/bad_qual_length.fastq",
            .stderr = "error: tests/data/synthetic/bad_qual_length.fastq: S005: " ++
                "sequence and quality lengths differ (record 0, line 4, offset 19)\n",
        },
        .{
            .path = "tests/data/synthetic/bad_quality_range.fastq",
            .stderr = "error: tests/data/synthetic/bad_quality_range.fastq: S006: " ++
                "quality byte must be ASCII 33 through 126 (record 0, line 4, offset 22)\n",
        },
    };
    for (exact_cases) |case| {
        try expectResult(
            try cli.run(allocator, &.{ "sample", "--count", "0", case.path }),
            1,
            "",
            case.stderr,
        );
    }
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "0",
            "tests/data/synthetic/bad_plus.fastq",
        }),
        1,
        "",
        "error: tests/data/synthetic/bad_plus.fastq: S001: " ++
            "plus line must start with '+' (record 0, line 3, offset 15)\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "0",
            "tests/data/synthetic/bad_quality_range.fastq",
        }),
        1,
        "",
        "error: tests/data/synthetic/bad_quality_range.fastq: S006: " ++
            "quality byte must be ASCII 33 through 126 (record 0, line 4, offset 22)\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "1",
            "--alphabet",
            "acgtn",
            "tests/data/synthetic/iupac_valid.fastq",
        }),
        1,
        "",
        "error: tests/data/synthetic/iupac_valid.fastq: S002: " ++
            "sequence byte is outside the selected alphabet (record 0, line 2, offset 11)\n",
    );
}

test "[cli] - [sample]: a later malformed record leaves the selected prefix" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const result = try cli.runWithStdin(
        arena.allocator(),
        &.{ "sample", "--fraction", "1", "-" },
        "@ok\nA\n+\n!\n@bad\nA\nx\n!\n",
        2,
    );

    try expectResult(
        result,
        1,
        "@ok\nA\n+\n!\n",
        "error: -: S001: plus line must start with '+' " ++
            "(record 1, line 3, offset 17)\n",
    );
}

test "[cli] - [sample]: piped input and output make progress together" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    var input: std.ArrayList(u8) = .empty;
    for (0..8192) |_| try input.appendSlice(allocator, "@r\nACGT\n+\n!!!!\n");
    const result = try cli.runWithStdin(
        allocator,
        &.{ "sample", "--fraction", "1", "-" },
        input.items,
        4096,
    );

    try expectResult(result, 0, input.items, "");
}

test "[cli] - [sample]: numeric grammar maps to exact exit classes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for ([_][]const u8{ "0", "1", "0.0", "0.000001", "0.999999", "1.000" }) |value| {
        try expectResult(
            try cli.run(allocator, &.{ "sample", "--fraction", value, "-" }),
            0,
            "",
            "",
        );
    }
    for ([_][]const u8{ "", "00", "0.", ".5", "1.", "1.1", "+0.5", "0e0", "NaN", " 0.5" }) |value| {
        try expectResult(
            try cli.run(allocator, &.{ "sample", "--fraction", value, EMPTY_PATH }),
            2,
            "",
            "error: invalid --fraction value\n",
        );
    }
    for ([_][]const u8{ "", "+1", "-1", "1.0", "1e2", " 1" }) |value| {
        try expectResult(
            try cli.run(allocator, &.{
                "sample",
                "--fraction",
                "0",
                "--seed",
                value,
                EMPTY_PATH,
            }),
            2,
            "",
            "error: invalid --seed value\n",
        );
    }
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "0",
            "--seed",
            "18446744073709551616",
            EMPTY_PATH,
        }),
        4,
        "",
        "error: --seed exceeds supported limit\n",
    );
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "0",
            "--seed",
            "18446744073709551615",
            EMPTY_PATH,
        }),
        0,
        "",
        "",
    );

    for ([_][]const u8{ "0", "1", "00011", "18446744073709551615" }) |value| {
        const expected = if ((try sampling.parseCount(value)) == 0) "" else EMPTY_RECORD;
        try expectResult(
            try cli.run(allocator, &.{ "sample", "--count", value, EMPTY_PATH }),
            0,
            expected,
            "",
        );
    }
    for ([_][]const u8{ "", "+1", "-1", "1.0", "1e2", " 1", "1 " }) |value| {
        try expectResult(
            try cli.run(allocator, &.{ "sample", "--count", value, EMPTY_PATH }),
            2,
            "",
            "error: invalid --count value\n",
        );
    }
    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--count",
            "18446744073709551616",
            EMPTY_PATH,
        }),
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
        .{ .args = &.{ "sample", "--fraction", "0" }, .stderr = "error: sample requires exactly one input\n" },
        .{ .args = &.{ "sample", "--count", "1" }, .stderr = "error: sample requires exactly one input\n" },
        .{ .args = &.{ "sample", "--fraction", "0", EMPTY_PATH, BASIC_PATH }, .stderr = "error: sample requires exactly one input\n" },
        .{ .args = &.{ "sample", "--fraction" }, .stderr = "error: --fraction requires a value\n" },
        .{ .args = &.{ "sample", "--count" }, .stderr = "error: --count requires a value\n" },
        .{
            .args = &.{ "sample", "--fraction", "0.5", "--count", "1", BASIC_PATH },
            .stderr = "error: --fraction and --count are mutually exclusive\n",
        },
        .{
            .args = &.{ "sample", "--count", "1", "-" },
            .stderr = "error: exact-count sampling requires a file path\n",
        },
        .{ .args = &.{ "sample", "--seed" }, .stderr = "error: --seed requires a value\n" },
        .{ .args = &.{ "sample", "--json", EMPTY_PATH }, .stderr = "error: unknown sample option: --json\n" },
        .{ .args = &.{ "sample", "--bogus", EMPTY_PATH }, .stderr = "error: unknown sample option: --bogus\n" },
    };
    for (cases) |case| {
        try expectResult(try cli.run(allocator, case.args), 2, "", case.stderr);
    }

    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "0", "--", "-missing" }),
        3,
        "",
        "error: -missing: file not found\n",
    );
}

test "[cli] - [sample]: line, input, and output failures retain their classes" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    try expectResult(
        try cli.run(allocator, &.{
            "sample",
            "--fraction",
            "1",
            "--max-line-bytes",
            "3",
            BASIC_PATH,
        }),
        4,
        "",
        "error: tests/data/synthetic/basic_valid.fastq: line length limit exceeded\n",
    );
    try expectResult(
        try cli.run(allocator, &.{ "sample", "--fraction", "1", "missing.fastq" }),
        3,
        "",
        "error: missing.fastq: file not found\n",
    );
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--fraction", "1", BASIC_PATH },
            "",
        ),
        3,
        "",
        "",
    );
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--fraction", "1", "-" },
            "@ok\nA\n+\n!\n@bad\nA\nx\n!\n",
        ),
        3,
        "",
        "error: -: S001: plus line must start with '+' " ++
            "(record 1, line 3, offset 17)\n",
    );
    try expectResult(
        try cli.runWithClosedStdin(
            allocator,
            &.{ "sample", "--fraction", "1", "-" },
        ),
        3,
        "",
        "error: -: I/O error\n",
    );
    try expectResult(
        try cli.runWithClosedStdin(
            allocator,
            &.{ "sample", "--count", "1", "-" },
        ),
        2,
        "",
        "error: exact-count sampling requires a file path\n",
    );
    try expectResult(
        try cli.runWithClosedStdout(
            allocator,
            &.{ "sample", "--count", "1", BASIC_PATH },
            "",
        ),
        3,
        "",
        "",
    );
}

fn expectResult(
    result: cli.CommandResult,
    exit_code: u8,
    stdout: []const u8,
    stderr: []const u8,
) !void {
    try std.testing.expectEqual(exit_code, result.exit_code);
    try std.testing.expectEqualStrings(stdout, result.stdout);
    try std.testing.expectEqualStrings(stderr, result.stderr);
}

fn tempPath(
    allocator: std.mem.Allocator,
    sub_path: []const u8,
    name: []const u8,
) ![]const u8 {
    return std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}", .{ sub_path, name });
}
