//! Installed `z-fastq count` results, stream output, options, and exit precedence.

const std = @import("std");
const zfastq = @import("z-fastq");
const cli = @import("utilities.zig");
const CommandResult = cli.CommandResult;
const runCli = cli.run;
const runCliWithStdin = cli.runWithStdin;
const runCliWithClosedStdout = cli.runWithClosedStdout;
const runCliWithClosedStdin = cli.runWithClosedStdin;
const GzipOptions = cli.GzipOptions;
const appendGzipMember = cli.appendGzipMember;
const fixturePath = cli.fixturePath;
const FIXTURE_DIR = cli.FIXTURE_DIR;
const OUTPUT_FAILURE = "error: standard output: I/O error\n";
const GZIP_OPTIONAL_FIELD_LIMIT = 64 * 1024;
const EXPECTED_USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count         Count records in plain or gzip FASTQ inputs
    \\  stats         Report aggregate FASTQ statistics
    \\  check         Validate FASTQ structure, sequence alphabet, and quality range
    \\  sample        Select records by deterministic probability or exact count
    \\  interleave    Validate and interleave paired FASTQ inputs
    \\  deinterleave  Validate and separate interleaved paired FASTQ input
    \\  process       Validate and trim paired FASTQ inputs into two plain or gzip outputs
    \\
    \\General options:
    \\  -h, --help           Show this help message
    \\  -V, --version        Print version
    \\
    \\Input options:
    \\  --max-line-bytes N   Override default line length limit
    \\
    \\Validation options:
    \\  --alphabet POLICY    Select iupac (default) or acgtn sequence symbols
    \\
    \\Machine output:
    \\  --json               Emit versioned JSON (stats and check only)
    \\
    \\Pair options:
    \\  --paired             Validate two inputs as paired reads
    \\  --interleaved        Validate consecutive records as paired reads
    \\  --pair-names POLICY  Select illumina (default) or exact pair names
    \\
    \\Sample options:
    \\  --fraction P         Use 0, 1, 0.DIGITS, or 1.ZEROES
    \\  --count K            Select exactly min(K, records or pairs) from paths
    \\  --seed S             Use an unsigned decimal u64 seed (default 11)
    \\
    \\Process options:
    \\  --gzip               Write both outputs as gzip
    \\  --max-pairs N        Stop after the first N pairs
    \\  --report PATH        Write a z-fastq/process-v1 JSON report (- for stdout)
    \\  --contaminants PATH  Name overrepresented sequences from a FastQC contaminant list
    \\  --pair-duplication   Report fastp's pair duplication rate (up to 1 GiB of memory)
    \\  --adapter1 SEQ       Trim mate 1 adapters and mate overlaps (with --adapter2)
    \\  --adapter2 SEQ       Trim mate 2 adapters and mate overlaps (with --adapter1)
    \\  --trim-poly-g N      Trim 3' poly-G runs, scanning at least N bases
    \\  --trim-tail-window W
    \\                       Trim the 3' end by a W-base quality window (with quality)
    \\  --trim-tail-quality Q
    \\                       Mean Phred quality the 3' window must reach (with window)
    \\  --min-length N       Drop pairs with a read shorter than N bases
    \\  --max-n N            Drop pairs with a read of more than N unknown (N) bases
    \\  --low-quality Q      Count bases below Phred Q as low quality (with percent)
    \\  --max-low-quality-percent P
    \\                       Drop pairs with a read over P% low-quality bases
    \\
    \\Count usage:
    \\  z-fastq count [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
    \\Stats usage:
    \\  z-fastq stats [--json] [--max-line-bytes N] <path|-> [<path|-> ...]
    \\
    \\Check usage:
    \\  z-fastq check [--json] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|-> [<path|-> ...]
    \\  z-fastq check --paired [--json] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <R1|-> <R2|->
    \\  z-fastq check --interleaved [--json] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|->
    \\
    \\Sample usage:
    \\  z-fastq sample --fraction P [--seed S] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|->
    \\  z-fastq sample --count K [--seed S] [--alphabet iupac|acgtn] [--max-line-bytes N] path
    \\  z-fastq sample --paired --fraction P [--seed S] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <R1|-> <R2|->
    \\  z-fastq sample --interleaved --fraction P [--seed S] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <path|->
    \\  z-fastq sample --paired --count K [--seed S] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] R1-path R2-path
    \\  z-fastq sample --interleaved --count K [--seed S] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] path
    \\
    \\Interleave usage:
    \\  z-fastq interleave [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] <R1|-> <R2|->
    \\
    \\Deinterleave usage:
    \\  z-fastq deinterleave [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] --out1 R1 --out2 R2 <path|->
    \\
    \\Process usage:
    \\  z-fastq process [--gzip] [--max-pairs N] [--report PATH] [--contaminants PATH] [--pair-duplication] [--adapter1 SEQ --adapter2 SEQ] [--trim-poly-g N] [--trim-tail-window W --trim-tail-quality Q] [--min-length N] [--max-n N] [--low-quality Q --max-low-quality-percent P] [--pair-names illumina|exact] [--alphabet iupac|acgtn] [--max-line-bytes N] --out1 R1 --out2 R2 <R1|-> <R2|->
    \\
;

const FixtureExpect = struct {
    path: []const u8,
    exit_code: u8,
    stdout: []const u8 = "",
    stderr: []const u8 = "",
    stdin_stderr: []const u8 = "",
};

const FIXTURES = [_]FixtureExpect{
    .{ .path = "basic_valid.fastq", .exit_code = 0, .stdout = "5\n" },
    .{ .path = "crlf.fastq", .exit_code = 0, .stdout = "2\n" },
    .{ .path = "missing_final_newline.fastq", .exit_code = 0, .stdout = "1\n" },
    .{
        .path = "bad_plus.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/bad_plus.fastq: S001: " ++
            "plus line must start with '+' (record 0, line 3, offset 15)\n",
        .stdin_stderr = "error: -: S001: " ++
            "plus line must start with '+' (record 0, line 3, offset 15)\n",
    },
    .{
        .path = "bad_qual_length.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/bad_qual_length.fastq: S005: " ++
            "sequence and quality lengths differ (record 0, line 4, offset 19)\n",
        .stdin_stderr = "error: -: S005: " ++
            "sequence and quality lengths differ (record 0, line 4, offset 19)\n",
    },
    .{
        .path = "bad_header.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/bad_header.fastq: S003: " ++
            "header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n",
        .stdin_stderr = "error: -: S003: " ++
            "header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n",
    },
    .{
        .path = "truncated_record.fastq",
        .exit_code = 1,
        .stderr = "error: tests/data/synthetic/truncated_record.fastq: S004: " ++
            "unexpected end of file in quality line (record 0, line 4, offset 17)\n",
        .stdin_stderr = "error: -: S004: " ++
            "unexpected end of file in quality line (record 0, line 4, offset 17)\n",
    },
};

fn runCountBytes(
    allocator: std.mem.Allocator,
    name: []const u8,
    bytes: []const u8,
) !CommandResult {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    {
        const file = try tmp.dir.createFile(io, name, .{});
        defer file.close(io);
        try std.Io.File.writeStreamingAll(file, io, bytes);
    }
    const path = try cli.tempPath(allocator, &tmp.sub_path, name);
    return runCount(allocator, path);
}

fn runCount(
    allocator: std.mem.Allocator,
    path: []const u8,
) !CommandResult {
    return runCli(allocator, &.{ "count", path });
}

test "[cli] - [count]: files and fragmented stdin produce exact fixture results" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for (FIXTURES) |fixture| {
        const path = try fixturePath(allocator, fixture.path);
        const file_result = try runCount(allocator, path);
        const data = try std.Io.Dir.cwd().readFileAlloc(
            std.testing.io,
            path,
            allocator,
            .limited(1024 * 1024),
        );
        const stdin_result = try runCliWithStdin(allocator, &.{ "count", "-" }, data, 1);

        try std.testing.expectEqual(fixture.exit_code, file_result.exit_code);
        try std.testing.expectEqualStrings(fixture.stdout, file_result.stdout);
        try std.testing.expectEqualStrings(fixture.stderr, file_result.stderr);
        try std.testing.expectEqual(fixture.exit_code, stdin_result.exit_code);
        try std.testing.expectEqualStrings(fixture.stdout, stdin_result.stdout);
        try std.testing.expectEqualStrings(fixture.stdin_stderr, stdin_result.stderr);
    }
}

test "[cli] - [count]: bytes select plain or chained gzip independently of suffix" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    var gzip: std.ArrayList(u8) = .empty;

    try appendGzipMember(allocator, &gzip, "", .{});
    try appendGzipMember(allocator, &gzip, "@a\r\nA\r", .{
        .extra = "xy",
        .name = "reads.fastq",
        .comment = "fixture",
        .header_crc = true,
    });
    try appendGzipMember(allocator, &gzip, "\n+\r\n!\r\n@b\r\nTT\n+name\r\n##\r\n@c\nA\n+\n!\n", .{});

    const file_result = try runCountBytes(allocator, "reads.bin", gzip.items);
    const plain_result = try runCountBytes(
        allocator,
        "plain.fastq.gz",
        "@plain\nA\n+\n!\n",
    );
    const stdin_result = try runCliWithStdin(
        allocator,
        &.{ "count", "-" },
        gzip.items,
        1,
    );

    for ([_]CommandResult{ file_result, stdin_result }) |result| {
        try std.testing.expectEqual(@as(u8, 0), result.exit_code);
        try std.testing.expectEqualStrings("3\n", result.stdout);
        try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    }
    try std.testing.expectEqual(@as(u8, 0), plain_result.exit_code);
    try std.testing.expectEqualStrings("1\n", plain_result.stdout);
    try std.testing.expectEqual(@as(usize, 0), plain_result.stderr.len);
}

test "[cli] - [count]: empty CRLF identifiers fail after dense warmup across gzip members" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const data = "@a\nA\n+\n!\n" ** 2 ++ "@\r\nA\n+\n!\n";
    var whole: std.ArrayList(u8) = .empty;
    var split: std.ArrayList(u8) = .empty;
    try appendGzipMember(allocator, &whole, data, .{});
    try appendGzipMember(allocator, &split, data[0..19], .{});
    try appendGzipMember(allocator, &split, data[19..], .{});

    for ([_][]const u8{ data, whole.items, split.items }) |input| {
        const result = try runCliWithStdin(allocator, &.{ "count", "-" }, input, input.len);

        try std.testing.expectEqual(@as(u8, 1), result.exit_code);
        try std.testing.expectEqualStrings("", result.stdout);
        try std.testing.expectEqualStrings(
            "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
                "(record 2, line 1, offset 18)\n",
            result.stderr,
        );
    }
}

test "[cli] - [count]: files and stdin select plain or gzip input and its exit class" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try cli.tempPath(allocator, &tmp.sub_path, "input");
    var valid: std.ArrayList(u8) = .empty;
    try appendGzipMember(allocator, &valid, "@r\nA\n+\n!\n", .{});
    var empty: std.ArrayList(u8) = .empty;
    try appendGzipMember(allocator, &empty, "", .{});
    var members: std.ArrayList(u8) = .empty;
    try members.appendSlice(allocator, valid.items);
    try appendGzipMember(allocator, &members, "@second\nTT\n+\n##\n", .{});

    const magic = try allocator.dupe(u8, valid.items);
    magic[0] ^= 1;
    const method = try allocator.dupe(u8, valid.items);
    method[2] = 0;
    const checksum = try allocator.dupe(u8, valid.items);
    checksum[cli.gzipCrcOffset(checksum.len)] ^= 1;
    const size = try allocator.dupe(u8, valid.items);
    size[cli.gzipSizeOffset(size.len)] ^= 1;
    const trailing = try std.mem.concat(allocator, u8, &.{ valid.items, "\x00\x00" });

    const plain_error = "S003: header line must start with '@' and contain a nonempty identifier " ++
        "(record 0, line 1, offset 0)\n";
    const cases = [_]struct {
        bytes: []const u8,
        exit_code: u8 = 3,
        stdout: []const u8 = "",
        message: []const u8 = "I/O error\n",
    }{
        .{ .bytes = "", .exit_code = 0, .stdout = "0\n", .message = "" },
        .{ .bytes = empty.items, .exit_code = 0, .stdout = "0\n", .message = "" },
        .{ .bytes = "@r\nA\n+\n!\n", .exit_code = 0, .stdout = "1\n", .message = "" },
        .{ .bytes = members.items, .exit_code = 0, .stdout = "2\n", .message = "" },
        .{
            .bytes = "@r",
            .exit_code = 1,
            .message = "S004: unexpected end of file in sequence line (record 0, line 2, offset 2)\n",
        },
        .{ .bytes = magic, .exit_code = 1, .message = plain_error },
        .{ .bytes = method },
        .{ .bytes = checksum },
        .{ .bytes = size },
        .{ .bytes = valid.items[0 .. valid.items.len - 1] },
        .{ .bytes = trailing },
    };
    for (cases, 0..) |case, index| {
        try tmp.dir.writeFile(io, .{ .sub_path = "input", .data = case.bytes });
        for ([_][]const u8{ path, "-" }, 0..) |label, transport| {
            errdefer std.debug.print("count input case {d}, input {s}\n", .{ index, label });
            const result = if (transport == 0)
                try runCount(allocator, path)
            else
                try runCliWithStdin(allocator, &.{ "count", "-" }, case.bytes, 1);
            const stderr = if (case.message.len == 0) "" else try std.fmt.allocPrint(
                allocator,
                "error: {s}: {s}",
                .{ label, case.message },
            );
            try cli.expectResult(result, case.exit_code, case.stdout, stderr);
        }
    }
}

test "[cli] - [count]: gzip optional name and comment fields enforce the 64 KiB limit" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const field_bytes = try allocator.alloc(u8, GZIP_OPTIONAL_FIELD_LIMIT);
    @memset(field_bytes, 'a');

    for ([_]struct {
        field: enum { name, comment },
        len: usize,
        exit_code: u8,
    }{
        .{ .field = .name, .len = GZIP_OPTIONAL_FIELD_LIMIT - 1, .exit_code = 0 },
        .{ .field = .name, .len = GZIP_OPTIONAL_FIELD_LIMIT, .exit_code = 3 },
        .{ .field = .comment, .len = GZIP_OPTIONAL_FIELD_LIMIT - 1, .exit_code = 0 },
        .{ .field = .comment, .len = GZIP_OPTIONAL_FIELD_LIMIT, .exit_code = 3 },
    }) |case| {
        const field = field_bytes[0..case.len];
        var member: std.ArrayList(u8) = .empty;
        try appendGzipMember(allocator, &member, "@r\nA\n+\n!\n", switch (case.field) {
            .name => .{ .name = field },
            .comment => .{ .comment = field },
        });
        const result = try runCliWithStdin(allocator, &.{ "count", "-" }, member.items, 37);

        if (case.exit_code == 0) {
            try cli.expectResult(result, 0, "1\n", "");
        } else {
            try cli.expectResult(result, 3, "", "error: -: I/O error\n");
        }
    }
}

test "[cli] - [count]: an over-limit stdin line exits with line-limit status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const limited = try runCliWithStdin(
        allocator,
        &.{ "count", "--max-line-bytes", "4", "-" },
        "@r\nAAAAA\n+\n!!!!!\n",
        2,
    );
    try cli.expectResult(limited, 4, "", "error: -: line length limit exceeded\n");
}

test "[cli] - [count]: stdin is explicit and duplicate use fails before input processing" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const implicit = try runCliWithStdin(allocator, &.{"count"}, "", 1);
    try std.testing.expectEqual(@as(u8, 2), implicit.exit_code);
    try std.testing.expectEqual(@as(usize, 0), implicit.stdout.len);
    try std.testing.expectEqualStrings(
        "error: count requires at least one input\n",
        implicit.stderr,
    );

    const duplicate = try cli.runWith(
        allocator,
        &.{
            "count",
            "tests/data/synthetic/does_not_exist.fastq",
            "-",
            "-",
        },
        .{ .stdin = .stalled },
    );
    try std.testing.expectEqual(@as(u8, 2), duplicate.exit_code);
    try std.testing.expectEqual(@as(usize, 0), duplicate.stdout.len);
    try std.testing.expectEqualStrings(
        "error: standard input may appear at most once\n",
        duplicate.stderr,
    );

    for ([_][]const []const u8{
        &.{ "count", "-", "/dev/stdin" },
        &.{ "count", "/dev/stdin", "/dev/stdin" },
        &.{ "stats", "/dev/stdin", "-" },
        &.{ "stats", "-", "-" },
        &.{ "stats", "--json", "-", "-" },
        &.{ "check", "--json", "-", "/dev/stdin" },
        &.{ "check", "missing.fastq", "-", "-" },
        &.{ "check", "--json", "missing.fastq", "-", "-" },
    }) |args| {
        try cli.expectResult(
            try runCliWithStdin(allocator, args, "@r\nA\n+\n!\n", 1),
            2,
            "",
            "error: standard input may appear at most once\n",
        );
    }
}

test "[cli] - [input]: a regular file on standard input is read from its start under each name" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "in.fastq", .data = "@r\nA\n+\n!\n" });
    const stdin_file = try tmp.dir.openFile(io, "in.fastq", .{});
    defer stdin_file.close(io);
    try cli.expectResult(
        try cli.runWithStdinFile(allocator, &.{ "count", "-", "/dev/stdin" }, stdin_file),
        0,
        "1\n1\n",
        "",
    );
}

test "[cli] - [arguments]: a double dash ends options before the path" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const after_double_dash = try runCliWithStdin(
        allocator,
        &.{ "count", "--", "-" },
        "@r\nA\n+\n!\n",
        1,
    );
    try std.testing.expectEqual(@as(u8, 0), after_double_dash.exit_code);
    try std.testing.expectEqualStrings("1\n", after_double_dash.stdout);
    try std.testing.expectEqual(@as(usize, 0), after_double_dash.stderr.len);
}

test "[cli] - [output]: stdout appended to an input is rejected before any output" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const contents = [_][]const u8{ "@r/1\nA\n+\n!\n", "@r/2\nA\n+\n!\n", "@r/1\nA\n+\n!\n@r/2\nA\n+\n!\n" };
    const names = [_][]const u8{ "r1.fastq", "r2.fastq", "pairs.fastq" };
    var paths: [3][]const u8 = undefined;
    for (names, &paths) |name, *path| path.* = try cli.tempPath(allocator, &tmp.sub_path, name);

    const Case = struct { args: []const []const u8, output: usize };
    const cases = [_]Case{
        .{ .args = &.{ "count", paths[0] }, .output = 0 },
        .{ .args = &.{ "count", paths[0], paths[1] }, .output = 1 },
        .{ .args = &.{ "stats", paths[0], paths[1] }, .output = 1 },
        .{ .args = &.{ "stats", "--json", paths[0], paths[1] }, .output = 1 },
        .{ .args = &.{ "check", paths[0] }, .output = 0 },
        .{ .args = &.{ "check", "--json", paths[0], paths[1] }, .output = 1 },
        .{ .args = &.{ "check", "--paired", "--json", paths[0], paths[1] }, .output = 1 },
        .{ .args = &.{ "check", "--interleaved", "--json", paths[2] }, .output = 2 },
    };
    for (cases) |case| {
        for (names, contents) |name, data| try tmp.dir.writeFile(io, .{ .sub_path = name, .data = data });
        const output = try cli.openAppendOutput(tmp.dir, names[case.output]);
        defer output.close(io);
        const diagnostic = try std.fmt.allocPrint(allocator, "error: {s}: input file is output file\n", .{paths[case.output]});
        try cli.expectResult(try cli.runWithStdoutFile(allocator, case.args, output, null), 2, "", diagnostic);
        const after = try tmp.dir.readFileAlloc(io, names[case.output], allocator, .limited(1024));
        try std.testing.expectEqualStrings(contents[case.output], after);
    }

    try tmp.dir.writeFile(io, .{ .sub_path = names[0], .data = contents[0] });
    const input = try tmp.dir.openFile(io, names[0], .{});
    defer input.close(io);
    const output = try cli.openAppendOutput(tmp.dir, names[0]);
    defer output.close(io);
    try cli.expectResult(try cli.runWithStdoutFile(allocator, &.{ "count", "-" }, output, input), 2, "", "error: -: input file is output file\n");
    try std.testing.expectEqualStrings(contents[0], try tmp.dir.readFileAlloc(io, names[0], allocator, .limited(1024)));
}

test "[cli] - [input]: a directory input is named before it is read" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const directory = "tests/data/synthetic";
    for ([_][]const []const u8{
        &.{ "count", directory },
        &.{ "stats", directory },
        &.{ "check", directory },
        &.{ "sample", "--fraction", "1", directory },
        &.{ "check", "--paired", directory, "tests/data/synthetic/basic_valid.fastq" },
    }) |args| {
        try cli.expectResult(try runCli(allocator, args), 3, "", "error: " ++ directory ++ ": input is a directory\n");
    }
}

test "[cli] - [count]: a missing input exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const result = try runCount(allocator, "tests/data/synthetic/does_not_exist.fastq");
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(
        "error: tests/data/synthetic/does_not_exist.fastq: file not found\n",
        result.stderr,
    );
}

test "[cli] - [count]: a long missing path is reported without truncation" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const path = try allocator.alloc(u8, 600);
    @memset(path, 'x');

    const result = try runCount(allocator, path);
    const expected = try std.fmt.allocPrint(
        allocator,
        "error: {s}: failed to open file\n",
        .{path},
    );

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqual(@as(usize, 0), result.stdout.len);
    try std.testing.expectEqualStrings(expected, result.stderr);
}

test "[cli] - [count]: all paths run and the highest exit class wins" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const result = try runCli(allocator, &.{
        "count",
        "tests/data/synthetic/missing_final_newline.fastq",
        "tests/data/synthetic/bad_header.fastq",
        "tests/data/synthetic/does_not_exist.fastq",
    });
    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqualStrings("1\n", result.stdout);
    try std.testing.expectEqualStrings(
        "error: tests/data/synthetic/bad_header.fastq: S003: " ++
            "header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n" ++
            "error: tests/data/synthetic/does_not_exist.fastq: file not found\n",
        result.stderr,
    );
}

test "[cli] - [count]: files around stdin retain output order and exit precedence" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const result = try runCliWithStdin(
        arena.allocator(),
        &.{
            "count",
            "tests/data/synthetic/missing_final_newline.fastq",
            "-",
            "tests/data/synthetic/does_not_exist.fastq",
            "tests/data/synthetic/basic_valid.fastq",
        },
        "@ description\nA\n+\n!\n",
        3,
    );

    try std.testing.expectEqual(@as(u8, 3), result.exit_code);
    try std.testing.expectEqualStrings("1\n5\n", result.stdout);
    try std.testing.expectEqualStrings(
        "error: -: S003: header line must start with '@' and contain a nonempty identifier " ++
            "(record 0, line 1, offset 0)\n" ++
            "error: tests/data/synthetic/does_not_exist.fastq: file not found\n",
        result.stderr,
    );
}

test "[cli] - [root]: version prints the package version in both spellings" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const version = try runCli(allocator, &.{"--version"});
    try cli.expectResult(version, 0, "z-fastq " ++ cli.EXPECTED_VERSION ++ "\n", "");

    const short_version = try runCli(allocator, &.{"-V"});
    try std.testing.expectEqual(@as(u8, 0), short_version.exit_code);
    try std.testing.expectEqualStrings(version.stdout, short_version.stdout);
    try std.testing.expectEqual(@as(usize, 0), short_version.stderr.len);
}

test "[cli] - [root]: help and usage failures are exact" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const help = try runCli(allocator, &.{"--help"});
    try std.testing.expectEqual(@as(u8, 0), help.exit_code);
    try std.testing.expectEqualStrings(EXPECTED_USAGE, help.stdout);
    try std.testing.expectEqual(@as(usize, 0), help.stderr.len);

    const short_help = try runCli(allocator, &.{"-h"});
    try std.testing.expectEqual(@as(u8, 0), short_help.exit_code);
    try std.testing.expectEqualStrings(help.stdout, short_help.stdout);
    try std.testing.expectEqual(@as(usize, 0), short_help.stderr.len);

    const command_help_cases = [_][2][]const u8{
        .{ "count", "--help" },
        .{ "count", "-h" },
        .{ "stats", "--help" },
        .{ "stats", "-h" },
        .{ "check", "--help" },
        .{ "check", "-h" },
        .{ "sample", "--help" },
        .{ "sample", "-h" },
        .{ "interleave", "--help" },
        .{ "interleave", "-h" },
        .{ "deinterleave", "--help" },
        .{ "deinterleave", "-h" },
    };
    for (command_help_cases) |args| {
        const result = try runCli(allocator, &args);
        try std.testing.expectEqual(@as(u8, 0), result.exit_code);
        try std.testing.expectEqualStrings(EXPECTED_USAGE, result.stdout);
        try std.testing.expectEqual(@as(usize, 0), result.stderr.len);
    }

    const missing_command = try runCli(allocator, &.{});
    try std.testing.expectEqual(@as(u8, 2), missing_command.exit_code);
    try std.testing.expectEqual(@as(usize, 0), missing_command.stdout.len);
    try std.testing.expectEqualStrings(EXPECTED_USAGE, missing_command.stderr);

    const unknown_command = try runCli(allocator, &.{"unknown"});
    try std.testing.expectEqual(@as(u8, 2), unknown_command.exit_code);
    try std.testing.expectEqual(@as(usize, 0), unknown_command.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unknown command: unknown\n" ++ EXPECTED_USAGE,
        unknown_command.stderr,
    );

    for ([_][]const u8{ "count", "stats" }) |command| {
        const without_input = try runCli(allocator, &.{command});
        const expected = try std.fmt.allocPrint(allocator, "error: {s} requires at least one input\n", .{command});
        try cli.expectResult(without_input, 2, "", expected);
    }
    const stats_json_without_input = try runCli(allocator, &.{ "stats", "--json" });
    try cli.expectResult(stats_json_without_input, 2, "", "error: stats requires at least one input\n");

    const missing_limit = try runCli(allocator, &.{ "count", "--max-line-bytes" });
    try std.testing.expectEqual(@as(u8, 2), missing_limit.exit_code);
    try std.testing.expectEqual(@as(usize, 0), missing_limit.stdout.len);
    try std.testing.expectEqualStrings(
        "error: --max-line-bytes requires a value\n",
        missing_limit.stderr,
    );

    const version_after_command = try runCli(allocator, &.{ "count", "--version" });
    try cli.expectResult(version_after_command, 2, "", "error: unknown count option: --version\n");
}

test "[cli] - [arguments]: line-limit grammar errors differ from overflow for every command" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const cases = [_]struct { value: []const u8, exit_code: u8 = 2 }{
        .{ .value = "+16" },
        .{ .value = "1_6" },
        .{ .value = "-0" },
        .{ .value = "-1" },
        .{ .value = "+0" },
        .{ .value = "" },
        .{ .value = "nope" },
        .{ .value = " 16" },
        .{ .value = "16 " },
        .{ .value = "1\t6" },
        .{ .value = "16\n" },
        .{ .value = "0x10" },
        .{ .value = "16.0" },
        .{ .value = "184467440737095516160x" },
        .{ .value = "18446744073709551616", .exit_code = 4 },
        .{ .value = "00018446744073709551616", .exit_code = 4 },
        .{ .value = "340282366920938463463374607431768211456", .exit_code = 4 },
    };
    for ([_][]const u8{ "count", "stats", "check", "sample", "interleave", "deinterleave" }) |command| {
        for (cases) |case| {
            const result = try runCli(allocator, &.{ command, "--max-line-bytes", case.value, "-" });
            try std.testing.expectEqual(case.exit_code, result.exit_code);
            try std.testing.expectEqualStrings("", result.stdout);
            try std.testing.expectEqualStrings(if (case.exit_code == 4)
                "error: --max-line-bytes exceeds supported limit\n"
            else
                "error: invalid --max-line-bytes value\n", result.stderr);
        }
    }
}

test "[cli] - [arguments]: a value option may be given only once" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for ([_]struct { args: []const []const u8, option: []const u8 }{
        .{ .args = &.{ "count", "--max-line-bytes", "8", "--max-line-bytes", "9", "-" }, .option = "--max-line-bytes" },
        .{ .args = &.{ "check", "--alphabet", "iupac", "--alphabet", "acgtn", "-" }, .option = "--alphabet" },
        .{ .args = &.{ "sample", "--fraction", "0.5", "--fraction", "0.1", "-" }, .option = "--fraction" },
        .{ .args = &.{ "sample", "--count", "1", "--count", "2", "-" }, .option = "--count" },
        .{ .args = &.{ "sample", "--fraction", "1", "--seed", "1", "--seed", "2", "-" }, .option = "--seed" },
        .{ .args = &.{ "interleave", "--pair-names", "exact", "--pair-names", "illumina", "a", "b" }, .option = "--pair-names" },
        .{ .args = &.{ "deinterleave", "--out1", "a", "--out1", "b", "--out2", "c", "-" }, .option = "--out1" },
        .{ .args = &.{ "deinterleave", "--out1", "a", "--out2", "b", "--out2", "c", "-" }, .option = "--out2" },
    }) |case| {
        const message = try std.fmt.allocPrint(allocator, "error: {s} may appear only once\n", .{case.option});
        try cli.expectResult(try runCli(allocator, case.args), 2, "", message);
    }
}

test "[cli] - [arguments]: usage failures name the command and exit with usage status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for ([_]struct { args: []const []const u8, stderr: []const u8 }{
        .{ .args = &.{ "count", "--json", "-" }, .stderr = "error: unknown count option: --json\n" },
        .{ .args = &.{ "count", "--bogus" }, .stderr = "error: unknown count option: --bogus\n" },
        .{ .args = &.{ "sample", "--fraction" }, .stderr = "error: --fraction requires a value\n" },
        .{ .args = &.{ "sample", "--fraction", "0" }, .stderr = "error: sample requires exactly one input\n" },
        .{ .args = &.{ "sample", "--json", "-" }, .stderr = "error: unknown sample option: --json\n" },
        .{ .args = &.{"check"}, .stderr = "error: check requires at least one input\n" },
        .{ .args = &.{ "check", "--json" }, .stderr = "error: check requires at least one input\n" },
        .{ .args = &.{ "check", "--alphabet", "dna", "-" }, .stderr = "error: --alphabet must be iupac or acgtn\n" },
        .{ .args = &.{ "check", "--paired", "missing.fastq", "other.fastq", "extra.fastq" }, .stderr = "error: check --paired requires exactly two inputs\n" },
        .{ .args = &.{ "check", "--interleaved" }, .stderr = "error: check --interleaved requires exactly one input\n" },
        .{ .args = &.{ "check", "--interleaved", "missing.fastq", "extra.fastq" }, .stderr = "error: check --interleaved requires exactly one input\n" },
        .{ .args = &.{ "check", "--paired", "--interleaved", "-" }, .stderr = "error: --paired and --interleaved are mutually exclusive\n" },
        .{ .args = &.{ "check", "--pair-names" }, .stderr = "error: --pair-names requires a value\n" },
        .{ .args = &.{ "check", "--pair-names", "prefix", "-" }, .stderr = "error: --pair-names must be illumina or exact\n" },
        .{ .args = &.{ "check", "--pair-names", "exact", "-" }, .stderr = "error: --pair-names requires --paired or --interleaved\n" },
        .{ .args = &.{ "check", "--paired", "-", "-" }, .stderr = "error: paired inputs may contain standard input at most once\n" },
        .{ .args = &.{ "interleave", "--pair-names" }, .stderr = "error: --pair-names requires a value\n" },
        .{ .args = &.{ "interleave", "-" }, .stderr = "error: interleave requires exactly two inputs\n" },
        .{ .args = &.{ "interleave", "-", "-" }, .stderr = "error: interleave inputs may contain standard input at most once\n" },
        .{ .args = &.{ "interleave", "--json", "a", "b" }, .stderr = "error: unknown interleave option: --json\n" },
        .{ .args = &.{ "deinterleave", "--out1" }, .stderr = "error: --out1 requires a value\n" },
        .{ .args = &.{ "deinterleave", "--out2", "b", "-" }, .stderr = "error: deinterleave requires --out1\n" },
        .{ .args = &.{ "deinterleave", "--out1", "a", "--out2", "a", "-" }, .stderr = "error: deinterleave output paths must differ\n" },
        .{ .args = &.{ "deinterleave", "--unknown", "-" }, .stderr = "error: unknown deinterleave option: --unknown\n" },
        .{ .args = &.{ "check", "--alphabet" }, .stderr = "error: --alphabet requires a value\n" },
        .{ .args = &.{ "check", "--bogus", "-" }, .stderr = "error: unknown check option: --bogus\n" },
        .{ .args = &.{ "check", "--paired", "missing.fastq" }, .stderr = "error: check --paired requires exactly two inputs\n" },
    }) |case| {
        const result = try cli.runWith(allocator, case.args, .{ .stdin = .stalled });
        try cli.expectResult(result, 2, "", case.stderr);
    }
}

test "[cli] - [arguments]: a double dash and escaped option bytes apply to every command" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for ([_][]const u8{ "sample", "interleave", "deinterleave", "check" }) |command| {
        const option = try std.fmt.allocPrint(allocator, "error: unknown {s} option: -bad\\x09\\\\\\x1B\n", .{command});
        const bad_option = try runCli(allocator, &.{ command, "-bad\t\\\x1b" });
        try cli.expectResult(bad_option, 2, "", option);
    }
    for ([_][]const []const u8{
        &.{ "sample", "--fraction", "1", "--", "-z-fastq-v002-missing/record.fastq" },
        &.{ "interleave", "--", "-z-fastq-v002-missing/record.fastq", "tests/data/synthetic/basic_valid.fastq" },
        &.{ "check", "--", "-z-fastq-v002-missing/record.fastq" },
        &.{ "count", "--", "-z-fastq-v002-missing/record.fastq" },
    }) |args| {
        try cli.expectResult(
            try runCli(allocator, args),
            3,
            "",
            "error: -z-fastq-v002-missing/record.fastq: file not found\n",
        );
    }
}

test "[cli] - [arguments]: decimal line limits preserve zero, leading zeros, and usize maximum" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    for ([_][]const u8{ "0", "000" }) |value| {
        const result = try runCli(allocator, &.{ "count", "--max-line-bytes", value, "-" });
        try std.testing.expectEqual(@as(u8, 0), result.exit_code);
        try std.testing.expectEqualStrings("0\n", result.stdout);
        try std.testing.expectEqualStrings("", result.stderr);
    }
    const record = "@r\nA\n+\n!\n";
    for ([_][]const u8{
        "2", "0002", "00016", "18446744073709551615", "00018446744073709551615", "0" ** 128 ++ "2",
    }) |value| {
        const result = try runCliWithStdin(allocator, &.{ "count", "--max-line-bytes", value, "-" }, record, 1);
        try std.testing.expectEqual(@as(u8, 0), result.exit_code);
        try std.testing.expectEqualStrings("1\n", result.stdout);
        try std.testing.expectEqualStrings("", result.stderr);
    }
    const limited = try runCliWithStdin(allocator, &.{ "count", "--max-line-bytes", "0001", "-" }, record, 1);
    try std.testing.expectEqual(@as(u8, 4), limited.exit_code);
    try std.testing.expectEqualStrings("", limited.stdout);
    try std.testing.expectEqualStrings("error: -: line length limit exceeded\n", limited.stderr);
}

test "[cli] - [output]: a closed stdout exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for ([_]struct {
        args: []const []const u8,
        stdin: []const u8 = "",
        exit_code: u8 = 3,
        stderr: []const u8 = OUTPUT_FAILURE,
    }{
        .{ .args = &.{"--help"} },
        .{ .args = &.{"--version"} },
        .{ .args = &.{ "count", "tests/data/synthetic/basic_valid.fastq" } },
        .{ .args = &.{ "count", "-" }, .stdin = "@r\nA\n+\n!\n" },
        .{ .args = &.{ "stats", "-" }, .stdin = "@r\nA\n+\n!\n" },
        .{ .args = &.{ "stats", "--json", "-" }, .stdin = "@r\nA\n+\n!\n" },
        .{ .args = &.{ "sample", "--fraction", "1", "-" }, .stdin = "@r\nA\n+\n!\n" },
        .{ .args = &.{ "check", "--json", "-" }, .stdin = "@r\nA\n+\n!\n" },
        .{
            .args = &.{
                "count",
                "--max-line-bytes",
                "16",
                "-",
                "tests/data/synthetic/missing_final_newline.fastq",
            },
            .stdin = "@r\nAAAAAAAAAAAAAAAAA\n+\n!!!!!!!!!!!!!!!!!\n",
            .exit_code = 4,
            .stderr = "error: -: line length limit exceeded\n" ++ OUTPUT_FAILURE,
        },
    }) |case| {
        const result = try runCliWithClosedStdout(allocator, case.args, case.stdin);
        try std.testing.expectEqual(case.exit_code, result.exit_code);
        try std.testing.expectEqualStrings(case.stderr, result.stderr);
    }
}

test "[cli] - [output]: a pipe closed by its reader exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for ([_][]const []const u8{
        &.{"--version"},
        &.{ "count", "tests/data/synthetic/basic_valid.fastq" },
        &.{ "stats", "tests/data/synthetic/basic_valid.fastq" },
        &.{ "stats", "--json", "tests/data/synthetic/basic_valid.fastq" },
        &.{ "check", "--json", "tests/data/synthetic/basic_valid.fastq" },
        &.{ "sample", "--fraction", "1", "tests/data/synthetic/basic_valid.fastq" },
        &.{ "sample", "--count", "1", "tests/data/synthetic/basic_valid.fastq" },
    }) |args| {
        var ends: [2]i32 = undefined;
        try std.testing.expectEqual(@as(usize, 0), std.os.linux.pipe2(&ends, .{ .CLOEXEC = true }));
        _ = std.os.linux.close(ends[0]);
        const write_end: std.Io.File = .{ .handle = ends[1], .flags = .{ .nonblocking = false } };
        defer _ = std.os.linux.close(ends[1]);
        try cli.expectResult(try cli.runWithStdoutFile(allocator, args, write_end, null), 3, "", OUTPUT_FAILURE);
    }
}

test "[cli] - [input]: a closed standard input exits with I/O status" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    for ([_][]const []const u8{
        &.{ "count", "-" },
        &.{ "stats", "-" },
        &.{ "sample", "--fraction", "1", "-" },
    }) |args| {
        const result = try runCliWithClosedStdin(allocator, args);
        try cli.expectResult(result, 3, "", "error: -: I/O error\n");
    }
}

test "[unit] - [test runner]: a stalled child is killed and reported as a timeout" {
    try std.testing.expectError(error.ChildProcessTimedOut, cli.runWith(
        std.testing.allocator,
        &.{ "count", "-" },
        .{ .stdin = .stalled, .timeout = .fromMilliseconds(100) },
    ));
}

test "[unit] - [test runner]: captured streams stop at the caller's limit and leave no child" {
    const allocator = std.testing.allocator;
    const limit = 64 * 1024;
    const record = "@record01\nA\n+\n!\n";
    const input = record ** (limit / record.len + 1);
    var missing_paths: [18][]const u8 = undefined;
    missing_paths[0] = "count";
    @memset(missing_paths[1..], "x" ** 4096);

    {
        const result = try cli.runWith(allocator, &.{ "sample", "--fraction", "1", "-" }, .{
            .stdin = .{ .data = .{ .bytes = input, .chunk_len = 4096 } },
            .output_limit = input.len,
        });
        defer allocator.free(result.stdout);
        defer allocator.free(result.stderr);
        try cli.expectResult(result, 0, input, "");
    }

    const cases = [_]struct {
        args: []const []const u8,
        stdin: []const u8 = "",
        stdout: cli.StdoutMode = .capture,
    }{
        .{ .args = &.{ "sample", "--fraction", "1", "-" }, .stdin = input },
        .{ .args = &missing_paths },
        .{ .args = &missing_paths, .stdout = .closed },
    };
    for (cases) |case| {
        const result = cli.runWith(allocator, case.args, .{
            .stdin = .{ .data = .{ .bytes = case.stdin, .chunk_len = 4096 } },
            .stdout = case.stdout,
            .output_limit = limit,
        });
        if (result) |unexpected| {
            defer allocator.free(unexpected.stdout);
            defer allocator.free(unexpected.stderr);
            return error.ExpectedOutputLimit;
        } else |err| {
            try std.testing.expectEqual(error.StreamTooLong, err);
        }

        var children: [128]u8 = undefined;
        try std.testing.expectEqualStrings("", try std.Io.Dir.cwd().readFile(
            std.testing.io,
            "/proc/thread-self/children",
            &children,
        ));
    }
}

test "[cli] - [diagnostics]: untrusted command, option, and path bytes use escaped ASCII" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const command = try runCli(allocator, &.{"bad\n\\\x1b"});
    try std.testing.expectEqual(@as(u8, 2), command.exit_code);
    try std.testing.expectEqual(@as(usize, 0), command.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unknown command: bad\\x0A\\\\\\x1B\n" ++ EXPECTED_USAGE,
        command.stderr,
    );

    const option = try runCli(allocator, &.{ "count", "-bad\t\\\x1b" });
    try std.testing.expectEqual(@as(u8, 2), option.exit_code);
    try std.testing.expectEqual(@as(usize, 0), option.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unknown count option: -bad\\x09\\\\\\x1B\n",
        option.stderr,
    );

    const path = try runCount(allocator, "unsafe\n\x1b.fastq");
    try std.testing.expectEqual(@as(u8, 3), path.exit_code);
    try std.testing.expectEqual(@as(usize, 0), path.stdout.len);
    try std.testing.expectEqualStrings(
        "error: unsafe\\x0A\\x1B.fastq: file not found\n",
        path.stderr,
    );
}

test "[edge] - [count]: a line one byte past the read buffer is accepted at the limit" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();

    const line_len: usize = zfastq.limits.COUNT_READ_BUFFER_BYTES + 1;
    const sequence = try allocator.alloc(u8, line_len);
    @memset(sequence, 'A');
    const quality = try allocator.alloc(u8, line_len);
    @memset(quality, '!');
    const data = try std.mem.concat(allocator, u8, &.{ "@long\n", sequence, "\n+\n", quality, "\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "long.fastq", .data = data });

    const path = try cli.tempPath(allocator, &tmp.sub_path, "long.fastq");
    try cli.expectResult(try runCount(allocator, path), 0, "1\n", "");

    const limited = try runCli(allocator, &.{
        "count",
        "--max-line-bytes",
        std.fmt.comptimePrint("{d}", .{zfastq.limits.COUNT_READ_BUFFER_BYTES}),
        path,
    });
    const limited_expected = try std.fmt.allocPrint(
        allocator,
        "error: {s}: line length limit exceeded\n",
        .{path},
    );
    try cli.expectResult(limited, 4, "", limited_expected);
}
