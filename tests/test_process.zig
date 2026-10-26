//! Installed `z-fastq process` behavior: paired validation, two outputs, the report, and exit status.

const std = @import("std");
const cli = @import("utilities.zig");

const R1 = "@p1/1 desc\nACGT\n+\nIIII\n@p2/1\nGGCC\n+p2\n#!#!\n";
const R2 = "@p1/2 desc\nTTTT\n+\nIIII\n@p2/2\nAANN\n+\n!!!!\n";
const R1_CRLF = "@p1/1 desc\r\nACGT\r\n+\r\nIIII\r\n@p2/1\r\nGGCC\r\n+p2\r\n#!#!\r\n";
const OUTPUT_READ_LIMIT = 4 * 1024 * 1024;
const READ_BUFFER_BYTES = 256 * 1024;

const Paths = struct {
    r1: []const u8,
    r2: []const u8,
    out1: []const u8,
    out2: []const u8,
    report: []const u8,

    fn init(allocator: std.mem.Allocator, tmp: *const std.testing.TmpDir) !Paths {
        return .{
            .r1 = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq"),
            .r2 = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq"),
            .out1 = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq"),
            .out2 = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq"),
            .report = try cli.tempPath(allocator, &tmp.sub_path, "report.json"),
        };
    }
};

fn writeInputs(tmp: *const std.testing.TmpDir, r1: []const u8, r2: []const u8) !void {
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "r1.fastq", .data = r1 });
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "r2.fastq", .data = r2 });
}

fn readPath(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    return std.Io.Dir.cwd().readFileAlloc(std.testing.io, path, allocator, .limited(OUTPUT_READ_LIMIT));
}

fn pathExists(path: []const u8) bool {
    std.Io.Dir.cwd().access(std.testing.io, path, .{}) catch return false;
    return true;
}

fn gunzipPath(allocator: std.mem.Allocator, path: []const u8) ![]u8 {
    const compressed = try readPath(allocator, path);
    var input: std.Io.Reader = .fixed(compressed);
    var decoded: std.Io.Writer.Allocating = .init(allocator);
    var decompress: std.compress.flate.Decompress = .init(&input, .gzip, &.{});
    _ = try decompress.reader.streamRemaining(&decoded.writer);
    return decoded.written();
}

fn parseReport(allocator: std.mem.Allocator, bytes: []const u8) !std.json.ObjectMap {
    const value = try std.json.parseFromSliceLeaky(std.json.Value, allocator, bytes, .{});
    try std.testing.expectEqualStrings("z-fastq/process-v1", value.object.get("schema").?.string);
    return value.object;
}

fn mate(report: std.json.ObjectMap, index: usize) std.json.ObjectMap {
    return report.get("mates").?.array.items[index].object;
}

test "[cli] - [process]: pairs pass through with LF endings and the report counts each mate" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const paths = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1_CRLF, R2);

    const result = try cli.run(allocator, &.{ "process", "--report", paths.report, "--out1", paths.out1, "--out2", paths.out2, paths.r1, paths.r2 });
    try cli.expectResult(result, 0, "", "");
    try std.testing.expectEqualStrings(R1, try readPath(allocator, paths.out1));
    try std.testing.expectEqualStrings(R2, try readPath(allocator, paths.out2));

    const report = try parseReport(allocator, try readPath(allocator, paths.report));
    try cli.expectJsonString(report.get("status"), "ok");
    try std.testing.expect(report.get("input_complete").?.bool);
    try cli.expectJsonInteger(report.get("pairs").?.object.get("read"), 2);
    try cli.expectJsonInteger(report.get("pairs").?.object.get("written"), 2);
    try cli.expectJsonInteger(mate(report, 0).get("bases"), 8);
    try cli.expectJsonInteger(mate(report, 0).get("c"), 3);
    try cli.expectJsonInteger(mate(report, 0).get("quality_sum"), 164);
    try cli.expectJsonInteger(mate(report, 1).get("n"), 2);
    try cli.expectJsonInteger(mate(report, 1).get("q20_bases"), 4);
    try std.testing.expectEqual(std.json.Value.null, report.get("error").?);
    try std.testing.expectEqual(std.json.Value.null, report.get("failed_input").?);
    try cli.expectJsonInteger(report.get("inputs").?.array.items[1].object.get("size"), R2.len);
}

test "[cli] - [process]: gzip outputs decode to the same records, and empty input gives empty members" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const r1 = try cli.tempPath(allocator, &tmp.sub_path, "r1.fastq");
    const r2 = try cli.tempPath(allocator, &tmp.sub_path, "r2.fastq");

    for ([_][2][]const u8{ .{ R1, R2 }, .{ "", "" } }, 0..) |inputs, case| {
        errdefer std.debug.print("case {d}\n", .{case});
        try writeInputs(&tmp, inputs[0], inputs[1]);
        const out1 = try cli.tempPath(allocator, &tmp.sub_path, try std.fmt.allocPrint(allocator, "{d}_1.fastq.gz", .{case}));
        const out2 = try cli.tempPath(allocator, &tmp.sub_path, try std.fmt.allocPrint(allocator, "{d}_2.fastq.gz", .{case}));
        const result = try cli.run(allocator, &.{ "process", "--gzip", "--out1", out1, "--out2", out2, r1, r2 });
        try cli.expectResult(result, 0, "", "");
        try std.testing.expectEqualStrings(inputs[0], try gunzipPath(allocator, out1));
        try std.testing.expectEqualStrings(inputs[1], try gunzipPath(allocator, out2));
    }
}

test "[cli] - [process]: reruns write identical outputs and reports" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const paths = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1, R2);

    const out1 = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq.gz");
    const out2 = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq.gz");
    var runs: [2][3][]u8 = undefined;
    for (&runs) |*files| {
        const run_paths = [3][]const u8{ out1, out2, paths.report };
        const result = try cli.run(allocator, &.{ "process", "--gzip", "--report", paths.report, "--out1", out1, "--out2", out2, paths.r1, paths.r2 });
        try cli.expectResult(result, 0, "", "");
        for (run_paths, files) |path, *bytes| {
            bytes.* = try readPath(allocator, path);
            try std.Io.Dir.cwd().deleteFile(std.testing.io, path);
        }
    }
    for (runs[0], runs[1]) |first, second| try std.testing.expectEqualSlices(u8, first, second);
}

test "[cli] - [process]: invalid invocations fail before any input is read or output created" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1, R2);
    const gz1 = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq.gz");
    const gz2 = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq.gz");

    const Case = struct { args: []const []const u8, exit_code: u8 = 2, stderr: []const u8 };
    const cases = [_]Case{
        .{ .args = &.{ "process", "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: process requires --out1\n" },
        .{ .args = &.{ "process", "--out1", p.out1, p.r1, p.r2 }, .stderr = "error: process requires --out2\n" },
        .{ .args = &.{ "process", "--out1", p.out1, "--out2", p.out2, p.r1 }, .stderr = "error: process requires exactly two inputs\n" },
        .{ .args = &.{ "process", "--out1", p.out1, "--out2", p.out2, "-", "-" }, .stderr = "error: process inputs may contain standard input at most once\n" },
        .{ .args = &.{ "process", "--out1", "-", "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: process output paths cannot be standard output\n" },
        .{ .args = &.{ "process", "--out1", gz1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: --out1 names a .gz file but --gzip is not given\n" },
        .{ .args = &.{ "process", "--out1", p.out1, "--out2", gz2, p.r1, p.r2 }, .stderr = "error: --out2 names a .gz file but --gzip is not given\n" },
        .{ .args = &.{ "process", "--out1", p.out1, "--out2", p.out1, p.r1, p.r2 }, .stderr = "error: process output paths must differ\n" },
        .{ .args = &.{ "process", "--report", p.out2, "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: process output paths must differ\n" },
        .{ .args = &.{ "process", "--report", p.report, "--report", p.report, p.r1, p.r2 }, .stderr = "error: --report may appear only once\n" },
        .{ .args = &.{ "process", "--max-pairs", "1x", p.r1, p.r2 }, .stderr = "error: invalid --max-pairs value\n" },
        .{ .args = &.{ "process", "--max-pairs", "18446744073709551616", p.r1, p.r2 }, .exit_code = 4, .stderr = "error: --max-pairs exceeds supported limit\n" },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        const result = try cli.runWith(allocator, case.args, .{ .stdin = .stalled });
        try cli.expectResult(result, case.exit_code, "", case.stderr);
        for ([_][]const u8{ p.out1, p.out2, gz1, gz2, p.report }) |path| try std.testing.expect(!pathExists(path));
    }
}

test "[cli] - [process]: an existing output fails and removes the outputs this run created" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1, R2);

    const Case = struct { out1: []const u8, out2: []const u8, report: []const u8, existing: []const u8 };
    const cases = [_]Case{
        .{ .out1 = p.r1, .out2 = p.out2, .report = p.report, .existing = p.r1 },
        .{ .out1 = p.out1, .out2 = p.r2, .report = p.report, .existing = p.r2 },
        .{ .out1 = p.out1, .out2 = p.out2, .report = p.r1, .existing = p.r1 },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        const result = try cli.run(allocator, &.{ "process", "--report", case.report, "--out1", case.out1, "--out2", case.out2, p.r1, p.r2 });
        try cli.expectResult(result, 3, "", try std.fmt.allocPrint(allocator, "error: {s}: output path already exists\n", .{case.existing}));
        try std.testing.expect(!pathExists(p.out1));
        try std.testing.expect(!pathExists(p.out2));
        try std.testing.expect(!pathExists(p.report));
        try std.testing.expectEqualStrings(R1, try readPath(allocator, p.r1));
        try std.testing.expectEqualStrings(R2, try readPath(allocator, p.r2));
    }
}

test "[cli] - [process]: a failing pair keeps the earlier pairs and the report names the failure" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const first_r1 = "@p1/1\nACGT\n+\nIIII\n";
    const first_r2 = "@p1/2\nTTTT\n+\nIIII\n";

    const Case = struct { r1: []const u8, r2: []const u8, code: []const u8, failed_input: ?[]const u8, stderr_first_line: []const u8 };
    const cases = [_]Case{
        .{ .r1 = first_r1 ++ "@p2/1\nA\n+\nI\n", .r2 = first_r2 ++ "@p3/2\nA\n+\nI\n", .code = "P001", .failed_input = null, .stderr_first_line = "P001: paired identifiers or mate markers do not match (pair 1)\n" },
        .{ .r1 = first_r1 ++ "@p2/1\nA\n+\nI\n", .r2 = first_r2, .code = "P002", .failed_input = null, .stderr_first_line = "P002: paired input is missing a mate (pair 1, remaining R1, last R1 record 1, last R2 record 0)\n" },
        .{ .r1 = first_r1 ++ "@p2/1\nA\n+\nI\n", .r2 = first_r2 ++ "@p2/2\nA\n+\nII\n", .code = "S005", .failed_input = "r2.fastq", .stderr_first_line = "S005: sequence and quality lengths differ (record 1, line 4, offset 28)\n" },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        try writeInputs(&tmp, case.r1, case.r2);
        const p = try Paths.init(allocator, &tmp);
        const out1 = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.out1, index });
        const out2 = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.out2, index });
        const report_path = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.report, index });
        const result = try cli.run(allocator, &.{ "process", "--report", report_path, "--out1", out1, "--out2", out2, p.r1, p.r2 });
        try std.testing.expectEqual(@as(u8, 1), result.exit_code);
        try std.testing.expect(std.mem.endsWith(u8, std.mem.sliceTo(result.stderr, '\n'), std.mem.trimEnd(u8, case.stderr_first_line, "\n")));
        try std.testing.expectEqualStrings(first_r1, try readPath(allocator, out1));
        try std.testing.expectEqualStrings(first_r2, try readPath(allocator, out2));

        const report = try parseReport(allocator, try readPath(allocator, report_path));
        try cli.expectJsonString(report.get("status"), "error");
        try std.testing.expect(!report.get("input_complete").?.bool);
        try cli.expectJsonInteger(report.get("pairs").?.object.get("written"), 1);
        try cli.expectJsonString(report.get("error").?.object.get("code"), case.code);
        if (case.failed_input) |name| {
            try std.testing.expect(std.mem.endsWith(u8, report.get("failed_input").?.string, name));
        } else {
            try std.testing.expectEqual(std.json.Value.null, report.get("failed_input").?);
        }
    }
}

test "[cli] - [process]: a pair cap stops after N pairs and reports the input incomplete" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const first_r1 = "@p1/1 desc\nACGT\n+\nIIII\n";
    const first_r2 = "@p1/2 desc\nTTTT\n+\nIIII\n";
    try writeInputs(&tmp, R1, R2 ++ "@lone/2\nA\n+\nI\n");

    const Case = struct { cap: []const u8, r1: []const u8, r2: []const u8, read: u64, complete: bool, exit_code: u8 };
    const cases = [_]Case{
        .{ .cap = "0", .r1 = "", .r2 = "", .read = 0, .complete = false, .exit_code = 0 },
        .{ .cap = "1", .r1 = first_r1, .r2 = first_r2, .read = 1, .complete = false, .exit_code = 0 },
        .{ .cap = "2", .r1 = R1, .r2 = R2, .read = 2, .complete = false, .exit_code = 0 },
        .{ .cap = "3", .r1 = R1, .r2 = R2, .read = 2, .complete = false, .exit_code = 1 },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("cap {s}\n", .{case.cap});
        const out1 = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.out1, index });
        const out2 = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.out2, index });
        const result = try cli.run(allocator, &.{ "process", "--max-pairs", case.cap, "--report", "-", "--out1", out1, "--out2", out2, p.r1, p.r2 });
        try std.testing.expectEqual(case.exit_code, result.exit_code);
        try std.testing.expectEqualStrings(case.r1, try readPath(allocator, out1));
        try std.testing.expectEqualStrings(case.r2, try readPath(allocator, out2));
        const report = try parseReport(allocator, result.stdout);
        try cli.expectJsonInteger(report.get("pairs").?.object.get("read"), case.read);
        try std.testing.expectEqual(case.complete, report.get("input_complete").?.bool);
        try cli.expectJsonInteger(report.get("options").?.object.get("max_pairs"), try std.fmt.parseInt(u64, case.cap, 10));
    }
}

test "[cli] - [process]: without a report the run writes only its outputs" {
    var tmp = std.testing.tmpDir(.{ .iterate = true });
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1, R2);

    const result = try cli.run(allocator, &.{ "process", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
    try cli.expectResult(result, 0, "", "");
    var entries: usize = 0;
    var iterator = tmp.dir.iterate();
    while (try iterator.next(std.testing.io)) |_| entries += 1;
    try std.testing.expectEqual(@as(usize, 4), entries);
}

test "[cli] - [process]: a report on standard output appended to an input is rejected before output" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1, R2);

    const stdout_file = try cli.openAppendOutput(tmp.dir, "r2.fastq");
    defer stdout_file.close(std.testing.io);
    const result = try cli.runWith(allocator, &.{ "process", "--report", "-", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .{ .stdout = .{ .file = stdout_file } });
    try cli.expectResult(result, 2, "", try std.fmt.allocPrint(allocator, "error: {s}: input file is output file\n", .{p.r2}));
    try std.testing.expectEqualStrings(R2, try readPath(allocator, p.r2));
    try std.testing.expect(!pathExists(p.out1));
}

test "[cli] - [process]: a record longer than the read buffer passes through unchanged" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const long_len = READ_BUFFER_BYTES + 1;
    const sequence = try allocator.alloc(u8, long_len);
    @memset(sequence, 'A');
    const quality = try allocator.alloc(u8, long_len);
    @memset(quality, 'I');
    const r1 = try std.fmt.allocPrint(allocator, "@short/1\nAC\n+\nII\n@long/1\n{s}\n+\n{s}\n", .{ sequence, quality });
    const r2 = "@short/2\nGT\n+\nII\n@long/2\nG\n+\nI\n";
    try writeInputs(&tmp, r1, r2);

    const result = try cli.run(allocator, &.{ "process", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
    try cli.expectResult(result, 0, "", "");
    try std.testing.expectEqualStrings(r1, try readPath(allocator, p.out1));
    try std.testing.expectEqualStrings(r2, try readPath(allocator, p.out2));
}
