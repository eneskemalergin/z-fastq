//! Installed `z-fastq process` behavior: paired validation, two outputs, the report, and exit status.

const std = @import("std");
const cli = @import("utilities.zig");

const R1 = "@p1/1 desc\nACGT\n+\nIIII\n@p2/1\nGGCC\n+p2\n#!#!\n";
const R2 = "@p1/2 desc\nTTTT\n+\nIIII\n@p2/2\nAANN\n+\n!!!!\n";
const R1_CRLF = "@p1/1 desc\r\nACGT\r\n+\r\nIIII\r\n@p2/1\r\nGGCC\r\n+p2\r\n#!#!\r\n";
const OUTPUT_READ_LIMIT = 4 * 1024 * 1024;
const READ_BUFFER_BYTES = 256 * 1024;
const RF1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const RF2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
const FF = "CTGTCTCTTATACACATCT";
const TRIM_STEPS = [_][]const u8{ "quality_tail", "poly_g", "overlap", "adapter" };

const Probe = struct { name: []const u8, rf: ?[2]usize, ff: ?[2]usize };

const PROBES = [_]Probe{
    .{ .name = "ov_rev_L45", .rf = .{ 45, 45 }, .ff = .{ 45, 45 } },
    .{ .name = "ov_rev_L31", .rf = .{ 31, 31 }, .ff = .{ 31, 31 } },
    .{ .name = "ov_rev_L30", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ov_rev_L29", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ov_fwd_offset0_L60", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ov_mm5_first50_L50", .rf = .{ 50, 50 }, .ff = .{ 50, 50 } },
    .{ .name = "ov_N5_first50_L50", .rf = .{ 50, 50 }, .ff = .{ 50, 50 } },
    .{ .name = "ov_mm6_first50_L50", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ov_N6_first50_L50", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ov_mm6_after50_L60", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ov_short_r1_mm4", .rf = .{ 24, 24 }, .ff = .{ 24, 24 } },
    .{ .name = "ov_short_r1_mm5", .rf = .{ 24, 60 }, .ff = .{ 24, 60 } },
    .{ .name = "ov_unequal_fwd_60_40", .rf = .{ 60, 40 }, .ff = .{ 60, 40 } },
    .{ .name = "ov_unequal_rev_40_60", .rf = .{ 40, 40 }, .ff = .{ 40, 40 } },
    .{ .name = "ad_full_after30", .rf = .{ 30, 60 }, .ff = .{ 30, 60 } },
    .{ .name = "ad_partial_k4", .rf = .{ 34, 60 }, .ff = .{ 34, 60 } },
    .{ .name = "ad_partial_k5", .rf = .{ 30, 60 }, .ff = .{ 30, 60 } },
    .{ .name = "ad_partial_k7_mm1", .rf = .{ 37, 60 }, .ff = .{ 37, 60 } },
    .{ .name = "ad_partial_k8_mm1", .rf = .{ 30, 60 }, .ff = .{ 30, 60 } },
    .{ .name = "ad_mm_at_limit", .rf = .{ 30, 60 }, .ff = .{ 30, 60 } },
    .{ .name = "ad_mm_limit_plus1", .rf = .{ 63, 60 }, .ff = .{ 49, 60 } },
    .{ .name = "ad_start_full", .rf = null, .ff = null },
    .{ .name = "ad_start_suffix4", .rf = null, .ff = null },
    .{ .name = "ad_internal_full", .rf = .{ 20, 60 }, .ff = .{ 20, 60 } },
    .{ .name = "ad_internal_lookalike", .rf = .{ 60, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ad_indel_internal", .rf = .{ 72, 60 }, .ff = .{ 58, 60 } },
    .{ .name = "ad_indel_start", .rf = .{ 41, 60 }, .ff = .{ 27, 60 } },
    .{ .name = "gap_start_del_cfill", .rf = .{ 41, 60 }, .ff = .{ 27, 60 } },
    .{ .name = "gap_start_ins_cfill", .rf = null, .ff = null },
    .{ .name = "gap_internal_del_cfill", .rf = .{ 72, 60 }, .ff = .{ 58, 60 } },
    .{ .name = "gap_internal_ins_cfill", .rf = .{ 74, 60 }, .ff = .{ 60, 60 } },
    .{ .name = "ad_N_in_adapter", .rf = .{ 30, 60 }, .ff = .{ 30, 60 } },
    .{ .name = "ad_r2_only", .rf = .{ 60, 30 }, .ff = .{ 60, 30 } },
    .{ .name = "ad_both_30", .rf = .{ 30, 30 }, .ff = .{ 30, 30 } },
    .{ .name = "ad_unequal_60_40", .rf = .{ 45, 40 }, .ff = .{ 45, 40 } },
    .{ .name = "dimer_both2", .rf = null, .ff = null },
    .{ .name = "dimer_both3", .rf = .{ 3, 3 }, .ff = .{ 3, 3 } },
    .{ .name = "dimer_r1_only_r2_60", .rf = .{ 2, 60 }, .ff = .{ 2, 60 } },
    .{ .name = "dimer_r1_trim_r2_short", .rf = null, .ff = null },
    .{ .name = "short_no_adapter_2_2", .rf = null, .ff = null },
    .{ .name = "qt_all_good_12", .rf = .{ 12, 12 }, .ff = .{ 12, 12 } },
    .{ .name = "qt_mean20_last_window", .rf = .{ 12, 12 }, .ff = .{ 12, 12 } },
    .{ .name = "qt_mean19_75_window", .rf = .{ 12, 12 }, .ff = .{ 8, 12 } },
    .{ .name = "qt_b20_window", .rf = .{ 14, 14 }, .ff = .{ 9, 14 } },
    .{ .name = "qt_b19_window", .rf = .{ 14, 14 }, .ff = .{ 8, 14 } },
    .{ .name = "qt_all_low_12", .rf = .{ 12, 12 }, .ff = .{ 1, 12 } },
    .{ .name = "qt_len4", .rf = .{ 4, 4 }, .ff = null },
    .{ .name = "qt_len5", .rf = .{ 5, 5 }, .ff = .{ 5, 5 } },
    .{ .name = "qt_trailing_N", .rf = .{ 22, 22 }, .ff = .{ 20, 22 } },
    .{ .name = "qt_middle_N", .rf = .{ 21, 21 }, .ff = .{ 21, 21 } },
    .{ .name = "qt_only_N", .rf = .{ 6, 6 }, .ff = null },
    .{ .name = "pg_g9_after_20", .rf = .{ 29, 30 }, .ff = .{ 20, 30 } },
    .{ .name = "pg_g10_after_20", .rf = .{ 30, 30 }, .ff = .{ 20, 30 } },
    .{ .name = "pg_g11_after_20", .rf = .{ 31, 30 }, .ff = .{ 20, 30 } },
    .{ .name = "pg_g9_only", .rf = .{ 9, 9 }, .ff = .{ 9, 9 } },
    .{ .name = "pg_g10_only", .rf = .{ 10, 10 }, .ff = null },
    .{ .name = "pg_run_mm_A", .rf = .{ 31, 31 }, .ff = .{ 20, 31 } },
    .{ .name = "pg_allG_then_A", .rf = .{ 13, 13 }, .ff = null },
    .{ .name = "pg_both_mates", .rf = .{ 30, 30 }, .ff = .{ 20, 20 } },
    .{ .name = "pg_r2_only", .rf = .{ 25, 30 }, .ff = .{ 25, 20 } },
    .{ .name = "pg_before_adapter", .rf = .{ 49, 39 }, .ff = .{ 30, 39 } },
};

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

fn mate(report: std.json.ObjectMap, population: []const u8, index: usize) std.json.ObjectMap {
    return report.get("mates").?.object.get(population).?.array.items[index].object;
}

fn trimStep(report: std.json.ObjectMap, index: usize, step: []const u8) std.json.Value {
    return report.get("trimming").?.array.items[index].object.get(step).?;
}

fn expectTrimStep(report: std.json.ObjectMap, index: usize, step: []const u8, reads: u64, bases: u64) !void {
    const counts = trimStep(report, index, step).object;
    try cli.expectJsonInteger(counts.get("reads"), reads);
    try cli.expectJsonInteger(counts.get("bases"), bases);
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
    try cli.expectJsonInteger(mate(report, "raw", 0).get("bases"), 8);
    try cli.expectJsonInteger(mate(report, "raw", 0).get("c"), 3);
    try cli.expectJsonInteger(mate(report, "raw", 0).get("quality_sum"), 164);
    try cli.expectJsonInteger(mate(report, "raw", 1).get("n"), 2);
    try cli.expectJsonInteger(mate(report, "raw", 1).get("q20_bases"), 4);
    try cli.expectJsonInteger(mate(report, "retained", 0).get("quality_sum"), 164);
    try cli.expectJsonInteger(mate(report, "retained", 1).get("n"), 2);
    for (0..2) |index| for (TRIM_STEPS) |step| {
        try std.testing.expectEqual(std.json.Value.null, trimStep(report, index, step));
    };
    try std.testing.expectEqual(std.json.Value.null, report.get("options").?.object.get("adapter1").?);
    try std.testing.expectEqual(std.json.Value.null, report.get("options").?.object.get("min_length").?);
    try std.testing.expectEqual(std.json.Value.null, report.get("filtering").?);
    try std.testing.expectEqual(std.json.Value.null, report.get("inserts").?);
    try cli.expectJsonInteger(report.get("pairs").?.object.get("dropped"), 0);
    try std.testing.expectEqual(std.json.Value.null, report.get("error").?);
    try std.testing.expectEqual(std.json.Value.null, report.get("failed_input").?);
    try cli.expectJsonInteger(report.get("inputs").?.array.items[1].object.get("size"), R2.len);

    const qc_mates = report.get("qc").?.object.get("raw").?.array.items;
    try std.testing.expectEqual(@as(usize, 2), qc_mates.len);
    try std.testing.expectEqual(@as(usize, 2), report.get("qc").?.object.get("retained").?.array.items.len);
    const read_quality = qc_mates[0].object.get("per_sequence_quality").?.array.items;
    try std.testing.expectEqual(@as(usize, 40), read_quality.len);
    try cli.expectJsonInteger(read_quality[0].object.get("quality"), 1);
    try cli.expectJsonInteger(read_quality[39].object.get("count"), 1);
    try std.testing.expectEqual(@as(usize, 101), qc_mates[1].object.get("per_sequence_gc").?.array.items.len);
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
        .{ .args = &.{ "process", "--adapter1", "ACGT", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: --adapter1 and --adapter2 must be given together\n" },
        .{ .args = &.{ "process", "--adapter1", "ACG", p.r1, p.r2 }, .stderr = "error: --adapter1 must be 4 to 128 bases of A, C, G, and T\n" },
        .{ .args = &.{ "process", "--adapter1", "A" ** 129, p.r1, p.r2 }, .stderr = "error: --adapter1 must be 4 to 128 bases of A, C, G, and T\n" },
        .{ .args = &.{ "process", "--adapter2", "ACGN", p.r1, p.r2 }, .stderr = "error: --adapter2 must be 4 to 128 bases of A, C, G, and T\n" },
        .{ .args = &.{ "process", "--adapter2", "acgt", p.r1, p.r2 }, .stderr = "error: --adapter2 must be 4 to 128 bases of A, C, G, and T\n" },
        .{ .args = &.{ "process", "--adapter1", "ACGT", "--adapter1", "ACGT", p.r1, p.r2 }, .stderr = "error: --adapter1 may appear only once\n" },
        .{ .args = &.{ "process", "--trim-tail-quality", "20", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: --trim-tail-window and --trim-tail-quality must be given together\n" },
        .{ .args = &.{ "process", "--trim-tail-window", "0", p.r1, p.r2 }, .stderr = "error: --trim-tail-window must be 1 to 1000\n" },
        .{ .args = &.{ "process", "--trim-tail-window", "1001", p.r1, p.r2 }, .stderr = "error: --trim-tail-window must be 1 to 1000\n" },
        .{ .args = &.{ "process", "--trim-tail-quality", "0", p.r1, p.r2 }, .stderr = "error: --trim-tail-quality must be 1 to 30\n" },
        .{ .args = &.{ "process", "--trim-tail-quality", "31", p.r1, p.r2 }, .stderr = "error: --trim-tail-quality must be 1 to 30\n" },
        .{ .args = &.{ "process", "--trim-poly-g", "0", p.r1, p.r2 }, .stderr = "error: --trim-poly-g must be 1 to 4294967295\n" },
        .{ .args = &.{ "process", "--trim-poly-g", "4294967296", p.r1, p.r2 }, .stderr = "error: --trim-poly-g must be 1 to 4294967295\n" },
        .{ .args = &.{ "process", "--trim-poly-g", "1x", p.r1, p.r2 }, .stderr = "error: --trim-poly-g must be 1 to 4294967295\n" },
        .{ .args = &.{ "process", "--min-length", "4294967296", p.r1, p.r2 }, .stderr = "error: --min-length must be 0 to 4294967295\n" },
        .{ .args = &.{ "process", "--min-length", "-1", p.r1, p.r2 }, .stderr = "error: --min-length must be 0 to 4294967295\n" },
        .{ .args = &.{ "process", "--max-n", "51", p.r1, p.r2 }, .stderr = "error: --max-n must be 0 to 50\n" },
        .{ .args = &.{ "process", "--low-quality", "94", "--max-low-quality-percent", "40", p.r1, p.r2 }, .stderr = "error: --low-quality must be 0 to 93\n" },
        .{ .args = &.{ "process", "--max-low-quality-percent", "101", p.r1, p.r2 }, .stderr = "error: --max-low-quality-percent must be 0 to 100\n" },
        .{ .args = &.{ "process", "--low-quality", "15", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: --low-quality and --max-low-quality-percent must be given together\n" },
        .{ .args = &.{ "process", "--max-n", "5", "--max-n", "5", p.r1, p.r2 }, .stderr = "error: --max-n may appear only once\n" },
        .{ .args = &.{ "process", "--contaminants", p.r1, "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: --contaminants requires --report\n" },
        .{ .args = &.{ "process", "--pair-duplication", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 }, .stderr = "error: --pair-duplication requires --report\n" },
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

test "[cli] - [process]: trimming writes the records fastp 1.1.0 writes on the probe pairs" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);

    const Family = struct { name: []const u8, options: []const []const u8 };
    const families = [_]Family{
        .{ .name = "rf", .options = &.{ "--adapter1", RF1, "--adapter2", RF2 } },
        .{ .name = "ff", .options = &.{ "--adapter1", FF, "--adapter2", FF, "--trim-poly-g", "10", "--trim-tail-window", "4", "--trim-tail-quality", "20" } },
    };
    for (families, 0..) |family, family_index| {
        const inputs = [2][]const u8{
            try std.fmt.allocPrint(allocator, "{s}/trim_{s}_1.fastq", .{ cli.FIXTURE_DIR, family.name }),
            try std.fmt.allocPrint(allocator, "{s}/trim_{s}_2.fastq", .{ cli.FIXTURE_DIR, family.name }),
        };
        const outputs = [2][]const u8{
            try std.fmt.allocPrint(allocator, "{s}.{s}", .{ p.out1, family.name }),
            try std.fmt.allocPrint(allocator, "{s}.{s}", .{ p.out2, family.name }),
        };
        const args = try std.mem.concat(allocator, []const u8, &.{
            &.{ "process", "--out1", outputs[0], "--out2", outputs[1] },
            family.options,
            &inputs,
        });
        try cli.expectResult(try cli.run(allocator, args), 0, "", "");
        for (inputs, outputs, 0..) |input, output, mate_index| {
            var source = std.mem.splitScalar(u8, try readPath(allocator, input), '\n');
            var written = std.mem.splitScalar(u8, try readPath(allocator, output), '\n');
            for (PROBES) |probe| {
                errdefer std.debug.print("{s} {s} mate {d}\n", .{ family.name, probe.name, mate_index + 1 });
                const lengths = (if (family_index == 0) probe.rf else probe.ff) orelse {
                    for (0..4) |_| _ = source.next().?;
                    continue;
                };
                const kept = lengths[mate_index];
                const header = source.next().?;
                try std.testing.expectEqualStrings(try std.fmt.allocPrint(allocator, "@{s}/{d}", .{ probe.name, mate_index + 1 }), header);
                try std.testing.expectEqualStrings(header, written.next().?);
                try std.testing.expectEqualStrings(source.next().?[0..kept], written.next().?);
                try std.testing.expectEqualStrings(source.next().?, written.next().?);
                try std.testing.expectEqualStrings(source.next().?[0..kept], written.next().?);
            }
            try std.testing.expectEqualStrings("", written.next().?);
            try std.testing.expectEqual(@as(?[]const u8, null), written.next());
        }
    }
}

test "[cli] - [process]: option range boundaries are accepted, and a pair with an emptied mate is dropped" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, R1, R2);

    const Case = struct { options: []const []const u8, r1: []const u8, r2: []const u8 };
    const cases = [_]Case{
        .{
            .options = &.{ "--adapter1", "ACGT", "--adapter2", "ACGT", "--trim-poly-g", "1", "--trim-tail-window", "1", "--trim-tail-quality", "1" },
            .r1 = "@p1/1 desc\nACGT\n+\nIIII\n@p2/1\nGGC\n+p2\n#!#\n",
            .r2 = "@p1/2 desc\nTTTT\n+\nIIII\n@p2/2\nA\n+\n!\n",
        },
        .{
            .options = &.{ "--adapter1", "A" ** 128, "--adapter2", "C" ** 128, "--trim-poly-g", "4294967295", "--trim-tail-window", "1000", "--trim-tail-quality", "30" },
            .r1 = "",
            .r2 = "",
        },
    };
    for (cases, 0..) |case, index| {
        errdefer std.debug.print("case {d}\n", .{index});
        const out1 = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.out1, index });
        const out2 = try std.fmt.allocPrint(allocator, "{s}.{d}", .{ p.out2, index });
        const args = try std.mem.concat(allocator, []const u8, &.{
            &.{ "process", "--report", "-", "--out1", out1, "--out2", out2 },
            case.options,
            &.{ p.r1, p.r2 },
        });
        const result = try cli.run(allocator, args);
        try std.testing.expectEqual(@as(u8, 0), result.exit_code);
        try std.testing.expectEqualStrings(case.r1, try readPath(allocator, out1));
        try std.testing.expectEqualStrings(case.r2, try readPath(allocator, out2));
        const report = try parseReport(allocator, result.stdout);
        const options = report.get("options").?.object;
        try cli.expectJsonString(options.get("adapter1"), case.options[1]);
        try cli.expectJsonString(options.get("adapter2"), case.options[3]);
        try cli.expectJsonInteger(options.get("trim_poly_g"), try std.fmt.parseInt(u64, case.options[5], 10));
        try cli.expectJsonInteger(options.get("trim_tail_window"), try std.fmt.parseInt(u64, case.options[7], 10));
        try cli.expectJsonInteger(options.get("trim_tail_quality"), try std.fmt.parseInt(u64, case.options[9], 10));
    }

    const report = try parseReport(allocator, (try cli.run(allocator, &.{
        "process", "--report", "-", "--out1", p.out1, "--out2", p.out2, "--trim-tail-window", "1000", "--trim-tail-quality", "30", p.r1, p.r2,
    })).stdout);
    try cli.expectJsonInteger(report.get("pairs").?.object.get("written"), 0);
    try cli.expectJsonInteger(report.get("pairs").?.object.get("dropped"), 2);
    try cli.expectJsonInteger(report.get("filtering").?.object.get("dropped_pairs").?.object.get("too_short"), 2);
    for (0..2) |index| {
        try expectTrimStep(report, index, "quality_tail", 2, 8);
        try cli.expectJsonInteger(mate(report, "raw", index).get("reads"), 2);
        try cli.expectJsonInteger(mate(report, "retained", index).get("reads"), 0);
        try cli.expectJsonInteger(mate(report, "retained", index).get("bases"), 0);
        try cli.expectJsonInteger(report.get("filtering").?.object.get("failed_reads").?.array.items[index].object.get("too_short"), 2);
        try std.testing.expectEqual(std.json.Value.null, trimStep(report, index, "poly_g"));
        try std.testing.expectEqual(std.json.Value.null, trimStep(report, index, "adapter"));
    }
}

test "[cli] - [process]: a trimmed record keeps its header and plus field and the report counts by step" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const insert = "ACCTTAATGGGCGTGGACGTCTCTTGCTTG";
    const r1 = "@a/1 x\r\n" ++ insert ++ FF ++ "\r\n+a\r\n" ++ "I" ** 49 ++ "\r\n@b/1\r\nACGTACGT\r\n+\r\nIIIIIIII\r\n";
    const r2 = "@a/2 x\nCGATTCAAATGACGGCAGCAGGCCGGGAGTCCCTGAGAGG\n+\n" ++ "I" ** 40 ++ "\n@b/2\nTTTT\n+\nIIII\n";
    try writeInputs(&tmp, r1, r2);

    const expected1 = "@a/1 x\n" ++ insert ++ "\n+a\n" ++ "I" ** 30 ++ "\n@b/1\nACGTACGT\n+\nIIIIIIII\n";
    const result = try cli.run(allocator, &.{ "process", "--report", "-", "--adapter1", FF, "--adapter2", FF, "--trim-poly-g", "10", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    try std.testing.expectEqualStrings(expected1, try readPath(allocator, p.out1));
    try std.testing.expectEqualStrings(r2, try readPath(allocator, p.out2));

    const gz1 = try cli.tempPath(allocator, &tmp.sub_path, "out1.fastq.gz");
    const gz2 = try cli.tempPath(allocator, &tmp.sub_path, "out2.fastq.gz");
    const gzip_result = try cli.run(allocator, &.{ "process", "--gzip", "--adapter1", FF, "--adapter2", FF, "--out1", gz1, "--out2", gz2, p.r1, p.r2 });
    try cli.expectResult(gzip_result, 0, "", "");
    try std.testing.expectEqualStrings(expected1, try gunzipPath(allocator, gz1));
    try std.testing.expectEqualStrings(r2, try gunzipPath(allocator, gz2));

    const report = try parseReport(allocator, result.stdout);
    try std.testing.expectEqual(std.json.Value.null, trimStep(report, 0, "quality_tail"));
    try expectTrimStep(report, 0, "poly_g", 0, 0);
    try expectTrimStep(report, 0, "overlap", 0, 0);
    try expectTrimStep(report, 0, "adapter", 1, FF.len);
    try expectTrimStep(report, 1, "adapter", 0, 0);
    try cli.expectJsonInteger(mate(report, "raw", 0).get("bases"), 57);
    try cli.expectJsonInteger(mate(report, "retained", 0).get("bases"), 38);
    try cli.expectJsonInteger(mate(report, "retained", 0).get("max_length"), 30);
    try cli.expectJsonInteger(mate(report, "retained", 1).get("bases"), 44);
    const lengths = report.get("qc").?.object.get("retained").?.array.items[0].object.get("length_distribution").?.array.items;
    try cli.expectJsonInteger(lengths[lengths.len - 1].object.get("count"), 1);
    try cli.expectJsonString(lengths[lengths.len - 1].object.get("length"), "30");
}

test "[cli] - [process]: filters drop whole pairs, count each reason, and keep the exit status" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const kept1 = "@ok/1\nACGTACGT\n+\nIIIIIIII\n";
    const kept2 = "@ok/2\nTTGGCCAA\n+\nIIIIIIII\n";
    const r1 = "@n/1\nANNACGTA\n+\nIIIIIIII\n" ++ kept1 ++ "@short/1\nACG\n+\nIII\n@low/1\nACGTACGT\n+\nIIIIIIII\n";
    const r2 = "@n/2\nACGTACGT\n+\nIIIIIIII\n" ++ kept2 ++ "@short/2\nACGTACGT\n+\nIIIIIIII\n@low/2\nACGTNACG\n+\nII!!!!II\n";
    try writeInputs(&tmp, r1, r2);

    const options = [_][]const u8{ "--min-length", "4", "--max-n", "1", "--low-quality", "15", "--max-low-quality-percent", "40" };
    const args = try std.mem.concat(allocator, []const u8, &.{ &.{ "process", "--report", "-", "--out1", p.out1, "--out2", p.out2 }, &options, &.{ p.r1, p.r2 } });
    const result = try cli.run(allocator, args);
    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    try std.testing.expectEqualStrings(kept1, try readPath(allocator, p.out1));
    try std.testing.expectEqualStrings(kept2, try readPath(allocator, p.out2));

    const report = try parseReport(allocator, result.stdout);
    const pairs = report.get("pairs").?.object;
    try cli.expectJsonInteger(pairs.get("read"), 4);
    try cli.expectJsonInteger(pairs.get("written"), 1);
    try cli.expectJsonInteger(pairs.get("dropped"), 3);
    const filtering = report.get("filtering").?.object;
    const dropped = filtering.get("dropped_pairs").?.object;
    try cli.expectJsonInteger(dropped.get("too_many_n"), 1);
    try cli.expectJsonInteger(dropped.get("too_short"), 1);
    try cli.expectJsonInteger(dropped.get("low_quality"), 1);
    try cli.expectJsonInteger(dropped.get("adapter_dimer"), 0);
    const failed = filtering.get("failed_reads").?.array.items;
    try cli.expectJsonInteger(failed[0].object.get("too_many_n"), 1);
    try cli.expectJsonInteger(failed[0].object.get("too_short"), 1);
    try cli.expectJsonInteger(failed[0].object.get("low_quality"), 0);
    try cli.expectJsonInteger(failed[1].object.get("low_quality"), 1);
    try cli.expectJsonInteger(failed[1].object.get("too_short"), 0);
    try cli.expectJsonInteger(mate(report, "raw", 0).get("reads"), 4);
    try cli.expectJsonInteger(mate(report, "retained", 0).get("reads"), 1);
    try cli.expectJsonInteger(mate(report, "retained", 1).get("bases"), 8);
    try cli.expectJsonInteger(mate(report, "retained", 1).get("n"), 0);
    try std.testing.expectEqual(std.json.Value.null, trimStep(report, 0, "adapter"));
    const options_json = report.get("options").?.object;
    try cli.expectJsonInteger(options_json.get("min_length"), 4);
    try cli.expectJsonInteger(options_json.get("max_n"), 1);
    try cli.expectJsonInteger(options_json.get("low_quality"), 15);
    try cli.expectJsonInteger(options_json.get("max_low_quality_percent"), 40);
    const retained_lengths = report.get("qc").?.object.get("retained").?.array.items[0].object.get("length_distribution").?.array.items;
    try std.testing.expectEqual(@as(usize, 1), retained_lengths.len);
    try cli.expectJsonInteger(retained_lengths[0].object.get("count"), 1);
}

test "[cli] - [process]: the report counts inserts from the mate overlap and adapter content per population" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const quality = "I" ** 50;
    const r1 = "@overlap/1\nACCTTAATGGGCGTGGACGTCTCTTGCTTGGAGGCTGGTGAAAGCGGCAC\n+\n" ++ quality ++ "\n@apart/1\nAGATCGGAAGAGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCC\n+\n" ++ quality ++ "\n";
    const r2 = "@overlap/2\nCACCAGCCTCCAAGCAAGAGACGTCCACGCCCATTAAGGTAAAGCGGCAC\n+\n" ++ quality ++ "\n@apart/2\nGGGGGGGGGGGGGGGGACACACACACACACACACACACACACACACACAC\n+\n" ++ quality ++ "\n";
    try writeInputs(&tmp, r1, r2);

    const result = try cli.run(allocator, &.{ "process", "--report", "-", "--adapter1", RF1, "--adapter2", RF2, "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    const report = try parseReport(allocator, result.stdout);
    const inserts = report.get("inserts").?.object;
    try cli.expectJsonInteger(inserts.get("pairs"), 2);
    try cli.expectJsonInteger(inserts.get("overlapping"), 1);
    try cli.expectJsonInteger(inserts.get("unknown"), 1);
    try cli.expectJsonInteger(inserts.get("peak"), 40);
    try cli.expectJsonInteger(inserts.get("median"), 40);
    const histogram = inserts.get("histogram").?.array.items;
    try std.testing.expectEqual(@as(usize, 512), histogram.len);
    try cli.expectJsonInteger(histogram[40], 1);

    const qc = report.get("qc").?.object;
    for ([_][]const u8{ "raw", "retained" }) |population| {
        const content = qc.get(population).?.array.items[0].object.get("adapter_content").?.object;
        try std.testing.expectEqual(@as(usize, 6), content.get("adapters").?.array.items.len);
        try cli.expectJsonString(content.get("adapters").?.array.items[0], "Illumina Universal Adapter");
        const rows = content.get("rows").?.array.items;
        try std.testing.expect(rows.len > 0);
        try std.testing.expectEqual(@as(?f64, 50), numberValue(rows[0].object.get("values").?.array.items[0]));
    }
    const raw_mate2 = qc.get("raw").?.array.items[1].object.get("adapter_content").?.object.get("rows").?.array.items;
    try std.testing.expectEqual(@as(?f64, 50), numberValue(raw_mate2[0].object.get("values").?.array.items[5]));
}

const SE_ADAPTER_1 = "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG";

test "[cli] - [process]: the report counts duplicated, overrepresented, and duplicate pairs" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const quality = "I" ** 33;
    const other = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAA";
    const mate2 = "TTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGG";
    const r1 = "@a/1\n" ++ SE_ADAPTER_1 ++ "\n+\n" ++ quality ++ "\n@b/1\n" ++ SE_ADAPTER_1 ++ "\n+\n" ++ quality ++ "\n@c/1\n" ++ other ++ "\n+\n" ++ quality ++ "\n";
    const r2 = "@a/2\n" ++ mate2 ++ "\n+\n" ++ quality ++ "\n@b/2\n" ++ mate2 ++ "\n+\n" ++ quality ++ "\n@c/2\n" ++ other ++ "\n+\n" ++ quality ++ "\n";
    try writeInputs(&tmp, r1, r2);
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "contaminants.txt", .data = "# list\nIllumina Single End Adapter 1\t\t" ++ SE_ADAPTER_1 ++ "\n" });
    const contaminants = try cli.tempPath(allocator, &tmp.sub_path, "contaminants.txt");

    const result = try cli.run(allocator, &.{ "process", "--report", "-", "--contaminants", contaminants, "--pair-duplication", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
    try std.testing.expectEqual(@as(u8, 0), result.exit_code);
    const report = try parseReport(allocator, result.stdout);
    try cli.expectJsonString(report.get("options").?.object.get("contaminants"), contaminants);
    try std.testing.expectEqual(std.json.Value{ .bool = true }, report.get("options").?.object.get("pair_duplication").?);
    const pairs = report.get("pair_duplication").?.object;
    try cli.expectJsonInteger(pairs.get("pairs"), 3);
    try cli.expectJsonInteger(pairs.get("duplicates"), 1);

    const mate1 = report.get("qc").?.object.get("raw").?.array.items[0].object;
    const duplication = mate1.get("duplication").?.object;
    try cli.expectJsonInteger(duplication.get("reads"), 3);
    try cli.expectJsonInteger(duplication.get("stored_sequences"), 2);
    try std.testing.expectEqual(@as(?f64, @as(f64, 2) / 3 * 100), numberValue(duplication.get("total_deduplicated_percentage").?));
    try std.testing.expectEqual(@as(usize, 16), duplication.get("levels").?.array.items.len);
    const rows = mate1.get("overrepresented").?.array.items;
    try std.testing.expectEqual(@as(usize, 2), rows.len);
    try cli.expectJsonString(rows[0].object.get("sequence"), SE_ADAPTER_1);
    try cli.expectJsonInteger(rows[0].object.get("count"), 2);
    const source = rows[0].object.get("possible_source").?.object;
    try cli.expectJsonString(source.get("name"), "Illumina Single End Adapter 1");
    try cli.expectJsonInteger(source.get("identity"), 100);
    try cli.expectJsonInteger(source.get("length"), 33);
    try std.testing.expectEqual(std.json.Value.null, rows[1].object.get("possible_source").?);

    const plain1 = try cli.tempPath(allocator, &tmp.sub_path, "plain1.fastq");
    const plain2 = try cli.tempPath(allocator, &tmp.sub_path, "plain2.fastq");
    const plain = try cli.run(allocator, &.{ "process", "--report", "-", "--out1", plain1, "--out2", plain2, p.r1, p.r2 });
    try std.testing.expectEqual(@as(u8, 0), plain.exit_code);
    const plain_report = try parseReport(allocator, plain.stdout);
    try std.testing.expectEqual(std.json.Value.null, plain_report.get("pair_duplication").?);
    const plain_rows = plain_report.get("qc").?.object.get("raw").?.array.items[0].object.get("overrepresented").?.array.items;
    try std.testing.expectEqual(std.json.Value.null, plain_rows[0].object.get("possible_source").?);
}

test "[cli] - [process]: an unreadable or invalid contaminant list fails before outputs exist" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    try writeInputs(&tmp, "@a/1\nACGT\n+\nIIII\n", "@a/2\nACGT\n+\nIIII\n");
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "bad.txt", .data = "good\tACGT\nbad\tACNT\n" });
    const Case = struct { file: []const u8, exit_code: u8, message: []const u8 };
    for ([_]Case{
        .{ .file = "bad.txt", .exit_code = 2, .message = "line 2: contaminant sequence has a base other than A, C, G or T" },
        .{ .file = "missing.txt", .exit_code = 3, .message = "file not found" },
    }) |case| {
        const path = try cli.tempPath(allocator, &tmp.sub_path, case.file);
        const result = try cli.run(allocator, &.{ "process", "--report", p.report, "--contaminants", path, "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
        try std.testing.expectEqual(case.exit_code, result.exit_code);
        try std.testing.expect(std.mem.endsWith(u8, std.mem.trimEnd(u8, result.stderr, "\n"), case.message));
        try std.testing.expect(!pathExists(p.out1));
        try std.testing.expect(!pathExists(p.out2));
        try std.testing.expect(!pathExists(p.report));
    }
}

fn numberValue(value: std.json.Value) ?f64 {
    return switch (value) {
        .integer => |integer| @floatFromInt(integer),
        .float => |float| float,
        else => null,
    };
}

test "[cli] - [process]: a malformed record after dropped pairs still fails the run" {
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const p = try Paths.init(allocator, &tmp);
    const r1 = "@a/1\nAC\n+\nII\n@b/1\nACGTACGT\n+\nIIIIIIII\n@c/1\nACGT\n+\nIII\n";
    const r2 = "@a/2\nACGTACGT\n+\nIIIIIIII\n@b/2\nACGTACGT\n+\nIIIIIIII\n@c/2\nACGT\n+\nIIII\n";
    try writeInputs(&tmp, r1, r2);

    const result = try cli.run(allocator, &.{ "process", "--min-length", "4", "--report", "-", "--out1", p.out1, "--out2", p.out2, p.r1, p.r2 });
    try std.testing.expectEqual(@as(u8, 1), result.exit_code);
    try std.testing.expectEqualStrings("@b/1\nACGTACGT\n+\nIIIIIIII\n", try readPath(allocator, p.out1));
    try std.testing.expectEqualStrings("@b/2\nACGTACGT\n+\nIIIIIIII\n", try readPath(allocator, p.out2));
    const report = try parseReport(allocator, result.stdout);
    try cli.expectJsonString(report.get("status"), "error");
    try cli.expectJsonString(report.get("error").?.object.get("code"), "S005");
    try cli.expectJsonInteger(report.get("pairs").?.object.get("dropped"), 1);
    try cli.expectJsonInteger(report.get("pairs").?.object.get("written"), 1);
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

    const with_adapter = try std.fmt.allocPrint(allocator, "@short/1\nAC\n+\nII\n@long/1\n{s}{s}\n+\n{s}{s}\n", .{ sequence, FF, quality, "I" ** FF.len });
    try writeInputs(&tmp, with_adapter, r2);
    const trimmed = try cli.run(allocator, &.{ "process", "--adapter1", FF, "--adapter2", FF, "--out1", try std.fmt.allocPrint(allocator, "{s}.trim", .{p.out1}), "--out2", try std.fmt.allocPrint(allocator, "{s}.trim", .{p.out2}), p.r1, p.r2 });
    try cli.expectResult(trimmed, 0, "", "");
    try std.testing.expectEqualStrings(r1, try readPath(allocator, try std.fmt.allocPrint(allocator, "{s}.trim", .{p.out1})));
    try std.testing.expectEqualStrings(r2, try readPath(allocator, try std.fmt.allocPrint(allocator, "{s}.trim", .{p.out2})));
}
