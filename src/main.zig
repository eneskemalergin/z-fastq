//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const build_options = @import("build_options");
const zfastq = @import("root.zig");
const fastq = @import("fastq.zig");
const io_layer = @import("io.zig");
const pairing = @import("pair.zig");
const sampling = @import("sample.zig");

const PAIR_DIAGNOSTIC_PREFIX_BYTES = 128;

const USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count         Count records in plain or gzip FASTQ inputs
    \\  stats         Report aggregate FASTQ statistics
    \\  check         Validate FASTQ structure, sequence alphabet, and quality range
    \\  sample        Select records by deterministic probability or exact count
    \\  interleave    Validate and interleave paired FASTQ inputs
    \\  deinterleave  Validate and separate interleaved paired FASTQ input
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
;

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;

    var args = std.process.Args.Iterator.initAllocator(init.minimal.args, gpa) catch {
        argumentError(io, "error: out of memory\n", 3);
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

    const command: Command = if (std.mem.eql(u8, cmd, "count"))
        .count
    else if (std.mem.eql(u8, cmd, "stats"))
        .stats
    else if (std.mem.eql(u8, cmd, "check"))
        .check
    else if (std.mem.eql(u8, cmd, "sample"))
        .sample
    else if (std.mem.eql(u8, cmd, "interleave"))
        .interleave
    else if (std.mem.eql(u8, cmd, "deinterleave"))
        .deinterleave
    else {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown command: ") catch {};
        writeEscaped(.stderr(), io, cmd);
        std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
        printUsageAndExit(io);
    };

    var max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES;
    var alphabet: zfastq.Alphabet = .iupac;
    var pair_mode: PairMode = .none;
    var pair_name_policy: pairing.NamePolicy = .illumina;
    var pair_names_set = false;
    var json_output = false;
    var fraction: ?sampling.Fraction = null;
    var sample_count: ?u64 = null;
    var sample_seed: u64 = 11;
    var output1: ?[]const u8 = null;
    var output2: ?[]const u8 = null;
    var positional = std.ArrayList([]const u8).empty;
    defer positional.deinit(gpa);

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--")) {
            while (args.next()) |path| {
                positional.append(gpa, path) catch {
                    argumentError(io, "error: out of memory\n", 3);
                };
            }
            break;
        }
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        }
        if (std.mem.eql(u8, arg, "--max-line-bytes")) {
            const value = args.next() orelse {
                argumentError(io, "error: --max-line-bytes requires a value\n", 2);
            };
            max_line_bytes = parseMaxLineBytes(value) catch |err| switch (err) {
                error.Overflow => {
                    argumentError(io, "error: --max-line-bytes exceeds supported limit\n", 4);
                },
                error.InvalidCharacter => {
                    argumentError(io, "error: invalid --max-line-bytes value\n", 2);
                },
            };
            continue;
        }
        if ((command == .stats or command == .check) and std.mem.eql(u8, arg, "--json")) {
            json_output = true;
            continue;
        }
        if ((command == .check or command == .sample or command == .interleave or
            command == .deinterleave) and
            std.mem.eql(u8, arg, "--alphabet"))
        {
            const value = args.next() orelse {
                argumentError(io, "error: --alphabet requires a value\n", 2);
            };
            alphabet = if (std.mem.eql(u8, value, "iupac"))
                .iupac
            else if (std.mem.eql(u8, value, "acgtn"))
                .acgtn
            else {
                argumentError(io, "error: --alphabet must be iupac or acgtn\n", 2);
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--fraction")) {
            const value = args.next() orelse {
                argumentError(io, "error: --fraction requires a value\n", 2);
            };
            fraction = sampling.Fraction.parse(value) catch {
                argumentError(io, "error: invalid --fraction value\n", 2);
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--count")) {
            const value = args.next() orelse {
                argumentError(io, "error: --count requires a value\n", 2);
            };
            sample_count = sampling.parseCount(value) catch |err| switch (err) {
                error.InvalidCount => {
                    argumentError(io, "error: invalid --count value\n", 2);
                },
                error.Overflow => {
                    argumentError(io, "error: --count exceeds supported limit\n", 4);
                },
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--seed")) {
            const value = args.next() orelse {
                argumentError(io, "error: --seed requires a value\n", 2);
            };
            sample_seed = sampling.parseSeed(value) catch |err| switch (err) {
                error.InvalidSeed => {
                    argumentError(io, "error: invalid --seed value\n", 2);
                },
                error.Overflow => {
                    argumentError(io, "error: --seed exceeds supported limit\n", 4);
                },
            };
            continue;
        }
        if ((command == .check or command == .sample) and
            std.mem.eql(u8, arg, "--paired"))
        {
            if (pair_mode == .interleaved) {
                argumentError(io, "error: --paired and --interleaved are mutually exclusive\n", 2);
            }
            pair_mode = .paired;
            continue;
        }
        if ((command == .check or command == .sample) and
            std.mem.eql(u8, arg, "--interleaved"))
        {
            if (pair_mode == .paired) {
                argumentError(io, "error: --paired and --interleaved are mutually exclusive\n", 2);
            }
            pair_mode = .interleaved;
            continue;
        }
        if ((command == .check or command == .sample or command == .interleave or
            command == .deinterleave) and
            std.mem.eql(u8, arg, "--pair-names"))
        {
            const value = args.next() orelse {
                argumentError(io, "error: --pair-names requires a value\n", 2);
            };
            pair_name_policy = if (std.mem.eql(u8, value, "illumina"))
                .illumina
            else if (std.mem.eql(u8, value, "exact"))
                .exact
            else {
                argumentError(io, "error: --pair-names must be illumina or exact\n", 2);
            };
            pair_names_set = true;
            continue;
        }
        if (command == .deinterleave and
            (std.mem.eql(u8, arg, "--out1") or std.mem.eql(u8, arg, "--out2")))
        {
            const value = args.next() orelse {
                const message = if (std.mem.eql(u8, arg, "--out1"))
                    "error: --out1 requires a value\n"
                else
                    "error: --out2 requires a value\n";
                argumentError(io, message, 2);
            };
            const destination = if (std.mem.eql(u8, arg, "--out1")) &output1 else &output2;
            if (destination.* != null) {
                const message = if (std.mem.eql(u8, arg, "--out1"))
                    "error: --out1 may appear only once\n"
                else
                    "error: --out2 may appear only once\n";
                argumentError(io, message, 2);
            }
            destination.* = value;
            continue;
        }
        if (arg.len > 1 and std.mem.startsWith(u8, arg, "-")) {
            std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown ") catch {};
            std.Io.File.writeStreamingAll(.stderr(), io, @tagName(command)) catch {};
            std.Io.File.writeStreamingAll(.stderr(), io, " option: ") catch {};
            writeEscaped(.stderr(), io, arg);
            argumentError(io, "\n", 2);
        }
        positional.append(gpa, arg) catch {
            argumentError(io, "error: out of memory\n", 3);
        };
    }

    const options = InputOptions{ .max_line_bytes = max_line_bytes };
    const code = switch (command) {
        .count => runCount(io, positional.items, options),
        .stats => runStats(io, gpa, positional.items, options, json_output),
        .check => if (pair_mode == .none)
            runSingleCheckCommand(io, positional.items, .{
                .max_line_bytes = max_line_bytes,
                .alphabet = alphabet,
            }, json_output, pair_names_set)
        else
            runPairedCheckCommand(io, gpa, positional.items, .{
                .max_line_bytes = max_line_bytes,
                .alphabet = alphabet,
                .pair_mode = pair_mode,
                .pair_name_policy = pair_name_policy,
            }, json_output),
        .sample => runSample(
            io,
            gpa,
            positional.items,
            .{
                .max_line_bytes = max_line_bytes,
                .alphabet = alphabet,
                .fraction = fraction,
                .count = sample_count,
                .seed = sample_seed,
                .pair_mode = pair_mode,
                .pair_name_policy = pair_name_policy,
                .pair_names_set = pair_names_set,
            },
        ),
        .interleave => runInterleave(io, gpa, positional.items, .{
            .max_line_bytes = max_line_bytes,
            .alphabet = alphabet,
            .pair_name_policy = pair_name_policy,
        }),
        .deinterleave => runDeinterleave(io, gpa, positional.items, output1, output2, .{
            .max_line_bytes = max_line_bytes,
            .alphabet = alphabet,
            .pair_name_policy = pair_name_policy,
        }),
    };
    std.process.exit(code);
}

const Command = enum { count, stats, check, sample, interleave, deinterleave };

fn parseMaxLineBytes(value: []const u8) std.fmt.ParseIntError!usize {
    if (value.len == 0) return error.InvalidCharacter;
    for (value) |byte| {
        if (byte < '0' or byte > '9') return error.InvalidCharacter;
    }
    return std.fmt.parseInt(usize, value, 10);
}

fn argumentError(io: std.Io, message: []const u8, exit_code: u8) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
    std.process.exit(exit_code);
}

fn printUsageAndExit(io: std.Io) noreturn {
    argumentError(io, USAGE, 2);
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

const InputOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
};

fn runCount(
    io: std.Io,
    inputs: []const []const u8,
    options: InputOptions,
) u8 {
    if (validateRecordCommandInputs(io, .count, inputs)) |exit_code| return exit_code;

    var exit_code: u8 = 0;
    for (inputs) |input| {
        const count = switch (countInput(io, input, options)) {
            .success => |count| count,
            .failure => |failure| {
                if (std.mem.eql(u8, failure.code, "out_of_memory")) {
                    std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                    return 3;
                }
                printCommandFailure(io, input, failure);
                exit_code = @max(exit_code, failure.exit_code);
                continue;
            },
        };
        printCount(io, count) catch return @max(exit_code, 3);
    }
    return exit_code;
}

const CountOutcome = union(enum) {
    success: u64,
    failure: CommandFailure,
};

fn countInput(io: std.Io, label: []const u8, options: InputOptions) CountOutcome {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, label, null)) |failure| return .{ .failure = failure };
    defer input.deinit(io);

    var scanner = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = options.max_line_bytes });
    var buffer: [io_layer.COUNT_DECOMPRESS_BUFFER_BYTES]u8 = undefined;
    const read_buffer = if (input.source == .plain)
        buffer[0..zfastq.limits.COUNT_READ_BUFFER_BYTES]
    else
        &buffer;
    while (input.readScannerChunk(read_buffer) catch return .{ .failure = IO_FAILURE }) |decoded| {
        _ = scanner.feed(decoded) catch |err| return .{ .failure = mapScanFailure(&scanner, err) };
    }
    scanner.finishEof() catch |err| return .{ .failure = mapScanFailure(&scanner, err) };
    return .{ .success = scanner.record_index };
}

fn mapScanFailure(scanner: *zfastq.count_scan.Scanner, err: zfastq.ReaderError) CommandFailure {
    return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => if (scanner.takeLastError()) |details| CommandFailure.lint(details) else .{
            .code = "format_error",
            .message = "",
            .exit_code = 1,
            .suppress_diagnostic = true,
        },
        error.Io => blk: {
            var failure = IO_FAILURE;
            failure.suppress_diagnostic = true;
            break :blk failure;
        },
        else => mapInputFailure(@errorCast(err)),
    };
}

fn printCount(io: std.Io, n: u64) !void {
    var buf: [32]u8 = undefined;
    const text = try std.fmt.bufPrint(&buf, "{d}\n", .{n});
    try std.Io.File.writeStreamingAll(.stdout(), io, text);
}

const CommandFailure = struct {
    code: []const u8,
    message: []const u8,
    exit_code: u8,
    record_index: ?u64 = null,
    byte_offset: ?u64 = null,
    line_in_record: ?u3 = null,
    suppress_diagnostic: bool = false,

    fn lint(details: zfastq.ParseError) CommandFailure {
        return .{
            .code = zfastq.codeTag(details.code),
            .message = details.message,
            .exit_code = 1,
            .record_index = details.record_index,
            .byte_offset = details.byte_offset,
            .line_in_record = details.line_in_record,
        };
    }

    fn plain(code: []const u8, message: []const u8, exit_code: u8) CommandFailure {
        return .{ .code = code, .message = message, .exit_code = exit_code };
    }
};

const OUT_OF_MEMORY = CommandFailure.plain("out_of_memory", "out of memory", 3);
const INPUT_LOCATION_LIMIT = CommandFailure.plain("arithmetic_limit", "input location exceeds supported limit", 4);
const RECORD_STAGING_LIMIT = CommandFailure.plain("arithmetic_limit", "record staging size exceeds supported limit", 4);
const IO_FAILURE = CommandFailure.plain("io_error", "I/O error", 3);
const LINE_LIMIT = CommandFailure.plain("line_limit", "line length limit exceeded", 4);

fn recordHasUnwritableEnding(record: zfastq.Record, canonical_span: ?[]const u8) bool {
    return canonical_span == null and fastq.recordHasTerminalCr(record);
}

fn writeCheckedRecord(
    writer: *zfastq.Writer,
    direct_writer: ?*std.Io.Writer,
    record: zfastq.Record,
    canonical_span: ?[]const u8,
) error{WriteFailed}!void {
    if (canonical_span) |span| {
        if (direct_writer) |output| {
            try output.writeAll(span);
        } else {
            try fastq.writeCanonicalRecordSpan(writer, span);
        }
    } else {
        try fastq.writeRecordFields(writer, record);
    }
}

fn unwritableRecordFailure() CommandFailure {
    return CommandFailure.plain(
        "unwritable_record",
        "record fields ending in CR cannot be written with LF endings",
        1,
    );
}

fn mapCurrentSemanticFailure(
    maybe_error: ?zfastq.SemanticError,
    reader: *const zfastq.Reader,
    record_index: u64,
) ?CommandFailure {
    const semantic_error = maybe_error orelse return null;
    const offsets = reader.currentRecordOffsets() orelse return CommandFailure.plain(
        "io_error",
        "record location is unavailable",
        3,
    );
    return mapSemanticFailure(semantic_error, offsets, record_index);
}

fn mapSemanticFailure(
    maybe_error: ?zfastq.SemanticError,
    offsets: zfastq.RecordOffsets,
    record_index: u64,
) ?CommandFailure {
    const semantic_error = maybe_error orelse return null;
    const field_offset = switch (semantic_error.field) {
        .sequence => offsets.sequence,
        .quality => offsets.quality,
    };
    const details = fastq.semanticParseError(semantic_error, record_index, field_offset) catch {
        return INPUT_LOCATION_LIMIT;
    };
    return CommandFailure.lint(details);
}

const StatsOutcome = union(enum) {
    success: zfastq.Stats,
    failure: CommandFailure,
};

fn validateRecordCommandInputs(
    io: std.Io,
    command: Command,
    inputs: []const []const u8,
) ?u8 {
    if (inputs.len == 0) {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
        std.Io.File.writeStreamingAll(.stderr(), io, @tagName(command)) catch {};
        std.Io.File.writeStreamingAll(.stderr(), io, " requires at least one input\n") catch {};
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
    return null;
}

// --- Check command ---

const PairMode = enum {
    none,
    paired,
    interleaved,
};

const CheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
};

const PairedCheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
    pair_mode: PairMode = .none,
    pair_name_policy: pairing.NamePolicy = .illumina,
};

const BoundedBytes = struct {
    prefix: [PAIR_DIAGNOSTIC_PREFIX_BYTES]u8 = undefined,
    prefix_len: u8,
    full_len: usize,

    fn init(value: []const u8) BoundedBytes {
        const prefix_len = @min(value.len, PAIR_DIAGNOSTIC_PREFIX_BYTES);
        var display: BoundedBytes = .{
            .prefix_len = @intCast(prefix_len),
            .full_len = value.len,
        };
        @memcpy(display.prefix[0..prefix_len], value[0..prefix_len]);
        return display;
    }

    fn bytes(self: *const BoundedBytes) []const u8 {
        return self.prefix[0..self.prefix_len];
    }

    fn truncated(self: *const BoundedBytes) bool {
        return self.full_len > @as(usize, self.prefix_len);
    }
};

const PairRecordDiagnostic = struct {
    record_index: u64,
    byte_offset: u64,
    first_token: BoundedBytes,
    normalized_id: BoundedBytes,
    mate_markers: u2,

    fn init(name: pairing.Name, record_index: u64, byte_offset: u64) PairRecordDiagnostic {
        return .{
            .record_index = record_index,
            .byte_offset = byte_offset,
            .first_token = .init(name.first_token),
            .normalized_id = .init(name.normalized_id),
            .mate_markers = name.mate_markers,
        };
    }
};

const StoredPairName = struct {
    normalized_id: std.ArrayList(u8) = .empty,
    first_token_len: usize = undefined,
    first_mate_marker: ?u2 = undefined,
    mate_markers: u2 = undefined,

    fn deinit(self: *StoredPairName, allocator: std.mem.Allocator) void {
        self.normalized_id.deinit(allocator);
    }

    fn store(self: *StoredPairName, allocator: std.mem.Allocator, parsed: pairing.Name) !void {
        try storeNormalizedId(allocator, &self.normalized_id, parsed.normalized_id);
        self.first_token_len = parsed.first_token.len;
        self.first_mate_marker = parsed.first_mate_marker;
        self.mate_markers = parsed.mate_markers;
    }

    fn name(self: *const StoredPairName) pairing.Name {
        return .{
            .first_token = self.normalized_id.items,
            .normalized_id = self.normalized_id.items,
            .mate_markers = self.mate_markers,
            .first_mate_marker = self.first_mate_marker,
        };
    }

    fn diagnostic(self: *const StoredPairName, record_index: u64, byte_offset: u64) PairRecordDiagnostic {
        var first_token = BoundedBytes.init(self.normalized_id.items);
        first_token.full_len = self.first_token_len;
        if (self.first_mate_marker) |marker| {
            if (first_token.prefix_len < PAIR_DIAGNOSTIC_PREFIX_BYTES) {
                first_token.prefix[first_token.prefix_len] = '/';
                first_token.prefix_len += 1;
            }
            if (first_token.prefix_len < PAIR_DIAGNOSTIC_PREFIX_BYTES) {
                first_token.prefix[first_token.prefix_len] = '0' + @as(u8, marker);
                first_token.prefix_len += 1;
            }
        }
        return .{
            .record_index = record_index,
            .byte_offset = byte_offset,
            .first_token = first_token,
            .normalized_id = .init(self.normalized_id.items),
            .mate_markers = self.mate_markers,
        };
    }
};

const PairNameFailure = struct {
    pair_index: u64,
    records: [2]PairRecordDiagnostic,
};

const PairCountFailure = struct {
    pair_index: u64,
    remaining_side: u1,
    record_indexes: [2]?u64,
};

const PairFailure = union(enum) {
    name_mismatch: PairNameFailure,
    count_mismatch: PairCountFailure,
};

const PairCommandFailure = union(enum) {
    command: struct {
        input_index: u1,
        details: CommandFailure,
    },
    pair: PairFailure,

    fn exitCode(self: PairCommandFailure) u8 {
        return switch (self) {
            .command => |failure| failure.details.exit_code,
            .pair => 1,
        };
    }
};

fn runSingleCheckCommand(
    io: std.Io,
    inputs: []const []const u8,
    options: CheckOptions,
    json_output: bool,
    pair_names_set: bool,
) u8 {
    if (pair_names_set) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: --pair-names requires --paired or --interleaved\n",
        ) catch {};
        return 2;
    }
    return runCheck(io, inputs, options, json_output);
}

fn runCheck(
    io: std.Io,
    inputs: []const []const u8,
    options: CheckOptions,
    json_output: bool,
) u8 {
    if (validateRecordCommandInputs(io, .check, inputs)) |exit_code| return exit_code;
    if (json_output) return runCheckJson(io, inputs, options);

    var exit_code: u8 = 0;
    for (inputs) |input| {
        const failure = checkInput(io, input, options);
        if (failure) |details| {
            printCommandFailure(io, input, details);
            exit_code = @max(exit_code, details.exit_code);
        }
    }
    return exit_code;
}

fn runPairedCheckCommand(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: PairedCheckOptions,
    json_output: bool,
) u8 {
    if (validatePairedCheckInputs(io, inputs, options.pair_mode)) |exit_code| {
        return exit_code;
    }
    if (json_output) return runPairedCheckJson(io, allocator, inputs, options);

    const failure = checkPairMode(io, allocator, inputs, options);
    if (failure) |details| printPairCommandFailure(io, inputs, options.pair_mode, details);
    return if (failure) |details| details.exitCode() else 0;
}

fn validatePairedCheckInputs(
    io: std.Io,
    inputs: []const []const u8,
    pair_mode: PairMode,
) ?u8 {
    const expected_inputs: usize = if (pair_mode == .paired) 2 else 1;
    if (inputs.len != expected_inputs) {
        const message = if (pair_mode == .paired)
            "error: check --paired requires exactly two inputs\n"
        else
            "error: check --interleaved requires exactly one input\n";
        std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
        return 2;
    }
    if (pair_mode == .paired and
        std.mem.eql(u8, inputs[0], "-") and
        std.mem.eql(u8, inputs[1], "-"))
    {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: paired inputs may contain standard input at most once\n",
        ) catch {};
        return 2;
    }
    return null;
}

fn runCheckJson(
    io: std.Io,
    inputs: []const []const u8,
    options: CheckOptions,
) u8 {
    var stdout_buffer: [16 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/check-v1") catch return 3;
    var exit_code: u8 = 0;
    for (inputs) |input| {
        const failure = checkInput(io, input, options);
        if (failure) |details| exit_code = @max(exit_code, details.exit_code);
        writeCheckJsonResult(&json, input, failure) catch return @max(exit_code, 3);
    }
    finishJsonDocument(&json) catch return @max(exit_code, 3);
    stdout_writer.interface.flush() catch return @max(exit_code, 3);
    return exit_code;
}

fn runPairedCheckJson(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: PairedCheckOptions,
) u8 {
    const failure = checkPairMode(io, allocator, inputs, options);
    if (failure) |details| {
        if (details.exitCode() == 2) {
            printPairCommandFailure(io, inputs, options.pair_mode, details);
            return 2;
        }
    }
    var stdout_buffer: [16 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/check-v1") catch return 3;
    const exit_code = if (failure) |details| details.exitCode() else 0;
    writePairedCheckJsonResult(&json, inputs, options.pair_mode, failure) catch
        return @max(exit_code, 3);
    finishJsonDocument(&json) catch return @max(exit_code, 3);
    stdout_writer.interface.flush() catch return @max(exit_code, 3);
    return exit_code;
}

fn checkPairMode(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: PairedCheckOptions,
) ?PairCommandFailure {
    return switch (options.pair_mode) {
        .none => unreachable,
        .paired => checkPaired(io, allocator, inputs, options),
        .interleaved => checkInterleaved(io, allocator, inputs[0], options),
    };
}

fn checkPaired(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: PairedCheckOptions,
) ?PairCommandFailure {
    var input1: RecordInput = undefined;
    var input2: RecordInput = undefined;
    if (initPairedRecordInputs(&input1, &input2, io, inputs, null)) |failure| return failure;
    defer input1.deinit(io);
    defer input2.deinit(io);

    if (comptime build_options.use_isa_l) {
        return checkPairedSources(
            allocator,
            input1.byteSource(),
            input2.byteSource(),
            options,
            null,
        );
    }
    return checkPairedSources(
        allocator,
        &input1,
        &input2,
        options,
        null,
    );
}

fn checkPairedSources(
    allocator: std.mem.Allocator,
    source1: anytype,
    source2: @TypeOf(source1),
    options: PairedCheckOptions,
    exact_selector: ?*sampling.ExactSelector,
) ?PairCommandFailure {
    const reader_options: fastq.Options = .{ .max_line_bytes = options.max_line_bytes };
    var reader1 = initSourceReader(allocator, source1, reader_options) catch
        return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader1.deinit();
    var reader2 = initSourceReader(allocator, source2, reader_options) catch
        return pairCommandFailure(1, OUT_OF_MEMORY);
    defer reader2.deinit();
    var validator1 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var validator2 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        const record1 = fastq.nextValidatedHeader(&reader1, &validator1) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader1, err));
        };
        const record2 = fastq.nextValidatedHeader(&reader2, &validator2) catch |err| {
            return pairCommandFailure(1, mapReaderFailure(&reader2, err));
        };

        if (record1 == null and record2 == null) return null;
        if (record1 == null or record2 == null) {
            const remaining_side: u1 = if (record1 != null) 0 else 1;
            const pair_index = if (record1 != null)
                reader1.recordIndex() - 1
            else
                reader2.recordIndex() - 1;
            return .{ .pair = .{ .count_mismatch = .{
                .pair_index = pair_index,
                .remaining_side = remaining_side,
                .record_indexes = .{
                    lastRecordIndex(&reader1),
                    lastRecordIndex(&reader2),
                },
            } } };
        }

        const record_index1 = reader1.recordIndex() - 1;
        const record_index2 = reader2.recordIndex() - 1;
        const offsets1 = reader1.currentRecordOffsets().?;
        const offsets2 = reader2.currentRecordOffsets().?;
        if (mapSemanticFailure(
            record1.?.semantic_error,
            offsets1,
            record_index1,
        )) |failure| {
            return pairCommandFailure(0, failure);
        }
        if (mapSemanticFailure(
            record2.?.semantic_error,
            offsets2,
            record_index2,
        )) |failure| {
            return pairCommandFailure(1, failure);
        }

        if (!pairing.headersMatch(
            record1.?.header,
            record2.?.header,
            options.pair_name_policy,
        )) {
            const name1 = pairing.parseName(record1.?.header, options.pair_name_policy);
            const name2 = pairing.parseName(record2.?.header, options.pair_name_policy);
            return .{ .pair = .{ .name_mismatch = .{
                .pair_index = record_index1,
                .records = .{
                    .init(name1, record_index1, offsets1.header),
                    .init(name2, record_index2, offsets2.header),
                },
            } } };
        }

        if (exact_selector) |selector| {
            selector.considerRecord(allocator, reader1.recordIndex()) catch |err| {
                return exactPairSelectionFailure(err);
            };
        }
    }
}

fn checkInterleaved(
    io: std.Io,
    allocator: std.mem.Allocator,
    input_label: []const u8,
    options: PairedCheckOptions,
) ?PairCommandFailure {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, input_label, null)) |failure| {
        return pairCommandFailure(0, failure);
    }
    defer input.deinit(io);

    return checkInterleavedSource(
        allocator,
        input.byteSource(),
        options,
        null,
    );
}

fn checkInterleavedSource(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    options: PairedCheckOptions,
    exact_selector: ?*sampling.ExactSelector,
) ?PairCommandFailure {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader.deinit();
    var stored_name: StoredPairName = .{};
    defer stored_name.deinit(allocator);
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        var record1 = fastq.nextValidatedHeader(&reader, &validator) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader, err));
        } orelse return null;
        const record_index1 = reader.recordIndex() - 1;
        const offsets1 = reader.currentRecordOffsets().?;
        const semantic1 = mapSemanticFailure(
            record1.semantic_error,
            offsets1,
            record_index1,
        );

        const paired_record2 = (if (semantic1 == null)
            fastq.nextPairedValidatedHeader(&reader, &record1.header, &validator)
        else
            fastq.nextBufferedValidatedRecord(&reader, &validator)) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader, err));
        };
        const both_headers_borrowed = paired_record2 != null;
        var semantic2: ?zfastq.SemanticError = null;
        const header2 = if (paired_record2) |validated2| header: {
            if (semantic1 == null) semantic2 = validated2.semantic_error;
            break :header validated2.record.header;
        } else header: {
            if (semantic1 == null) {
                const name1 = pairing.parseName(record1.header, options.pair_name_policy);
                stored_name.store(allocator, name1) catch return pairCommandFailure(0, OUT_OF_MEMORY);
            }
            if (semantic1) |failure| {
                const record2 = fastq.nextRecordWithoutId(&reader) catch |err| {
                    return pairCommandFailure(0, mapReaderFailure(&reader, err));
                } orelse return missingInterleavedMateFailure(record_index1);
                _ = record2;
                return pairCommandFailure(0, failure);
            }
            const record2 = fastq.nextValidatedHeader(&reader, &validator) catch |err| {
                return pairCommandFailure(0, mapReaderFailure(&reader, err));
            } orelse return missingInterleavedMateFailure(record_index1);
            semantic2 = record2.semantic_error;
            break :header record2.header;
        };

        if (semantic1) |failure| {
            return pairCommandFailure(0, failure);
        }

        const record_index2 = reader.recordIndex() - 1;
        const offsets2 = reader.currentRecordOffsets().?;
        if (mapSemanticFailure(
            semantic2,
            offsets2,
            record_index2,
        )) |failure| {
            return pairCommandFailure(0, failure);
        }

        const names_match = if (both_headers_borrowed)
            pairing.headersMatch(record1.header, header2, options.pair_name_policy)
        else blk: {
            const name1 = stored_name.name();
            const name2 = pairing.parseName(header2, options.pair_name_policy);
            break :blk pairing.namesMatch(name1, name2);
        };
        if (!names_match) {
            const name2 = pairing.parseName(header2, options.pair_name_policy);
            return .{ .pair = .{ .name_mismatch = .{
                .pair_index = record_index1 / 2,
                .records = .{
                    if (both_headers_borrowed)
                        .init(
                            pairing.parseName(record1.header, options.pair_name_policy),
                            record_index1,
                            offsets1.header,
                        )
                    else
                        stored_name.diagnostic(record_index1, offsets1.header),
                    .init(name2, record_index2, offsets2.header),
                },
            } } };
        }

        if (exact_selector) |selector| {
            selector.considerRecord(allocator, reader.recordIndex() / 2) catch |err| {
                return exactPairSelectionFailure(err);
            };
        }
    }
}

fn exactPairSelectionFailure(err: sampling.ReservoirError) PairCommandFailure {
    return pairCommandFailure(0, exactSelectionFailure(err));
}

fn storeNormalizedId(
    allocator: std.mem.Allocator,
    storage: *std.ArrayList(u8),
    normalized_id: []const u8,
) std.mem.Allocator.Error!void {
    storage.clearRetainingCapacity();
    try storage.ensureTotalCapacityPrecise(allocator, normalized_id.len);
    storage.appendSliceAssumeCapacity(normalized_id);
}

fn exactSelectionFailure(err: sampling.ReservoirError) CommandFailure {
    return switch (err) {
        error.OutOfMemory => OUT_OF_MEMORY,
        error.Overflow => CommandFailure.plain(
            "arithmetic_limit",
            "sample index storage exceeds supported limit",
            4,
        ),
    };
}

fn initRecordInput(
    input: *RecordInput,
    io: std.Io,
    label: []const u8,
    output_identity: ?FileIdentity,
) ?CommandFailure {
    const file = openRecordFile(io, label) catch |err| return inputOpenFailure(err);
    const owns_file = !std.mem.eql(u8, label, "-");
    var transferred = false;
    defer if (!transferred and owns_file) file.close(io);
    if (outputFileFailure(io, file, output_identity)) |failure| return failure;
    input.init(io, file, owns_file) catch
        return IO_FAILURE;
    transferred = true;
    return null;
}

fn openRecordFile(io: std.Io, label: []const u8) std.Io.File.OpenError!std.Io.File {
    if (std.mem.eql(u8, label, "-")) return .stdin();
    return std.Io.Dir.cwd().openFile(io, label, .{});
}

fn inputOpenFailure(err: (std.Io.File.OpenError || std.posix.OpenError)) CommandFailure {
    return CommandFailure.plain("io_error", if (err == error.FileNotFound)
        "file not found"
    else
        "failed to open file", 3);
}

fn initPairedRecordInputs(
    input1: *RecordInput,
    input2: *RecordInput,
    io: std.Io,
    inputs: []const []const u8,
    output_identity: ?FileIdentity,
) ?PairCommandFailure {
    for (inputs, 0..) |label, index| {
        if (std.mem.eql(u8, label, "-")) {
            // Opening a path could reuse descriptor 0 if stdin is closed.
            _ = fileIdentity(io, .stdin()) catch
                return pairCommandFailure(@intCast(index), IO_FAILURE);
        }
    }
    var transferred = false;
    const owns1 = !std.mem.eql(u8, inputs[0], "-");
    const file1 = openRecordFile(io, inputs[0]) catch |err|
        return pairCommandFailure(0, inputOpenFailure(err));
    defer if (!transferred and owns1) file1.close(io);

    const owns2 = !std.mem.eql(u8, inputs[1], "-");
    var file2 = openNonblockingRecordFile(inputs[1]) catch |err|
        return pairCommandFailure(1, inputOpenFailure(err));
    defer if (!transferred and owns2) file2.close(io);

    if (pairedFilesFailure(io, file1, file2, output_identity)) |failure| return failure;
    input1.init(io, file1, owns1) catch return pairCommandFailure(0, IO_FAILURE);
    if (owns2) {
        prepareRecordFile(io, &file2) catch return pairCommandFailure(1, IO_FAILURE);
    }
    input2.init(io, file2, owns2) catch return pairCommandFailure(1, IO_FAILURE);
    transferred = true;
    return null;
}

fn openNonblockingRecordFile(label: []const u8) std.posix.OpenError!std.Io.File {
    if (std.mem.eql(u8, label, "-")) return .stdin();
    // R2's FIFO writer may wait for R1 to drain before opening R2.
    const handle = try std.posix.openat(std.Io.Dir.cwd().handle, label, .{
        .ACCMODE = .RDONLY,
        .NONBLOCK = true,
        .CLOEXEC = true,
        .NOCTTY = true,
    }, 0);
    return .{ .handle = handle, .flags = .{ .nonblocking = true } };
}

fn prepareRecordFile(io: std.Io, file: *std.Io.File) (std.Io.Cancelable || error{ Io, BadFileDescriptor })!void {
    const linux = std.os.linux;
    const stat = try statRecordFile(io, file.*);
    if (stat.mode & linux.S.IFMT == linux.S.IFIFO) {
        // A nonblocking FIFO open can precede its writer; reading now would report early EOF.
        var poll_fd = [1]linux.pollfd{.{ .fd = file.handle, .events = linux.POLL.IN, .revents = 0 }};
        while (true) switch (linux.errno(linux.poll(&poll_fd, 1, -1))) {
            .SUCCESS => break,
            .INTR => try io.checkCancel(),
            else => return error.Io,
        };
        if (poll_fd[0].revents & (linux.POLL.ERR | linux.POLL.NVAL) != 0) return error.Io;
    }
    while (true) switch (linux.errno(linux.fcntl(file.handle, linux.F.SETFL, 0))) {
        .SUCCESS => break,
        .INTR => try io.checkCancel(),
        else => return error.Io,
    };
    file.flags.nonblocking = false;
}

const FileIdentity = struct {
    device_major: u32,
    device_minor: u32,
    inode: u64,
};

fn fileIdentity(io: std.Io, file: std.Io.File) (std.Io.Cancelable || error{ Io, BadFileDescriptor })!FileIdentity {
    // std.Io.File.Stat does not retain the device containing the inode.
    const stat = try statRecordFile(io, file);
    return .{ .device_major = stat.dev_major, .device_minor = stat.dev_minor, .inode = stat.ino };
}

fn statRecordFile(io: std.Io, file: std.Io.File) (std.Io.Cancelable || error{ Io, BadFileDescriptor })!std.os.linux.Statx {
    const linux = std.os.linux;
    var stat: linux.Statx = undefined;
    while (true) switch (linux.errno(linux.statx(file.handle, "", linux.AT.EMPTY_PATH, .{ .INO = true, .TYPE = true }, &stat))) {
        .SUCCESS => break,
        .INTR => try io.checkCancel(),
        .BADF => return error.BadFileDescriptor,
        else => return error.Io,
    };
    if (!stat.mask.INO or !stat.mask.TYPE) return error.Io;
    return stat;
}

fn stdoutFileIdentity(io: std.Io) (std.Io.Cancelable || error{Io})!?FileIdentity {
    const stat = statRecordFile(io, .stdout()) catch |err| switch (err) {
        // Inspect before opening inputs, which could reuse a closed descriptor 1.
        error.BadFileDescriptor => return null,
        else => |other| return other,
    };
    if (stat.mode & std.os.linux.S.IFMT != std.os.linux.S.IFREG) return null;
    return .{ .device_major = stat.dev_major, .device_minor = stat.dev_minor, .inode = stat.ino };
}

fn outputFileFailure(io: std.Io, file: std.Io.File, output_identity: ?FileIdentity) ?CommandFailure {
    const output = output_identity orelse return null;
    const input = fileIdentity(io, file) catch return inputInspectionFailure(io, output_identity);
    if (std.meta.eql(input, output)) {
        return outputAliasFailure(io, &.{input});
    }
    return null;
}

fn canReportInputFailure(io: std.Io, known_inputs: ?[]const FileIdentity) bool {
    const stat = statRecordFile(io, .stderr()) catch return false;
    // Without input identities, a regular stderr file could be an input.
    const inputs = known_inputs orelse return stat.mode & std.os.linux.S.IFMT != std.os.linux.S.IFREG;
    const stderr_identity: FileIdentity = .{
        .device_major = stat.dev_major,
        .device_minor = stat.dev_minor,
        .inode = stat.ino,
    };
    for (inputs) |input| {
        if (std.meta.eql(input, stderr_identity)) return false;
    }
    return true;
}

fn outputAliasFailure(io: std.Io, inputs: []const FileIdentity) CommandFailure {
    var failure = CommandFailure.plain("same_output", "input file is output file", 2);
    failure.suppress_diagnostic = !canReportInputFailure(io, inputs);
    return failure;
}

fn inputInspectionFailure(io: std.Io, output_identity: ?FileIdentity) CommandFailure {
    var failure = CommandFailure.plain("io_error", "failed to inspect file", 3);
    if (output_identity != null) failure.suppress_diagnostic = !canReportInputFailure(io, null);
    return failure;
}

fn pairedFilesFailure(io: std.Io, file1: std.Io.File, file2: std.Io.File, output_identity: ?FileIdentity) ?PairCommandFailure {
    const first = fileIdentity(io, file1) catch
        return pairCommandFailure(0, inputInspectionFailure(io, output_identity));
    const second = fileIdentity(io, file2) catch
        return pairCommandFailure(1, inputInspectionFailure(io, output_identity));
    if (std.meta.eql(first, second)) {
        var failure = pairCommandFailure(1, CommandFailure.plain("same_input", "paired inputs refer to the same file", 2));
        if (output_identity) |output| {
            if (std.meta.eql(first, output)) {
                failure.command.details.suppress_diagnostic = !canReportInputFailure(io, &.{first});
            }
        }
        return failure;
    }
    if (output_identity) |output| {
        for ([_]FileIdentity{ first, second }, 0..) |input, index| {
            if (std.meta.eql(input, output)) {
                return pairCommandFailure(@intCast(index), outputAliasFailure(io, &.{ first, second }));
            }
        }
    }
    return null;
}

fn pairCommandFailure(input_index: u1, details: CommandFailure) PairCommandFailure {
    return .{ .command = .{ .input_index = input_index, .details = details } };
}

fn lastRecordIndex(reader: *const zfastq.Reader) ?u64 {
    return if (reader.recordIndex() == 0) null else reader.recordIndex() - 1;
}

fn checkInput(
    io: std.Io,
    label: []const u8,
    options: CheckOptions,
) ?CommandFailure {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, label, null)) |failure| return failure;
    defer input.deinit(io);
    return checkRecordInput(&input, options);
}

fn checkRecordInput(
    input: *RecordInput,
    options: CheckOptions,
) ?CommandFailure {
    var scanner = fastq.CheckScanner.init(
        .{ .max_line_bytes = options.max_line_bytes },
        .{ .alphabet = options.alphabet },
    );
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const chunk = input.readScannerChunk(&buf) catch
            return IO_FAILURE;
        const decoded = chunk orelse break;
        _ = scanner.feed(decoded) catch |err| {
            return mapCheckScannerFailure(&scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapCheckScannerFailure(&scanner, err);
    };
    return null;
}

fn mapCheckScannerFailure(
    scanner: *fastq.CheckScanner,
    err: fastq.CheckScannerError,
) CommandFailure {
    return switch (err) {
        error.Format => validationFailure(scanner.takeLastError()),
        else => mapInputFailure(@errorCast(err)),
    };
}

// --- Stats command ---

fn runStats(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: InputOptions,
    json_output: bool,
) u8 {
    if (validateRecordCommandInputs(io, .stats, inputs)) |exit_code| return exit_code;
    if (json_output) return runStatsJson(io, allocator, inputs, options);

    var exit_code: u8 = 0;
    var printed_block = false;
    for (inputs) |input| {
        const outcome = statsInput(io, allocator, input, options);
        const stats = switch (outcome) {
            .success => |stats| stats,
            .failure => |failure| {
                printCommandFailure(io, input, failure);
                exit_code = @max(exit_code, failure.exit_code);
                continue;
            },
        };

        printStats(io, input, stats.result(), printed_block) catch
            return @max(exit_code, 3);
        printed_block = true;
    }
    return exit_code;
}

fn runStatsJson(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: InputOptions,
) u8 {
    var stdout_buffer: [16 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/stats-v1") catch return 3;
    var exit_code: u8 = 0;
    for (inputs) |input| {
        const outcome = statsInput(io, allocator, input, options);
        switch (outcome) {
            .success => {},
            .failure => |failure| exit_code = @max(exit_code, failure.exit_code),
        }
        writeStatsJsonResult(&json, input, outcome) catch return @max(exit_code, 3);
    }
    finishJsonDocument(&json) catch return @max(exit_code, 3);
    stdout_writer.interface.flush() catch return @max(exit_code, 3);
    return exit_code;
}

fn statsInput(
    io: std.Io,
    allocator: std.mem.Allocator,
    label: []const u8,
    options: InputOptions,
) StatsOutcome {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, label, null)) |failure| return .{ .failure = failure };
    defer input.deinit(io);
    return collectStats(allocator, input.byteSource(), options);
}

const RecordInput = struct {
    file: std.Io.File,
    owns_file: bool,
    read_buffer: [64 * 1024]u8,
    file_reader: std.Io.File.Reader,
    source: union(enum) {
        plain: io_layer.PlainFileSource,
        gzip: io_layer.GzipSource,
    },

    fn init(
        self: *RecordInput,
        io: std.Io,
        file: std.Io.File,
        owns_file: bool,
    ) error{Io}!void {
        self.file = file;
        self.owns_file = owns_file;
        self.file_reader = file.readerStreaming(io, &self.read_buffer);
        const prefix = self.file_reader.interface.peek(2) catch |err| switch (err) {
            error.EndOfStream => null,
            error.ReadFailed => return error.Io,
        };
        if (prefix) |bytes| {
            if (std.mem.eql(u8, bytes, &.{ 0x1f, 0x8b })) {
                self.source = .{ .gzip = io_layer.GzipSource.init(&self.file_reader.interface) };
                return;
            }
        }
        self.source = .{ .plain = io_layer.PlainFileSource.init(&self.file_reader) };
    }

    fn deinit(self: *RecordInput, io: std.Io) void {
        if (self.owns_file) self.file.close(io);
        self.* = undefined;
    }

    fn byteSource(self: *RecordInput) zfastq.io.ByteSource {
        return switch (self.source) {
            .plain => |*source| source.byteSource(),
            .gzip => |*source| source.byteSource(),
        };
    }

    fn initReader(
        self: *RecordInput,
        allocator: std.mem.Allocator,
        options: fastq.Options,
    ) !zfastq.Reader {
        return switch (self.source) {
            .plain => |*source| zfastq.Reader.init(allocator, source.byteSource(), options),
            .gzip => |*source| if (build_options.use_isa_l)
                zfastq.Reader.init(allocator, source.byteSource(), options)
            else
                fastq.initBorrowedGzipReader(allocator, source, options),
        };
    }

    fn readScannerChunk(self: *RecordInput, buffer: []u8) error{Io}!?[]const u8 {
        return switch (self.source) {
            .plain => |*source| plain: {
                const byte_source = source.byteSource();
                const count = byte_source.read(buffer) catch return error.Io;
                break :plain if (count == 0) null else buffer[0..count];
            },
            .gzip => |*source| io_layer.readGzipChunk(source, buffer) catch return error.Io,
        };
    }
};

fn initSourceReader(
    allocator: std.mem.Allocator,
    source: anytype,
    options: fastq.Options,
) !zfastq.Reader {
    return if (comptime @TypeOf(source) == zfastq.io.ByteSource)
        zfastq.Reader.init(allocator, source, options)
    else
        source.initReader(allocator, options);
}

fn collectStats(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    options: InputOptions,
) StatsOutcome {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return .{ .failure = OUT_OF_MEMORY };
    defer reader.deinit();

    var stats: zfastq.Stats = .{};
    while (fastq.nextPredictedPayload(&reader) catch |err| {
        return .{ .failure = mapReaderFailure(&reader, err) };
    }) |predicted| {
        var payload = predicted.payload;
        var checkpoint = predicted.checkpoint;
        while (true) {
            stats.addRecord(.{
                .header = "",
                .id = "",
                .sequence = payload.sequence,
                .plus = "",
                .quality = payload.quality,
            }) catch |err| {
                if (err == error.S006InvalidQuality) {
                    if (checkpoint) |saved| {
                        payload = saved.reread(&reader) catch |reader_error| {
                            return .{ .failure = mapReaderFailure(&reader, reader_error) };
                        };
                        checkpoint = null;
                        continue;
                    }
                }
                switch (err) {
                    error.S006InvalidQuality => {
                        const quality_error = stats.takeLastQualityError() orelse {
                            return .{ .failure = CommandFailure.plain(
                                "io_error",
                                "quality validation failed without details",
                                3,
                            ) };
                        };
                        const offsets = reader.currentRecordOffsets() orelse {
                            return .{ .failure = CommandFailure.plain(
                                "io_error",
                                "record location is unavailable",
                                3,
                            ) };
                        };
                        const details = fastq.semanticParseError(
                            fastq.semanticQualityError(quality_error.byte_index),
                            reader.recordIndex() - 1,
                            offsets.quality,
                        ) catch return .{ .failure = CommandFailure.plain(
                            "arithmetic_limit",
                            "statistics arithmetic limit exceeded",
                            4,
                        ) };
                        return .{ .failure = CommandFailure.lint(details) };
                    },
                    error.S005LengthMismatch => @panic("Reader returned unequal sequence and quality lengths"),
                    error.Overflow => return .{ .failure = CommandFailure.plain(
                        "arithmetic_limit",
                        "statistics arithmetic limit exceeded",
                        4,
                    ) },
                }
            };
            break;
        }
    }
    return .{ .success = stats };
}

fn mapReaderFailure(reader: *zfastq.Reader, err: zfastq.ReaderError) CommandFailure {
    return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => validationFailure(reader.takeLastError()),
        else => mapInputFailure(@errorCast(err)),
    };
}

fn validationFailure(details: ?zfastq.ParseError) CommandFailure {
    return if (details) |failure|
        CommandFailure.lint(failure)
    else
        CommandFailure.plain("io_error", "validation failed without details", 3);
}

fn mapInputFailure(err: error{ LineTooLong, ArithmeticLimit, OutOfMemory, Io }) CommandFailure {
    return switch (err) {
        error.LineTooLong => LINE_LIMIT,
        error.ArithmeticLimit => INPUT_LOCATION_LIMIT,
        error.OutOfMemory => OUT_OF_MEMORY,
        error.Io => IO_FAILURE,
    };
}

// --- Sample command ---

const SAMPLE_SCAN_BUFFER_BYTES = 64 * 1024;

const SampleOptions = struct {
    max_line_bytes: usize,
    alphabet: zfastq.Alphabet,
    fraction: ?sampling.Fraction,
    count: ?u64,
    seed: u64,
    pair_mode: PairMode = .none,
    pair_name_policy: pairing.NamePolicy = .illumina,
    pair_names_set: bool = false,
    output_identity: ?FileIdentity = null,
};

const SampleMode = union(enum) {
    fraction: sampling.Fraction,
    count: u64,
};

const ExactOutputCursor = struct {
    select_all: bool,
    indexes: sampling.ExactIndexes,
    selected_cursor: usize = 0,
    unit_count: u64 = 0,

    fn nextSelected(self: ExactOutputCursor) ?u64 {
        if (self.select_all or self.selected_cursor == self.indexes.len()) return null;
        return self.indexes.at(self.selected_cursor);
    }

    fn selectNextCached(self: *ExactOutputCursor, next_selected: *?u64) bool {
        if (self.select_all) return true;
        const selected = next_selected.* orelse return false;
        if (selected != self.unit_count + 1) return false;

        self.selected_cursor += 1;
        next_selected.* = self.nextSelected();
        return true;
    }

    fn completeUnit(self: *ExactOutputCursor) void {
        self.unit_count += 1;
    }

    fn selectionComplete(self: ExactOutputCursor) bool {
        return self.select_all or self.selected_cursor == self.indexes.len();
    }
};

const PairOutputSelector = union(enum) {
    all,
    fraction: *sampling.Selector,

    fn selectPair(self: *PairOutputSelector) bool {
        return switch (self.*) {
            .all => true,
            .fraction => |selector| selector.selectRecord(),
        };
    }
};

fn runSample(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    parsed_options: SampleOptions,
) u8 {
    var options = parsed_options;
    const mode: SampleMode = if (options.fraction) |fraction| blk: {
        if (options.count != null) {
            std.Io.File.writeStreamingAll(
                .stderr(),
                io,
                "error: --fraction and --count are mutually exclusive\n",
            ) catch {};
            return 2;
        }
        break :blk .{ .fraction = fraction };
    } else if (options.count) |count|
        .{ .count = count }
    else {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: sample requires --fraction P or --count K\n",
        ) catch {};
        return 2;
    };
    if (options.pair_mode != .none) {
        return runPairSampleCommand(io, allocator, inputs, mode, options);
    }
    if (options.pair_names_set) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: --pair-names requires --paired or --interleaved\n",
        ) catch {};
        return 2;
    }
    if (inputs.len != 1) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: sample requires exactly one input\n",
        ) catch {};
        return 2;
    }
    const input = inputs[0];
    if (mode == .count and std.mem.eql(u8, input, "-")) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: exact-count sampling requires a file path\n",
        ) catch {};
        return 2;
    }

    options.output_identity = stdoutFileIdentity(io) catch {
        if (canReportInputFailure(io, null)) printOutputFailure(io);
        return 3;
    };
    var stdout_buffer: [64 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var sink_adapter = io_layer.WriterSink.init(&stdout_writer.interface);
    var writer = zfastq.Writer.init(sink_adapter.byteSink());
    const failure = (switch (mode) {
        .fraction => |fraction| blk: {
            var selector = sampling.Selector.init(fraction, options.seed);
            break :blk sampleFractionInput(
                io,
                allocator,
                input,
                &writer,
                &selector,
                options,
            );
        },
        .count => |count| sampleExactFile(
            io,
            allocator,
            input,
            &writer,
            count,
            options,
        ),
    }) catch {
        writer.flush() catch {};
        return 3;
    };

    var output_exit: u8 = 0;
    writer.flush() catch {
        output_exit = 3;
    };
    if (failure) |details| {
        printCommandFailure(io, input, details);
        return @max(output_exit, details.exit_code);
    }
    return output_exit;
}

fn runPairSampleCommand(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    mode: SampleMode,
    parsed_options: SampleOptions,
) u8 {
    var options = parsed_options;
    const expected_inputs: usize = if (options.pair_mode == .paired) 2 else 1;
    if (inputs.len != expected_inputs) {
        const message = if (options.pair_mode == .paired)
            "error: sample --paired requires exactly two inputs\n"
        else
            "error: sample --interleaved requires exactly one input\n";
        std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
        return 2;
    }
    switch (mode) {
        .count => {
            for (inputs) |input| {
                if (!std.mem.eql(u8, input, "-")) continue;
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: paired exact-count sampling requires file paths\n",
                ) catch {};
                return 2;
            }
        },
        .fraction => if (options.pair_mode == .paired and
            std.mem.eql(u8, inputs[0], "-") and
            std.mem.eql(u8, inputs[1], "-"))
        {
            std.Io.File.writeStreamingAll(
                .stderr(),
                io,
                "error: paired sample inputs may contain standard input at most once\n",
            ) catch {};
            return 2;
        },
    }

    const needs_staging = options.pair_mode == .interleaved and switch (mode) {
        .fraction => |fraction| fraction != .none,
        .count => |count| count != 0,
    };
    const staging_limit = if (needs_staging)
        deinterleaveStagingLimit(options.max_line_bytes) catch {
            if (canReportInputFailure(io, null)) {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: sample record staging size exceeds supported limit\n",
                ) catch {};
            }
            return 4;
        }
    else
        0;
    options.output_identity = stdoutFileIdentity(io) catch {
        if (canReportInputFailure(io, null)) printOutputFailure(io);
        return 3;
    };
    var stdout_buffer: [64 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var sink_adapter = io_layer.WriterSink.init(&stdout_writer.interface);
    var writer = zfastq.Writer.init(sink_adapter.byteSink());
    const failure = (switch (mode) {
        .fraction => |fraction| blk: {
            var selector = sampling.Selector.init(fraction, options.seed);
            var output_selector: PairOutputSelector = .{ .fraction = &selector };
            break :blk switch (options.pair_mode) {
                .none => unreachable,
                .paired => interleaveInputs(io, allocator, inputs, &writer, null, .{
                    .max_line_bytes = options.max_line_bytes,
                    .alphabet = options.alphabet,
                    .pair_name_policy = options.pair_name_policy,
                    .selection = output_selector,
                    .output_identity = options.output_identity,
                }),
                .interleaved => sampleInterleavedFractionInput(
                    io,
                    allocator,
                    inputs[0],
                    &writer,
                    &output_selector,
                    staging_limit,
                    options,
                ),
            };
        },
        .count => |count| sampleExactPairs(
            io,
            allocator,
            inputs,
            &writer,
            count,
            staging_limit,
            options,
        ),
    }) catch {
        writer.flush() catch {};
        return 3;
    };

    var output_exit: u8 = 0;
    writer.flush() catch {
        output_exit = 3;
    };
    if (failure) |details| {
        printPairCommandFailure(io, inputs, options.pair_mode, details);
        return @max(output_exit, details.exitCode());
    }
    return output_exit;
}

fn sampleInterleavedFractionInput(
    io: std.Io,
    allocator: std.mem.Allocator,
    input_label: []const u8,
    writer: *zfastq.Writer,
    selection: *PairOutputSelector,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, input_label, options.output_identity)) |failure| {
        return pairCommandFailure(0, failure);
    }
    defer input.deinit(io);

    return sampleInterleavedSource(
        allocator,
        input.byteSource(),
        writer,
        selection,
        staging_limit,
        options,
    );
}

const InterleavedFirstRecordStorage = enum {
    unused,
    reader,
    retained,
    staged,
};

const RefillSpanningMate = struct {
    record: ?fastq.ValidatedRecord,
    first_storage: InterleavedFirstRecordStorage,
};

fn missingInterleavedMateFailure(record_index1: u64) PairCommandFailure {
    return .{ .pair = .{ .count_mismatch = .{
        .pair_index = record_index1 / 2,
        .remaining_side = 0,
        .record_indexes = .{ record_index1, if (record_index1 == 0) null else record_index1 - 1 },
    } } };
}

fn mapPreservedMateFailure(
    reader: *zfastq.Reader,
    err: (zfastq.ReaderError || error{RecordStagingLimit}),
) PairCommandFailure {
    return pairCommandFailure(0, switch (err) {
        error.RecordStagingLimit => RECORD_STAGING_LIMIT,
        else => mapReaderFailure(reader, @errorCast(err)),
    });
}

fn nextAfterPreservingInterleavedMate1(
    allocator: std.mem.Allocator,
    reader: *zfastq.Reader,
    retained: *fastq.RetainedRecordStorage,
    staging: *std.ArrayList(u8),
    record1: zfastq.Record,
    canonical_span1: ?[]const u8,
    staging_limit: usize,
    validator: *fastq.AdaptiveRecordValidator,
) (zfastq.ReaderError || error{RecordStagingLimit})!RefillSpanningMate {
    if (fastq.retainFallbackRecordStorage(reader, retained, record1)) {
        const transferred_record2 = try fastq.nextBufferedAfterFallbackTransfer(
            reader,
            validator,
        );
        if (transferred_record2) |record2| return .{
            .record = record2,
            .first_storage = .retained,
        };
        // Reader's four line limits cover CLI staging limits, including markers.
        // Internal callers can supply a smaller limit and still need its error.
        if (staging_limit < 4 or (staging_limit - 4) / 4 < reader.options.max_line_bytes) {
            if (try canonicalRecordSize(record1, canonical_span1) > staging_limit) {
                return error.RecordStagingLimit;
            }
        }
        return .{
            .record = try fastq.nextFallbackValidatedRecord(reader, validator),
            .first_storage = .retained,
        };
    }

    try stageCanonicalRecord(
        allocator,
        staging,
        record1,
        canonical_span1,
        staging_limit,
    );
    return .{
        .record = try fastq.nextValidatedRecord(reader, validator),
        .first_storage = .staged,
    };
}

fn writePreservedInterleavedMate1(
    writer: *zfastq.Writer,
    reader: *zfastq.Reader,
    retained: *fastq.RetainedRecordStorage,
    staged: []const u8,
    record: zfastq.Record,
    canonical_span: ?[]const u8,
    storage: InterleavedFirstRecordStorage,
) error{WriteFailed}!void {
    switch (storage) {
        .unused => unreachable,
        .reader, .retained => try writeCheckedRecord(writer, null, record, canonical_span),
        .staged => try fastq.writeCanonicalRecordSpan(writer, staged),
    }
    if (storage == .retained) {
        fastq.restoreFallbackRecordStorage(reader, retained);
    }
}

fn sampleInterleavedSource(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    writer: *zfastq.Writer,
    selection: *PairOutputSelector,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader.deinit();
    var stored_name: StoredPairName = .{};
    defer stored_name.deinit(allocator);
    var staged_record: std.ArrayList(u8) = .empty;
    defer staged_record.deinit(allocator);
    var retained_record_storage: fastq.RetainedRecordStorage = .{};
    defer retained_record_storage.deinit(allocator);
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        var validated1 = fastq.nextValidatedRecord(&reader, &validator) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader, err));
        } orelse return null;
        const record_index1 = reader.recordIndex() - 1;
        const offsets1 = reader.currentRecordOffsets().?;
        const semantic1 = mapSemanticFailure(
            validated1.semantic_error,
            offsets1,
            record_index1,
        );

        const selected = semantic1 == null and selection.selectPair();
        const unwritable1 = selected and recordHasUnwritableEnding(validated1.record, validated1.canonical_span);

        var record1_storage: InterleavedFirstRecordStorage = .unused;
        const paired_record2 = (if (selected)
            fastq.nextPairedValidatedRecord(&reader, &validated1, &validator)
        else if (semantic1 == null)
            fastq.nextPairedValidatedHeader(&reader, &validated1.record.header, &validator)
        else
            fastq.nextBufferedValidatedRecord(&reader, &validator)) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader, err));
        };
        const record1 = validated1.record;
        const canonical_span1 = validated1.canonical_span;
        const both_headers_borrowed = paired_record2 != null;
        const validated2 = (if (paired_record2) |complete_record2| buffered: {
            if (selected) record1_storage = .reader;
            break :buffered complete_record2;
        } else record: {
            if (semantic1 == null) {
                const name1 = pairing.parseName(record1.header, options.pair_name_policy);
                stored_name.store(allocator, name1) catch return pairCommandFailure(0, OUT_OF_MEMORY);
            }
            if (!selected) break :record fastq.nextValidatedRecord(
                &reader,
                &validator,
            ) catch |err| {
                return pairCommandFailure(0, mapReaderFailure(&reader, err));
            };

            const preserved = nextAfterPreservingInterleavedMate1(
                allocator,
                &reader,
                &retained_record_storage,
                &staged_record,
                record1,
                canonical_span1,
                staging_limit,
                &validator,
            ) catch |err| return mapPreservedMateFailure(&reader, err);
            record1_storage = preserved.first_storage;
            break :record preserved.record;
        }) orelse return missingInterleavedMateFailure(record_index1);

        const record2 = validated2.record;
        const canonical_span2 = validated2.canonical_span;
        if (semantic1) |failure| {
            return pairCommandFailure(0, failure);
        }

        const record_index2 = reader.recordIndex() - 1;
        const offsets2 = reader.currentRecordOffsets().?;
        if (mapSemanticFailure(
            validated2.semantic_error,
            offsets2,
            record_index2,
        )) |failure| {
            return pairCommandFailure(0, failure);
        }

        if (both_headers_borrowed) {
            if (!pairing.headersMatch(
                record1.header,
                record2.header,
                options.pair_name_policy,
            )) {
                const name1 = pairing.parseName(record1.header, options.pair_name_policy);
                const name2 = pairing.parseName(record2.header, options.pair_name_policy);
                return .{ .pair = .{ .name_mismatch = .{
                    .pair_index = record_index1 / 2,
                    .records = .{
                        .init(name1, record_index1, offsets1.header),
                        .init(name2, record_index2, offsets2.header),
                    },
                } } };
            }
        } else {
            const name1 = stored_name.name();
            const name2 = pairing.parseName(record2.header, options.pair_name_policy);
            if (!pairing.namesMatch(name1, name2)) {
                return .{ .pair = .{ .name_mismatch = .{
                    .pair_index = record_index1 / 2,
                    .records = .{
                        stored_name.diagnostic(record_index1, offsets1.header),
                        .init(name2, record_index2, offsets2.header),
                    },
                } } };
            }
        }
        if (!selected) continue;
        if (unwritable1 or recordHasUnwritableEnding(record2, canonical_span2)) {
            return pairCommandFailure(0, unwritableRecordFailure());
        }

        try writePreservedInterleavedMate1(
            writer,
            &reader,
            &retained_record_storage,
            staged_record.items,
            record1,
            canonical_span1,
            record1_storage,
        );
        try writeCheckedRecord(writer, null, record2, canonical_span2);
    }
}

fn sampleFractionInput(
    io: std.Io,
    allocator: std.mem.Allocator,
    label: []const u8,
    writer: *zfastq.Writer,
    selector: *sampling.Selector,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, label, options.output_identity)) |failure| return failure;
    defer input.deinit(io);
    if (comptime build_options.use_isa_l) {
        if (selector.* == .none) {
            return checkRecordInput(&input, .{
                .max_line_bytes = options.max_line_bytes,
                .alphabet = options.alphabet,
            });
        }
        return sampleFractionSource(allocator, input.byteSource(), writer, selector, options);
    }
    return sampleFractionSource(allocator, &input, writer, selector, options);
}

fn sampleFractionSource(
    allocator: std.mem.Allocator,
    source: anytype,
    writer: *zfastq.Writer,
    selector: *sampling.Selector,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    const reader_options: fastq.Options = .{ .max_line_bytes = options.max_line_bytes };
    var reader = initSourceReader(allocator, source, reader_options) catch
        return OUT_OF_MEMORY;
    defer reader.deinit();
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) switch (nextValidatedSampleRecord(&reader, &validator)) {
        .done => return null,
        .failure => |failure| return failure,
        .record => |validated| {
            if (!selector.selectRecord()) continue;
            if (recordHasUnwritableEnding(validated.record, validated.canonical_span)) {
                return unwritableRecordFailure();
            }
            try writeCheckedRecord(writer, null, validated.record, validated.canonical_span);
        },
    };
}

const ValidatedSampleRecord = struct {
    record: zfastq.Record,
    canonical_span: ?[]const u8,
};

const SampleRecordOutcome = union(enum) {
    done,
    record: ValidatedSampleRecord,
    failure: CommandFailure,
};

fn nextValidatedSampleRecord(
    reader: *zfastq.Reader,
    validator: *fastq.AdaptiveRecordValidator,
) SampleRecordOutcome {
    const validated = fastq.nextValidatedRecord(reader, validator) catch |err| {
        return .{ .failure = mapReaderFailure(reader, err) };
    } orelse return .done;
    if (mapCurrentSemanticFailure(
        validated.semantic_error,
        reader,
        reader.recordIndex() - 1,
    )) |failure| {
        return .{ .failure = failure };
    }
    return .{ .record = .{
        .record = validated.record,
        .canonical_span = validated.canonical_span,
    } };
}

const FileSnapshot = struct {
    inode: std.Io.File.INode,
    size: u64,
    mtime_nanoseconds: i96,
};

const ExactFirstPass = union(enum) {
    success: struct {
        snapshot: FileSnapshot,
    },
    failure: CommandFailure,
};

fn sampleExactFile(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    writer: *zfastq.Writer,
    count: u64,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    var selector = sampling.ExactSelector.init(count, options.seed);
    defer selector.deinit(allocator);

    const first = sampleExactFirstPass(io, allocator, path, &selector, options);
    const completed = switch (first) {
        .failure => |failure| return failure,
        .success => |success| success,
    };
    return switch (selector.finish()) {
        .none => null,
        .all => sampleExactSecondPass(
            true,
            io,
            allocator,
            path,
            writer,
            .empty,
            completed.snapshot,
            selector.record_count,
            options,
        ),
        .indexes => |indexes| sampleExactSecondPass(
            false,
            io,
            allocator,
            path,
            writer,
            indexes,
            completed.snapshot,
            selector.record_count,
            options,
        ),
    };
}

fn sampleExactFirstPass(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    selector: *sampling.ExactSelector,
    options: SampleOptions,
) ExactFirstPass {
    var input: RecordInput = undefined;
    const initial = switch (initExactInput(&input, io, path, null, options.output_identity)) {
        .failure => |failure| return .{ .failure = failure },
        .success => |snapshot| snapshot,
    };
    defer input.deinit(io);

    var failure: ?CommandFailure = null;
    var scanner = fastq.CheckScanner.init(
        .{ .max_line_bytes = options.max_line_bytes },
        .{ .alphabet = options.alphabet },
    );
    var buf: [SAMPLE_SCAN_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const chunk = input.readScannerChunk(&buf) catch {
            failure = IO_FAILURE;
            break;
        };
        const decoded = chunk orelse break;
        _ = scanner.feed(decoded) catch |err| {
            failure = considerExactRecords(allocator, selector, scanner.record_index);
            if (failure == null) failure = mapCheckScannerFailure(&scanner, err);
            break;
        };
        failure = considerExactRecords(allocator, selector, scanner.record_index);
        if (failure != null) break;
    }
    if (failure == null) {
        scanner.finishEof() catch |err| {
            failure = mapCheckScannerFailure(&scanner, err);
        };
        if (failure == null) {
            failure = considerExactRecords(allocator, selector, scanner.record_index);
        }
    }

    if (exactInputSnapshotFailure(&input, io, initial)) |snapshot_failure| {
        return .{ .failure = snapshot_failure };
    }
    if (failure) |details| return .{ .failure = details };
    return .{ .success = .{
        .snapshot = initial,
    } };
}

fn considerExactRecords(
    allocator: std.mem.Allocator,
    selector: *sampling.ExactSelector,
    completed_records: u64,
) ?CommandFailure {
    selector.considerRecordsThrough(allocator, completed_records) catch |err|
        return exactSelectionFailure(err);
    return null;
}

fn sampleExactSecondPass(
    comptime select_all: bool,
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    writer: *zfastq.Writer,
    selected: sampling.ExactIndexes,
    expected_snapshot: FileSnapshot,
    expected_count: u64,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    var input: RecordInput = undefined;
    switch (initExactInput(&input, io, path, expected_snapshot, options.output_identity)) {
        .failure => |failure| return failure,
        .success => {},
    }
    defer input.deinit(io);

    var reader = zfastq.Reader.init(
        allocator,
        input.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return OUT_OF_MEMORY;
    defer reader.deinit();
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    var cursor = ExactOutputCursor{
        .select_all = select_all,
        .indexes = selected,
    };
    var next_selected = cursor.nextSelected();
    var failure: ?CommandFailure = null;
    while (cursor.unit_count < expected_count) {
        if (cursor.selectNextCached(&next_selected)) {
            const validated = fastq.nextValidatedRecord(&reader, &validator) catch |err| {
                failure = mapReaderFailure(&reader, err);
                break;
            } orelse break;
            if (validated.semantic_error != null) {
                failure = inputChangedFailure();
                break;
            }
            const record = validated.record;
            if (recordHasUnwritableEnding(record, validated.canonical_span)) {
                failure = unwritableRecordFailure();
                break;
            }
            try writeCheckedRecord(writer, null, record, validated.canonical_span);
        } else {
            const advanced = reader.advance() catch |err| {
                failure = mapReaderFailure(&reader, err);
                break;
            };
            if (!advanced) break;
        }
        cursor.completeUnit();
    }
    if (failure == null and cursor.unit_count == expected_count) {
        if (reader.advance() catch |err| extra: {
            failure = mapReaderFailure(&reader, err);
            break :extra false;
        }) {
            failure = inputChangedFailure();
        }
    }

    if (exactInputSnapshotFailure(&input, io, expected_snapshot)) |snapshot_failure| {
        return snapshot_failure;
    }
    if (failure) |details| return details;
    if (cursor.unit_count != expected_count or !cursor.selectionComplete()) {
        return inputChangedFailure();
    }
    return null;
}

fn fileSnapshot(
    file: std.Io.File,
    io: std.Io,
) (std.Io.File.StatError || error{NotRegularFile})!FileSnapshot {
    const stat = try file.stat(io);
    if (stat.kind != .file) return error.NotRegularFile;
    return .{
        .inode = stat.inode,
        .size = stat.size,
        .mtime_nanoseconds = stat.mtime.nanoseconds,
    };
}

fn sameFileSnapshot(left: FileSnapshot, right: FileSnapshot) bool {
    return left.inode == right.inode and
        left.size == right.size and
        left.mtime_nanoseconds == right.mtime_nanoseconds;
}

fn inputChangedFailure() CommandFailure {
    return CommandFailure.plain(
        "input_changed",
        "input changed during exact sampling",
        3,
    );
}

const ExactPairSnapshots = union(enum) {
    paired: [2]FileSnapshot,
    interleaved: FileSnapshot,
};

const ExactPairFirstPass = union(enum) {
    success: struct {
        snapshots: ExactPairSnapshots,
    },
    failure: PairCommandFailure,
};

const ExactInput = union(enum) {
    success: FileSnapshot,
    failure: CommandFailure,
};

fn sampleExactPairs(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    writer: *zfastq.Writer,
    count: u64,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    var selector = sampling.ExactSelector.init(count, options.seed);
    defer selector.deinit(allocator);

    const first = sampleExactPairFirstPass(io, allocator, inputs, &selector, options);
    const completed = switch (first) {
        .failure => |failure| return failure,
        .success => |success| success,
    };
    return switch (selector.finish()) {
        .none => null,
        .all => sampleExactPairSecondPass(
            true,
            io,
            allocator,
            inputs,
            writer,
            .empty,
            completed.snapshots,
            selector.record_count,
            staging_limit,
            options,
        ),
        .indexes => |indexes| sampleExactPairSecondPass(
            false,
            io,
            allocator,
            inputs,
            writer,
            indexes,
            completed.snapshots,
            selector.record_count,
            staging_limit,
            options,
        ),
    };
}

fn sampleExactPairFirstPass(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    selector: *sampling.ExactSelector,
    options: SampleOptions,
) ExactPairFirstPass {
    const pair_options = PairedCheckOptions{
        .max_line_bytes = options.max_line_bytes,
        .alphabet = options.alphabet,
        .pair_mode = options.pair_mode,
        .pair_name_policy = options.pair_name_policy,
    };
    return switch (options.pair_mode) {
        .none => unreachable,
        .paired => sampleExactPairedFirstPass(
            io,
            allocator,
            inputs,
            selector,
            pair_options,
            options.output_identity,
        ),
        .interleaved => sampleExactInterleavedFirstPass(
            io,
            allocator,
            inputs[0],
            selector,
            pair_options,
            options.output_identity,
        ),
    };
}

fn sampleExactPairedFirstPass(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    selector: *sampling.ExactSelector,
    options: PairedCheckOptions,
    output_identity: ?FileIdentity,
) ExactPairFirstPass {
    var input1: RecordInput = undefined;
    var input2: RecordInput = undefined;
    const snapshots = switch (initExactPairedInputs(&input1, &input2, io, inputs, null, output_identity)) {
        .failure => |failure| return .{ .failure = failure },
        .success => |snapshots| snapshots,
    };
    defer input1.deinit(io);
    defer input2.deinit(io);

    const failure = checkPairedSources(
        allocator,
        input1.byteSource(),
        input2.byteSource(),
        options,
        selector,
    );
    if (exactPairSnapshotFailure(&input1, io, snapshots[0], 0)) |changed| {
        return .{ .failure = changed };
    }
    if (exactPairSnapshotFailure(&input2, io, snapshots[1], 1)) |changed| {
        return .{ .failure = changed };
    }
    if (failure) |details| return .{ .failure = details };
    return .{ .success = .{
        .snapshots = .{ .paired = snapshots },
    } };
}

fn sampleExactInterleavedFirstPass(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    selector: *sampling.ExactSelector,
    options: PairedCheckOptions,
    output_identity: ?FileIdentity,
) ExactPairFirstPass {
    var input: RecordInput = undefined;
    const snapshot = switch (initExactInput(&input, io, path, null, output_identity)) {
        .failure => |failure| return .{ .failure = pairCommandFailure(0, failure) },
        .success => |captured| captured,
    };
    defer input.deinit(io);

    const failure = checkInterleavedSource(
        allocator,
        input.byteSource(),
        options,
        selector,
    );
    if (exactPairSnapshotFailure(&input, io, snapshot, 0)) |changed| {
        return .{ .failure = changed };
    }
    if (failure) |details| return .{ .failure = details };
    return .{ .success = .{
        .snapshots = .{ .interleaved = snapshot },
    } };
}

fn sampleExactPairSecondPass(
    select_all: bool,
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    writer: *zfastq.Writer,
    indexes: sampling.ExactIndexes,
    snapshots: ExactPairSnapshots,
    expected_count: u64,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    return switch (snapshots) {
        .paired => |paired_snapshots| sampleExactPairedSecondPass(
            select_all,
            io,
            allocator,
            inputs,
            writer,
            indexes,
            paired_snapshots,
            expected_count,
            options,
        ),
        .interleaved => |snapshot| sampleExactInterleavedSecondPass(
            select_all,
            io,
            allocator,
            inputs[0],
            writer,
            indexes,
            snapshot,
            expected_count,
            staging_limit,
            options,
        ),
    };
}

fn sampleExactPairedSecondPass(
    select_all: bool,
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    writer: *zfastq.Writer,
    indexes: sampling.ExactIndexes,
    snapshots: [2]FileSnapshot,
    expected_count: u64,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    var input1: RecordInput = undefined;
    var input2: RecordInput = undefined;
    switch (initExactPairedInputs(&input1, &input2, io, inputs, snapshots, options.output_identity)) {
        .failure => |failure| return failure,
        .success => {},
    }
    defer input1.deinit(io);
    defer input2.deinit(io);

    var reader1 = zfastq.Reader.init(
        allocator,
        input1.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader1.deinit();
    var reader2 = zfastq.Reader.init(
        allocator,
        input2.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return pairCommandFailure(1, OUT_OF_MEMORY);
    defer reader2.deinit();
    var validator1 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var validator2 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    var cursor = ExactOutputCursor{
        .select_all = select_all,
        .indexes = indexes,
    };
    var next_selected = cursor.nextSelected();
    var failure: ?PairCommandFailure = null;
    while (cursor.unit_count < expected_count) {
        const selected = cursor.selectNextCached(&next_selected);
        var canonical_span1: ?[]const u8 = null;
        var record1: ?zfastq.Record = null;
        var semantic1: ?zfastq.SemanticError = null;
        const got1 = if (selected) selected_record: {
            const validated = fastq.nextValidatedRecord(&reader1, &validator1) catch |err| failed: {
                failure = pairCommandFailure(0, mapReaderFailure(&reader1, err));
                break :failed null;
            };
            if (validated) |value| {
                record1 = value.record;
                canonical_span1 = value.canonical_span;
                semantic1 = value.semantic_error;
            }
            break :selected_record record1 != null;
        } else reader1.advance() catch |err| failed: {
            failure = pairCommandFailure(0, mapReaderFailure(&reader1, err));
            break :failed false;
        };
        if (failure != null) break;

        var canonical_span2: ?[]const u8 = null;
        var record2: ?zfastq.Record = null;
        var semantic2: ?zfastq.SemanticError = null;
        const got2 = if (selected) selected_record: {
            const validated = fastq.nextValidatedRecord(&reader2, &validator2) catch |err| failed: {
                failure = pairCommandFailure(1, mapReaderFailure(&reader2, err));
                break :failed null;
            };
            if (validated) |value| {
                record2 = value.record;
                canonical_span2 = value.canonical_span;
                semantic2 = value.semantic_error;
            }
            break :selected_record record2 != null;
        } else reader2.advance() catch |err| failed: {
            failure = pairCommandFailure(1, mapReaderFailure(&reader2, err));
            break :failed false;
        };
        if (failure != null) break;

        if (!got1 or !got2) {
            failure = if (got1 or got2)
                inputChangedPairFailure(if (got1) 0 else 1)
            else
                inputChangedPairFailure(0);
            break;
        }
        if (selected) {
            if (semantic1 != null) {
                failure = inputChangedPairFailure(0);
                break;
            }
            if (semantic2 != null) {
                failure = inputChangedPairFailure(1);
                break;
            }
            if (!pairing.headersMatch(record1.?.header, record2.?.header, options.pair_name_policy)) {
                failure = inputChangedPairFailure(0);
                break;
            }
            if (recordHasUnwritableEnding(record1.?, canonical_span1)) {
                failure = pairCommandFailure(0, unwritableRecordFailure());
                break;
            }
            if (recordHasUnwritableEnding(record2.?, canonical_span2)) {
                failure = pairCommandFailure(1, unwritableRecordFailure());
                break;
            }
            try writeCheckedRecord(writer, null, record1.?, canonical_span1);
            try writeCheckedRecord(writer, null, record2.?, canonical_span2);
        }
        cursor.completeUnit();
    }

    if (failure == null and cursor.unit_count == expected_count) {
        const extra1 = reader1.advance() catch |err| failed: {
            failure = pairCommandFailure(0, mapReaderFailure(&reader1, err));
            break :failed false;
        };
        if (failure == null) {
            const extra2 = reader2.advance() catch |err| failed: {
                failure = pairCommandFailure(1, mapReaderFailure(&reader2, err));
                break :failed false;
            };
            if (failure == null and (extra1 or extra2)) {
                failure = inputChangedPairFailure(if (extra1 != extra2)
                    (if (extra1) 0 else 1)
                else
                    0);
            }
        }
    }
    if (exactPairSnapshotFailure(&input1, io, snapshots[0], 0)) |changed| return changed;
    if (exactPairSnapshotFailure(&input2, io, snapshots[1], 1)) |changed| return changed;
    if (failure) |details| return details;
    if (cursor.unit_count != expected_count or !cursor.selectionComplete()) {
        return inputChangedPairFailure(0);
    }
    return null;
}

fn sampleExactInterleavedSecondPass(
    select_all: bool,
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    writer: *zfastq.Writer,
    indexes: sampling.ExactIndexes,
    snapshot: FileSnapshot,
    expected_count: u64,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    var input: RecordInput = undefined;
    switch (initExactInput(&input, io, path, snapshot, options.output_identity)) {
        .failure => |failure| return pairCommandFailure(0, failure),
        .success => {},
    }
    defer input.deinit(io);

    var reader = zfastq.Reader.init(
        allocator,
        input.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader.deinit();
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var staged_record: std.ArrayList(u8) = .empty;
    defer staged_record.deinit(allocator);
    var retained_record_storage: fastq.RetainedRecordStorage = .{};
    defer retained_record_storage.deinit(allocator);

    var cursor = ExactOutputCursor{
        .select_all = select_all,
        .indexes = indexes,
    };
    var next_selected = cursor.nextSelected();
    var failure: ?PairCommandFailure = null;
    while (cursor.unit_count < expected_count) {
        const selected = cursor.selectNextCached(&next_selected);
        if (!selected) {
            const mate1 = reader.advance() catch |err| failed: {
                failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
                break :failed false;
            };
            if (failure != null or !mate1) break;
            const mate2 = reader.advance() catch |err| failed: {
                failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
                break :failed false;
            };
            if (failure != null or !mate2) break;
            cursor.completeUnit();
            continue;
        }

        var validated1 = fastq.nextValidatedRecord(&reader, &validator) catch |err| failed: {
            failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
            break :failed null;
        } orelse break;

        const invalid1 = validated1.semantic_error != null;
        const header1_len = validated1.record.header.len;
        const unwritable1 = recordHasUnwritableEnding(validated1.record, validated1.canonical_span);
        var record1_storage: InterleavedFirstRecordStorage = .reader;
        var paired_record2 = fastq.nextBufferedValidatedRecord(
            &reader,
            &validator,
        ) catch |err| failed: {
            failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
            break :failed null;
        };
        // HACK: Preserve retry diagnostics until the exact-output bug is fixed.
        if (paired_record2 == null and failure == null) {
            paired_record2 = fastq.nextPairedValidatedRecord(
                &reader,
                &validated1,
                &validator,
            ) catch |err| {
                failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
                break;
            };
        }
        const record1 = validated1.record;
        const canonical_span1 = validated1.canonical_span;
        const validated2 = paired_record2 orelse record: {
            const preserved = nextAfterPreservingInterleavedMate1(
                allocator,
                &reader,
                &retained_record_storage,
                &staged_record,
                record1,
                canonical_span1,
                staging_limit,
                &validator,
            ) catch |err| failed: {
                failure = mapPreservedMateFailure(&reader, err);
                break :failed null;
            };
            if (failure != null) break :record null;
            record1_storage = preserved.?.first_storage;
            break :record preserved.?.record;
        } orelse break;

        const record2 = validated2.record;
        const canonical_span2 = validated2.canonical_span;
        if (invalid1 or validated2.semantic_error != null) {
            failure = inputChangedPairFailure(0);
            break;
        }
        const header1 = switch (record1_storage) {
            .unused => unreachable,
            .reader, .retained => record1.header,
            .staged => staged_record.items[1 .. 1 + header1_len],
        };
        if (!pairing.headersMatch(header1, record2.header, options.pair_name_policy)) {
            failure = inputChangedPairFailure(0);
            break;
        }
        if (unwritable1 or recordHasUnwritableEnding(record2, canonical_span2)) {
            failure = pairCommandFailure(0, unwritableRecordFailure());
            break;
        }
        try writePreservedInterleavedMate1(
            writer,
            &reader,
            &retained_record_storage,
            staged_record.items,
            record1,
            canonical_span1,
            record1_storage,
        );
        try writeCheckedRecord(writer, null, record2, canonical_span2);
        cursor.completeUnit();
    }

    if (failure == null and cursor.unit_count == expected_count) {
        const extra_mate1 = reader.advance() catch |err| failed: {
            failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
            break :failed false;
        };
        if (failure == null and extra_mate1) {
            _ = reader.advance() catch |err| failed: {
                failure = pairCommandFailure(0, mapReaderFailure(&reader, err));
                break :failed false;
            };
            if (failure == null) failure = inputChangedPairFailure(0);
        }
    }
    if (exactPairSnapshotFailure(&input, io, snapshot, 0)) |changed| return changed;
    if (failure) |details| return details;
    if (cursor.unit_count != expected_count or !cursor.selectionComplete()) {
        return inputChangedPairFailure(0);
    }
    return null;
}

fn initExactInput(
    input: *RecordInput,
    io: std.Io,
    path: []const u8,
    expected_snapshot: ?FileSnapshot,
    output_identity: ?FileIdentity,
) ExactInput {
    const opened = switch (openExactInput(io, path, expected_snapshot)) {
        .failure => |failure| return .{ .failure = failure },
        .success => |opened| opened,
    };
    var transferred = false;
    defer if (!transferred) opened.file.close(io);
    if (outputFileFailure(io, opened.file, output_identity)) |failure| {
        if (expected_snapshot != null and failure.exit_code == 2) {
            var changed = inputChangedFailure();
            changed.suppress_diagnostic = failure.suppress_diagnostic;
            return .{ .failure = changed };
        }
        return .{ .failure = failure };
    }
    input.init(io, opened.file, true) catch
        return .{ .failure = IO_FAILURE };
    transferred = true;
    return .{ .success = opened.snapshot };
}

const ExactFile = union(enum) {
    success: struct { file: std.Io.File, snapshot: FileSnapshot },
    failure: CommandFailure,
};

fn openExactInput(io: std.Io, path: []const u8, expected_snapshot: ?FileSnapshot) ExactFile {
    // A FIFO must reach the descriptor type check without waiting for a writer.
    const handle = std.posix.openat(std.Io.Dir.cwd().handle, path, .{
        .ACCMODE = .RDONLY,
        .NONBLOCK = true,
        .CLOEXEC = true,
        .NOCTTY = true,
    }, 0) catch |err| {
        return .{ .failure = if (expected_snapshot != null and err == error.FileNotFound)
            inputChangedFailure()
        else if (expected_snapshot != null)
            CommandFailure.plain("io_error", "failed to reopen file", 3)
        else if (err == error.FileNotFound)
            CommandFailure.plain("io_error", "file not found", 3)
        else
            CommandFailure.plain("io_error", "failed to open file", 3) };
    };
    const file: std.Io.File = .{ .handle = handle, .flags = .{ .nonblocking = true } };
    const snapshot = fileSnapshot(file, io) catch |err| {
        file.close(io);
        return .{ .failure = if (expected_snapshot != null and err == error.NotRegularFile)
            inputChangedFailure()
        else if (err == error.NotRegularFile)
            CommandFailure.plain(
                "io_error",
                "exact-count sampling requires a regular file",
                3,
            )
        else
            CommandFailure.plain("io_error", "failed to inspect file", 3) };
    };
    if (expected_snapshot) |expected| {
        if (!sameFileSnapshot(expected, snapshot)) {
            file.close(io);
            return .{ .failure = inputChangedFailure() };
        }
    }
    return .{ .success = .{ .file = file, .snapshot = snapshot } };
}

const ExactPairedInputs = union(enum) {
    success: [2]FileSnapshot,
    failure: PairCommandFailure,
};

fn initExactPairedInputs(
    input1: *RecordInput,
    input2: *RecordInput,
    io: std.Io,
    inputs: []const []const u8,
    expected_snapshots: ?[2]FileSnapshot,
    output_identity: ?FileIdentity,
) ExactPairedInputs {
    var transferred = false;
    const first = switch (openExactInput(io, inputs[0], if (expected_snapshots) |expected| expected[0] else null)) {
        .failure => |failure| return .{ .failure = pairCommandFailure(0, failure) },
        .success => |opened| opened,
    };
    defer if (!transferred) first.file.close(io);
    const second = switch (openExactInput(io, inputs[1], if (expected_snapshots) |expected| expected[1] else null)) {
        .failure => |failure| return .{ .failure = pairCommandFailure(1, failure) },
        .success => |opened| opened,
    };
    defer if (!transferred) second.file.close(io);

    if (pairedFilesFailure(io, first.file, second.file, output_identity)) |failure| {
        if (expected_snapshots != null and failure.exitCode() == 2) {
            var changed = inputChangedPairFailure(failure.command.input_index);
            changed.command.details.suppress_diagnostic = failure.command.details.suppress_diagnostic;
            return .{ .failure = changed };
        }
        return .{ .failure = failure };
    }
    input1.init(io, first.file, true) catch
        return .{ .failure = pairCommandFailure(0, IO_FAILURE) };
    input2.init(io, second.file, true) catch
        return .{ .failure = pairCommandFailure(1, IO_FAILURE) };
    transferred = true;
    return .{ .success = .{ first.snapshot, second.snapshot } };
}

fn exactPairSnapshotFailure(
    input: *RecordInput,
    io: std.Io,
    expected: FileSnapshot,
    input_index: u1,
) ?PairCommandFailure {
    const failure = exactInputSnapshotFailure(input, io, expected) orelse return null;
    return pairCommandFailure(input_index, failure);
}

fn exactInputSnapshotFailure(
    input: *RecordInput,
    io: std.Io,
    expected: FileSnapshot,
) ?CommandFailure {
    const final = fileSnapshot(input.file, io) catch {
        return CommandFailure.plain("io_error", "failed to inspect file", 3);
    };
    if (!sameFileSnapshot(expected, final)) return inputChangedFailure();
    return null;
}

fn inputChangedPairFailure(input_index: u1) PairCommandFailure {
    return pairCommandFailure(input_index, inputChangedFailure());
}

// --- Interleave command ---

const InterleaveOptions = struct {
    max_line_bytes: usize,
    alphabet: zfastq.Alphabet,
    pair_name_policy: pairing.NamePolicy,
    selection: PairOutputSelector = .all,
    output_identity: ?FileIdentity = null,
};

fn runInterleave(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    parsed_options: InterleaveOptions,
) u8 {
    if (validateInterleaveInputs(io, inputs)) |exit_code| return exit_code;

    var options = parsed_options;
    options.output_identity = stdoutFileIdentity(io) catch {
        if (canReportInputFailure(io, null)) printOutputFailure(io);
        return 3;
    };
    var stdout_buffer: [64 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var sink_adapter = io_layer.WriterSink.init(&stdout_writer.interface);
    var writer = zfastq.Writer.init(sink_adapter.byteSink());

    const failure = interleaveInputs(
        io,
        allocator,
        inputs,
        &writer,
        &stdout_writer.interface,
        options,
    ) catch {
        writer.flush() catch {};
        printOutputFailure(io);
        return 3;
    };
    if (failure) |details| {
        printPairCommandFailure(io, inputs, .paired, details);
    }

    var exit_code: u8 = if (failure) |details| details.exitCode() else 0;
    writer.flush() catch {
        printOutputFailure(io);
        exit_code = @max(exit_code, 3);
    };
    return exit_code;
}

fn validateInterleaveInputs(io: std.Io, inputs: []const []const u8) ?u8 {
    if (inputs.len != 2) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: interleave requires exactly two inputs\n",
        ) catch {};
        return 2;
    }
    if (std.mem.eql(u8, inputs[0], "-") and std.mem.eql(u8, inputs[1], "-")) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: interleave inputs may contain standard input at most once\n",
        ) catch {};
        return 2;
    }
    return null;
}

fn interleaveInputs(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    writer: *zfastq.Writer,
    direct_writer: ?*std.Io.Writer,
    options: InterleaveOptions,
) error{WriteFailed}!?PairCommandFailure {
    var input1: RecordInput = undefined;
    var input2: RecordInput = undefined;
    if (initPairedRecordInputs(&input1, &input2, io, inputs, options.output_identity)) |failure| return failure;
    defer input1.deinit(io);
    defer input2.deinit(io);

    var selection = options.selection;
    if (comptime build_options.use_isa_l) {
        return interleaveSources(
            allocator,
            input1.byteSource(),
            input2.byteSource(),
            writer,
            direct_writer,
            &selection,
            options,
        );
    }
    return interleaveSources(
        allocator,
        &input1,
        &input2,
        writer,
        direct_writer,
        &selection,
        options,
    );
}

fn interleaveSources(
    allocator: std.mem.Allocator,
    source1: anytype,
    source2: @TypeOf(source1),
    writer: *zfastq.Writer,
    direct_writer: ?*std.Io.Writer,
    selection: *PairOutputSelector,
    options: InterleaveOptions,
) error{WriteFailed}!?PairCommandFailure {
    const reader_options: fastq.Options = .{ .max_line_bytes = options.max_line_bytes };
    var reader1 = initSourceReader(allocator, source1, reader_options) catch
        return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader1.deinit();
    var reader2 = initSourceReader(allocator, source2, reader_options) catch
        return pairCommandFailure(1, OUT_OF_MEMORY);
    defer reader2.deinit();
    var validator1 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var validator2 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        const validated1 = fastq.nextValidatedRecord(&reader1, &validator1) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader1, err));
        };
        const validated2 = fastq.nextValidatedRecord(&reader2, &validator2) catch |err| {
            return pairCommandFailure(1, mapReaderFailure(&reader2, err));
        };
        const record1 = if (validated1) |value| value.record else null;
        const record2 = if (validated2) |value| value.record else null;
        const canonical_span1 = if (validated1) |value| value.canonical_span else null;
        const canonical_span2 = if (validated2) |value| value.canonical_span else null;

        if (record1 == null and record2 == null) return null;
        if (record1 == null or record2 == null) {
            const remaining_side: u1 = if (record1 != null) 0 else 1;
            const pair_index = if (record1 != null)
                reader1.recordIndex() - 1
            else
                reader2.recordIndex() - 1;
            return .{ .pair = .{ .count_mismatch = .{
                .pair_index = pair_index,
                .remaining_side = remaining_side,
                .record_indexes = .{
                    lastRecordIndex(&reader1),
                    lastRecordIndex(&reader2),
                },
            } } };
        }

        const record_index1 = reader1.recordIndex() - 1;
        const record_index2 = reader2.recordIndex() - 1;
        const offsets1 = reader1.currentRecordOffsets().?;
        const offsets2 = reader2.currentRecordOffsets().?;
        if (mapSemanticFailure(
            validated1.?.semantic_error,
            offsets1,
            record_index1,
        )) |failure| {
            return pairCommandFailure(0, failure);
        }
        if (mapSemanticFailure(
            validated2.?.semantic_error,
            offsets2,
            record_index2,
        )) |failure| {
            return pairCommandFailure(1, failure);
        }

        if (!pairing.headersMatch(
            record1.?.header,
            record2.?.header,
            options.pair_name_policy,
        )) {
            const name1 = pairing.parseName(record1.?.header, options.pair_name_policy);
            const name2 = pairing.parseName(record2.?.header, options.pair_name_policy);
            return .{ .pair = .{ .name_mismatch = .{
                .pair_index = record_index1,
                .records = .{
                    .init(name1, record_index1, offsets1.header),
                    .init(name2, record_index2, offsets2.header),
                },
            } } };
        }

        if (!selection.selectPair()) continue;
        if (recordHasUnwritableEnding(record1.?, canonical_span1)) {
            return pairCommandFailure(0, unwritableRecordFailure());
        }
        if (recordHasUnwritableEnding(record2.?, canonical_span2)) {
            return pairCommandFailure(1, unwritableRecordFailure());
        }

        try writeCheckedRecord(writer, direct_writer, record1.?, canonical_span1);
        try writeCheckedRecord(writer, direct_writer, record2.?, canonical_span2);
    }
}

fn printOutputFailure(io: std.Io) void {
    std.Io.File.writeStreamingAll(
        .stderr(),
        io,
        "error: standard output: I/O error\n",
    ) catch {};
}

// --- Deinterleave command ---

const DeinterleaveOptions = struct {
    max_line_bytes: usize,
    alphabet: zfastq.Alphabet,
    pair_name_policy: pairing.NamePolicy,
};

const DeinterleaveWriteError = error{
    Output1WriteFailed,
    Output2WriteFailed,
};

const DeinterleaveOutput = struct {
    path: []const u8,
    file: ?std.Io.File = null,
    write_buffer: [64 * 1024]u8 = undefined,
    sink: io_layer.FileSink = undefined,
    writer: zfastq.Writer = undefined,

    fn init(path: []const u8) DeinterleaveOutput {
        return .{ .path = path };
    }

    fn create(self: *DeinterleaveOutput, io: std.Io) ?CommandFailure {
        const file = std.Io.Dir.cwd().createFile(io, self.path, .{ .exclusive = true }) catch |err| {
            return if (err == error.PathAlreadyExists)
                CommandFailure.plain("io_error", "output path already exists", 3)
            else
                CommandFailure.plain("io_error", "failed to create output", 3);
        };
        self.file = file;
        self.sink = io_layer.FileSink.init(io, file, &self.write_buffer);
        self.writer = zfastq.Writer.init(self.sink.byteSink());
        return null;
    }

    fn close(self: *DeinterleaveOutput, io: std.Io) void {
        if (self.file) |file| file.close(io);
        self.file = null;
    }
};

fn runDeinterleave(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    output1_path: ?[]const u8,
    output2_path: ?[]const u8,
    options: DeinterleaveOptions,
) u8 {
    const paths = validateDeinterleaveArguments(io, inputs, output1_path, output2_path) orelse
        return 2;
    const staging_limit = deinterleaveStagingLimit(options.max_line_bytes) catch {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: deinterleave record staging size exceeds supported limit\n",
        ) catch {};
        return 4;
    };

    const input_label = inputs[0];
    const owns_input = !std.mem.eql(u8, input_label, "-");
    const input_file = openRecordFile(io, input_label) catch |err| {
        const failure = inputOpenFailure(err);
        printCommandFailure(io, input_label, failure);
        return failure.exit_code;
    };
    defer if (owns_input) input_file.close(io);

    var output1 = DeinterleaveOutput.init(paths[0]);
    defer output1.close(io);
    var output2 = DeinterleaveOutput.init(paths[1]);
    defer output2.close(io);
    if (output1.create(io)) |failure| {
        printCommandFailure(io, output1.path, failure);
        return failure.exit_code;
    }
    if (output2.create(io)) |failure| {
        printCommandFailure(io, output2.path, failure);
        return failure.exit_code;
    }

    var input: RecordInput = undefined;
    input.init(io, input_file, false) catch {
        const failure = IO_FAILURE;
        printCommandFailure(io, input_label, failure);
        return failure.exit_code;
    };
    defer input.deinit(io);

    const failure = deinterleaveSource(
        allocator,
        input.byteSource(),
        &output1.writer,
        &output2.writer,
        staging_limit,
        options,
    ) catch |err| {
        const failed_output = if (err == error.Output1WriteFailed) &output1 else &output2;
        printPathError(io, failed_output.path, "I/O error");
        return 3;
    };
    if (failure) |details| {
        printPairCommandFailure(io, inputs, .interleaved, details);
    }

    var exit_code: u8 = if (failure) |details| details.exitCode() else 0;
    flushDeinterleaveWriters(&output1.writer, &output2.writer) catch |err| {
        const failed_output = if (err == error.Output1WriteFailed) &output1 else &output2;
        printPathError(io, failed_output.path, "I/O error");
        exit_code = @max(exit_code, 3);
    };
    return exit_code;
}

fn validateDeinterleaveArguments(
    io: std.Io,
    inputs: []const []const u8,
    output1_path: ?[]const u8,
    output2_path: ?[]const u8,
) ?[2][]const u8 {
    if (inputs.len != 1) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: deinterleave requires exactly one input\n",
        ) catch {};
        return null;
    }
    const path1 = output1_path orelse {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: deinterleave requires --out1\n") catch {};
        return null;
    };
    const path2 = output2_path orelse {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: deinterleave requires --out2\n") catch {};
        return null;
    };
    if (std.mem.eql(u8, path1, "-") or std.mem.eql(u8, path2, "-")) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: deinterleave output paths cannot be standard output\n",
        ) catch {};
        return null;
    }
    if (std.mem.eql(u8, path1, path2)) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: deinterleave output paths must differ\n",
        ) catch {};
        return null;
    }
    return .{ path1, path2 };
}

fn deinterleaveStagingLimit(max_line_bytes: usize) error{ArithmeticLimit}!usize {
    const fields = std.math.mul(usize, max_line_bytes, 4) catch
        return error.ArithmeticLimit;
    return std.math.add(usize, fields, 4) catch error.ArithmeticLimit;
}

fn deinterleaveSource(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    writer1: *zfastq.Writer,
    writer2: *zfastq.Writer,
    staging_limit: usize,
    options: DeinterleaveOptions,
) DeinterleaveWriteError!?PairCommandFailure {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return pairCommandFailure(0, OUT_OF_MEMORY);
    defer reader.deinit();
    var staged_record: std.ArrayList(u8) = .empty;
    defer staged_record.deinit(allocator);
    var retained_record_storage: fastq.RetainedRecordStorage = .{};
    defer retained_record_storage.deinit(allocator);
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        var validated1 = fastq.nextValidatedRecord(&reader, &validator) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader, err));
        } orelse return null;
        const record_index1 = reader.recordIndex() - 1;
        const offsets1 = reader.currentRecordOffsets().?;
        const header1_len = validated1.record.header.len;
        const semantic1 = mapSemanticFailure(
            validated1.semantic_error,
            offsets1,
            record_index1,
        );

        const unwritable1 = recordHasUnwritableEnding(validated1.record, validated1.canonical_span);
        var record1_storage: InterleavedFirstRecordStorage = .reader;
        const paired_record2 = (if (semantic1 == null)
            fastq.nextPairedValidatedRecord(&reader, &validated1, &validator)
        else
            fastq.nextBufferedValidatedRecord(&reader, &validator)) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader, err));
        };
        const record1 = validated1.record;
        const canonical_span1 = validated1.canonical_span;
        const validated2 = paired_record2 orelse record: {
            @branchHint(.cold);
            if (semantic1 != null) {
                break :record fastq.nextValidatedRecord(
                    &reader,
                    &validator,
                ) catch |err| {
                    return pairCommandFailure(0, mapReaderFailure(&reader, err));
                } orelse return missingInterleavedMateFailure(record_index1);
            }

            const preserved = nextAfterPreservingInterleavedMate1(
                allocator,
                &reader,
                &retained_record_storage,
                &staged_record,
                record1,
                canonical_span1,
                staging_limit,
                &validator,
            ) catch |err| return mapPreservedMateFailure(&reader, err);
            record1_storage = preserved.first_storage;
            break :record preserved.record orelse return missingInterleavedMateFailure(record_index1);
        };

        const record2 = validated2.record;
        const canonical_span2 = validated2.canonical_span;
        if (semantic1) |details| {
            return pairCommandFailure(0, details);
        }

        const record_index2 = reader.recordIndex() - 1;
        const offsets2 = reader.currentRecordOffsets().?;
        if (mapSemanticFailure(
            validated2.semantic_error,
            offsets2,
            record_index2,
        )) |details| {
            return pairCommandFailure(0, details);
        }

        const header1 = switch (record1_storage) {
            .unused => unreachable,
            .reader, .retained => record1.header,
            .staged => staged_record.items[1 .. 1 + header1_len],
        };
        if (!pairing.headersMatch(header1, record2.header, options.pair_name_policy)) {
            const name1 = pairing.parseName(header1, options.pair_name_policy);
            const name2 = pairing.parseName(record2.header, options.pair_name_policy);
            return .{ .pair = .{ .name_mismatch = .{
                .pair_index = record_index1 / 2,
                .records = .{
                    .init(name1, record_index1, offsets1.header),
                    .init(name2, record_index2, offsets2.header),
                },
            } } };
        }

        if (unwritable1 or recordHasUnwritableEnding(record2, canonical_span2)) {
            return pairCommandFailure(0, unwritableRecordFailure());
        }
        writePreservedInterleavedMate1(
            writer1,
            &reader,
            &retained_record_storage,
            staged_record.items,
            record1,
            canonical_span1,
            record1_storage,
        ) catch return error.Output1WriteFailed;
        writeCheckedRecord(writer2, null, record2, canonical_span2) catch
            return error.Output2WriteFailed;
    }
}

fn canonicalRecordSize(
    record: zfastq.Record,
    canonical_span: ?[]const u8,
) error{RecordStagingLimit}!usize {
    return if (canonical_span) |span|
        span.len
    else blk: {
        var total: usize = 6;
        inline for (&.{ record.header, record.sequence, record.plus, record.quality }) |field| {
            total = std.math.add(usize, total, field.len) catch
                return error.RecordStagingLimit;
        }
        break :blk total;
    };
}

fn stageCanonicalRecord(
    allocator: std.mem.Allocator,
    staging: *std.ArrayList(u8),
    record: zfastq.Record,
    canonical_span: ?[]const u8,
    staging_limit: usize,
) error{ OutOfMemory, RecordStagingLimit }!void {
    const required = try canonicalRecordSize(record, canonical_span);
    if (required > staging_limit) return error.RecordStagingLimit;

    if (staging.capacity > required and
        staging.capacity - required > zfastq.limits.DEFAULT_READER_BUFFER_BYTES)
    {
        staging.clearAndFree(allocator);
    } else {
        staging.clearRetainingCapacity();
    }
    staging.ensureTotalCapacityPrecise(allocator, required) catch return error.OutOfMemory;
    if (canonical_span) |span| {
        staging.appendSliceAssumeCapacity(span);
        return;
    }
    staging.appendAssumeCapacity('@');
    staging.appendSliceAssumeCapacity(record.header);
    staging.appendAssumeCapacity('\n');
    staging.appendSliceAssumeCapacity(record.sequence);
    staging.appendSliceAssumeCapacity("\n+");
    staging.appendSliceAssumeCapacity(record.plus);
    staging.appendAssumeCapacity('\n');
    staging.appendSliceAssumeCapacity(record.quality);
    staging.appendAssumeCapacity('\n');
}

fn flushDeinterleaveWriters(
    writer1: *zfastq.Writer,
    writer2: *zfastq.Writer,
) DeinterleaveWriteError!void {
    writer1.flush() catch return error.Output1WriteFailed;
    writer2.flush() catch return error.Output2WriteFailed;
}

// --- Machine output ---

fn beginJsonDocument(json: *std.json.Stringify, schema: []const u8) !void {
    try json.beginObject();
    try json.objectField("schema");
    try json.write(schema);
    try json.objectField("tool");
    try json.beginObject();
    try json.objectField("name");
    try json.write("z-fastq");
    try json.objectField("version");
    try json.write(zfastq.VERSION);
    try json.endObject();
    try json.objectField("byte_strings");
    try json.write("escaped-bytes-v1");
    try json.objectField("results");
    try json.beginArray();
}

fn finishJsonDocument(json: *std.json.Stringify) !void {
    try json.endArray();
    try json.endObject();
    try json.writer.writeByte('\n');
}

fn writeCheckJsonResult(
    json: *std.json.Stringify,
    input: []const u8,
    failure: ?CommandFailure,
) !void {
    try json.beginObject();
    try json.objectField("input");
    try writeEscapedJsonString(json, input);
    try json.objectField("status");
    try json.write(if (failure == null) "ok" else "error");
    if (failure) |details| try writeJsonFailure(json, details);
    try json.endObject();
}

fn writePairedCheckJsonResult(
    json: *std.json.Stringify,
    inputs: []const []const u8,
    pair_mode: PairMode,
    failure: ?PairCommandFailure,
) !void {
    try json.beginObject();
    if (pair_mode == .paired) {
        try json.objectField("inputs");
        try json.beginArray();
        for (inputs) |input| try writeEscapedJsonString(json, input);
        try json.endArray();
    } else {
        try json.objectField("input");
        try writeEscapedJsonString(json, inputs[0]);
    }
    try json.objectField("status");
    try json.write(if (failure == null) "ok" else "error");
    if (failure) |details| switch (details) {
        .command => |command_failure| {
            if (pair_mode == .paired) {
                try json.objectField("failed_input");
                try writeEscapedJsonString(json, inputs[command_failure.input_index]);
            }
            try writeJsonFailure(json, command_failure.details);
        },
        .pair => |pair_failure| try writePairJsonFailure(json, pair_failure),
    };
    try json.endObject();
}

fn writePairJsonFailure(json: *std.json.Stringify, failure: PairFailure) !void {
    try json.objectField("error");
    try json.beginObject();
    switch (failure) {
        .name_mismatch => |details| {
            try json.objectField("code");
            try json.write("P001");
            try json.objectField("message");
            try json.write("paired identifiers or mate markers do not match");
            try json.objectField("pair_index");
            try json.write(details.pair_index);
            try json.objectField("record_indexes");
            try json.beginArray();
            for (details.records) |record| try json.write(record.record_index);
            try json.endArray();
            try json.objectField("byte_offsets");
            try json.beginArray();
            for (details.records) |record| try json.write(record.byte_offset);
            try json.endArray();
            try json.objectField("first_tokens");
            try json.beginArray();
            for (details.records) |record| try writeBoundedJsonBytes(json, &record.first_token);
            try json.endArray();
            try json.objectField("normalized_ids");
            try json.beginArray();
            for (details.records) |record| try writeBoundedJsonBytes(json, &record.normalized_id);
            try json.endArray();
            try json.objectField("mate_markers");
            try json.beginArray();
            for (details.records) |record| try writeMateMarkersJson(json, record.mate_markers);
            try json.endArray();
        },
        .count_mismatch => |details| {
            try json.objectField("code");
            try json.write("P002");
            try json.objectField("message");
            try json.write("paired input is missing a mate");
            try json.objectField("pair_index");
            try json.write(details.pair_index);
            try json.objectField("remaining_side");
            try json.write(if (details.remaining_side == 0) "R1" else "R2");
            try json.objectField("record_indexes");
            try json.beginArray();
            for (details.record_indexes) |record_index| try json.write(record_index);
            try json.endArray();
        },
    }
    try json.endObject();
}

fn writeBoundedJsonBytes(json: *std.json.Stringify, display: *const BoundedBytes) !void {
    try json.beginObject();
    try json.objectField("prefix");
    try writeEscapedJsonString(json, display.bytes());
    try json.objectField("length");
    try json.write(display.full_len);
    try json.objectField("truncated");
    try json.write(display.truncated());
    try json.endObject();
}

fn writeMateMarkersJson(json: *std.json.Stringify, markers: u2) !void {
    try json.beginArray();
    if (markers & 0b01 != 0) try json.write(@as(u8, 1));
    if (markers & 0b10 != 0) try json.write(@as(u8, 2));
    try json.endArray();
}

fn writeStatsJsonResult(
    json: *std.json.Stringify,
    input: []const u8,
    outcome: StatsOutcome,
) !void {
    try json.beginObject();
    try json.objectField("input");
    try writeEscapedJsonString(json, input);
    switch (outcome) {
        .failure => |failure| {
            try json.objectField("status");
            try json.write("error");
            try writeJsonFailure(json, failure);
        },
        .success => |stats| {
            const result = stats.result();
            try json.objectField("status");
            try json.write("ok");
            try json.objectField("reads");
            try json.write(result.reads);
            try json.objectField("bases");
            try json.write(result.bases);
            try json.objectField("min_length");
            try json.write(result.min_length);
            try json.objectField("max_length");
            try json.write(result.max_length);
            try json.objectField("mean_length");
            try json.write(result.mean_length);
            try json.objectField("a");
            try json.write(result.a);
            try json.objectField("c");
            try json.write(result.c);
            try json.objectField("g");
            try json.write(result.g);
            try json.objectField("t");
            try json.write(result.t);
            try json.objectField("n");
            try json.write(result.n);
            try json.objectField("other_bases");
            try json.write(result.other_bases);
            try json.objectField("gc_fraction");
            try json.write(result.gc_fraction);
            try json.objectField("quality_sum");
            try json.write(result.quality_sum);
            try json.objectField("mean_quality");
            try json.write(result.mean_quality);
            try json.objectField("q20_bases");
            try json.write(result.q20_bases);
            try json.objectField("q20_fraction");
            try json.write(result.q20_fraction);
            try json.objectField("q30_bases");
            try json.write(result.q30_bases);
            try json.objectField("q30_fraction");
            try json.write(result.q30_fraction);
        },
    }
    try json.endObject();
}

fn writeJsonFailure(json: *std.json.Stringify, failure: CommandFailure) !void {
    try json.objectField("error");
    try json.beginObject();
    try json.objectField("code");
    try json.write(failure.code);
    try json.objectField("message");
    try json.write(failure.message);
    try json.objectField("record_index");
    try json.write(failure.record_index);
    try json.objectField("byte_offset");
    try json.write(failure.byte_offset);
    try json.objectField("line_in_record");
    try json.write(failure.line_in_record);
    try json.endObject();
}

fn writeEscapedJsonString(json: *std.json.Stringify, bytes: []const u8) !void {
    try json.beginWriteRaw();
    try json.writer.writeByte('"');
    try writeEscapedBytes(json.writer, bytes, true);
    try json.writer.writeByte('"');
    json.endWriteRaw();
}

fn printStats(
    io: std.Io,
    label: []const u8,
    result: zfastq.StatsResult,
    separator: bool,
) !void {
    if (separator) try std.Io.File.writeStreamingAll(.stdout(), io, "\n");
    try std.Io.File.writeStreamingAll(.stdout(), io, "input: ");
    try writeEscapedAll(.stdout(), io, label);
    try std.Io.File.writeStreamingAll(.stdout(), io, "\n");
    var buffer: [1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &buffer);
    const output = &stdout_writer.interface;
    try writeUnsignedField(output, "reads", result.reads);
    try writeUnsignedField(output, "bases", result.bases);
    try writeOptionalUnsignedField(output, "min_length", result.min_length);
    try writeOptionalUnsignedField(output, "max_length", result.max_length);
    try writeRatioField(output, "mean_length", result.bases, result.reads);
    try writeUnsignedField(output, "a", result.a);
    try writeUnsignedField(output, "c", result.c);
    try writeUnsignedField(output, "g", result.g);
    try writeUnsignedField(output, "t", result.t);
    try writeUnsignedField(output, "n", result.n);
    try writeUnsignedField(output, "other_bases", result.other_bases);
    try writeRatioField(
        output,
        "gc_fraction",
        @as(u128, result.g) + result.c,
        @as(u128, result.a) + result.c + result.g + result.t,
    );
    try writeUnsignedField(output, "quality_sum", result.quality_sum);
    try writeRatioField(output, "mean_quality", result.quality_sum, result.bases);
    try writeUnsignedField(output, "q20_bases", result.q20_bases);
    try writeRatioField(output, "q20_fraction", result.q20_bases, result.bases);
    try writeUnsignedField(output, "q30_bases", result.q30_bases);
    try writeRatioField(output, "q30_fraction", result.q30_bases, result.bases);
    try output.flush();
}

fn writeUnsignedField(output: *std.Io.Writer, name: []const u8, value: u64) !void {
    try output.print("{s}: {d}\n", .{ name, value });
}

fn writeOptionalUnsignedField(output: *std.Io.Writer, name: []const u8, value: ?u64) !void {
    if (value) |number| return writeUnsignedField(output, name, number);
    try output.print("{s}: -\n", .{name});
}

fn writeRatioField(
    output: *std.Io.Writer,
    name: []const u8,
    numerator: u128,
    denominator: u128,
) !void {
    var buf: [96]u8 = undefined;
    if (denominator == 0) {
        const line = try std.fmt.bufPrint(&buf, "{s}: -\n", .{name});
        return output.writeAll(line);
    }
    const scale = 1_000_000;
    const rounded = (numerator * scale + denominator / 2) / denominator;
    const line = try std.fmt.bufPrint(
        &buf,
        "{s}: {d}.{d:0>6}\n",
        .{ name, rounded / scale, rounded % scale },
    );
    try output.writeAll(line);
}

fn printPathError(io: std.Io, path: []const u8, message: []const u8) void {
    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
}

fn printCommandFailure(io: std.Io, path: []const u8, failure: CommandFailure) void {
    if (failure.suppress_diagnostic) return;
    if (failure.record_index == null or
        failure.byte_offset == null or
        failure.line_in_record == null)
    {
        return printPathError(io, path, failure.message);
    }

    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, failure.code) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, failure.message) catch {};

    var buf: [96]u8 = undefined;
    const suffix = std.fmt.bufPrint(
        &buf,
        " (record {d}, line {d}, offset {d})\n",
        .{
            failure.record_index.?,
            failure.line_in_record.?,
            failure.byte_offset.?,
        },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, suffix) catch {};
}

fn printPairCommandFailure(
    io: std.Io,
    inputs: []const []const u8,
    pair_mode: PairMode,
    failure: PairCommandFailure,
) void {
    switch (failure) {
        .command => |details| printCommandFailure(
            io,
            inputs[details.input_index],
            details.details,
        ),
        .pair => |details| printPairFailure(io, inputs, pair_mode, details),
    }
}

fn printPairFailure(
    io: std.Io,
    inputs: []const []const u8,
    pair_mode: PairMode,
    failure: PairFailure,
) void {
    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, inputs[0]);
    if (pair_mode == .paired) {
        std.Io.File.writeStreamingAll(.stderr(), io, " + ") catch {};
        writeEscaped(.stderr(), io, inputs[1]);
    }
    switch (failure) {
        .name_mismatch => |details| {
            std.Io.File.writeStreamingAll(
                .stderr(),
                io,
                ": P001: paired identifiers or mate markers do not match",
            ) catch {};
            var buf: [64]u8 = undefined;
            const pair = std.fmt.bufPrint(&buf, " (pair {d})\n", .{details.pair_index}) catch
                return;
            std.Io.File.writeStreamingAll(.stderr(), io, pair) catch {};
            for (details.records, 0..) |record, index| {
                writePairRecordDiagnostic(
                    io,
                    inputs[if (pair_mode == .paired) index else 0],
                    index,
                    &record,
                );
            }
        },
        .count_mismatch => |details| {
            std.Io.File.writeStreamingAll(
                .stderr(),
                io,
                ": P002: paired input is missing a mate",
            ) catch {};
            var buf: [96]u8 = undefined;
            const prefix = std.fmt.bufPrint(
                &buf,
                " (pair {d}, remaining R{d}, last R1 record ",
                .{ details.pair_index, @as(u8, details.remaining_side) + 1 },
            ) catch return;
            std.Io.File.writeStreamingAll(.stderr(), io, prefix) catch {};
            writeOptionalIndex(io, details.record_indexes[0]);
            std.Io.File.writeStreamingAll(.stderr(), io, ", last R2 record ") catch {};
            writeOptionalIndex(io, details.record_indexes[1]);
            std.Io.File.writeStreamingAll(.stderr(), io, ")\n") catch {};
        },
    }
}

fn writePairRecordDiagnostic(
    io: std.Io,
    input: []const u8,
    side_index: usize,
    record: *const PairRecordDiagnostic,
) void {
    var buf: [96]u8 = undefined;
    const prefix = std.fmt.bufPrint(
        &buf,
        "  R{d}: input=",
        .{side_index + 1},
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, prefix) catch {};
    writeEscaped(.stderr(), io, input);
    const location = std.fmt.bufPrint(
        &buf,
        ", record={d}, offset={d}, first_token=",
        .{ record.record_index, record.byte_offset },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, location) catch {};
    writeBoundedHuman(io, &record.first_token);
    std.Io.File.writeStreamingAll(.stderr(), io, ", normalized_id=") catch {};
    writeBoundedHuman(io, &record.normalized_id);
    std.Io.File.writeStreamingAll(.stderr(), io, ", mate_markers=") catch {};
    writeMateMarkersHuman(io, record.mate_markers);
    std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
}

fn writeBoundedHuman(io: std.Io, display: *const BoundedBytes) void {
    writeEscaped(.stderr(), io, display.bytes());
    var buf: [64]u8 = undefined;
    const suffix = std.fmt.bufPrint(
        &buf,
        " [length={d}, truncated={s}]",
        .{ display.full_len, if (display.truncated()) "true" else "false" },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, suffix) catch {};
}

fn writeMateMarkersHuman(io: std.Io, markers: u2) void {
    if (markers == 0) {
        std.Io.File.writeStreamingAll(.stderr(), io, "none") catch {};
        return;
    }
    if (markers & 0b01 != 0) {
        std.Io.File.writeStreamingAll(.stderr(), io, "1") catch {};
    }
    if (markers == 0b11) {
        std.Io.File.writeStreamingAll(.stderr(), io, ",") catch {};
    }
    if (markers & 0b10 != 0) {
        std.Io.File.writeStreamingAll(.stderr(), io, "2") catch {};
    }
}

fn writeOptionalIndex(io: std.Io, record_index: ?u64) void {
    if (record_index) |value| {
        var buf: [32]u8 = undefined;
        const text = std.fmt.bufPrint(&buf, "{d}", .{value}) catch return;
        std.Io.File.writeStreamingAll(.stderr(), io, text) catch {};
    } else {
        std.Io.File.writeStreamingAll(.stderr(), io, "none") catch {};
    }
}

const HEX = "0123456789ABCDEF";

fn writeEscaped(file: std.Io.File, io: std.Io, bytes: []const u8) void {
    writeEscapedAll(file, io, bytes) catch {};
}

fn writeEscapedAll(file: std.Io.File, io: std.Io, bytes: []const u8) !void {
    var output = file.writerStreaming(io, &.{});
    try writeEscapedBytes(&output.interface, bytes, false);
    try output.interface.flush();
}

fn writeEscapedBytes(output: *std.Io.Writer, bytes: []const u8, comptime json_string: bool) !void {
    var run_start: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte >= 0x20 and byte <= 0x7e and byte != '\\') continue;

        try writeDisplayRun(output, bytes[run_start..index], json_string);
        if (byte == '\\') {
            try writeDisplayRun(output, "\\\\", json_string);
        } else {
            const escaped = [4]u8{ '\\', 'x', HEX[byte >> 4], HEX[byte & 0x0f] };
            try writeDisplayRun(output, &escaped, json_string);
        }
        run_start = index + 1;
    }
    try writeDisplayRun(output, bytes[run_start..], json_string);
}

fn writeDisplayRun(output: *std.Io.Writer, bytes: []const u8, comptime json_string: bool) !void {
    if (json_string) {
        try std.json.Stringify.encodeJsonStringChars(bytes, .{}, output);
    } else {
        try output.writeAll(bytes);
    }
}

test "[unit] - [interleaved staging]: releases oversized slack before a smaller mate" {
    var staging: std.ArrayList(u8) = .empty;
    defer staging.deinit(std.testing.allocator);
    const old_capacity = zfastq.limits.DEFAULT_READER_BUFFER_BYTES * 2;
    try staging.ensureTotalCapacityPrecise(std.testing.allocator, old_capacity);

    try stageCanonicalRecord(
        std.testing.allocator,
        &staging,
        .{
            .header = "pair/1",
            .id = "pair/1",
            .sequence = "A",
            .plus = "",
            .quality = "!",
        },
        null,
        old_capacity,
    );

    try std.testing.expectEqualStrings("@pair/1\nA\n+\n!\n", staging.items);
    try std.testing.expect(staging.capacity < old_capacity);
}

test "[unit] - [machine output]: handled non-lint errors share one exact shape" {
    const cases = [_]struct {
        code: []const u8,
        message: []const u8,
        exit_code: u8,
    }{
        .{ .code = "io_error", .message = "I/O error", .exit_code = 3 },
        .{ .code = "line_limit", .message = "line length limit exceeded", .exit_code = 4 },
        .{ .code = "arithmetic_limit", .message = "arithmetic limit exceeded", .exit_code = 4 },
        .{ .code = "out_of_memory", .message = "out of memory", .exit_code = 3 },
    };

    for (cases) |case| {
        var storage: [256]u8 = undefined;
        var writer = std.Io.Writer.fixed(&storage);
        var json: std.json.Stringify = .{ .writer = &writer };
        try json.beginObject();
        try writeJsonFailure(
            &json,
            CommandFailure.plain(case.code, case.message, case.exit_code),
        );
        try json.endObject();

        var parsed = try std.json.parseFromSlice(
            std.json.Value,
            std.testing.allocator,
            writer.buffered(),
            .{},
        );
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqual(@as(usize, 1), object.count());
        const failure = object.get("error").?.object;
        try std.testing.expectEqual(@as(usize, 5), failure.count());
        try std.testing.expectEqualStrings(case.code, failure.get("code").?.string);
        try std.testing.expectEqualStrings(case.message, failure.get("message").?.string);
        try std.testing.expect(failure.get("record_index").? == .null);
        try std.testing.expect(failure.get("byte_offset").? == .null);
        try std.testing.expect(failure.get("line_in_record").? == .null);
    }
}

test "[unit] - [structural arithmetic]: Reader and count preserve the CLI limit class" {
    var source = zfastq.io.plain.SliceSource.init("");
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();

    const reader_failure = mapReaderFailure(&reader, error.ArithmeticLimit);
    try std.testing.expectEqualStrings("arithmetic_limit", reader_failure.code);
    try std.testing.expectEqualStrings(
        "input location exceeds supported limit",
        reader_failure.message,
    );
    try std.testing.expectEqual(@as(u8, 4), reader_failure.exit_code);

    var scanner = zfastq.count_scan.Scanner.init(.{});
    try std.testing.expectEqualDeep(
        reader_failure,
        mapScanFailure(&scanner, error.ArithmeticLimit),
    );
    var check_scanner = fastq.CheckScanner.init(.{}, .{});
    try std.testing.expectEqualDeep(
        reader_failure,
        mapCheckScannerFailure(&check_scanner, error.ArithmeticLimit),
    );
}

test "[edge] - [semantic diagnostics]: field locations cross 32 bits and reject overflow" {
    const cases = [_]struct {
        semantic_error: zfastq.SemanticError,
        code: []const u8,
        line: u3,
    }{
        .{
            .semantic_error = .{
                .code = .s002_invalid_sequence_alphabet,
                .message = "sequence byte is outside the selected alphabet",
                .field = .sequence,
                .byte_index = 1,
            },
            .code = "S002",
            .line = 2,
        },
        .{
            .semantic_error = .{
                .code = .s006_invalid_quality_range,
                .message = "quality byte must be ASCII 33 through 126",
                .field = .quality,
                .byte_index = 1,
            },
            .code = "S006",
            .line = 4,
        },
    };
    const maximum = std.math.maxInt(u64);
    for (cases) |case| {
        for ([_]u64{ (1 << 32) - 1, maximum - 1, maximum }) |field_offset| {
            var offsets: zfastq.RecordOffsets = .{
                .header = 0,
                .sequence = 7,
                .plus = 9,
                .quality = 11,
            };
            switch (case.semantic_error.field) {
                .sequence => offsets.sequence = field_offset,
                .quality => offsets.quality = field_offset,
            }
            const actual = mapSemanticFailure(case.semantic_error, offsets, maximum).?;
            const expected: CommandFailure = if (field_offset == maximum) .{
                .code = "arithmetic_limit",
                .message = "input location exceeds supported limit",
                .exit_code = 4,
            } else .{
                .code = case.code,
                .message = case.semantic_error.message,
                .exit_code = 1,
                .record_index = maximum,
                .byte_offset = field_offset + 1,
                .line_in_record = case.line,
            };
            try std.testing.expectEqualDeep(expected, actual);
        }
    }
}

test "[edge] - [stats-json]: preserves the maximum u64 counter" {
    var storage: [1024]u8 = undefined;
    var writer = std.Io.Writer.fixed(&storage);
    var json: std.json.Stringify = .{ .writer = &writer };
    try writeStatsJsonResult(&json, "-", .{ .success = .{
        .reads = std.math.maxInt(u64),
    } });

    var parsed = try std.json.parseFromSlice(
        std.json.Value,
        std.testing.allocator,
        writer.buffered(),
        .{ .parse_numbers = false },
    );
    defer parsed.deinit();
    const result = parsed.value.object.get("reads").?;
    try std.testing.expectEqualStrings("18446744073709551615", result.number_string);
}

test "[failure] - [stats command]: reader allocation failure becomes a handled result" {
    var source = zfastq.io.plain.SliceSource.init("");
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 0,
    });

    const outcome = collectStats(failing.allocator(), source.byteSource(), .{});
    const failure = switch (outcome) {
        .success => return error.ExpectedFailure,
        .failure => |failure| failure,
    };
    try std.testing.expect(failing.has_induced_failure);
    try std.testing.expectEqualStrings("out_of_memory", failure.code);
    try std.testing.expectEqualStrings("out of memory", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expect(failure.record_index == null);
    try std.testing.expect(failure.byte_offset == null);
    try std.testing.expect(failure.line_in_record == null);
}

test "[integration] - [stats command]: handled failures preserve input order and exit precedence" {
    const Capture = struct {
        dir: std.Io.Dir,
        stdout: std.Io.Writer,
        stderr: std.Io.Writer,
        opened: usize = 0,
        closed: usize = 0,

        fn openFile(
            ctx: ?*anyopaque,
            _: std.Io.Dir,
            path: []const u8,
            options: std.Io.Dir.OpenFileOptions,
        ) std.Io.File.OpenError!std.Io.File {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            const file = try self.dir.openFile(std.testing.io, path, options);
            self.opened += 1;
            return file;
        }

        fn closeFiles(ctx: ?*anyopaque, files: []const std.Io.File) void {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            for (files) |file| file.close(std.testing.io);
            self.closed += files.len;
        }

        fn operate(ctx: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            if (operation == .file_write_streaming) {
                const write = operation.file_write_streaming;
                const writer = if (write.file.handle == std.Io.File.stdout().handle)
                    &self.stdout
                else if (write.file.handle == std.Io.File.stderr().handle)
                    &self.stderr
                else
                    return std.testing.io.operate(operation);
                return .{ .file_write_streaming = writer.writeSplatHeader(
                    write.header,
                    write.data,
                    write.splat,
                ) catch error.NoSpaceLeft };
            }
            return std.testing.io.operate(operation);
        }
    };

    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const window = io_layer.DEFAULT_READER_BUFFER_BYTES;
    const long_line = try std.testing.allocator.alloc(u8, window + 2);
    defer std.testing.allocator.free(long_line);
    try tmp.dir.writeFile(io, .{ .sub_path = "small.fastq", .data = "@r\nA\n+\n!\n" });
    {
        const file = try tmp.dir.createFile(io, "large.fastq", .{});
        defer file.close(io);
        try file.writeStreamingAll(io, "@r\nA\n+\n!\n@long\n");
        @memset(long_line, 'A');
        try file.writeStreamingAll(io, long_line[0 .. window + 1]);
        try file.writeStreamingAll(io, "\n+\n");
        @memset(long_line, '!');
        try file.writeStreamingAll(io, long_line[0 .. window + 1]);
        try file.writeStreamingAll(io, "\n");
    }
    @memset(long_line, 'x');
    long_line[0] = '@';
    try tmp.dir.writeFile(io, .{ .sub_path = "limit.fastq", .data = long_line });

    // Both Reader buffers fit; growing a retained field fails without exhausting the host.
    const storage = try std.testing.allocator.alloc(u8, 2 * window);
    defer std.testing.allocator.free(storage);
    const block = "input: small.fastq\n" ++
        "reads: 1\nbases: 1\nmin_length: 1\nmax_length: 1\nmean_length: 1.000000\n" ++
        "a: 1\nc: 0\ng: 0\nt: 0\nn: 0\nother_bases: 0\ngc_fraction: 0.000000\n" ++
        "quality_sum: 0\nmean_quality: 0.000000\nq20_bases: 0\nq20_fraction: 0.000000\n" ++
        "q30_bases: 0\nq30_fraction: 0.000000\n";
    const partial = block ++ "\ninput: small.fastq\n";
    const cases = [_]struct {
        inputs: []const []const u8,
        stdout: []const u8,
        stderr: []const u8,
        status: u8,
        opened: usize,
        stdout_capacity: usize = 1024,
    }{
        .{
            .inputs = &.{ "large.fastq", "small.fastq" },
            .stdout = block,
            .stderr = "error: large.fastq: out of memory\n",
            .status = 3,
            .opened = 2,
        },
        .{
            .inputs = &.{ "limit.fastq", "large.fastq", "small.fastq" },
            .stdout = block,
            .stderr = "error: limit.fastq: line length limit exceeded\n" ++
                "error: large.fastq: out of memory\n",
            .status = 4,
            .opened = 3,
        },
        .{
            .inputs = &.{ "small.fastq", "large.fastq", "small.fastq" },
            .stdout = block ++ "\n" ++ block,
            .stderr = "error: large.fastq: out of memory\n",
            .status = 3,
            .opened = 3,
        },
        .{
            .inputs = &.{ "small.fastq", "limit.fastq", "small.fastq", "small.fastq" },
            .stdout = partial,
            .stderr = "error: limit.fastq: line length limit exceeded\n",
            .status = 4,
            .opened = 3,
            .stdout_capacity = partial.len,
        },
    };
    var vtable = std.Io.failing.vtable.*;
    vtable.dirOpenFile = Capture.openFile;
    vtable.fileClose = Capture.closeFiles;
    vtable.operate = Capture.operate;
    for (cases) |case| {
        var stdout_buffer: [1024]u8 = undefined;
        var stderr_buffer: [256]u8 = undefined;
        var capture: Capture = .{
            .dir = tmp.dir,
            .stdout = .fixed(stdout_buffer[0..case.stdout_capacity]),
            .stderr = .fixed(&stderr_buffer),
        };
        var bounded = std.heap.FixedBufferAllocator.init(storage);

        const status = runStats(
            .{ .userdata = &capture, .vtable = &vtable },
            bounded.allocator(),
            case.inputs,
            .{ .max_line_bytes = window + 1 },
            false,
        );

        try std.testing.expectEqual(case.status, status);
        try std.testing.expectEqualStrings(case.stdout, capture.stdout.buffered());
        try std.testing.expectEqualStrings(case.stderr, capture.stderr.buffered());
        try std.testing.expectEqual(case.opened, capture.opened);
        try std.testing.expectEqual(capture.opened, capture.closed);
        try std.testing.expectEqual(@as(usize, 0), bounded.end_index);
    }
}

fn snapshotTestFile(io: std.Io, path: []const u8) !FileSnapshot {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    return fileSnapshot(file, io);
}

test "[unit] - [exact sample]: every retained file-change signal is compared" {
    const baseline = FileSnapshot{
        .inode = 17,
        .size = 23,
        .mtime_nanoseconds = 29,
    };
    try std.testing.expect(sameFileSnapshot(baseline, baseline));

    var changed = baseline;
    changed.inode += 1;
    try std.testing.expect(!sameFileSnapshot(baseline, changed));
    changed = baseline;
    changed.size += 1;
    try std.testing.expect(!sameFileSnapshot(baseline, changed));
    changed = baseline;
    changed.mtime_nanoseconds += 1;
    try std.testing.expect(!sameFileSnapshot(baseline, changed));
}

test "[failure] - [exact sample]: selection allocation failure precedes later block format error" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{
        .sub_path = "records.fastq",
        .data = "@one\nA\n+\n!\n" ++
            "@two\nC\n+\n#\n" ++
            "@bad\nG\nx\n$\n",
    });
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(
        &path_buffer,
        ".zig-cache/tmp/{s}/records.fastq",
        .{tmp.sub_path},
    );
    var selector = sampling.ExactSelector.init(1, 11);
    defer selector.deinit(std.testing.allocator);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 0,
    });

    const result = sampleExactFirstPass(io, failing.allocator(), path, &selector, .{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 1,
        .seed = 11,
    });
    const failure = switch (result) {
        .success => return error.ExpectedFailure,
        .failure => |failure| failure,
    };

    try std.testing.expect(failing.has_induced_failure);
    try std.testing.expectEqualStrings("out_of_memory", failure.code);
    try std.testing.expectEqual(@as(u64, 2), selector.record_count);
}

test "[failure] - [exact sample]: final record-count change keeps a valid prefix" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const name = "record.fastq";
    {
        const file = try tmp.dir.createFile(io, name, .{});
        defer file.close(io);
        try std.Io.File.writeStreamingAll(file, io, "@one\nA\n+\n!\n");
    }
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(
        &path_buffer,
        ".zig-cache/tmp/{s}/{s}",
        .{ tmp.sub_path, name },
    );
    const snapshot = try snapshotTestFile(io, path);

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactSecondPass(
        false,
        io,
        std.testing.allocator,
        path,
        &writer,
        .{ .low_words = &.{1}, .middle_bytes = &.{0} },
        snapshot,
        2,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
        },
    )).?;

    try std.testing.expectEqualStrings("input_changed", failure.code);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expectEqualStrings("@one\nA\n+\n!\n", sink.written());
}

test "[failure] - [exact sample]: changed metadata stops the second pass before output" {
    const io = std.testing.io;
    const path = "tests/data/synthetic/basic_valid.fastq";
    var snapshot = try snapshotTestFile(io, path);
    snapshot.size += 1;

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactSecondPass(
        false,
        io,
        std.testing.allocator,
        path,
        &writer,
        .{ .low_words = &.{1}, .middle_bytes = &.{0} },
        snapshot,
        5,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
        },
    )).?;

    try std.testing.expectEqualStrings("input_changed", failure.code);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}

test "[failure] - [exact sample]: a FIFO replacement reports input changed on reopen" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const name = "record.fastq";
    try tmp.dir.writeFile(io, .{ .sub_path = name, .data = "@one\nA\n+\n!\n" });
    const snapshot = snapshot: {
        const file = try tmp.dir.openFile(io, name, .{});
        defer file.close(io);
        break :snapshot try fileSnapshot(file, io);
    };
    try tmp.dir.deleteFile(io, name);
    try std.testing.expectEqual(.SUCCESS, std.os.linux.errno(std.os.linux.mknodat(
        tmp.dir.handle,
        name,
        std.os.linux.S.IFIFO | 0o600,
        0,
    )));
    // A connected FIFO keeps a blocking-open regression from hanging this in-process test.
    const keeper: std.Io.File = .{
        .handle = try std.posix.openat(tmp.dir.handle, name, .{
            .ACCMODE = .RDWR,
            .NONBLOCK = true,
            .CLOEXEC = true,
        }, 0),
        .flags = .{ .nonblocking = true },
    };
    defer keeper.close(io);
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(
        &path_buffer,
        ".zig-cache/tmp/{s}/{s}",
        .{ tmp.sub_path, name },
    );
    var input: RecordInput = undefined;
    const failure = switch (initExactInput(&input, io, path, snapshot, null)) {
        .failure => |failure| failure,
        .success => {
            input.deinit(io);
            return error.ExpectedFailure;
        },
    };
    try std.testing.expectEqualStrings("input_changed", failure.code);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
}

fn writeExactTestInput(
    io: std.Io,
    path: []const u8,
    bytes: []const u8,
    gzip: bool,
    snapshot: ?FileSnapshot,
) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var buffer: [4096]u8 = undefined;
    var output = file.writer(io, &buffer);
    if (gzip) {
        const len = std.math.cast(u16, bytes.len) orelse return error.TestFixtureTooLarge;
        try output.interface.writeAll(&.{ 0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255, 1 });
        try output.interface.writeInt(u16, len, .little);
        try output.interface.writeInt(u16, ~len, .little);
        try output.interface.writeAll(bytes);
        try output.interface.writeInt(u32, std.hash.Crc32.hash(bytes), .little);
        try output.interface.writeInt(u32, len, .little);
    } else {
        try output.interface.writeAll(bytes);
    }
    try output.interface.flush();
    if (snapshot) |expected| {
        try file.setTimestamps(io, .{
            .modify_timestamp = .{ .new = .{ .nanoseconds = expected.mtime_nanoseconds } },
        });
        try std.testing.expect(sameFileSnapshot(expected, try fileSnapshot(file, io)));
    }
}

test "[integration] - [exact sample]: selected records are revalidated after metadata-preserving changes" {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/input", .{tmp.sub_path});
    defer allocator.free(path);
    const prefix = "@a\nC\n+\n#\n";
    const cases = [_]struct {
        before: []const u8 = "@b\nA\n+\n!\n",
        after: []const u8,
        alphabet: zfastq.Alphabet = .iupac,
        valid_output: ?[]const u8 = null,
    }{
        .{ .after = "@b\n.\n+\n!\n" },
        .{ .after = "@b\nA\n+\n \n" },
        .{ .after = "@b\nA\n+\n\x7f\n" },
        .{ .after = "@b\nR\n+\n!\n", .alphabet = .acgtn },
        .{ .after = "@b\nR\n+\n!\n", .valid_output = "@b\nR\n+\n!\n" },
        .{ .before = "@b\r\nA\r\n+\r\n!\r\n", .after = "@b\r\n.\r\n+\r\n!\r\n" },
        .{
            .before = "@b\r\nA\r\n+\r\n!\r\n",
            .after = "@b\r\nR\r\n+\r\n~\r\n",
            .valid_output = "@b\nR\n+\n~\n",
        },
    };
    for ([_]bool{ false, true }) |gzip| {
        inline for (.{ .all, .selected, .unselected }) |selection| {
            for (cases, 0..) |case, case_index| {
                errdefer std.debug.print("gzip={} selection={s} case={d}\n", .{ gzip, @tagName(selection), case_index });
                var input_buffer: [128]u8 = undefined;
                const before = try std.fmt.bufPrint(&input_buffer, "{s}{s}", .{ prefix, case.before });
                try writeExactTestInput(io, path, before, gzip, null);
                const options = SampleOptions{
                    .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
                    .alphabet = case.alphabet,
                    .fraction = null,
                    .count = if (selection == .all) 2 else 1,
                    .seed = 11,
                };
                var selector = sampling.ExactSelector.init(options.count.?, options.seed);
                defer selector.deinit(allocator);
                const first = sampleExactFirstPass(io, allocator, path, &selector, options);
                try std.testing.expect(first == .success);
                try std.testing.expectEqual(@as(u64, 2), selector.record_count);
                const snapshot = first.success.snapshot;
                const after = try std.fmt.bufPrint(&input_buffer, "{s}{s}", .{ prefix, case.after });
                try writeExactTestInput(io, path, after, gzip, snapshot);

                var output: [128]u8 = undefined;
                var sink = io_layer.SliceSink.init(&output);
                var writer = zfastq.Writer.init(sink.byteSink());
                const failure = try sampleExactSecondPass(
                    selection == .all,
                    io,
                    allocator,
                    path,
                    &writer,
                    if (selection == .all) .empty else .{
                        .low_words = &.{if (selection == .selected) 2 else 1},
                        .middle_bytes = &.{0},
                    },
                    snapshot,
                    selector.record_count,
                    options,
                );
                if (selection != .unselected and case.valid_output == null) {
                    try std.testing.expect(failure != null);
                    try std.testing.expectEqualStrings("input_changed", failure.?.code);
                    try std.testing.expectEqualStrings("input changed during exact sampling", failure.?.message);
                    try std.testing.expectEqual(@as(u8, 3), failure.?.exit_code);
                } else {
                    try std.testing.expect(failure == null);
                }
                var expected: [128]u8 = undefined;
                try std.testing.expectEqualStrings(try std.fmt.bufPrint(&expected, "{s}{s}", .{
                    if (selection == .selected) "" else prefix,
                    if (selection == .unselected) "" else case.valid_output orelse "",
                }), sink.written());
            }
        }
    }
}

test "[integration] - [paired exact sample]: selected pairs are revalidated after metadata-preserving changes" {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path1 = try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/r1", .{tmp.sub_path});
    defer allocator.free(path1);
    const path2 = try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/r2", .{tmp.sub_path});
    defer allocator.free(path2);
    const prefix1 = "@a\nC\n+\n#\n";
    const prefix2 = "@a\nG\n+\n$\n";
    const prefix = prefix1 ++ prefix2;
    const cases = [_]struct {
        before1: []const u8 = "@b/1\nA\n+\n!\n",
        before2: []const u8 = "@b/2\nT\n+\n!\n",
        after1: []const u8 = "@b/1\nA\n+\n!\n",
        after2: []const u8 = "@b/2\nT\n+\n!\n",
        alphabet: zfastq.Alphabet = .iupac,
        policy: pairing.NamePolicy = .illumina,
        failed_input: ?u1 = 0,
        code: []const u8 = "input_changed",
        exit_code: u8 = 3,
    }{
        .{ .after1 = "@b/1\n.\n+\n!\n" },
        .{ .after2 = "@b/2\n.\n+\n!\n", .failed_input = 1 },
        .{ .after1 = "@b/1\nA\n+\n \n" },
        .{ .after1 = "@b/1\nR\n+\n \n" },
        .{ .after2 = "@b/2\nT\n+\n\x7f\n", .failed_input = 1 },
        .{ .after2 = "@c/2\nT\n+\n!\n" },
        .{ .after2 = "@b/1\nT\n+\n!\n" },
        .{ .before2 = "@b/1\nT\n+\n!\n", .policy = .exact },
        .{
            .before1 = "@b/1 1:a\nA\n+\n!\n",
            .after1 = "@b/1 1:a\nA\n+\n!\n",
            .before2 = "@b/1 2:a\nT\n+\n!\n",
            .after2 = "@b/1 2:b\nT\n+\n!\n",
            .policy = .exact,
            .failed_input = null,
        },
        .{ .after2 = "@b/2\nR\n+\n!\n", .alphabet = .acgtn, .failed_input = 1 },
        .{ .after2 = "@b/2\nR\n+\n!\n", .failed_input = null },
        .{
            .before1 = "@b/1\r\nA\r\n+\r\n!\r\n",
            .after1 = "@b/1\r\n.\r\n+\r\n!\r\n",
        },
        .{
            .after1 = "@b/1\n.\n+\n!\n",
            .after2 = "@c/2\nT\nx\n!\n",
            .failed_input = 1,
            .code = "S001",
            .exit_code = 1,
        },
        .{
            .after1 = "@b/1\nR\n+\n \n",
            .after2 = "@c/2\nT\nx\n!\n",
            .failed_input = 1,
            .code = "S001",
            .exit_code = 1,
        },
    };
    for ([_]PairMode{ .paired, .interleaved }) |mode| {
        const inputs: []const []const u8 = if (mode == .paired) &.{ path1, path2 } else &.{path1};
        for ([_]bool{ false, true }) |gzip| {
            inline for (.{ .all, .selected, .unselected }) |selection| {
                for (cases, 0..) |case, case_index| {
                    errdefer std.debug.print("mode={s} gzip={} selection={s} case={d}\n", .{ @tagName(mode), gzip, @tagName(selection), case_index });
                    const options = SampleOptions{
                        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
                        .alphabet = case.alphabet,
                        .fraction = null,
                        .count = if (selection == .all) 2 else 1,
                        .seed = 11,
                        .pair_mode = mode,
                        .pair_name_policy = case.policy,
                    };
                    var buffer1: [256]u8 = undefined;
                    var buffer2: [256]u8 = undefined;
                    const before1 = if (mode == .paired)
                        try std.fmt.bufPrint(&buffer1, "{s}{s}", .{ prefix1, case.before1 })
                    else
                        try std.fmt.bufPrint(&buffer1, "{s}{s}{s}", .{ prefix, case.before1, case.before2 });
                    try writeExactTestInput(io, path1, before1, gzip, null);
                    if (mode == .paired) {
                        const before2 = try std.fmt.bufPrint(&buffer2, "{s}{s}", .{ prefix2, case.before2 });
                        try writeExactTestInput(io, path2, before2, gzip, null);
                    }
                    var selector = sampling.ExactSelector.init(options.count.?, options.seed);
                    defer selector.deinit(allocator);
                    const first = sampleExactPairFirstPass(io, allocator, inputs, &selector, options);
                    try std.testing.expect(first == .success);
                    try std.testing.expectEqual(@as(u64, 2), selector.record_count);
                    const snapshots = first.success.snapshots;
                    const after1 = if (mode == .paired)
                        try std.fmt.bufPrint(&buffer1, "{s}{s}", .{ prefix1, case.after1 })
                    else
                        try std.fmt.bufPrint(&buffer1, "{s}{s}{s}", .{ prefix, case.after1, case.after2 });
                    try writeExactTestInput(io, path1, after1, gzip, if (mode == .paired)
                        snapshots.paired[0]
                    else
                        snapshots.interleaved);
                    if (mode == .paired) {
                        const after2 = try std.fmt.bufPrint(&buffer2, "{s}{s}", .{ prefix2, case.after2 });
                        try writeExactTestInput(io, path2, after2, gzip, snapshots.paired[1]);
                    }

                    var output: [256]u8 = undefined;
                    var sink = io_layer.SliceSink.init(&output);
                    var writer = zfastq.Writer.init(sink.byteSink());
                    const failure = try sampleExactPairSecondPass(
                        selection == .all,
                        io,
                        allocator,
                        inputs,
                        &writer,
                        if (selection == .all) .empty else .{
                            .low_words = &.{if (selection == .selected) 2 else 1},
                            .middle_bytes = &.{0},
                        },
                        snapshots,
                        selector.record_count,
                        try deinterleaveStagingLimit(options.max_line_bytes),
                        options,
                    );
                    const reject = case.failed_input != null and
                        (selection != .unselected or std.mem.eql(u8, case.code, "S001"));
                    if (reject) {
                        try std.testing.expect(failure != null);
                        try std.testing.expect(failure.? == .command);
                        const command = failure.?.command;
                        try std.testing.expectEqual(if (mode == .paired) case.failed_input.? else 0, command.input_index);
                        try std.testing.expectEqualStrings(case.code, command.details.code);
                        try std.testing.expectEqual(case.exit_code, command.details.exit_code);
                    } else {
                        try std.testing.expect(failure == null);
                    }
                    var expected: [256]u8 = undefined;
                    const emit_changed = !reject and selection != .unselected;
                    try std.testing.expectEqualStrings(try std.fmt.bufPrint(&expected, "{s}{s}{s}", .{
                        if (selection == .selected) "" else prefix,
                        if (emit_changed) case.after1 else "",
                        if (emit_changed) case.after2 else "",
                    }), sink.written());
                }
            }
        }
    }
}

test "[integration] - [interleaved exact sample]: revalidation survives mate storage changes and output failures" {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/pairs", .{tmp.sub_path});
    defer allocator.free(path);
    const small_prefix = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n";
    const options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 2,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    const staging_limit = try deinterleaveStagingLimit(options.max_line_bytes);
    const Change = enum { valid, alphabet1, quality1, alphabet2, quality2, name, structure, write };
    for ([_]InterleavedFirstRecordStorage{ .reader, .staged, .retained }) |storage| {
        var input: std.ArrayList(u8) = .empty;
        defer input.deinit(allocator);
        if (storage == .reader) {
            // Put the selected pair across the transport refill while both mates stay small.
            try input.appendSlice(allocator, "@a/1 ");
            try input.appendNTimes(allocator, 'x', zfastq.limits.DEFAULT_READER_BUFFER_BYTES - 300 - small_prefix.len - 1);
            try input.appendSlice(allocator, small_prefix[4..]);
        } else try input.appendSlice(allocator, small_prefix);
        const prefix_len = input.items.len;
        var header_starts: [2]usize = undefined;
        var sequence_starts: [2]usize = undefined;
        var plus_starts: [2]usize = undefined;
        var quality_starts: [2]usize = undefined;
        for (0..2) |mate| {
            const field_len: usize = if (storage == .retained or (storage == .staged and mate == 1))
                zfastq.limits.DEFAULT_READER_BUFFER_BYTES
            else
                8;
            try input.append(allocator, '@');
            header_starts[mate] = input.items.len;
            try input.appendNTimes(allocator, 'b', 192);
            try input.appendSlice(allocator, if (mate == 0) "/1\n" else "/2\n");
            sequence_starts[mate] = input.items.len;
            try input.appendNTimes(allocator, if (mate == 0) 'A' else 'T', field_len);
            try input.append(allocator, '\n');
            plus_starts[mate] = input.items.len;
            try input.appendSlice(allocator, "+\n");
            quality_starts[mate] = input.items.len;
            try input.appendNTimes(allocator, '!', field_len);
            try input.append(allocator, '\n');
        }
        const prefix = input.items[0..prefix_len];
        {
            var source = io_layer.SliceSource.init(input.items);
            var reader = try zfastq.Reader.init(allocator, source.byteSource(), .{});
            defer reader.deinit();
            try std.testing.expect(try reader.advance());
            try std.testing.expect(try reader.advance());
            var retained: fastq.RetainedRecordStorage = .{};
            defer retained.deinit(allocator);
            var staged: std.ArrayList(u8) = .empty;
            defer staged.deinit(allocator);
            var validator = fastq.AdaptiveRecordValidator.init(.{});
            var first_record = (try fastq.nextValidatedRecord(&reader, &validator)).?;
            if (storage == .reader) {
                try std.testing.expect((try fastq.nextBufferedValidatedRecord(&reader, &validator)) == null);
            }
            if (try fastq.nextPairedValidatedRecord(&reader, &first_record, &validator)) |_| {
                try std.testing.expectEqual(.reader, storage);
            } else {
                const next = try nextAfterPreservingInterleavedMate1(
                    allocator,
                    &reader,
                    &retained,
                    &staged,
                    first_record.record,
                    first_record.canonical_span,
                    staging_limit,
                    &validator,
                );
                try std.testing.expectEqual(storage, next.first_storage);
                try std.testing.expect(next.record != null);
            }
        }
        try writeExactTestInput(io, path, input.items, false, null);
        var selector = sampling.ExactSelector.init(2, options.seed);
        defer selector.deinit(allocator);
        const first = sampleExactPairFirstPass(io, allocator, &.{path}, &selector, options);
        try std.testing.expect(first == .success);
        try std.testing.expectEqual(@as(u64, 2), selector.record_count);
        const snapshot = first.success.snapshots.interleaved;
        const changed = try allocator.dupe(u8, input.items);
        defer allocator.free(changed);
        const output = try allocator.alloc(u8, input.items.len);
        defer allocator.free(output);
        for ([_]bool{ false, true }) |select_all| {
            for (std.enums.values(Change)) |change| {
                errdefer std.debug.print("storage={s} select_all={} change={s}\n", .{ @tagName(storage), select_all, @tagName(change) });
                @memcpy(changed, input.items);
                switch (change) {
                    .valid, .write => {},
                    .alphabet1 => changed[sequence_starts[0]] = '.',
                    .quality1 => changed[quality_starts[0]] = ' ',
                    .alphabet2 => changed[sequence_starts[1]] = '.',
                    .quality2 => changed[quality_starts[1]] = 0x7f,
                    .name => changed[header_starts[1] + 180] = 'c',
                    .structure => {
                        changed[sequence_starts[0]] = '.';
                        changed[plus_starts[1]] = 'x';
                    },
                }
                try writeExactTestInput(io, path, changed, false, snapshot);
                const expected_prefix = if (select_all) prefix else "";
                var sink = io_layer.SliceSink.init(if (change == .write)
                    output[0..expected_prefix.len]
                else
                    output);
                var writer = zfastq.Writer.init(sink.byteSink());
                const result = sampleExactInterleavedSecondPass(
                    select_all,
                    io,
                    allocator,
                    path,
                    &writer,
                    if (select_all) .empty else .{ .low_words = &.{2}, .middle_bytes = &.{0} },
                    snapshot,
                    2,
                    staging_limit,
                    options,
                );
                if (change == .write) {
                    try std.testing.expectError(error.WriteFailed, result);
                } else if (change == .valid) {
                    try std.testing.expect(try result == null);
                    try std.testing.expectEqualStrings(if (select_all) input.items else input.items[prefix.len..], sink.written());
                    continue;
                } else {
                    const failure = try result;
                    try std.testing.expect(failure != null);
                    try std.testing.expect(failure.? == .command);
                    const details = failure.?.command.details;
                    try std.testing.expectEqualStrings(if (change == .structure) "S001" else "input_changed", details.code);
                    try std.testing.expectEqual(@as(u8, if (change == .structure) 1 else 3), details.exit_code);
                }
                try std.testing.expectEqualStrings(expected_prefix, sink.written());
            }
        }
    }
}

test "[integration] - [interleaved exact sample]: preserves buffered mate failure diagnostics" {
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/pairs", .{tmp.sub_path});
    defer allocator.free(path);
    const original = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n@b/1\nA\n+\n!\n@b/2\nT\n+\n#\n";
    const changed = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n@b/1\nA\n+\n!\n@b/2\nT\nx\n#\n";
    try writeExactTestInput(io, path, original, false, null);
    const snapshot = try snapshotTestFile(io, path);
    try writeExactTestInput(io, path, changed, false, snapshot);
    var output: [original.len]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const options = SampleOptions{
        .max_line_bytes = 64,
        .alphabet = .iupac,
        .fraction = null,
        .count = 2,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    const failure = (try sampleExactInterleavedSecondPass(
        true,
        io,
        allocator,
        path,
        &writer,
        .empty,
        snapshot,
        2,
        260,
        options,
    )).?;
    try std.testing.expect(failure == .command);
    const details = failure.command.details;
    try std.testing.expectEqualStrings("S001", details.code);
    try std.testing.expectEqual(@as(u8, 1), details.exit_code);
    try std.testing.expectEqual(@as(?u64, 3), details.record_index);
    try std.testing.expectEqual(@as(?u3, 3), details.line_in_record);
    // The existing buffered-error path retries at quality and replaces offset 40.
    try std.testing.expectEqual(@as(?u64, 42), details.byte_offset);
    try std.testing.expectEqualStrings(original[0..22], sink.written());
}

test "[failure] - [paired exact sample]: each input snapshot is checked independently" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@a/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@a/2\nT\n+\n#\n" });
    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n",
    });

    var path1_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path1 = try std.fmt.bufPrint(
        &path1_buffer,
        ".zig-cache/tmp/{s}/r1.fastq",
        .{tmp.sub_path},
    );
    var path2_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path2 = try std.fmt.bufPrint(
        &path2_buffer,
        ".zig-cache/tmp/{s}/r2.fastq",
        .{tmp.sub_path},
    );
    var pairs_path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const pairs_path = try std.fmt.bufPrint(
        &pairs_path_buffer,
        ".zig-cache/tmp/{s}/pairs.fastq",
        .{tmp.sub_path},
    );
    const inputs = [_][]const u8{ path1, path2 };

    const snapshot1 = try snapshotTestFile(io, path1);
    const snapshot2 = try snapshotTestFile(io, path2);
    const pair_snapshot = try snapshotTestFile(io, pairs_path);
    const options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 1,
        .seed = 11,
        .pair_mode = .paired,
    };

    for (0..2) |changed_input| {
        var snapshots = [2]FileSnapshot{ snapshot1, snapshot2 };
        snapshots[changed_input].size += 1;
        var output: [64]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        const failure = (try sampleExactPairedSecondPass(
            true,
            io,
            std.testing.allocator,
            &inputs,
            &writer,
            .empty,
            snapshots,
            1,
            options,
        )).?;
        const command = switch (failure) {
            .command => |command| command,
            .pair => return error.ExpectedCommandFailure,
        };
        try std.testing.expectEqual(@as(u1, @intCast(changed_input)), command.input_index);
        try std.testing.expectEqualStrings("input_changed", command.details.code);
        try std.testing.expectEqual(@as(usize, 0), sink.written().len);
    }

    var changed_pair_snapshot = pair_snapshot;
    changed_pair_snapshot.size += 1;
    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactInterleavedSecondPass(
        true,
        io,
        std.testing.allocator,
        pairs_path,
        &writer,
        .empty,
        changed_pair_snapshot,
        1,
        try deinterleaveStagingLimit(zfastq.limits.DEFAULT_MAX_LINE_BYTES),
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
            .pair_mode = .interleaved,
        },
    )).?;
    const command = switch (failure) {
        .command => |command| command,
        .pair => return error.ExpectedCommandFailure,
    };
    try std.testing.expectEqualStrings("input_changed", command.details.code);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}

test "[integration] - [paired exact sample]: the output pass checks structure and count" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var path1_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path1 = try std.fmt.bufPrint(
        &path1_buffer,
        ".zig-cache/tmp/{s}/r1.fastq",
        .{tmp.sub_path},
    );
    var path2_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path2 = try std.fmt.bufPrint(
        &path2_buffer,
        ".zig-cache/tmp/{s}/r2.fastq",
        .{tmp.sub_path},
    );
    var pairs_path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const pairs_path = try std.fmt.bufPrint(
        &pairs_path_buffer,
        ".zig-cache/tmp/{s}/pairs.fastq",
        .{tmp.sub_path},
    );
    const inputs = [_][]const u8{ path1, path2 };
    const first_pair = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n";
    const paired_options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 1,
        .seed = 11,
        .pair_mode = .paired,
    };

    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@a/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@a/2\nT\n+\n#\n" });
    const snapshot1 = try snapshotTestFile(io, path1);
    const snapshot2 = try snapshotTestFile(io, path2);

    var paired_output: [64]u8 = undefined;
    var paired_sink = io_layer.SliceSink.init(&paired_output);
    var paired_writer = zfastq.Writer.init(paired_sink.byteSink());
    const count_failure = (try sampleExactPairedSecondPass(
        true,
        io,
        std.testing.allocator,
        &inputs,
        &paired_writer,
        .empty,
        .{ snapshot1, snapshot2 },
        2,
        paired_options,
    )).?;
    const count_command = switch (count_failure) {
        .command => |command| command,
        .pair => return error.ExpectedCommandFailure,
    };
    try std.testing.expectEqualStrings("input_changed", count_command.details.code);
    try std.testing.expectEqualStrings(first_pair, paired_sink.written());

    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = first_pair ++ "@odd/1\nC\n+\n$\n",
    });
    const pair_snapshot = try snapshotTestFile(io, pairs_path);
    var interleaved_output: [64]u8 = undefined;
    var interleaved_sink = io_layer.SliceSink.init(&interleaved_output);
    var interleaved_writer = zfastq.Writer.init(interleaved_sink.byteSink());
    const odd_failure = (try sampleExactInterleavedSecondPass(
        true,
        io,
        std.testing.allocator,
        pairs_path,
        &interleaved_writer,
        .empty,
        pair_snapshot,
        1,
        try deinterleaveStagingLimit(zfastq.limits.DEFAULT_MAX_LINE_BYTES),
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
            .pair_mode = .interleaved,
        },
    )).?;
    const odd_command = switch (odd_failure) {
        .command => |command| command,
        .pair => return error.ExpectedCommandFailure,
    };
    try std.testing.expectEqualStrings("input_changed", odd_command.details.code);
    try std.testing.expectEqualStrings(first_pair, interleaved_sink.written());

    const structural_r1 = "@bad/1\nA\n+\n!\n";
    const structural_r2 = "@bad/2\nT\nx\n#\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = structural_r1 });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = structural_r2 });
    const structural_snapshot1 = try snapshotTestFile(io, path1);
    const structural_snapshot2 = try snapshotTestFile(io, path2);
    var structural_output: [64]u8 = undefined;
    var structural_sink = io_layer.SliceSink.init(&structural_output);
    var structural_writer = zfastq.Writer.init(structural_sink.byteSink());
    const structural_failure = (try sampleExactPairedSecondPass(
        true,
        io,
        std.testing.allocator,
        &inputs,
        &structural_writer,
        .empty,
        .{ structural_snapshot1, structural_snapshot2 },
        1,
        paired_options,
    )).?;
    const structural_command = switch (structural_failure) {
        .command => |command| command,
        .pair => return error.ExpectedCommandFailure,
    };
    try std.testing.expectEqual(@as(u1, 1), structural_command.input_index);
    try std.testing.expectEqualStrings("S001", structural_command.details.code);
    try std.testing.expectEqual(@as(usize, 0), structural_sink.written().len);
}

test "[property] - [paired identifiers]: storage reuses its largest allocation and releases failures" {
    const bytes = try std.testing.allocator.alloc(u8, 256 * 1024 + 1);
    defer std.testing.allocator.free(bytes);
    @memset(bytes, 'r');
    bytes[bytes.len - 1] = 'z';
    const Exercise = struct {
        fn run(allocator: std.mem.Allocator, name: []const u8) !void {
            var tracking = std.testing.FailingAllocator.init(allocator, .{});
            {
                var storage: std.ArrayList(u8) = .empty;
                defer storage.deinit(tracking.allocator());
                try storeNormalizedId(tracking.allocator(), &storage, name[0..3]);
                try storeNormalizedId(tracking.allocator(), &storage, name);
                const allocated = tracking.allocated_bytes;
                for ([_]usize{ 7, name.len, 0, name.len }) |len| {
                    try storeNormalizedId(tracking.allocator(), &storage, name[0..len]);
                    try std.testing.expectEqualStrings(name[0..len], storage.items);
                    try std.testing.expectEqual(name.len, storage.capacity);
                    try std.testing.expectEqual(allocated, tracking.allocated_bytes);
                }
            }
            try std.testing.expectEqual(tracking.allocated_bytes, tracking.freed_bytes);
        }
    };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Exercise.run, .{bytes});
}

fn exerciseExactAllocations(allocator: std.mem.Allocator) !void {
    const r1 = "@a/1\nA\n+\n!\n@b/1\nC\n+\n#\n@c/1\nG\n+\n$\n";
    const r2 = "@a/2\nT\n+\n!\n@b/2\nG\n+\n#\n@c/2\nC\n+\n$\n";
    const interleaved =
        "@a/1\nA\n+\n!\n@a/2\nT\n+\n!\n" ++
        "@b/1\nC\n+\n#\n@b/2\nG\n+\n#\n" ++
        "@c/1\nG\n+\n$\n@c/2\nC\n+\n$\n";
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const names = [_][]const u8{ "r1.fastq", "r2.fastq", "interleaved.fastq" };
    var paths: [3][]const u8 = undefined;
    var path_count: usize = 0;
    defer for (paths[0..path_count]) |path| std.testing.allocator.free(path);
    for (names, [_][]const u8{ r1, r2, interleaved }, 0..) |name, bytes, index| {
        try tmp.dir.writeFile(io, .{ .sub_path = name, .data = bytes });
        paths[index] = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/{s}", .{ tmp.sub_path, name });
        path_count += 1;
    }

    const descriptors = try openDescriptorCount();
    for ([_]PairMode{ .none, .paired, .interleaved }) |mode| {
        for ([_]u64{ 0, 1, 2, 3, 4 }) |count| {
            var output: [128]u8 = undefined;
            var sink = io_layer.SliceSink.init(&output);
            var writer = zfastq.Writer.init(sink.byteSink());
            const options = SampleOptions{
                .max_line_bytes = 8192,
                .alphabet = .iupac,
                .fraction = null,
                .count = count,
                .seed = 11,
                .pair_mode = mode,
            };
            if (mode == .none) {
                if (try sampleExactFile(io, allocator, paths[0], &writer, count, options)) |failure| {
                    try std.testing.expectEqual(descriptors, try openDescriptorCount());
                    return pairAllocationFailure(.{ .command = .{ .input_index = 0, .details = failure } });
                }
            } else {
                const inputs = if (mode == .paired) paths[0..2] else paths[2..3];
                if (try sampleExactPairs(io, allocator, inputs, &writer, count, try deinterleaveStagingLimit(options.max_line_bytes), options)) |failure| {
                    try std.testing.expectEqual(descriptors, try openDescriptorCount());
                    return pairAllocationFailure(failure);
                }
            }
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
            const copies: usize = if (mode == .none) 1 else 2;
            try std.testing.expectEqual(@as(usize, @intCast(@min(count, 3))) * copies * (r1.len / 3), sink.written().len);
            if (count >= 3) try std.testing.expectEqualStrings(if (mode == .none) r1 else interleaved, sink.written());
        }
    }
}

fn openDescriptorCount() !usize {
    const io = std.testing.io;
    var directory = try std.Io.Dir.openDirAbsolute(io, "/proc/self/fd", .{ .iterate = true });
    defer directory.close(io);
    var entries = directory.iterate();
    var count: usize = 0;
    while (try entries.next(io)) |_| count += 1;
    return count;
}

test "[integration] - [paired inputs]: identity follows open descriptors after path replacement" {
    const Capture = struct {
        dir: std.Io.Dir,
        aliases: bool,
        opened: usize = 0,
        closed: usize = 0,
        early_reads: usize = 0,

        fn openFile(ctx: ?*anyopaque, _: std.Io.Dir, path: []const u8, options: std.Io.Dir.OpenFileOptions) std.Io.File.OpenError!std.Io.File {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            const io = std.testing.io;
            const file = try std.Io.Dir.cwd().openFile(io, path, options);
            errdefer file.close(io);
            self.opened += 1;
            if (self.opened == 1) {
                self.dir.rename("r1", self.dir, "saved", io) catch return error.Unexpected;
                if (self.aliases) {
                    self.dir.writeFile(io, .{ .sub_path = "r1", .data = "replacement" }) catch return error.Unexpected;
                } else {
                    self.dir.hardLink("r2", self.dir, "r1", io, .{}) catch return error.Unexpected;
                }
            }
            return file;
        }

        fn closeFiles(ctx: ?*anyopaque, files: []const std.Io.File) void {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            for (files) |file| file.close(std.testing.io);
            self.closed += files.len;
        }

        fn operate(ctx: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            if (operation == .file_read_streaming and self.aliases) self.early_reads += 1;
            return std.testing.io.operate(operation);
        }
    };

    const r1 = "@same\nA\n+\n!\n";
    const r2 = "@same\nT\n+\n#\n";
    for ([_]bool{ false, true }) |aliases| {
        var tmp = std.testing.tmpDir(.{});
        defer tmp.cleanup();
        try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "r1", .data = r1 });
        if (aliases) {
            try tmp.dir.hardLink("r1", tmp.dir, "r2", std.testing.io, .{});
        } else {
            try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "r2", .data = r2 });
        }
        var capture: Capture = .{ .dir = tmp.dir, .aliases = aliases };
        var vtable = std.Io.failing.vtable.*;
        vtable.dirOpenFile = Capture.openFile;
        vtable.fileClose = Capture.closeFiles;
        vtable.operate = Capture.operate;
        var output: [128]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        const descriptors = try openDescriptorCount();
        const r1_path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/r1", .{tmp.sub_path});
        defer std.testing.allocator.free(r1_path);
        const r2_path = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}/r2", .{tmp.sub_path});
        defer std.testing.allocator.free(r2_path);

        const failure = try interleaveInputs(
            .{ .userdata = &capture, .vtable = &vtable },
            std.testing.allocator,
            &.{ r1_path, r2_path },
            &writer,
            null,
            .{ .max_line_bytes = 8192, .alphabet = .iupac, .pair_name_policy = .illumina },
        );
        try writer.flush();

        if (aliases) {
            try std.testing.expectEqual(@as(u8, 2), failure.?.exitCode());
            try std.testing.expectEqualStrings("same_input", failure.?.command.details.code);
        } else {
            try std.testing.expect(failure == null);
        }
        try std.testing.expectEqualStrings(if (aliases) "" else r1 ++ r2, sink.written());
        try std.testing.expectEqual(@as(usize, 0), capture.early_reads);
        try std.testing.expectEqual(@as(usize, 1), capture.opened);
        try std.testing.expectEqual(@as(usize, 2), capture.closed);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
    }
}

test "[integration] - [output aliases]: rejection closes inputs on both exact passes" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const allocator = arena.allocator();
    const names = [_][]const u8{ "r1", "r2" };
    var paths: [2][]const u8 = undefined;
    var snapshots: [2]FileSnapshot = undefined;
    var identities: [2]FileIdentity = undefined;
    for (names, 0..) |name, index| {
        try tmp.dir.writeFile(io, .{ .sub_path = name, .data = "@same\nAC\n+\nII\n" });
        paths[index] = try std.fmt.allocPrint(allocator, ".zig-cache/tmp/{s}/{s}", .{ tmp.sub_path, name });
        const file = try tmp.dir.openFile(io, name, .{});
        defer file.close(io);
        snapshots[index] = try fileSnapshot(file, io);
        identities[index] = try fileIdentity(io, file);
    }
    const descriptors = try openDescriptorCount();
    for (0..8) |_| {
        var input1: RecordInput = undefined;
        var input2: RecordInput = undefined;
        const stream_failure = initRecordInput(&input1, io, paths[0], identities[0]).?;
        try std.testing.expectEqualStrings("same_output", stream_failure.code);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
        for ([_]?FileSnapshot{ null, snapshots[0] }) |snapshot| {
            const failure = initExactInput(&input1, io, paths[0], snapshot, identities[0]).failure;
            try std.testing.expectEqualStrings(if (snapshot == null) "same_output" else "input_changed", failure.code);
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
        }

        var output: [128]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        var options: SampleOptions = .{
            .max_line_bytes = 8192,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
            .output_identity = identities[0],
        };
        const single_failure = (try sampleExactSecondPass(true, io, allocator, paths[0], &writer, .empty, snapshots[0], 1, options)).?;
        try std.testing.expectEqualStrings("input_changed", single_failure.code);
        try std.testing.expectEqual(@as(u8, 3), single_failure.exit_code);
        const interleaved_failure = (try sampleExactInterleavedSecondPass(true, io, allocator, paths[0], &writer, .empty, snapshots[0], 1, try deinterleaveStagingLimit(options.max_line_bytes), options)).?;
        try std.testing.expectEqualStrings("input_changed", interleaved_failure.command.details.code);
        try std.testing.expectEqual(@as(u8, 3), interleaved_failure.exitCode());
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
        for (identities, 0..) |identity, side| {
            options.output_identity = identity;
            const paired_stream = initPairedRecordInputs(&input1, &input2, io, &paths, identity).?;
            try std.testing.expectEqualStrings("same_output", paired_stream.command.details.code);
            try std.testing.expectEqual(@as(u1, @intCast(side)), paired_stream.command.input_index);
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
            const paired_failure = (try sampleExactPairedSecondPass(true, io, allocator, &paths, &writer, .empty, snapshots, 1, options)).?;
            try std.testing.expectEqualStrings("input_changed", paired_failure.command.details.code);
            try std.testing.expectEqual(@as(u8, 3), paired_failure.exitCode());
            try std.testing.expectEqual(@as(u1, @intCast(side)), paired_failure.command.input_index);
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
        }
        try writer.flush();
        try std.testing.expectEqual(@as(usize, 0), sink.written().len);
    }
}

test "[integration] - [input resources]: repeated failures close owned files and preserve borrowed files" {
    const QuietStderr = struct {
        fn operate(_: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            if (operation == .file_write_streaming) {
                const write = operation.file_write_streaming;
                if (write.file.handle == std.Io.File.stderr().handle) {
                    var discarded = write.header.len;
                    for (write.data) |data| discarded += data.len * write.splat;
                    return .{ .file_write_streaming = discarded };
                }
            }
            return std.testing.io.operate(operation);
        }
    };

    const io = std.testing.io;
    var quiet_vtable = io.vtable.*;
    quiet_vtable.operate = QuietStderr.operate;
    const quiet_io = std.Io{ .userdata = io.userdata, .vtable = &quiet_vtable };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "r1", .data = "@a/1\nA\n+\n!\n" });
    try tmp.dir.createDir(io, "directory", .default_dir);
    const base = try std.fmt.allocPrint(std.testing.allocator, ".zig-cache/tmp/{s}", .{tmp.sub_path});
    defer std.testing.allocator.free(base);
    const path = try std.fmt.allocPrint(std.testing.allocator, "{s}/r1", .{base});
    defer std.testing.allocator.free(path);
    const directory = try std.fmt.allocPrint(std.testing.allocator, "{s}/directory", .{base});
    defer std.testing.allocator.free(directory);
    const missing = try std.fmt.allocPrint(std.testing.allocator, "{s}/missing", .{base});
    defer std.testing.allocator.free(missing);
    const output_path = try std.fmt.allocPrint(std.testing.allocator, "{s}/out1", .{base});
    defer std.testing.allocator.free(output_path);
    var changed = try snapshotTestFile(io, path);
    changed.size += 1;

    const descriptors = try openDescriptorCount();
    for (0..8) |_| {
        var input: RecordInput = undefined;
        try std.testing.expect(initRecordInput(&input, io, path, null) == null);
        input.deinit(io);
        try std.testing.expectEqualStrings("io_error", initRecordInput(&input, io, directory, null).?.code);
        try std.testing.expectEqualStrings("io_error", initExactInput(&input, io, directory, null, null).failure.code);
        try std.testing.expectEqualStrings("input_changed", initExactInput(&input, io, path, changed, null).failure.code);
        try std.testing.expectEqualStrings("input_changed", initExactInput(&input, io, missing, changed, null).failure.code);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());

        const pairs = [_][2][]const u8{
            .{ path, path },
            .{ path, missing },
            .{ missing, path },
            .{ directory, path },
            .{ path, directory },
        };
        for (pairs, 0..) |paths, index| {
            var input1: RecordInput = undefined;
            var input2: RecordInput = undefined;
            const expected = if (index == 0) "same_input" else "io_error";
            const stream_failure = initPairedRecordInputs(&input1, &input2, io, &paths, null).?;
            try std.testing.expectEqualStrings(expected, stream_failure.command.details.code);
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
            const exact_failure = initExactPairedInputs(&input1, &input2, io, &paths, null, null).failure;
            try std.testing.expectEqualStrings(expected, exact_failure.command.details.code);
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
        }

        {
            const borrowed = try tmp.dir.openFile(io, "r1", .{});
            defer borrowed.close(io);
            try input.init(io, borrowed, false);
            input.deinit(io);
            try std.testing.expectEqual(@as(u64, "@a/1\nA\n+\n!\n".len), (try borrowed.stat(io)).size);
            const invalid: std.Io.File = .{ .handle = -1, .flags = .{ .nonblocking = false } };
            for ([_][2]std.Io.File{ .{ invalid, borrowed }, .{ borrowed, invalid } }, 0..) |files, side| {
                const failure = pairedFilesFailure(io, files[0], files[1], null).?;
                try std.testing.expectEqual(@as(u8, 3), failure.exitCode());
                try std.testing.expectEqual(@as(u1, @intCast(side)), failure.command.input_index);
                try std.testing.expectEqualStrings("failed to inspect file", failure.command.details.message);
            }
        }
        var output: [64]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        const inputs = [_][]const u8{ path, missing };
        const failure = (try interleaveInputs(io, std.testing.allocator, &inputs, &writer, null, .{
            .max_line_bytes = 8192,
            .alphabet = .iupac,
            .pair_name_policy = .illumina,
        })).?;
        try std.testing.expectEqual(@as(u1, 1), failure.command.input_index);
        try std.testing.expectEqualStrings("io_error", failure.command.details.code);
        try std.testing.expectEqual(@as(usize, 0), sink.written().len);
        try std.testing.expectEqual(@as(u8, 3), runDeinterleave(quiet_io, std.testing.allocator, &.{path}, output_path, path, .{
            .max_line_bytes = 8192,
            .alphabet = .iupac,
            .pair_name_policy = .illumina,
        }));
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
        try std.testing.expectEqual(@as(u64, 0), (try snapshotTestFile(io, output_path)).size);
        try tmp.dir.deleteFile(io, "out1");
    }
}

fn pairAllocationFailure(failure: PairCommandFailure) error{ OutOfMemory, UnexpectedFailure } {
    return switch (failure) {
        .command => |details| if (std.mem.eql(u8, details.details.code, "out_of_memory"))
            error.OutOfMemory
        else
            error.UnexpectedFailure,
        .pair => error.UnexpectedFailure,
    };
}

test "[failure] - [exact sample]: allocation failures in both passes release memory and descriptors" {
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseExactAllocations,
        .{},
    );
}

test "[edge] - [paired exact sample]: target above pair count keeps indexes implicit" {
    var source1 = io_layer.SliceSource.init("@a/1\nA\n+\n!\n@b/1\nC\n+\n#\n");
    var source2 = io_layer.SliceSource.init("@a/2\nT\n+\n!\n@b/2\nG\n+\n#\n");
    var selector = sampling.ExactSelector.init(3, 11);
    defer selector.deinit(std.testing.allocator);

    const failure = checkPairedSources(
        std.testing.allocator,
        source1.byteSource(),
        source2.byteSource(),
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .pair_mode = .paired,
            .pair_name_policy = .illumina,
        },
        &selector,
    );

    try std.testing.expect(failure == null);
    try std.testing.expect(selector.finish() == .all);
    try std.testing.expectEqual(@as(usize, 0), selector.indexes.capacity);
}

test "[failure] - [paired exact sample]: a read failure follows complete earlier pairs" {
    const FailingSource = struct {
        data: []const u8,
        position: usize = 0,
        fail_at: usize,

        fn byteSource(self: *@This()) zfastq.io.ByteSource {
            return .{ .ctx = self, .vtable = &vtable };
        }

        const vtable = zfastq.io.ByteSource.VTable{ .read = read };

        fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
            const self: *@This() = @ptrCast(@alignCast(ctx));
            if (dest.len == 0) return 0;
            if (self.position == self.fail_at) return error.ReadFailed;
            const end = @min(self.fail_at, self.position + dest.len);
            const amount = end - self.position;
            @memcpy(dest[0..amount], self.data[self.position..end]);
            self.position = end;
            return amount;
        }
    };
    const first1 = "@a/1\nA\n+\n!\n";
    const first2 = "@a/2\nT\n+\n!\n";
    var source1 = FailingSource{
        .data = first1 ++ "@b/1\nC\n+\n#\n",
        .fail_at = first1.len,
    };
    var source2 = io_layer.SliceSource.init(first2 ++ "@b/2\nG\n+\n#\n");
    var selector = sampling.ExactSelector.init(1, 11);
    defer selector.deinit(std.testing.allocator);

    const failure = checkPairedSources(
        std.testing.allocator,
        source1.byteSource(),
        source2.byteSource(),
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .pair_mode = .paired,
            .pair_name_policy = .illumina,
        },
        &selector,
    ).?;

    const command = switch (failure) {
        .command => |command| command,
        .pair => return error.ExpectedCommandFailure,
    };
    try std.testing.expectEqual(@as(u1, 0), command.input_index);
    try std.testing.expectEqualStrings("io_error", command.details.code);
    try std.testing.expectEqual(@as(u64, 1), selector.record_count);
}

test "[property] - [interleaved fraction sample]: source chunks preserve selected pairs" {
    const ChunkedSource = struct {
        data: []const u8,
        chunk_len: usize,
        position: usize = 0,

        fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
            const self: *@This() = @ptrCast(@alignCast(ctx));
            const amount = @min(dest.len, self.chunk_len, self.data.len - self.position);
            @memcpy(dest[0..amount], self.data[self.position..][0..amount]);
            self.position += amount;
            return amount;
        }
    };
    const data = "@a/1 left\r\nAC\r\n+left\r\n!#\r\n@a/2 right\nGT\n+right\n$%\n" ++
        "@b/1\nG\n+\n&\n@b/2\nC\n+\n'\n" ++
        "@c/1\nTTA\n+one\n()*\n@c/2\nAAT\n+two\n+,-";
    // Seed 11 at one half selects pair indexes 0 and 2 in the frozen CLI vectors.
    const expected = "@a/1 left\nAC\n+left\n!#\n@a/2 right\nGT\n+right\n$%\n" ++
        "@c/1\nTTA\n+one\n()*\n@c/2\nAAT\n+two\n+,-\n";
    for (1..data.len + 1) |chunk_len| {
        var source: ChunkedSource = .{ .data = data, .chunk_len = chunk_len };
        var output: [expected.len]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        var selector = sampling.Selector.init(.{ .probability = 0.5 }, 11);
        var selection: PairOutputSelector = .{ .fraction = &selector };

        const failure = try sampleInterleavedSource(
            std.testing.allocator,
            .{ .ctx = &source, .vtable = &.{ .read = ChunkedSource.read } },
            &writer,
            &selection,
            1024,
            .{
                .max_line_bytes = 128,
                .alphabet = .iupac,
                .fraction = .{ .probability = 0.5 },
                .count = null,
                .seed = 11,
                .pair_mode = .interleaved,
            },
        );
        try writer.flush();

        try std.testing.expect(failure == null);
        try std.testing.expectEqual(data.len, source.position);
        try std.testing.expectEqualStrings(expected, sink.written());
    }
}

test "[edge] - [single fraction sample]: fraction zero avoids allocation in the default backend" {
    if (!build_options.use_isa_l) return error.SkipZigTest;
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(std.testing.allocator);
    try input.appendSlice(std.testing.allocator, "@large ");
    try input.appendNTimes(std.testing.allocator, 'x', io_layer.DEFAULT_READER_BUFFER_BYTES + 1);
    const header_len = input.items.len;
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(
        &path_buffer,
        ".zig-cache/tmp/{s}/input.fastq",
        .{tmp.sub_path},
    );

    for ([_][]const u8{ "\nA\n+\n!\n", "\nA\n+\n!!\n" }, 0..) |suffix, case_index| {
        input.shrinkRetainingCapacity(header_len);
        try input.appendSlice(std.testing.allocator, suffix);
        try tmp.dir.writeFile(io, .{ .sub_path = "input.fastq", .data = input.items });
        var output: [1]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        var selector = sampling.Selector.init(.none, 11);
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
            .fail_index = 0,
        });

        const failure = try sampleFractionInput(io, failing.allocator(), path, &writer, &selector, .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = .none,
            .count = null,
            .seed = 11,
        });
        try std.testing.expect(!failing.has_induced_failure);
        try std.testing.expectEqual(@as(usize, 0), sink.written().len);
        if (case_index == 0) {
            try std.testing.expect(failure == null);
        } else {
            const details = failure orelse return error.ExpectedFailure;
            try std.testing.expectEqualStrings("S005", details.code);
            try std.testing.expectEqualStrings("sequence and quality lengths differ", details.message);
            try std.testing.expectEqual(@as(u8, 1), details.exit_code);
            try std.testing.expectEqual(@as(?u64, 0), details.record_index);
            try std.testing.expectEqual(@as(?u64, header_len + 5), details.byte_offset);
            try std.testing.expectEqual(@as(?u3, 4), details.line_in_record);
        }
    }
}

test "[edge] - [paired fraction sample]: fraction zero keeps buffered mate one borrowed" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var input_bytes: std.ArrayList(u8) = .empty;
    defer input_bytes.deinit(std.testing.allocator);
    try input_bytes.appendSlice(std.testing.allocator, "@large/1 ");
    try input_bytes.appendNTimes(std.testing.allocator, 'x', 128 * 1024);
    try input_bytes.appendSlice(
        std.testing.allocator,
        "\nA\n+\n!\n@large/2\nT\n+\n#\n",
    );
    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = input_bytes.items,
    });
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(
        &path_buffer,
        ".zig-cache/tmp/{s}/pairs.fastq",
        .{tmp.sub_path},
    );
    var output: [1]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    var selector = sampling.Selector.init(.none, 11);
    var output_selector: PairOutputSelector = .{ .fraction = &selector };
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 1,
    });

    const failure = try sampleInterleavedFractionInput(
        io,
        failing.allocator(),
        path,
        &writer,
        &output_selector,
        0,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = .none,
            .count = null,
            .seed = 11,
            .pair_mode = .interleaved,
        },
    );
    try std.testing.expect(failure == null);
    try std.testing.expect(!failing.has_induced_failure);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}

const InterleavedTestSource = struct {
    data: []const u8,
    first_chunk: usize,
    position: usize = 0,

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *@This() = @ptrCast(@alignCast(ctx));
        const size = if (self.position == 0) self.first_chunk else 3;
        const n = @min(size, @min(dest.len, self.data.len - self.position));
        @memcpy(dest[0..n], self.data[self.position..][0..n]);
        self.position += n;
        return n;
    }
};

test "[property] - [interleaved input]: small pairs need no preservation allocations" {
    const expected1 = "@pair/1\nR\n+left\n!\n";
    const expected2 = "@pair/2\nT\n+right\n#\n";
    for ([_][2][]const u8{
        .{ expected1, expected2 },
        .{ "@pair/1\r\nR\r\n+left\r\n!\r\n", "@pair/2\r\nT\r\n+right\r\n#\r\n" },
        .{ "@pair/1\r\nR\n+left\r\n!\n", "@pair/2\nT\r\n+right\n#\r\n" },
    }) |mates| {
        const input = try std.mem.concat(std.testing.allocator, u8, &mates);
        defer std.testing.allocator.free(input);
        for ([_]usize{ mates[0].len, mates[0].len + 7 }) |first_chunk| {
            for (0..4) |command| {
                var source: InterleavedTestSource = .{ .data = input, .first_chunk = first_chunk };
                const bytes: zfastq.io.ByteSource = .{ .ctx = &source, .vtable = &.{ .read = InterleavedTestSource.read } };
                var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 2 });
                var out1: [128]u8 = undefined;
                var out2: [128]u8 = undefined;
                var sink1 = io_layer.SliceSink.init(&out1);
                var sink2 = io_layer.SliceSink.init(&out2);
                var writer1 = zfastq.Writer.init(sink1.byteSink());
                var writer2 = zfastq.Writer.init(sink2.byteSink());
                var selector = sampling.Selector.init(if (command == 1) .none else .all, 11);
                var selection: PairOutputSelector = .{ .fraction = &selector };
                const failure = switch (command) {
                    0 => checkInterleavedSource(failing.allocator(), bytes, .{
                        .max_line_bytes = 64,
                        .alphabet = .iupac,
                        .pair_mode = .interleaved,
                        .pair_name_policy = .illumina,
                    }, null),
                    1, 2 => try sampleInterleavedSource(failing.allocator(), bytes, &writer1, &selection, 260, .{
                        .max_line_bytes = 64,
                        .alphabet = .iupac,
                        .fraction = if (command == 1) .none else .all,
                        .count = null,
                        .seed = 11,
                        .pair_mode = .interleaved,
                        .pair_name_policy = .illumina,
                    }),
                    3 => try deinterleaveSource(failing.allocator(), bytes, &writer1, &writer2, 260, .{
                        .max_line_bytes = 64,
                        .alphabet = .iupac,
                        .pair_name_policy = .illumina,
                    }),
                    else => unreachable,
                };

                try std.testing.expect(failure == null);
                try std.testing.expect(!failing.has_induced_failure);
                try std.testing.expectEqualStrings(switch (command) {
                    0, 1 => "",
                    2 => expected1 ++ expected2,
                    3 => expected1,
                    else => unreachable,
                }, sink1.written());
                try std.testing.expectEqualStrings(if (command == 3) expected2 else "", sink2.written());
            }
        }
    }
}

test "[failure] - [deinterleave]: preserves CRLF output on a failed write" {
    const first = "@pair/1\r\nR\r\n+left\r\n!\r\n";
    const second = "@pair/2\r\nT\r\n+right\r\n#\r\n";
    var source: InterleavedTestSource = .{ .data = first ++ second, .first_chunk = first.len };
    var output: ["@pair/1\nR\n+left\n!\n".len - 1]u8 = undefined;
    var sink1 = io_layer.SliceSink.init(&output);
    var sink2 = io_layer.SliceSink.init(&.{});
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());

    try std.testing.expectError(error.Output1WriteFailed, deinterleaveSource(
        std.testing.allocator,
        .{ .ctx = &source, .vtable = &.{ .read = InterleavedTestSource.read } },
        &writer1,
        &writer2,
        260,
        .{ .max_line_bytes = 64, .alphabet = .iupac, .pair_name_policy = .illumina },
    ));
    // Canonical staging submits the complete mate to this all-or-error sink.
    try std.testing.expectEqualStrings("", sink1.written());
    try std.testing.expectEqualStrings("", sink2.written());
}

const DeinterleaveTestSink = struct {
    buffer: [256]u8 = undefined,
    length: usize = 0,
    fail_write: bool = false,
    fail_flush: bool = false,
    flush_count: usize = 0,

    fn byteSink(self: *DeinterleaveTestSink) zfastq.io.ByteSink {
        return .{ .ctx = self, .vtable = &vtable };
    }

    const vtable = zfastq.io.ByteSink.VTable{
        .write = write,
        .flush = flush,
    };

    fn write(ctx: *anyopaque, bytes: []const u8) error{WriteFailed}!void {
        const self: *DeinterleaveTestSink = @ptrCast(@alignCast(ctx));
        if (self.fail_write) return error.WriteFailed;
        const end = std.math.add(usize, self.length, bytes.len) catch
            return error.WriteFailed;
        if (end > self.buffer.len) return error.WriteFailed;
        @memcpy(self.buffer[self.length..end], bytes);
        self.length = end;
    }

    fn flush(ctx: *anyopaque) error{WriteFailed}!void {
        const self: *DeinterleaveTestSink = @ptrCast(@alignCast(ctx));
        self.flush_count += 1;
        if (self.fail_flush) return error.WriteFailed;
    }
};

test "[failure] - [deinterleave]: identifies both output write positions" {
    const input = "@pair/1\nA\n+left\n!\n@pair/2\nT\n+right\n#\n";
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };
    const staging_limit = try deinterleaveStagingLimit(options.max_line_bytes);

    {
        var source = io_layer.SliceSource.init(input);
        var sink1 = DeinterleaveTestSink{ .fail_write = true };
        var sink2 = DeinterleaveTestSink{};
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());

        try std.testing.expectError(
            error.Output1WriteFailed,
            deinterleaveSource(
                std.testing.allocator,
                source.byteSource(),
                &writer1,
                &writer2,
                staging_limit,
                options,
            ),
        );
        try std.testing.expectEqual(@as(usize, 0), sink2.length);
    }

    {
        var source = io_layer.SliceSource.init(input);
        var sink1 = DeinterleaveTestSink{};
        var sink2 = DeinterleaveTestSink{ .fail_write = true };
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());

        try std.testing.expectError(
            error.Output2WriteFailed,
            deinterleaveSource(
                std.testing.allocator,
                source.byteSource(),
                &writer1,
                &writer2,
                staging_limit,
                options,
            ),
        );
        try std.testing.expectEqualStrings("@pair/1\nA\n+left\n!\n", sink1.buffer[0..sink1.length]);
    }
}

test "[failure] - [deinterleave]: flushes outputs in order and stops after failure" {
    var sink1 = DeinterleaveTestSink{ .fail_flush = true };
    var sink2 = DeinterleaveTestSink{};
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());

    try std.testing.expectError(
        error.Output1WriteFailed,
        flushDeinterleaveWriters(&writer1, &writer2),
    );
    try std.testing.expectEqual(@as(usize, 1), sink1.flush_count);
    try std.testing.expectEqual(@as(usize, 0), sink2.flush_count);

    sink1.fail_flush = false;
    sink2.fail_flush = true;
    try std.testing.expectError(
        error.Output2WriteFailed,
        flushDeinterleaveWriters(&writer1, &writer2),
    );
    try std.testing.expectEqual(@as(usize, 2), sink1.flush_count);
    try std.testing.expectEqual(@as(usize, 1), sink2.flush_count);
}

test "[integration] - [deinterleave]: output flush failures preserve diagnostics and exit precedence" {
    const Capture = struct {
        dir: std.Io.Dir,
        stderr: std.Io.Writer,
        fail_write: usize,
        writes: usize = 0,
        closed: usize = 0,

        fn openFile(
            ctx: ?*anyopaque,
            _: std.Io.Dir,
            path: []const u8,
            options: std.Io.Dir.OpenFileOptions,
        ) std.Io.File.OpenError!std.Io.File {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            return self.dir.openFile(std.testing.io, path, options);
        }

        fn createFile(
            ctx: ?*anyopaque,
            _: std.Io.Dir,
            path: []const u8,
            options: std.Io.Dir.CreateFileOptions,
        ) std.Io.File.OpenError!std.Io.File {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            return self.dir.createFile(std.testing.io, path, options);
        }

        fn closeFiles(ctx: ?*anyopaque, files: []const std.Io.File) void {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            for (files) |file| file.close(std.testing.io);
            self.closed += files.len;
        }

        fn readPositional(
            _: ?*anyopaque,
            file: std.Io.File,
            data: []const []u8,
            offset: u64,
        ) std.Io.File.ReadPositionalError!usize {
            const io = std.testing.io;
            return io.vtable.fileReadPositional(io.userdata, file, data, offset);
        }

        fn writePositional(
            ctx: ?*anyopaque,
            file: std.Io.File,
            header: []const u8,
            data: []const []const u8,
            splat: usize,
            offset: u64,
        ) std.Io.File.WritePositionalError!usize {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            self.writes += 1;
            if (self.writes == self.fail_write) return error.NoSpaceLeft;
            const io = std.testing.io;
            return io.vtable.fileWritePositional(io.userdata, file, header, data, splat, offset);
        }

        fn operate(ctx: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
            const self: *@This() = @ptrCast(@alignCast(ctx.?));
            if (operation == .file_write_streaming) {
                const write = operation.file_write_streaming;
                if (write.file.handle == std.Io.File.stderr().handle) {
                    return .{ .file_write_streaming = self.stderr.writeSplatHeader(
                        write.header,
                        write.data,
                        write.splat,
                    ) catch error.NoSpaceLeft };
                }
            }
            return std.testing.io.operate(operation);
        }
    };

    const r1 = "@ok/1\nA\n+\n!\n";
    const r2 = "@ok/2\nT\n+\n#\n";
    const cases = [_]struct { tail: []const u8, exit_code: u8, stderr: []const u8 }{
        .{ .tail = "", .exit_code = 0, .stderr = "" },
        .{
            .tail = "@odd/1\nA\n+\n!\n",
            .exit_code = 1,
            .stderr = "error: input: P002: paired input is missing a mate " ++
                "(pair 1, remaining R1, last R1 record 2, last R2 record 1)\n",
        },
        .{
            .tail = "@too-long-header\nA\n+\n!\n",
            .exit_code = 4,
            .stderr = "error: input: line length limit exceeded\n",
        },
    };
    var vtable = std.Io.failing.vtable.*;
    vtable.dirOpenFile = Capture.openFile;
    vtable.dirCreateFile = Capture.createFile;
    vtable.fileClose = Capture.closeFiles;
    vtable.fileReadPositional = Capture.readPositional;
    vtable.fileWritePositional = Capture.writePositional;
    vtable.operate = Capture.operate;
    for (cases) |case| {
        for (0..3) |fail_write| {
            var tmp = std.testing.tmpDir(.{});
            defer tmp.cleanup();
            var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
            defer arena.deinit();
            const allocator = arena.allocator();
            const bytes = try std.mem.concat(allocator, u8, &.{ r1, r2, case.tail });
            try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "input", .data = bytes });
            var stderr_buffer: [512]u8 = undefined;
            var capture: Capture = .{
                .dir = tmp.dir,
                .stderr = .fixed(&stderr_buffer),
                .fail_write = fail_write,
            };

            const status = runDeinterleave(
                .{ .userdata = &capture, .vtable = &vtable },
                std.testing.allocator,
                &.{"input"},
                "out1",
                "out2",
                .{ .max_line_bytes = 8, .alphabet = .iupac, .pair_name_policy = .illumina },
            );

            const expected_status = if (fail_write == 0) case.exit_code else @max(case.exit_code, 3);
            const expected_stderr = if (fail_write == 0) case.stderr else try std.fmt.allocPrint(
                allocator,
                "{s}error: out{d}: I/O error\n",
                .{ case.stderr, fail_write },
            );
            try std.testing.expectEqual(expected_status, status);
            try std.testing.expectEqualStrings(expected_stderr, capture.stderr.buffered());
            try std.testing.expectEqual(@as(usize, if (fail_write == 1) 1 else 2), capture.writes);
            try std.testing.expectEqual(@as(usize, 3), capture.closed);
            const output1 = try tmp.dir.readFileAlloc(std.testing.io, "out1", allocator, .limited(64));
            const output2 = try tmp.dir.readFileAlloc(std.testing.io, "out2", allocator, .limited(64));
            try std.testing.expectEqualStrings(if (fail_write == 1) "" else r1, output1);
            try std.testing.expectEqualStrings(if (fail_write == 0) r2 else "", output2);
        }
    }
}

test "[failure] - [deinterleave]: staging allocation failure emits no output" {
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(std.testing.allocator);
    try input.appendSlice(std.testing.allocator, "@pair/1\nA\n+\n!\n@pair/2\n");
    try input.appendNTimes(
        std.testing.allocator,
        'T',
        zfastq.limits.DEFAULT_READER_BUFFER_BYTES,
    );
    try input.appendSlice(std.testing.allocator, "\n+\n");
    try input.appendNTimes(
        std.testing.allocator,
        '#',
        zfastq.limits.DEFAULT_READER_BUFFER_BYTES,
    );
    try input.append(std.testing.allocator, '\n');
    var source = io_layer.SliceSource.init(input.items);
    var sink1 = DeinterleaveTestSink{};
    var sink2 = DeinterleaveTestSink{};
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 2,
    });
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };

    const failure = (try deinterleaveSource(
        failing.allocator(),
        source.byteSource(),
        &writer1,
        &writer2,
        try deinterleaveStagingLimit(options.max_line_bytes),
        options,
    )).?;

    try std.testing.expect(failing.has_induced_failure);
    try std.testing.expectEqual(@as(u8, 3), failure.exitCode());
    const command_failure = switch (failure) {
        .command => |command| command.details,
        .pair => return error.ExpectedCommandFailure,
    };
    try std.testing.expectEqualStrings("out_of_memory", command_failure.code);
    try std.testing.expectEqual(@as(usize, 0), sink1.length);
    try std.testing.expectEqual(@as(usize, 0), sink2.length);
    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseRetainedPairAllocations,
        .{input.items},
    );
}

test "[failure] - [deinterleave]: both fallback owners survive allocation and output failures" {
    const first_field_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES - 7;
    const second_field_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES;
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(std.testing.allocator);
    try input.appendSlice(std.testing.allocator, "@pair/1\n");
    try input.appendNTimes(std.testing.allocator, 'A', first_field_len);
    try input.appendSlice(std.testing.allocator, "\n+\n");
    try input.appendNTimes(std.testing.allocator, '!', first_field_len);
    try input.appendSlice(std.testing.allocator, "\n@pair/2\n");
    try input.appendNTimes(std.testing.allocator, 'T', second_field_len);
    try input.appendSlice(std.testing.allocator, "\n+\n");
    try input.appendNTimes(std.testing.allocator, '#', second_field_len);
    try input.append(std.testing.allocator, '\n');

    {
        var source = io_layer.SliceSource.init(input.items);
        var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
        defer reader.deinit();
        var retained: fastq.RetainedRecordStorage = .{};
        defer retained.deinit(std.testing.allocator);
        var staging: std.ArrayList(u8) = .empty;
        defer staging.deinit(std.testing.allocator);
        var span: ?[]const u8 = null;
        var validator = fastq.AdaptiveRecordValidator.init(.{});
        const record1 = (try fastq.nextWithoutId(&reader, &span)).?;
        try std.testing.expect((try fastq.nextBufferedWithoutId(&reader, &span)) == null);
        const next = try nextAfterPreservingInterleavedMate1(
            std.testing.allocator,
            &reader,
            &retained,
            &staging,
            record1,
            null,
            try deinterleaveStagingLimit(reader.options.max_line_bytes),
            &validator,
        );
        try std.testing.expect(next.first_storage == .retained);
        try std.testing.expectEqual(@as(usize, 0), staging.capacity);
        try std.testing.expect(std.mem.allEqual(u8, record1.sequence, 'A'));
        try std.testing.expectEqual(second_field_len, next.record.?.record.sequence.len);
        try std.testing.expect(std.mem.allEqual(u8, next.record.?.record.sequence, 'T'));
        try std.testing.expect(next.record.?.semantic_error == null);
    }

    try std.testing.checkAllAllocationFailures(
        std.testing.allocator,
        exerciseRetainedPairAllocations,
        .{input.items},
    );

    const output = try std.testing.allocator.alloc(u8, input.items.len);
    defer std.testing.allocator.free(output);
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };

    for (0..3) |failure_mode| {
        var source = io_layer.SliceSource.init(input.items);
        var sink1 = io_layer.SliceSink.init(if (failure_mode == 0) output[0..0] else output);
        var sink2 = io_layer.SliceSink.init(output[0..0]);
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());
        const result = deinterleaveSource(
            std.testing.allocator,
            source.byteSource(),
            &writer1,
            &writer2,
            if (failure_mode == 2) 1 else try deinterleaveStagingLimit(options.max_line_bytes),
            options,
        );
        if (failure_mode < 2) {
            try std.testing.expectError(
                if (failure_mode == 0) error.Output1WriteFailed else error.Output2WriteFailed,
                result,
            );
            if (failure_mode == 0) try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
            if (failure_mode == 1) {
                try std.testing.expectEqualStrings(input.items[0 .. 2 * first_field_len + 12], sink1.written());
            }
        } else {
            const failure = (try result).?;
            try std.testing.expectEqualStrings("arithmetic_limit", failure.command.details.code);
            try std.testing.expectEqual(@as(usize, 0), sink1.written().len);
            try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
        }
    }

    const split = std.mem.find(u8, input.items, "@pair/2").?;
    for (0..5) |failure_mode| {
        if (failure_mode == 1) input.items[input.items.len - 2] = ' ';
        if (failure_mode == 2) input.items[split + 1] = 'x';
        if (failure_mode == 4) input.items["@pair/1\n".len] = '.';
        const end = switch (failure_mode) {
            0, 4 => input.items.len - 2,
            3 => split,
            else => input.items.len,
        };
        var source = io_layer.SliceSource.init(input.items[0..end]);
        var sink1 = DeinterleaveTestSink{};
        var sink2 = DeinterleaveTestSink{};
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());
        const failure = (try deinterleaveSource(
            std.testing.allocator,
            source.byteSource(),
            &writer1,
            &writer2,
            try deinterleaveStagingLimit(options.max_line_bytes),
            options,
        )).?;
        switch (failure_mode) {
            0, 4 => try std.testing.expectEqualStrings("S005", failure.command.details.code),
            1 => try std.testing.expectEqualStrings("S006", failure.command.details.code),
            2 => try std.testing.expect(failure.pair == .name_mismatch),
            3 => try std.testing.expect(failure.pair == .count_mismatch),
            else => unreachable,
        }
        try std.testing.expectEqual(@as(usize, 0), sink1.length);
        try std.testing.expectEqual(@as(usize, 0), sink2.length);
        input.items[input.items.len - 2] = '#';
        input.items[split + 1] = 'p';
        input.items["@pair/1\n".len] = 'A';
    }
}

fn exerciseRetainedPairAllocations(allocator: std.mem.Allocator, input: []const u8) !void {
    const output = try std.testing.allocator.alloc(u8, input.len);
    defer std.testing.allocator.free(output);
    const split = std.mem.find(u8, input, "@pair/2").?;
    var source = io_layer.SliceSource.init(input);
    var sink1 = io_layer.SliceSink.init(output[0..split]);
    var sink2 = io_layer.SliceSink.init(output[split..]);
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };
    if (try deinterleaveSource(
        allocator,
        source.byteSource(),
        &writer1,
        &writer2,
        split,
        options,
    )) |failure| {
        try std.testing.expectEqual(@as(usize, 0), sink1.written().len);
        try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
        return pairAllocationFailure(failure);
    }
    try std.testing.expectEqualStrings(input[0..split], sink1.written());
    try std.testing.expectEqualStrings(input[split..], sink2.written());
}

test "[property] - [byte display]: human and JSON encode all byte values" {
    var bytes: [256]u8 = undefined;
    for (&bytes, 0..) |*byte, index| byte.* = @intCast(index);
    const expected = "\\x00\\x01\\x02\\x03\\x04\\x05\\x06\\x07\\x08\\x09\\x0A\\x0B\\x0C\\x0D\\x0E\\x0F\\x10\\x11" ++
        "\\x12\\x13\\x14\\x15\\x16\\x17\\x18\\x19\\x1A\\x1B\\x1C\\x1D\\x1E\\x1F !\"#$%&'()*+,-./" ++
        "0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\\\]^_`abcdefghijklmnopqrstuv" ++
        "wxyz{|}~\\x7F\\x80\\x81\\x82\\x83\\x84\\x85\\x86\\x87\\x88\\x89\\x8A\\x8B\\x8C\\x8D\\x8E" ++
        "\\x8F\\x90\\x91\\x92\\x93\\x94\\x95\\x96\\x97\\x98\\x99\\x9A\\x9B\\x9C\\x9D\\x9E\\x9F\\xA0" ++
        "\\xA1\\xA2\\xA3\\xA4\\xA5\\xA6\\xA7\\xA8\\xA9\\xAA\\xAB\\xAC\\xAD\\xAE\\xAF\\xB0\\xB1\\xB2" ++
        "\\xB3\\xB4\\xB5\\xB6\\xB7\\xB8\\xB9\\xBA\\xBB\\xBC\\xBD\\xBE\\xBF\\xC0\\xC1\\xC2\\xC3\\xC4" ++
        "\\xC5\\xC6\\xC7\\xC8\\xC9\\xCA\\xCB\\xCC\\xCD\\xCE\\xCF\\xD0\\xD1\\xD2\\xD3\\xD4\\xD5\\xD6" ++
        "\\xD7\\xD8\\xD9\\xDA\\xDB\\xDC\\xDD\\xDE\\xDF\\xE0\\xE1\\xE2\\xE3\\xE4\\xE5\\xE6\\xE7\\xE8" ++
        "\\xE9\\xEA\\xEB\\xEC\\xED\\xEE\\xEF\\xF0\\xF1\\xF2\\xF3\\xF4\\xF5\\xF6\\xF7\\xF8\\xF9\\xFA" ++
        "\\xFB\\xFC\\xFD\\xFE\\xFF";

    var human_buffer: [1024]u8 = undefined;
    var human = std.Io.Writer.fixed(&human_buffer);
    try writeEscapedBytes(&human, &bytes, false);
    try std.testing.expectEqualStrings(expected, human.buffered());

    var json_buffer: [2048]u8 = undefined;
    var output = std.Io.Writer.fixed(&json_buffer);
    var json: std.json.Stringify = .{ .writer = &output };
    try writeEscapedJsonString(&json, &bytes);
    var parsed = try std.json.parseFromSlice([]const u8, std.testing.allocator, output.buffered(), .{});
    defer parsed.deinit();
    try std.testing.expectEqualStrings(expected, parsed.value);
}
