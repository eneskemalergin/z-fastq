//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("root.zig");
const fastq = @import("fastq.zig");
const io_layer = @import("io.zig");
const pairing = @import("pair.zig");
const sampling = @import("sample.zig");
const statistics = @import("stats.zig");

// z-fastq never opens a network connection; this drops std.Io's DNS and socket code.
pub const std_options: std.Options = .{ .networking = false };

// --- Arguments ---

const Command = enum { count, stats, check, sample, interleave, deinterleave };

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

fn argumentError(io: std.Io, message: []const u8, exit_code: u8) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
    std.process.exit(exit_code);
}

fn printUsageAndExit(io: std.Io) noreturn {
    argumentError(io, USAGE, 2);
}

fn printHelpAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, USAGE) catch std.process.exit(stdoutFailed(io, 0));
    std.process.exit(0);
}

/// Value options seen so far: each may be given once, so a later value never
/// silently replaces an earlier one.
const GivenOptions = struct {
    max_line_bytes: bool = false,
    alphabet: bool = false,
    fraction: bool = false,
    count: bool = false,
    seed: bool = false,
    pair_names: bool = false,
    out1: bool = false,
    out2: bool = false,
};

fn onlyOnce(io: std.Io, given: *bool, comptime option: []const u8) void {
    if (given.*) argumentError(io, "error: " ++ option ++ " may appear only once\n", 2);
    given.* = true;
}

fn printVersionAndExit(io: std.Io) noreturn {
    var buf: [64]u8 = undefined;
    const line = std.fmt.bufPrint(&buf, "z-fastq {s}\n", .{zfastq.VERSION}) catch "z-fastq\n";
    std.Io.File.writeStreamingAll(.stdout(), io, line) catch std.process.exit(stdoutFailed(io, 0));
    std.process.exit(0);
}

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

    const output_identity = checkedStdoutIdentity(io) catch return 3;
    return inputAliasFailure(io, inputs, output_identity);
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;
    // A write to a closed pipe then fails with EPIPE, which ends the command
    // with the documented I/O status instead of a signal.
    std.posix.sigaction(.PIPE, &.{
        .handler = .{ .handler = std.posix.SIG.IGN },
        .mask = std.posix.sigemptyset(),
        .flags = 0,
    }, null);

    var args = std.process.Args.Iterator.initAllocator(init.minimal.args, gpa) catch {
        argumentError(io, OUT_OF_MEMORY_LINE, 3);
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
    var given: GivenOptions = .{};
    var positional = std.ArrayList([]const u8).empty;
    defer positional.deinit(gpa);

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--")) {
            while (args.next()) |path| {
                positional.append(gpa, path) catch {
                    argumentError(io, OUT_OF_MEMORY_LINE, 3);
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
            onlyOnce(io, &given.max_line_bytes, "--max-line-bytes");
            max_line_bytes = sampling.parseDecimal(usize, value) catch |err| switch (err) {
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
            onlyOnce(io, &given.alphabet, "--alphabet");
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
            onlyOnce(io, &given.fraction, "--fraction");
            fraction = sampling.Fraction.parse(value) catch {
                argumentError(io, "error: invalid --fraction value\n", 2);
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--count")) {
            const value = args.next() orelse {
                argumentError(io, "error: --count requires a value\n", 2);
            };
            onlyOnce(io, &given.count, "--count");
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
            onlyOnce(io, &given.seed, "--seed");
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
            onlyOnce(io, &given.pair_names, "--pair-names");
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
            if (std.mem.eql(u8, arg, "--out1")) {
                onlyOnce(io, &given.out1, "--out1");
                output1 = value;
            } else {
                onlyOnce(io, &given.out2, "--out2");
                output2 = value;
            }
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
            argumentError(io, OUT_OF_MEMORY_LINE, 3);
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

// --- Command failures and diagnostics ---

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

const STATISTICS_LIMIT = CommandFailure.plain("arithmetic_limit", "statistics arithmetic limit exceeded", 4);

const INSPECT_FAILURE = CommandFailure.plain("io_error", "failed to inspect file", 3);

const RECORD_LOCATION_UNAVAILABLE = CommandFailure.plain("io_error", "record location is unavailable", 3);

/// For failures reported before a command's input labels exist.
const OUT_OF_MEMORY_LINE = "error: " ++ OUT_OF_MEMORY.message ++ "\n";

/// Output buffers: a JSON document is flushed when it outgrows its buffer;
/// FASTQ output is written in whole pages from its buffer.
const JSON_OUTPUT_BUFFER_BYTES = 16 * 1024;
const RECORD_OUTPUT_BUFFER_BYTES = 64 * 1024;

const PAIR_DIAGNOSTIC_PREFIX_BYTES = 128;

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

    /// Stable validation code, shared by text and JSON output.
    fn code(self: PairFailure) []const u8 {
        return switch (self) {
            .name_mismatch => "P001",
            .count_mismatch => "P002",
        };
    }

    fn message(self: PairFailure) []const u8 {
        return switch (self) {
            .name_mismatch => "paired identifiers or mate markers do not match",
            .count_mismatch => "paired input is missing a mate",
        };
    }
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

fn mapCurrentSemanticFailure(
    maybe_error: ?zfastq.SemanticError,
    reader: *const zfastq.Reader,
    record_index: u64,
) ?CommandFailure {
    const semantic_error = maybe_error orelse return null;
    const offsets = reader.currentRecordOffsets() orelse return RECORD_LOCATION_UNAVAILABLE;
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

fn pairCommandFailure(input_index: u1, details: CommandFailure) PairCommandFailure {
    return .{ .command = .{ .input_index = input_index, .details = details } };
}

fn unwritableRecordFailure() CommandFailure {
    return CommandFailure.plain(
        "unwritable_record",
        "record fields ending in CR cannot be written with LF endings",
        1,
    );
}

/// In interleaved input, mate 1 of pair `p` is record `2p`.
fn interleavedPairIndex(mate1_record_index: u64) u64 {
    return mate1_record_index / 2;
}

fn missingInterleavedMateFailure(record_index1: u64) PairCommandFailure {
    return .{ .pair = .{ .count_mismatch = .{
        .pair_index = interleavedPairIndex(record_index1),
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
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, failure.code()) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, failure.message()) catch {};
    switch (failure) {
        .name_mismatch => |details| {
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

fn printOutputFailure(io: std.Io) void {
    std.Io.File.writeStreamingAll(
        .stderr(),
        io,
        "error: standard output: " ++ IO_FAILURE.message ++ "\n",
    ) catch {};
}

/// Reports a failed write to stdout and raises the exit status to the I/O class.
fn stdoutFailed(io: std.Io, exit_code: u8) u8 {
    printOutputFailure(io);
    return @max(exit_code, 3);
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

// --- Input opening ---

const InputOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
};

const FileIdentity = struct {
    device_major: u32,
    device_minor: u32,
    inode: u64,

    // std.Io.File.Stat does not retain the device containing the inode, so identities come from statx.
    fn of(stat: std.os.linux.Statx) FileIdentity {
        return .{ .device_major = stat.dev_major, .device_minor = stat.dev_minor, .inode = stat.ino };
    }
};

const GZIP_INPUT_BUFFER_BYTES = 64 * 1024;

const RecordInput = struct {
    file: std.Io.File,
    owns_file: bool,
    read_buffer: [zfastq.limits.DEFAULT_READER_BUFFER_BYTES]u8,
    file_reader: std.Io.File.Reader,
    kind: enum { plain, gzip },
    // Kept apart from `kind`: assigning null or a payload-free union tag zeroes this whole field.
    gzip: io_layer.GzipSource,

    fn init(
        self: *RecordInput,
        io: std.Io,
        file: std.Io.File,
        owns_file: bool,
    ) error{Io}!void {
        self.file = file;
        self.owns_file = owns_file;
        self.file_reader = file.readerStreaming(io, self.read_buffer[0..GZIP_INPUT_BUFFER_BYTES]);
        const prefix = self.file_reader.interface.peek(2) catch |err| switch (err) {
            error.EndOfStream => null,
            error.ReadFailed => return error.Io,
        };
        if (prefix) |bytes| {
            if (std.mem.eql(u8, bytes, &.{ 0x1f, 0x8b })) {
                self.gzip = io_layer.GzipSource.init(&self.file_reader.interface);
                self.kind = .gzip;
                return;
            }
        }
        const reader = &self.file_reader.interface;
        const sniff_filled = reader.bufferedLen() == reader.buffer.len;
        // The sniff buffer is a prefix of `read_buffer`, so widening keeps the sniffed bytes.
        reader.buffer = &self.read_buffer;
        // A full sniff read can end inside a record; topping up keeps it out of spill storage.
        if (sniff_filled) {
            reader.fillMore() catch |err| switch (err) {
                error.EndOfStream => {},
                error.ReadFailed => return error.Io,
            };
        }
        self.kind = .plain;
    }

    fn deinit(self: *RecordInput, io: std.Io) void {
        if (self.owns_file) self.file.close(io);
        self.* = undefined;
    }

    fn stream(self: *RecordInput) *std.Io.Reader {
        return switch (self.kind) {
            .plain => &self.file_reader.interface,
            .gzip => io_layer.gzipReader(&self.gzip),
        };
    }

    fn readScannerChunk(self: *RecordInput) error{Io}!?[]const u8 {
        return io_layer.readChunk(self.stream()) catch error.Io;
    }
};

fn initRecordInput(
    input: *RecordInput,
    io: std.Io,
    label: []const u8,
    output_identity: ?FileIdentity,
) ?CommandFailure {
    const file = openRecordFile(io, label) catch |err| return inputOpenFailure(err);
    const owns_file = !namesStdin(label);
    var transferred = false;
    defer if (!transferred and owns_file) file.close(io);
    if (inputTypeFailure(io, file)) |failure| return failure;
    if (outputFileFailure(io, file, output_identity)) |failure| return failure;
    input.init(io, file, owns_file) catch
        return IO_FAILURE;
    transferred = true;
    return null;
}

/// A directory opens like a file but fails on its first read; name it instead.
fn inputTypeFailure(io: std.Io, file: std.Io.File) ?CommandFailure {
    const stat = statRecordFile(io, file) catch return null;
    if (stat.mode & std.os.linux.S.IFMT == std.os.linux.S.IFDIR) {
        return CommandFailure.plain("io_error", "input is a directory", 3);
    }
    return null;
}

/// The CLI names standard input `-`.
fn namesStdin(label: []const u8) bool {
    return std.mem.eql(u8, label, "-");
}

/// Standard input can be read only once, so it may be named only once.
fn stdinRepeats(inputs: []const []const u8) bool {
    var seen = false;
    for (inputs) |label| {
        if (!namesStdin(label)) continue;
        if (seen) return true;
        seen = true;
    }
    return false;
}

fn openRecordFile(io: std.Io, label: []const u8) std.Io.File.OpenError!std.Io.File {
    if (namesStdin(label)) return .stdin();
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
        if (namesStdin(label)) {
            // Opening a path could reuse descriptor 0 if stdin is closed.
            _ = fileIdentity(io, .stdin()) catch
                return pairCommandFailure(@intCast(index), IO_FAILURE);
        }
    }
    var transferred = false;
    const owns1 = !namesStdin(inputs[0]);
    const file1 = openRecordFile(io, inputs[0]) catch |err|
        return pairCommandFailure(0, inputOpenFailure(err));
    defer if (!transferred and owns1) file1.close(io);

    const owns2 = !namesStdin(inputs[1]);
    var file2 = openNonblockingRecordFile(inputs[1]) catch |err|
        return pairCommandFailure(1, inputOpenFailure(err));
    defer if (!transferred and owns2) file2.close(io);

    if (inputTypeFailure(io, file1)) |failure| return pairCommandFailure(0, failure);
    if (inputTypeFailure(io, file2)) |failure| return pairCommandFailure(1, failure);
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
    if (namesStdin(label)) return .stdin();
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

fn fileIdentity(io: std.Io, file: std.Io.File) (std.Io.Cancelable || error{ Io, BadFileDescriptor })!FileIdentity {
    const stat = try statRecordFile(io, file);
    return .of(stat);
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
    return .of(stat);
}

/// `stdoutFileIdentity`, with a failure to inspect stdout already reported (exit 3).
fn checkedStdoutIdentity(io: std.Io) error{Reported}!?FileIdentity {
    return stdoutFileIdentity(io) catch {
        if (canReportInputFailure(io, null)) printOutputFailure(io);
        return error.Reported;
    };
}

/// The identity of the file named `path`, without opening it; null when it
/// cannot be inspected, which the command reports when it reaches the input.
fn pathIdentity(io: std.Io, path: []const u8) ?FileIdentity {
    const linux = std.os.linux;
    const path_z = std.posix.toPosixPath(path) catch return null;
    var stat: linux.Statx = undefined;
    while (true) switch (linux.errno(linux.statx(linux.AT.FDCWD, &path_z, 0, .{ .INO = true, .TYPE = true }, &stat))) {
        .SUCCESS => break,
        .INTR => io.checkCancel() catch return null,
        else => return null,
    };
    if (!stat.mask.INO or !stat.mask.TYPE) return null;
    return .of(stat);
}

/// Rejects, before any input is read or output written, an input that is
/// stdout's file, and standard input named more than once under any name
/// while it is a pipe or FIFO, which can be read only once. A regular file on
/// standard input may also be named by path: each name reads it from the
/// start.
fn inputAliasFailure(io: std.Io, inputs: []const []const u8, output_identity: ?FileIdentity) ?u8 {
    const stdin_stat = statRecordFile(io, .stdin()) catch null;
    const stdin_identity: ?FileIdentity = if (stdin_stat) |stat| .of(stat) else null;
    const stdin_is_stream = if (stdin_stat) |stat| stat.mode & std.os.linux.S.IFMT != std.os.linux.S.IFREG else false;
    var stdin_uses: usize = 0;
    for (inputs) |label| {
        const identity = if (namesStdin(label)) stdin_identity else pathIdentity(io, label);
        const reads_stdin = namesStdin(label) or
            (stdin_is_stream and identity != null and std.meta.eql(identity.?, stdin_identity.?));
        if (reads_stdin) {
            stdin_uses += 1;
            if (stdin_uses > 1) {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: standard input may appear at most once\n",
                ) catch {};
                return 2;
            }
        }
        const output = output_identity orelse continue;
        if (identity) |input| if (std.meta.eql(input, output)) {
            printCommandFailure(io, label, outputAliasFailure(io, &.{input}));
            return 2;
        };
    }
    return null;
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
    const stderr_identity: FileIdentity = .of(stat);
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
    var failure = INSPECT_FAILURE;
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

// --- Record output ---

fn recordHasUnwritableEnding(record: zfastq.Record, canonical_span: ?[]const u8) bool {
    return canonical_span == null and fastq.recordHasTerminalCr(record);
}

fn writeCheckedRecord(
    writer: *zfastq.Writer,
    record: zfastq.Record,
    canonical_span: ?[]const u8,
) error{WriteFailed}!void {
    if (canonical_span) |span| {
        try fastq.writeCanonicalRecordSpan(writer, span);
    } else {
        try fastq.writeRecordUnchecked(writer, record);
    }
}

/// A span is writable unless a field still ends in CR once its line ending
/// is removed, which canonical LF output cannot represent.
fn spanIsWritable(span: fastq.RecordSpan) bool {
    return span.canonical or !fastq.recordHasTerminalCr(span.record());
}

fn writeRecordSpan(writer: *zfastq.Writer, span: fastq.RecordSpan) error{WriteFailed}!void {
    if (span.canonical) return fastq.writeCanonicalRecordSpan(writer, span.bytes);
    return fastq.writeRecordUnchecked(writer, span.record());
}

fn headerNameMismatch(
    header1: []const u8,
    header2: []const u8,
    policy: pairing.NamePolicy,
    pair_index: u64,
    record_indexes: [2]u64,
    header_offsets: [2]u64,
) PairCommandFailure {
    return .{ .pair = .{ .name_mismatch = .{
        .pair_index = pair_index,
        .records = .{
            .init(pairing.parseName(header1, policy), record_indexes[0], header_offsets[0]),
            .init(pairing.parseName(header2, policy), record_indexes[1], header_offsets[1]),
        },
    } } };
}

/// Four lines at the line limit, each with its LF: the largest canonical
/// record. It saturates, so every command accepts the same line limits; no
/// record reaches the size of the address space.
fn recordStagingLimit(max_line_bytes: usize) usize {
    const line = std.math.add(usize, max_line_bytes, 1) catch return std.math.maxInt(usize);
    return std.math.mul(usize, line, 4) catch std.math.maxInt(usize);
}

fn canonicalRecordSize(
    record: zfastq.Record,
    canonical_span: ?[]const u8,
) error{RecordStagingLimit}!usize {
    if (canonical_span) |span| return span.len;
    return fastq.canonicalRecordLen(record) orelse error.RecordStagingLimit;
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
    _ = fastq.copyCanonicalRecord(staging.addManyAsSliceAssumeCapacity(required), record);
}

// --- Pair names and retained mates ---

const PairMode = enum {
    none,
    paired,
    interleaved,
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

fn storeNormalizedId(
    allocator: std.mem.Allocator,
    storage: *std.ArrayList(u8),
    normalized_id: []const u8,
) std.mem.Allocator.Error!void {
    storage.clearRetainingCapacity();
    try storage.ensureTotalCapacityPrecise(allocator, normalized_id.len);
    storage.appendSliceAssumeCapacity(normalized_id);
}

fn lastRecordIndex(reader: *const zfastq.Reader) ?u64 {
    return if (reader.recordIndex() == 0) null else reader.recordIndex() - 1;
}

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

const InterleavedFirstRecordStorage = enum {
    reader,
    retained,
    staged,
};

const RefillSpanningMate = struct {
    record: ?fastq.ValidatedRecord,
    first_storage: InterleavedFirstRecordStorage,
};

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
        const transferred_record2 = try fastq.nextValidatedAfterFallbackTransfer(
            reader,
            validator,
        );
        if (transferred_record2) |record2| return .{
            .record = record2,
            .first_storage = .retained,
        };
        // Reader's four line limits cover CLI staging limits, including markers.
        // Internal callers can supply a smaller limit and still need its error.
        if (staging_limit < recordStagingLimit(reader.options.max_line_bytes)) {
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
        .reader, .retained => try writeCheckedRecord(writer, record, canonical_span),
        .staged => try fastq.writeCanonicalRecordSpan(writer, staged),
    }
    if (storage == .retained) {
        fastq.restoreFallbackRecordStorage(reader, retained);
    }
}

// --- Count command ---

const CountOutcome = union(enum) {
    success: u64,
    failure: CommandFailure,
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
                printCommandFailure(io, input, failure);
                exit_code = @max(exit_code, failure.exit_code);
                continue;
            },
        };
        printCount(io, count) catch return stdoutFailed(io, exit_code);
    }
    return exit_code;
}

fn countInput(io: std.Io, label: []const u8, options: InputOptions) CountOutcome {
    var input: RecordInput = undefined;
    if (initRecordInput(&input, io, label, null)) |failure| return .{ .failure = failure };
    defer input.deinit(io);

    var scanner = zfastq.count_scan.Scanner.init(.{ .max_line_bytes = options.max_line_bytes });
    while (input.readScannerChunk() catch return .{ .failure = IO_FAILURE }) |decoded| {
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

// --- Stats command ---

const StatsOutcome = union(enum) {
    success: zfastq.Stats,
    failure: CommandFailure,
};

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

        printStats(io, input, &stats, printed_block) catch
            return stdoutFailed(io, exit_code);
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
    var stdout_buffer: [JSON_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/stats-v1") catch return stdoutFailed(io, 0);
    var exit_code: u8 = 0;
    for (inputs) |input| {
        const outcome = statsInput(io, allocator, input, options);
        switch (outcome) {
            .success => {},
            .failure => |failure| exit_code = @max(exit_code, failure.exit_code),
        }
        writeStatsJsonResult(&json, input, outcome) catch return stdoutFailed(io, exit_code);
    }
    finishJsonDocument(&json) catch return stdoutFailed(io, exit_code);
    stdout_writer.interface.flush() catch return stdoutFailed(io, exit_code);
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
    return collectStats(allocator, input.stream(), options);
}

fn collectStats(
    allocator: std.mem.Allocator,
    source: *std.Io.Reader,
    options: InputOptions,
) StatsOutcome {
    var reader = fastq.initStreamReader(allocator, source, .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();

    var stats: zfastq.Stats = .{};
    // A record larger than what is buffered streams its totals here.
    var streamed: statistics.StreamedTotals = .{};
    while (true) {
        const next = fastq.nextPayloadOrStreamed(&reader, statistics.StreamedTotals, &streamed) catch |err| {
            return .{ .failure = mapReaderFailure(&reader, err) };
        } orelse break;
        switch (next) {
            .streamed => {
                stats.addStreamed(&streamed) catch |err| {
                    return .{ .failure = statsFailure(&stats, &reader, err) };
                };
                streamed = .{};
            },
            .payload => |predicted| {
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
                        return .{ .failure = statsFailure(&stats, &reader, err) };
                    };
                    break;
                }
            },
        }
    }
    return .{ .success = stats };
}

fn statsFailure(stats: *zfastq.Stats, reader: *zfastq.Reader, err: statistics.StatsError) CommandFailure {
    switch (err) {
        error.S006InvalidQuality => {
            const quality_error = stats.takeLastQualityError() orelse {
                return CommandFailure.plain(
                    "io_error",
                    "quality validation failed without details",
                    3,
                );
            };
            const offsets = reader.currentRecordOffsets() orelse {
                return RECORD_LOCATION_UNAVAILABLE;
            };
            const details = fastq.semanticParseError(
                fastq.semanticQualityError(quality_error.byte_index),
                reader.recordIndex() - 1,
                offsets.quality,
            ) catch return INPUT_LOCATION_LIMIT;
            return CommandFailure.lint(details);
        },
        // The Reader rejects unequal lengths as S005 before a record reaches
        // the statistics; a panic keeps a broken invariant a defined crash.
        error.S005LengthMismatch => @panic("Reader returned unequal sequence and quality lengths"),
        error.Overflow => return STATISTICS_LIMIT,
    }
}

// --- Check command ---

const CheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
};

const PairedCheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
    pair_mode: PairMode = .none,
    pair_name_policy: pairing.NamePolicy = .illumina,
    /// Stdout's file, which no input may be (JSON output would land in it).
    output_identity: ?FileIdentity = null,
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
    var checked_options = options;
    checked_options.output_identity = checkedStdoutIdentity(io) catch return 3;
    if (json_output) return runPairedCheckJson(io, allocator, inputs, checked_options);

    const failure = checkPairMode(io, allocator, inputs, checked_options);
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
    if (pair_mode == .paired and stdinRepeats(inputs)) {
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
    var stdout_buffer: [JSON_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/check-v1") catch return stdoutFailed(io, 0);
    var exit_code: u8 = 0;
    for (inputs) |input| {
        const failure = checkInput(io, input, options);
        if (failure) |details| exit_code = @max(exit_code, details.exit_code);
        writeCheckJsonResult(&json, input, failure) catch return stdoutFailed(io, exit_code);
    }
    finishJsonDocument(&json) catch return stdoutFailed(io, exit_code);
    stdout_writer.interface.flush() catch return stdoutFailed(io, exit_code);
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
    var stdout_buffer: [JSON_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/check-v1") catch return stdoutFailed(io, 0);
    const exit_code = if (failure) |details| details.exitCode() else 0;
    writePairedCheckJsonResult(&json, inputs, options.pair_mode, failure) catch
        return stdoutFailed(io, exit_code);
    finishJsonDocument(&json) catch return stdoutFailed(io, exit_code);
    stdout_writer.interface.flush() catch return stdoutFailed(io, exit_code);
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
    if (initPairedRecordInputs(&input1, &input2, io, inputs, options.output_identity)) |failure| return failure;
    defer input1.deinit(io);
    defer input2.deinit(io);

    return checkPairedSources(
        allocator,
        input1.stream(),
        input2.stream(),
        options,
        null,
    );
}

const PairStep = union(enum) {
    /// Both inputs ended together.
    end,
    pair,
    failure: PairCommandFailure,
};

/// The verdict on the mates just read from two files (null at the end of an
/// input): count mismatch, then mate 1's and mate 2's semantic failures, then
/// their names.
inline fn verifyPair(
    reader1: *const fastq.Reader,
    reader2: *const fastq.Reader,
    mate1: ?fastq.ValidatedHeader,
    mate2: ?fastq.ValidatedHeader,
    policy: pairing.NamePolicy,
) PairStep {
    if (mate1 == null and mate2 == null) return .end;
    if (mate1 == null or mate2 == null) {
        const remaining_side: u1 = if (mate1 != null) 0 else 1;
        const pair_index = if (mate1 != null)
            reader1.recordIndex() - 1
        else
            reader2.recordIndex() - 1;
        return .{ .failure = .{ .pair = .{ .count_mismatch = .{
            .pair_index = pair_index,
            .remaining_side = remaining_side,
            .record_indexes = .{
                lastRecordIndex(reader1),
                lastRecordIndex(reader2),
            },
        } } } };
    }

    const record_index1 = reader1.recordIndex() - 1;
    const record_index2 = reader2.recordIndex() - 1;
    const offsets1 = reader1.currentRecordOffsets().?;
    const offsets2 = reader2.currentRecordOffsets().?;
    if (mapSemanticFailure(mate1.?.semantic_error, offsets1, record_index1)) |failure| {
        return .{ .failure = pairCommandFailure(0, failure) };
    }
    if (mapSemanticFailure(mate2.?.semantic_error, offsets2, record_index2)) |failure| {
        return .{ .failure = pairCommandFailure(1, failure) };
    }
    if (!pairing.headersMatch(mate1.?.header, mate2.?.header, policy)) {
        return .{ .failure = headerNameMismatch(
            mate1.?.header,
            mate2.?.header,
            policy,
            record_index1,
            .{ record_index1, record_index2 },
            .{ offsets1.header, offsets2.header },
        ) };
    }
    return .pair;
}

fn checkPairedSources(
    allocator: std.mem.Allocator,
    source1: *std.Io.Reader,
    source2: *std.Io.Reader,
    options: PairedCheckOptions,
    exact_selector: ?*sampling.ExactSelector,
) ?PairCommandFailure {
    const reader_options: fastq.ReaderOptions = .{ .max_line_bytes = options.max_line_bytes };
    var reader1 = fastq.initStreamReader(allocator, source1, reader_options);
    defer reader1.deinit();
    var reader2 = fastq.initStreamReader(allocator, source2, reader_options);
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
        switch (verifyPair(&reader1, &reader2, record1, record2, options.pair_name_policy)) {
            .end => return null,
            .failure => |failure| return failure,
            .pair => {},
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
    if (initRecordInput(&input, io, input_label, options.output_identity)) |failure| {
        return pairCommandFailure(0, failure);
    }
    defer input.deinit(io);

    return checkInterleavedSource(
        allocator,
        input.stream(),
        options,
        null,
    );
}

fn checkInterleavedSource(
    allocator: std.mem.Allocator,
    source: *std.Io.Reader,
    options: PairedCheckOptions,
    exact_selector: ?*sampling.ExactSelector,
) ?PairCommandFailure {
    var reader = fastq.initStreamReader(allocator, source, .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();
    var stored_name: StoredPairName = .{};
    defer stored_name.deinit(allocator);
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        switch (verifyInterleavedPair(allocator, &reader, &validator, &stored_name, options.pair_name_policy, null)) {
            .end => return null,
            .failure => |failure| return failure,
            .pair => {},
        }
        if (exact_selector) |selector| {
            selector.considerRecord(allocator, reader.recordIndex() / 2) catch |err| {
                return exactPairSelectionFailure(err);
            };
        }
    }
}

/// Reads and verifies one interleaved pair without keeping either record:
/// mate 1's semantic failure (after mate 2 proves it exists), mate 2's, then
/// their names. `first` is mate 1 when the caller has already read it.
inline fn verifyInterleavedPair(
    allocator: std.mem.Allocator,
    reader: *fastq.Reader,
    validator: *fastq.AdaptiveRecordValidator,
    stored_name: *StoredPairName,
    policy: pairing.NamePolicy,
    first: ?fastq.ValidatedHeader,
) PairStep {
    var record1 = first orelse (fastq.nextValidatedHeader(reader, validator) catch |err| {
        return .{ .failure = pairCommandFailure(0, mapReaderFailure(reader, err)) };
    }) orelse return .end;
    const record_index1 = reader.recordIndex() - 1;
    const offsets1 = reader.currentRecordOffsets().?;
    const semantic1 = mapSemanticFailure(
        record1.semantic_error,
        offsets1,
        record_index1,
    );

    const paired_record2 = (if (semantic1 == null)
        fastq.nextValidatedMatePreservingHeader(reader, &record1.header, validator)
    else
        fastq.nextBufferedValidatedRecord(reader, validator)) catch |err| {
        return .{ .failure = pairCommandFailure(0, mapReaderFailure(reader, err)) };
    };
    const both_headers_borrowed = paired_record2 != null;
    var semantic2: ?zfastq.SemanticError = null;
    const header2 = if (paired_record2) |validated2| header: {
        if (semantic1 == null) semantic2 = validated2.semantic_error;
        break :header validated2.record.header;
    } else header: {
        if (semantic1) |failure| {
            // Mate 2 must exist, and its structural error comes first. The
            // call stays out of line: this path ends the command.
            const exists = @call(.never_inline, fastq.Reader.advance, .{reader}) catch |err| {
                return .{ .failure = pairCommandFailure(0, mapReaderFailure(reader, err)) };
            };
            if (!exists) return .{ .failure = missingInterleavedMateFailure(record_index1) };
            return .{ .failure = pairCommandFailure(0, failure) };
        }
        const name1 = pairing.parseName(record1.header, policy);
        stored_name.store(allocator, name1) catch return .{ .failure = pairCommandFailure(0, OUT_OF_MEMORY) };
        const record2 = fastq.nextValidatedHeader(reader, validator) catch |err| {
            return .{ .failure = pairCommandFailure(0, mapReaderFailure(reader, err)) };
        } orelse return .{ .failure = missingInterleavedMateFailure(record_index1) };
        semantic2 = record2.semantic_error;
        break :header record2.header;
    };

    if (semantic1) |failure| {
        return .{ .failure = pairCommandFailure(0, failure) };
    }

    const record_index2 = reader.recordIndex() - 1;
    const offsets2 = reader.currentRecordOffsets().?;
    if (mapSemanticFailure(
        semantic2,
        offsets2,
        record_index2,
    )) |failure| {
        return .{ .failure = pairCommandFailure(0, failure) };
    }

    const names_match = if (both_headers_borrowed)
        pairing.headersMatch(record1.header, header2, policy)
    else blk: {
        const name1 = stored_name.name();
        const name2 = pairing.parseName(header2, policy);
        break :blk pairing.namesMatch(name1, name2);
    };
    if (!names_match) {
        const name2 = pairing.parseName(header2, policy);
        return .{ .failure = .{ .pair = .{ .name_mismatch = .{
            .pair_index = interleavedPairIndex(record_index1),
            .records = .{
                if (both_headers_borrowed)
                    .init(
                        pairing.parseName(record1.header, policy),
                        record_index1,
                        offsets1.header,
                    )
                else
                    stored_name.diagnostic(record_index1, offsets1.header),
                .init(name2, record_index2, offsets2.header),
            },
        } } } };
    }
    return .pair;
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
    while (true) {
        const chunk = input.readScannerChunk() catch
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

// --- Sample command ---

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

/// One unit (a record or a pair) of an exact-count second pass; a failed
/// unit leaves its failure in the pass.
const ExactUnit = enum { done, end, failed };

/// The second pass of exact-count sampling: units in input order, selected
/// ones written and the rest skipped; then no unit may follow and the inputs
/// must be unchanged. `pass` reads its units with `write` and `skip`, checks
/// for one more with `extra`, compares snapshots with `changed`, and names
/// the failure for a changed input with `inputChanged`; failures go to
/// `pass.failure`.
inline fn runExactSecondPass(
    comptime Failure: type,
    cursor: *ExactOutputCursor,
    expected_count: u64,
    pass: anytype,
) error{WriteFailed}!?Failure {
    var next_selected = cursor.nextSelected();
    while (cursor.unit_count < expected_count) {
        const unit = if (cursor.selectNextCached(&next_selected)) try pass.write() else pass.skip();
        if (unit != .done) break;
        cursor.completeUnit();
    }
    if (pass.failure.* == null and cursor.unit_count == expected_count) pass.extra();
    if (pass.changed()) |changed| return changed;
    if (pass.failure.*) |details| return details;
    if (cursor.unit_count != expected_count or !cursor.selectionComplete()) return pass.inputChanged();
    return null;
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

fn exactPairSelectionFailure(err: sampling.ReservoirError) PairCommandFailure {
    return pairCommandFailure(0, exactSelectionFailure(err));
}

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
    if (mode == .count and namesStdin(input)) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: exact-count sampling requires a file path\n",
        ) catch {};
        return 2;
    }

    options.output_identity = checkedStdoutIdentity(io) catch return 3;
    var stdout_buffer: [RECORD_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var stdout_sink = io_layer.FileSink.initStreaming(io, .stdout(), &stdout_buffer);
    var writer = zfastq.Writer.init(stdout_sink.byteSink());
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
        return stdoutFailed(io, 0);
    };

    // Messages follow detection order: a record's failure before the final flush's.
    const flushed = if (writer.flush()) true else |_| false;
    if (failure) |details| printCommandFailure(io, input, details);
    const output_exit: u8 = if (flushed) 0 else stdoutFailed(io, 0);
    return if (failure) |details| @max(output_exit, details.exit_code) else output_exit;
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
                if (!namesStdin(input)) continue;
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: paired exact-count sampling requires file paths\n",
                ) catch {};
                return 2;
            }
        },
        .fraction => if (options.pair_mode == .paired and stdinRepeats(inputs)) {
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
    const staging_limit = if (needs_staging) recordStagingLimit(options.max_line_bytes) else 0;
    options.output_identity = checkedStdoutIdentity(io) catch return 3;
    var stdout_buffer: [RECORD_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var stdout_sink = io_layer.FileSink.initStreaming(io, .stdout(), &stdout_buffer);
    var writer = zfastq.Writer.init(stdout_sink.byteSink());
    const failure = (switch (mode) {
        .fraction => |fraction| blk: {
            var selector = sampling.Selector.init(fraction, options.seed);
            var output_selector: PairOutputSelector = .{ .fraction = &selector };
            break :blk switch (options.pair_mode) {
                .none => unreachable,
                .paired => interleaveInputs(io, allocator, inputs, &writer, .{
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
        return stdoutFailed(io, 0);
    };

    // Messages follow detection order: a pair's failure before the final flush's.
    const flushed = if (writer.flush()) true else |_| false;
    if (failure) |details| printPairCommandFailure(io, inputs, options.pair_mode, details);
    const output_exit: u8 = if (flushed) 0 else stdoutFailed(io, 0);
    return if (failure) |details| @max(output_exit, details.exitCode()) else output_exit;
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
        input.stream(),
        writer,
        selection,
        staging_limit,
        options,
    );
}

fn sampleInterleavedSource(
    allocator: std.mem.Allocator,
    source: *std.Io.Reader,
    writer: *zfastq.Writer,
    selection: *PairOutputSelector,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    var reader = fastq.initStreamReader(allocator, source, .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();
    var stored_name: StoredPairName = .{};
    defer stored_name.deinit(allocator);
    var mate1_storage: PreservedMate1 = .{};
    defer mate1_storage.deinit(allocator);
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        var span1: ?fastq.RecordSpan = null;
        if (splitSpanPairs(
            &reader,
            &validator,
            options.pair_name_policy,
            selection,
            .{ writer, writer },
            &span1,
        ) catch return error.WriteFailed) |failure| return failure;
        // Drawn before the pair is read: an unselected pair is validated
        // without being kept.
        if (!selection.selectPair()) {
            const first: ?fastq.ValidatedHeader = if (span1) |mate1| .{ .header = mate1.header(), .semantic_error = null } else null;
            switch (verifyInterleavedPair(allocator, &reader, &validator, &stored_name, options.pair_name_policy, first)) {
                .end => return null,
                .failure => |failure| return failure,
                .pair => continue,
            }
        }
        const pair = switch (readInterleavedPair(
            allocator,
            &reader,
            &validator,
            &stored_name,
            &mate1_storage,
            staging_limit,
            options.pair_name_policy,
            span1,
        )) {
            .end => return null,
            .failure => |failure| return failure,
            .pair => |pair| pair,
        };
        try mate1_storage.write(writer, &reader, pair);
        try writeCheckedRecord(writer, pair.record2, pair.canonical_span2);
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
    if (selector.* == .none) {
        return checkRecordInput(&input, .{
            .max_line_bytes = options.max_line_bytes,
            .alphabet = options.alphabet,
        });
    }
    return sampleFractionSource(allocator, input.stream(), writer, selector, options);
}

fn sampleFractionSource(
    allocator: std.mem.Allocator,
    source: *std.Io.Reader,
    writer: *zfastq.Writer,
    selector: *sampling.Selector,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    var reader = fastq.initStreamReader(allocator, source, .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        // One draw per record, made before reading it: an unselected record
        // is validated without being copied.
        const selected = selector.selectRecord();
        if (fastq.nextRecordSpan(&reader, &validator)) |span| {
            if (!selected) continue;
            if (!spanIsWritable(span)) return unwritableRecordFailure();
            try writeRecordSpan(writer, span);
            continue;
        }
        if (!selected) {
            const header = fastq.nextValidatedHeader(&reader, &validator) catch |err| {
                return mapReaderFailure(&reader, err);
            } orelse return null;
            if (mapCurrentSemanticFailure(
                header.semantic_error,
                &reader,
                reader.recordIndex() - 1,
            )) |failure| return failure;
            continue;
        }
        var next = fastq.nextValidatedRecord(&reader, &validator) catch |err| {
            return mapReaderFailure(&reader, err);
        };
        // A pointer capture reads the record where the reader built it.
        const validated = if (next) |*validated| validated else return null;
        if (mapCurrentSemanticFailure(
            validated.semantic_error,
            &reader,
            reader.recordIndex() - 1,
        )) |failure| return failure;
        if (recordHasUnwritableEnding(validated.record, validated.canonical_span)) {
            return unwritableRecordFailure();
        }
        try writeCheckedRecord(writer, validated.record, validated.canonical_span);
    }
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
    // Both passes reuse this storage, so their read buffers touch the same stack pages.
    var input: RecordInput = undefined;

    const first = sampleExactFirstPass(io, allocator, &input, path, &selector, options);
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
            &input,
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
            &input,
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
    input: *RecordInput,
    path: []const u8,
    selector: *sampling.ExactSelector,
    options: SampleOptions,
) ExactFirstPass {
    const initial = switch (initExactInput(input, io, path, null, options.output_identity)) {
        .failure => |failure| return .{ .failure = failure },
        .success => |snapshot| snapshot,
    };
    defer input.deinit(io);

    var failure: ?CommandFailure = null;
    var scanner = fastq.CheckScanner.init(
        .{ .max_line_bytes = options.max_line_bytes },
        .{ .alphabet = options.alphabet },
    );
    while (true) {
        const chunk = input.readScannerChunk() catch {
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

    if (exactInputSnapshotFailure(input, io, initial)) |snapshot_failure| {
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
    input: *RecordInput,
    path: []const u8,
    writer: *zfastq.Writer,
    selected: sampling.ExactIndexes,
    expected_snapshot: FileSnapshot,
    expected_count: u64,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    switch (initExactInput(input, io, path, expected_snapshot, options.output_identity)) {
        .failure => |failure| return failure,
        .success => {},
    }
    defer input.deinit(io);

    var reader = fastq.initStreamReader(allocator, input.stream(), .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    var cursor = ExactOutputCursor{
        .select_all = select_all,
        .indexes = selected,
    };
    var failure: ?CommandFailure = null;
    const Pass = struct {
        reader: *zfastq.Reader,
        validator: *fastq.AdaptiveRecordValidator,
        writer: *zfastq.Writer,
        input: *RecordInput,
        io: std.Io,
        snapshot: FileSnapshot,
        failure: *?CommandFailure,

        fn fail(self: @This(), details: CommandFailure) ExactUnit {
            self.failure.* = details;
            return .failed;
        }

        inline fn write(self: @This()) error{WriteFailed}!ExactUnit {
            const validated = fastq.nextValidatedRecord(self.reader, self.validator) catch |err| {
                return self.fail(mapReaderFailure(self.reader, err));
            } orelse return .end;
            if (validated.semantic_error != null) return self.fail(inputChangedFailure());
            if (recordHasUnwritableEnding(validated.record, validated.canonical_span)) {
                return self.fail(unwritableRecordFailure());
            }
            try writeCheckedRecord(self.writer, validated.record, validated.canonical_span);
            return .done;
        }

        inline fn skip(self: @This()) ExactUnit {
            const advanced = self.reader.advance() catch |err| return self.fail(mapReaderFailure(self.reader, err));
            return if (advanced) .done else .end;
        }

        fn extra(self: @This()) void {
            const advanced = self.reader.advance() catch |err| {
                _ = self.fail(mapReaderFailure(self.reader, err));
                return;
            };
            if (advanced) _ = self.fail(inputChangedFailure());
        }

        fn changed(self: @This()) ?CommandFailure {
            return exactInputSnapshotFailure(self.input, self.io, self.snapshot);
        }

        fn inputChanged(_: @This()) CommandFailure {
            return inputChangedFailure();
        }
    };
    return runExactSecondPass(CommandFailure, &cursor, expected_count, Pass{
        .reader = &reader,
        .validator = &validator,
        .writer = writer,
        .input = input,
        .io = io,
        .snapshot = expected_snapshot,
        .failure = &failure,
    });
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
    // Both passes reuse this storage, so their read buffers touch the same stack pages.
    var record_inputs: [2]RecordInput = undefined;

    const first = sampleExactPairFirstPass(
        io,
        allocator,
        &record_inputs,
        inputs,
        &selector,
        options,
    );
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
            &record_inputs,
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
            &record_inputs,
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
    record_inputs: *[2]RecordInput,
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
            record_inputs,
            inputs,
            selector,
            pair_options,
            options.output_identity,
        ),
        .interleaved => sampleExactInterleavedFirstPass(
            io,
            allocator,
            &record_inputs[0],
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
    record_inputs: *[2]RecordInput,
    inputs: []const []const u8,
    selector: *sampling.ExactSelector,
    options: PairedCheckOptions,
    output_identity: ?FileIdentity,
) ExactPairFirstPass {
    const input1 = &record_inputs[0];
    const input2 = &record_inputs[1];
    const snapshots = switch (initExactPairedInputs(
        input1,
        input2,
        io,
        inputs,
        null,
        output_identity,
    )) {
        .failure => |failure| return .{ .failure = failure },
        .success => |snapshots| snapshots,
    };
    defer input1.deinit(io);
    defer input2.deinit(io);

    const failure = checkPairedSources(
        allocator,
        input1.stream(),
        input2.stream(),
        options,
        selector,
    );
    if (exactPairSnapshotFailure(input1, io, snapshots[0], 0)) |changed| {
        return .{ .failure = changed };
    }
    if (exactPairSnapshotFailure(input2, io, snapshots[1], 1)) |changed| {
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
    input: *RecordInput,
    path: []const u8,
    selector: *sampling.ExactSelector,
    options: PairedCheckOptions,
    output_identity: ?FileIdentity,
) ExactPairFirstPass {
    const snapshot = switch (initExactInput(input, io, path, null, output_identity)) {
        .failure => |failure| return .{ .failure = pairCommandFailure(0, failure) },
        .success => |captured| captured,
    };
    defer input.deinit(io);

    const failure = checkInterleavedSource(
        allocator,
        input.stream(),
        options,
        selector,
    );
    if (exactPairSnapshotFailure(input, io, snapshot, 0)) |changed| {
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
    record_inputs: *[2]RecordInput,
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
            record_inputs,
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
            &record_inputs[0],
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
    record_inputs: *[2]RecordInput,
    inputs: []const []const u8,
    writer: *zfastq.Writer,
    indexes: sampling.ExactIndexes,
    snapshots: [2]FileSnapshot,
    expected_count: u64,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    const input1 = &record_inputs[0];
    const input2 = &record_inputs[1];
    switch (initExactPairedInputs(input1, input2, io, inputs, snapshots, options.output_identity)) {
        .failure => |failure| return failure,
        .success => {},
    }
    defer input1.deinit(io);
    defer input2.deinit(io);

    var reader1 = fastq.initStreamReader(allocator, input1.stream(), .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader1.deinit();
    var reader2 = fastq.initStreamReader(allocator, input2.stream(), .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader2.deinit();
    var validator1 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var validator2 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    var cursor = ExactOutputCursor{
        .select_all = select_all,
        .indexes = indexes,
    };
    var failure: ?PairCommandFailure = null;
    const Pass = struct {
        readers: [2]*zfastq.Reader,
        validators: [2]*fastq.AdaptiveRecordValidator,
        writer: *zfastq.Writer,
        inputs: [2]*RecordInput,
        io: std.Io,
        snapshots: [2]FileSnapshot,
        policy: pairing.NamePolicy,
        failure: *?PairCommandFailure,

        fn fail(self: @This(), details: PairCommandFailure) ExactUnit {
            self.failure.* = details;
            return .failed;
        }

        fn write(self: @This()) error{WriteFailed}!ExactUnit {
            var mates: [2]?fastq.ValidatedRecord = undefined;
            for (&mates, 0..) |*mate, side| {
                mate.* = fastq.nextValidatedRecord(self.readers[side], self.validators[side]) catch |err| {
                    return self.fail(pairCommandFailure(@intCast(side), mapReaderFailure(self.readers[side], err)));
                };
            }
            if (mates[0] == null or mates[1] == null) return self.fail(missingMate(mates[0] != null, mates[1] != null));
            const mate1 = mates[0].?;
            const mate2 = mates[1].?;
            if (mate1.semantic_error != null) return self.fail(inputChangedPairFailure(0));
            if (mate2.semantic_error != null) return self.fail(inputChangedPairFailure(1));
            if (!pairing.headersMatch(mate1.record.header, mate2.record.header, self.policy)) {
                return self.fail(inputChangedPairFailure(0));
            }
            if (recordHasUnwritableEnding(mate1.record, mate1.canonical_span)) {
                return self.fail(pairCommandFailure(0, unwritableRecordFailure()));
            }
            if (recordHasUnwritableEnding(mate2.record, mate2.canonical_span)) {
                return self.fail(pairCommandFailure(1, unwritableRecordFailure()));
            }
            try writeCheckedRecord(self.writer, mate1.record, mate1.canonical_span);
            try writeCheckedRecord(self.writer, mate2.record, mate2.canonical_span);
            return .done;
        }

        fn skip(self: @This()) ExactUnit {
            var got: [2]bool = undefined;
            for (&got, 0..) |*mate, side| {
                mate.* = self.readers[side].advance() catch |err| {
                    return self.fail(pairCommandFailure(@intCast(side), mapReaderFailure(self.readers[side], err)));
                };
            }
            if (!got[0] or !got[1]) return self.fail(missingMate(got[0], got[1]));
            return .done;
        }

        /// Blames the input that still has a record, or R1 when neither has.
        fn missingMate(got1: bool, got2: bool) PairCommandFailure {
            return inputChangedPairFailure(if (got1 != got2 and got2) 1 else 0);
        }

        fn extra(self: @This()) void {
            var extra_mates: [2]bool = undefined;
            for (&extra_mates, 0..) |*mate, side| {
                mate.* = self.readers[side].advance() catch |err| {
                    _ = self.fail(pairCommandFailure(@intCast(side), mapReaderFailure(self.readers[side], err)));
                    return;
                };
            }
            if (extra_mates[0] or extra_mates[1]) {
                _ = self.fail(inputChangedPairFailure(if (extra_mates[0] != extra_mates[1] and extra_mates[1]) 1 else 0));
            }
        }

        fn changed(self: @This()) ?PairCommandFailure {
            for (self.inputs, self.snapshots, 0..) |input, snapshot, side| {
                if (exactPairSnapshotFailure(input, self.io, snapshot, @intCast(side))) |changed_input| return changed_input;
            }
            return null;
        }

        fn inputChanged(_: @This()) PairCommandFailure {
            return inputChangedPairFailure(0);
        }
    };
    return runExactSecondPass(PairCommandFailure, &cursor, expected_count, Pass{
        .readers = .{ &reader1, &reader2 },
        .validators = .{ &validator1, &validator2 },
        .writer = writer,
        .inputs = .{ input1, input2 },
        .io = io,
        .snapshots = snapshots,
        .policy = options.pair_name_policy,
        .failure = &failure,
    });
}

fn sampleExactInterleavedSecondPass(
    select_all: bool,
    io: std.Io,
    allocator: std.mem.Allocator,
    input: *RecordInput,
    path: []const u8,
    writer: *zfastq.Writer,
    indexes: sampling.ExactIndexes,
    snapshot: FileSnapshot,
    expected_count: u64,
    staging_limit: usize,
    options: SampleOptions,
) error{WriteFailed}!?PairCommandFailure {
    switch (initExactInput(input, io, path, snapshot, options.output_identity)) {
        .failure => |failure| return pairCommandFailure(0, failure),
        .success => {},
    }
    defer input.deinit(io);

    var reader = fastq.initStreamReader(allocator, input.stream(), .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var mate1_storage: PreservedMate1 = .{};
    defer mate1_storage.deinit(allocator);

    var cursor = ExactOutputCursor{
        .select_all = select_all,
        .indexes = indexes,
    };
    var failure: ?PairCommandFailure = null;
    const Pass = struct {
        allocator: std.mem.Allocator,
        reader: *zfastq.Reader,
        validator: *fastq.AdaptiveRecordValidator,
        mate1_storage: *PreservedMate1,
        staging_limit: usize,
        writer: *zfastq.Writer,
        input: *RecordInput,
        io: std.Io,
        snapshot: FileSnapshot,
        policy: pairing.NamePolicy,
        failure: *?PairCommandFailure,

        fn fail(self: @This(), details: PairCommandFailure) ExactUnit {
            self.failure.* = details;
            return .failed;
        }

        fn write(self: @This()) error{WriteFailed}!ExactUnit {
            var validated1 = fastq.nextValidatedRecord(self.reader, self.validator) catch |err| {
                return self.fail(pairCommandFailure(0, mapReaderFailure(self.reader, err)));
            } orelse return .end;
            const invalid1 = validated1.semantic_error != null;
            const header1_len = validated1.record.header.len;
            const unwritable1 = recordHasUnwritableEnding(validated1.record, validated1.canonical_span);
            var record1_storage: InterleavedFirstRecordStorage = .reader;
            const paired_record2 = fastq.nextValidatedMatePreservingRecord(self.reader, &validated1, self.validator) catch |err| {
                return self.fail(pairCommandFailure(0, mapReaderFailure(self.reader, err)));
            };
            const record1 = validated1.record;
            const canonical_span1 = validated1.canonical_span;
            const validated2 = paired_record2 orelse record: {
                const preserved = nextAfterPreservingInterleavedMate1(
                    self.allocator,
                    self.reader,
                    &self.mate1_storage.retained,
                    &self.mate1_storage.staged_record,
                    record1,
                    canonical_span1,
                    self.staging_limit,
                    self.validator,
                ) catch |err| return self.fail(mapPreservedMateFailure(self.reader, err));
                record1_storage = preserved.first_storage;
                break :record preserved.record orelse return .end;
            };
            const record2 = validated2.record;
            const canonical_span2 = validated2.canonical_span;
            if (invalid1 or validated2.semantic_error != null) return self.fail(inputChangedPairFailure(0));
            const header1 = switch (record1_storage) {
                .reader, .retained => record1.header,
                .staged => self.mate1_storage.staged_record.items[1 .. 1 + header1_len],
            };
            if (!pairing.headersMatch(header1, record2.header, self.policy)) return self.fail(inputChangedPairFailure(0));
            if (unwritable1 or recordHasUnwritableEnding(record2, canonical_span2)) {
                return self.fail(pairCommandFailure(0, unwritableRecordFailure()));
            }
            try self.mate1_storage.write(self.writer, self.reader, .{
                .record1 = record1,
                .canonical_span1 = canonical_span1,
                .record1_storage = record1_storage,
                .record2 = record2,
                .canonical_span2 = canonical_span2,
            });
            try writeCheckedRecord(self.writer, record2, canonical_span2);
            return .done;
        }

        fn skip(self: @This()) ExactUnit {
            for (0..2) |_| {
                const advanced = self.reader.advance() catch |err| {
                    return self.fail(pairCommandFailure(0, mapReaderFailure(self.reader, err)));
                };
                if (!advanced) return .end;
            }
            return .done;
        }

        fn extra(self: @This()) void {
            const extra_mate1 = self.reader.advance() catch |err| {
                _ = self.fail(pairCommandFailure(0, mapReaderFailure(self.reader, err)));
                return;
            };
            if (!extra_mate1) return;
            _ = self.reader.advance() catch |err| {
                _ = self.fail(pairCommandFailure(0, mapReaderFailure(self.reader, err)));
                return;
            };
            _ = self.fail(inputChangedPairFailure(0));
        }

        fn changed(self: @This()) ?PairCommandFailure {
            return exactPairSnapshotFailure(self.input, self.io, self.snapshot, 0);
        }

        fn inputChanged(_: @This()) PairCommandFailure {
            return inputChangedPairFailure(0);
        }
    };
    return runExactSecondPass(PairCommandFailure, &cursor, expected_count, Pass{
        .allocator = allocator,
        .reader = &reader,
        .validator = &validator,
        .mate1_storage = &mate1_storage,
        .staging_limit = staging_limit,
        .writer = writer,
        .input = input,
        .io = io,
        .snapshot = snapshot,
        .policy = options.pair_name_policy,
        .failure = &failure,
    });
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
        else
            inputOpenFailure(err) };
    };
    const file: std.Io.File = .{ .handle = handle, .flags = .{ .nonblocking = true } };
    const snapshot = fileSnapshot(file, io) catch |err| {
        file.close(io);
        return .{
            .failure = if (expected_snapshot != null and err == error.NotRegularFile)
                inputChangedFailure()
            else if (err == error.NotRegularFile)
                // An invalid request, as `-` is: exact-count sampling reads its input twice.
                CommandFailure.plain(
                    "usage",
                    "exact-count sampling requires a regular file",
                    2,
                )
            else
                INSPECT_FAILURE,
        };
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
        return INSPECT_FAILURE;
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
    options.output_identity = checkedStdoutIdentity(io) catch return 3;
    var stdout_buffer: [RECORD_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var stdout_sink = io_layer.FileSink.initStreaming(io, .stdout(), &stdout_buffer);
    var writer = zfastq.Writer.init(stdout_sink.byteSink());

    const failure = interleaveInputs(
        io,
        allocator,
        inputs,
        &writer,
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

/// Writes pairs straight from both readers' buffers until a mate needs the
/// general path; `pending1` returns mate 1 when only it was taken. Kept out
/// of line so the per-pair loop has its own small frame.
noinline fn joinSpanPairs(
    writer: *zfastq.Writer,
    reader1: *zfastq.Reader,
    reader2: *zfastq.Reader,
    validator1: *fastq.AdaptiveRecordValidator,
    validator2: *fastq.AdaptiveRecordValidator,
    selection: *PairOutputSelector,
    policy: pairing.NamePolicy,
    pending1: *?fastq.RecordSpan,
) error{WriteFailed}!?PairCommandFailure {
    while (fastq.nextRecordSpan(reader1, validator1)) |mate1| {
        const mate2 = fastq.nextRecordSpan(reader2, validator2) orelse {
            pending1.* = mate1;
            return null;
        };
        if (!pairing.headersMatch(mate1.header(), mate2.header(), policy)) {
            const record_index1 = reader1.recordIndex() - 1;
            return headerNameMismatch(
                mate1.header(),
                mate2.header(),
                policy,
                record_index1,
                .{ record_index1, reader2.recordIndex() - 1 },
                .{ reader1.currentRecordOffsets().?.header, reader2.currentRecordOffsets().?.header },
            );
        }
        if (!selection.selectPair()) continue;
        if (!spanIsWritable(mate1)) return pairCommandFailure(0, unwritableRecordFailure());
        if (!spanIsWritable(mate2)) return pairCommandFailure(1, unwritableRecordFailure());
        try writeRecordSpan(writer, mate1);
        try writeRecordSpan(writer, mate2);
    }
    return null;
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
    if (stdinRepeats(inputs)) {
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
    options: InterleaveOptions,
) error{WriteFailed}!?PairCommandFailure {
    var input1: RecordInput = undefined;
    var input2: RecordInput = undefined;
    if (initPairedRecordInputs(&input1, &input2, io, inputs, options.output_identity)) |failure| return failure;
    defer input1.deinit(io);
    defer input2.deinit(io);

    var selection = options.selection;
    return interleaveSources(
        allocator,
        input1.stream(),
        input2.stream(),
        writer,
        &selection,
        options,
    );
}

fn interleaveSources(
    allocator: std.mem.Allocator,
    source1: *std.Io.Reader,
    source2: *std.Io.Reader,
    writer: *zfastq.Writer,
    selection: *PairOutputSelector,
    options: InterleaveOptions,
) error{WriteFailed}!?PairCommandFailure {
    const reader_options: fastq.ReaderOptions = .{ .max_line_bytes = options.max_line_bytes };
    var reader1 = fastq.initStreamReader(allocator, source1, reader_options);
    defer reader1.deinit();
    var reader2 = fastq.initStreamReader(allocator, source2, reader_options);
    defer reader2.deinit();
    var validator1 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var validator2 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });

    while (true) {
        var span1: ?fastq.RecordSpan = null;
        if (try joinSpanPairs(
            writer,
            &reader1,
            &reader2,
            &validator1,
            &validator2,
            selection,
            options.pair_name_policy,
            &span1,
        )) |failure| return failure;
        // Drawn before the pair is read: an unselected pair is validated
        // without being kept.
        if (!selection.selectPair()) {
            const header1 = if (span1) |mate1| fastq.ValidatedHeader{ .header = mate1.header(), .semantic_error = null } else fastq.nextValidatedHeader(&reader1, &validator1) catch |err| {
                return pairCommandFailure(0, mapReaderFailure(&reader1, err));
            };
            const header2 = fastq.nextValidatedHeader(&reader2, &validator2) catch |err| {
                return pairCommandFailure(1, mapReaderFailure(&reader2, err));
            };
            switch (verifyPair(&reader1, &reader2, header1, header2, options.pair_name_policy)) {
                .end => return null,
                .failure => |failure| return failure,
                .pair => continue,
            }
        }
        var next1 = if (span1) |mate1| fastq.validatedSpan(mate1) else fastq.nextValidatedRecord(&reader1, &validator1) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader1, err));
        };
        var next2 = fastq.nextValidatedRecord(&reader2, &validator2) catch |err| {
            return pairCommandFailure(1, mapReaderFailure(&reader2, err));
        };
        // Pointer captures read both records where the readers built them.
        const maybe1: ?*const fastq.ValidatedRecord = if (next1) |*value| value else null;
        const maybe2: ?*const fastq.ValidatedRecord = if (next2) |*value| value else null;
        switch (verifyPair(
            &reader1,
            &reader2,
            if (maybe1) |mate1| mate1.validatedHeader() else null,
            if (maybe2) |mate2| mate2.validatedHeader() else null,
            options.pair_name_policy,
        )) {
            .end => return null,
            .failure => |failure| return failure,
            .pair => {},
        }
        const record1 = &maybe1.?.record;
        const record2 = &maybe2.?.record;
        const canonical_span1 = maybe1.?.canonical_span;
        const canonical_span2 = maybe2.?.canonical_span;
        if (recordHasUnwritableEnding(record1.*, canonical_span1)) {
            return pairCommandFailure(0, unwritableRecordFailure());
        }
        if (recordHasUnwritableEnding(record2.*, canonical_span2)) {
            return pairCommandFailure(1, unwritableRecordFailure());
        }

        try writeCheckedRecord(writer, record1.*, canonical_span1);
        try writeCheckedRecord(writer, record2.*, canonical_span2);
    }
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
    write_buffer: [RECORD_OUTPUT_BUFFER_BYTES]u8 = undefined,
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
    const staging_limit = recordStagingLimit(options.max_line_bytes);

    const input_label = inputs[0];
    const owns_input = !namesStdin(input_label);
    const input_file = openRecordFile(io, input_label) catch |err| {
        const failure = inputOpenFailure(err);
        printCommandFailure(io, input_label, failure);
        return failure.exit_code;
    };
    defer if (owns_input) input_file.close(io);
    // Checked before any output exists, so a failure here creates nothing.
    if (inputTypeFailure(io, input_file)) |failure| {
        printCommandFailure(io, input_label, failure);
        return failure.exit_code;
    }

    var output1 = DeinterleaveOutput.init(paths[0]);
    defer output1.close(io);
    var output2 = DeinterleaveOutput.init(paths[1]);
    defer output2.close(io);
    if (output1.create(io)) |failure| {
        printCommandFailure(io, output1.path, failure);
        return failure.exit_code;
    }
    if (output2.create(io)) |failure| {
        // Output 1 was created empty a moment ago; leaving it would block a rerun.
        output1.close(io);
        std.Io.Dir.cwd().deleteFile(io, output1.path) catch {};
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
        input.stream(),
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
    if (namesStdin(path1) or namesStdin(path2)) {
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

fn deinterleaveSource(
    allocator: std.mem.Allocator,
    source: *std.Io.Reader,
    writer1: *zfastq.Writer,
    writer2: *zfastq.Writer,
    staging_limit: usize,
    options: DeinterleaveOptions,
) DeinterleaveWriteError!?PairCommandFailure {
    var reader = fastq.initStreamReader(allocator, source, .{
        .max_line_bytes = options.max_line_bytes,
    });
    defer reader.deinit();
    var stored_name: StoredPairName = .{};
    defer stored_name.deinit(allocator);
    var mate1_storage: PreservedMate1 = .{};
    defer mate1_storage.deinit(allocator);
    var validator = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var every_pair: PairOutputSelector = .all;

    while (true) {
        var span1: ?fastq.RecordSpan = null;
        if (try splitSpanPairs(
            &reader,
            &validator,
            options.pair_name_policy,
            &every_pair,
            .{ writer1, writer2 },
            &span1,
        )) |failure| return failure;
        const pair = switch (readInterleavedPair(
            allocator,
            &reader,
            &validator,
            &stored_name,
            &mate1_storage,
            staging_limit,
            options.pair_name_policy,
            span1,
        )) {
            .end => return null,
            .failure => |failure| return failure,
            .pair => |pair| pair,
        };
        mate1_storage.write(writer1, &reader, pair) catch return error.Output1WriteFailed;
        writeCheckedRecord(writer2, pair.record2, pair.canonical_span2) catch
            return error.Output2WriteFailed;
    }
}

/// Where mate 1 of an interleaved pair lives while mate 2 is read, when the
/// reader's buffer cannot keep it.
const PreservedMate1 = struct {
    staged_record: std.ArrayList(u8) = .empty,
    retained: fastq.RetainedRecordStorage = .{},

    fn deinit(self: *PreservedMate1, allocator: std.mem.Allocator) void {
        self.staged_record.deinit(allocator);
        self.retained.deinit(allocator);
    }

    fn write(self: *PreservedMate1, writer: *zfastq.Writer, reader: *fastq.Reader, pair: InterleavedPair) error{WriteFailed}!void {
        return writePreservedInterleavedMate1(
            writer,
            reader,
            &self.retained,
            self.staged_record.items,
            pair.record1,
            pair.canonical_span1,
            pair.record1_storage,
        );
    }
};

const InterleavedPair = struct {
    record1: fastq.Record,
    canonical_span1: ?[]const u8,
    record1_storage: InterleavedFirstRecordStorage,
    record2: fastq.Record,
    canonical_span2: ?[]const u8,
};

const InterleavedPairStep = union(enum) {
    end,
    pair: InterleavedPair,
    failure: PairCommandFailure,
};

/// Reads one interleaved pair for writing. Mate 1 stays alive through mate
/// 2's read in the reader or in `mate1_storage`; both mates are validated,
/// their names compared, and their endings checked writable. A pair whose
/// mate 1 fails is reported as `verifyInterleavedPair` reports it. `span1`
/// is mate 1 when the span loop has already taken it.
fn readInterleavedPair(
    allocator: std.mem.Allocator,
    reader: *fastq.Reader,
    validator: *fastq.AdaptiveRecordValidator,
    stored_name: *StoredPairName,
    mate1_storage: *PreservedMate1,
    staging_limit: usize,
    policy: pairing.NamePolicy,
    span1: ?fastq.RecordSpan,
) InterleavedPairStep {
    var validated1 = if (span1) |mate1| fastq.validatedSpan(mate1) else fastq.nextValidatedRecord(reader, validator) catch |err| {
        return .{ .failure = pairCommandFailure(0, mapReaderFailure(reader, err)) };
    } orelse return .end;
    if (validated1.semantic_error != null) {
        return switch (verifyInterleavedPair(allocator, reader, validator, stored_name, policy, validated1.validatedHeader())) {
            .failure => |failure| .{ .failure = failure },
            // A mate 1 that fails semantic validation fails its pair.
            .end, .pair => unreachable,
        };
    }
    const record_index1 = reader.recordIndex() - 1;
    const offsets1 = reader.currentRecordOffsets().?;
    const header1_len = validated1.record.header.len;

    var record1_storage: InterleavedFirstRecordStorage = .reader;
    const paired_record2 = fastq.nextValidatedMatePreservingRecord(reader, &validated1, validator) catch |err| {
        return .{ .failure = pairCommandFailure(0, mapReaderFailure(reader, err)) };
    };
    const record1 = validated1.record;
    const canonical_span1 = validated1.canonical_span;
    const validated2 = paired_record2 orelse record: {
        @branchHint(.cold);
        const preserved = nextAfterPreservingInterleavedMate1(
            allocator,
            reader,
            &mate1_storage.retained,
            &mate1_storage.staged_record,
            record1,
            canonical_span1,
            staging_limit,
            validator,
        ) catch |err| return .{ .failure = mapPreservedMateFailure(reader, err) };
        record1_storage = preserved.first_storage;
        break :record preserved.record orelse return .{ .failure = missingInterleavedMateFailure(record_index1) };
    };

    const record_index2 = reader.recordIndex() - 1;
    const offsets2 = reader.currentRecordOffsets().?;
    if (mapSemanticFailure(
        validated2.semantic_error,
        offsets2,
        record_index2,
    )) |details| {
        return .{ .failure = pairCommandFailure(0, details) };
    }

    const header1 = switch (record1_storage) {
        .reader, .retained => record1.header,
        .staged => mate1_storage.staged_record.items[1 .. 1 + header1_len],
    };
    if (!pairing.headersMatch(header1, validated2.record.header, policy)) {
        return .{ .failure = headerNameMismatch(
            header1,
            validated2.record.header,
            policy,
            interleavedPairIndex(record_index1),
            .{ record_index1, record_index2 },
            .{ offsets1.header, offsets2.header },
        ) };
    }

    if (recordHasUnwritableEnding(record1, canonical_span1) or
        recordHasUnwritableEnding(validated2.record, validated2.canonical_span))
    {
        return .{ .failure = pairCommandFailure(0, unwritableRecordFailure()) };
    }
    return .{ .pair = .{
        .record1 = record1,
        .canonical_span1 = canonical_span1,
        .record1_storage = record1_storage,
        .record2 = validated2.record,
        .canonical_span2 = validated2.canonical_span,
    } };
}

/// Takes interleaved pairs straight from the reader's buffer until a mate
/// needs the general path; `pending1` returns mate 1 when only it was taken.
/// Mate 1 goes to `writers[0]` and mate 2 to `writers[1]`. The selection
/// draws once per pair, after both mates are taken and their names match, so
/// a pair handed to the general path is drawn there. Kept out of line so the
/// per-pair loop has its own small frame.
noinline fn splitSpanPairs(
    reader: *zfastq.Reader,
    validator: *fastq.AdaptiveRecordValidator,
    policy: pairing.NamePolicy,
    selection: *PairOutputSelector,
    writers: [2]*zfastq.Writer,
    pending1: *?fastq.RecordSpan,
) DeinterleaveWriteError!?PairCommandFailure {
    while (fastq.nextRecordSpan(reader, validator)) |mate1| {
        const offset1 = reader.currentRecordOffsets().?.header;
        const mate2 = fastq.nextRecordSpan(reader, validator) orelse {
            pending1.* = mate1;
            return null;
        };
        if (!pairing.headersMatch(mate1.header(), mate2.header(), policy)) {
            const record_index2 = reader.recordIndex() - 1;
            return headerNameMismatch(
                mate1.header(),
                mate2.header(),
                policy,
                interleavedPairIndex(record_index2 - 1),
                .{ record_index2 - 1, record_index2 },
                .{ offset1, reader.currentRecordOffsets().?.header },
            );
        }
        if (!selection.selectPair()) continue;
        if (!spanIsWritable(mate1) or !spanIsWritable(mate2)) {
            return pairCommandFailure(0, unwritableRecordFailure());
        }
        writeRecordSpan(writers[0], mate1) catch return error.Output1WriteFailed;
        writeRecordSpan(writers[1], mate2) catch return error.Output2WriteFailed;
    }
    return null;
}

fn flushDeinterleaveWriters(
    writer1: *zfastq.Writer,
    writer2: *zfastq.Writer,
) DeinterleaveWriteError!void {
    writer1.flush() catch return error.Output1WriteFailed;
    writer2.flush() catch return error.Output2WriteFailed;
}

// --- JSON output ---

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
    try json.objectField("code");
    try json.write(failure.code());
    try json.objectField("message");
    try json.write(failure.message());
    switch (failure) {
        .name_mismatch => |details| {
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
            inline for (@typeInfo(zfastq.StatsResult).@"struct".fields) |field| {
                try json.objectField(field.name);
                try json.write(@field(result, field.name));
            }
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

// --- Human output ---

fn printCount(io: std.Io, n: u64) !void {
    var buf: [32]u8 = undefined;
    const text = try std.fmt.bufPrint(&buf, "{d}\n", .{n});
    try std.Io.File.writeStreamingAll(.stdout(), io, text);
}

fn printStats(
    io: std.Io,
    label: []const u8,
    stats: *const zfastq.Stats,
    separator: bool,
) !void {
    if (separator) try std.Io.File.writeStreamingAll(.stdout(), io, "\n");
    try std.Io.File.writeStreamingAll(.stdout(), io, "input: ");
    try writeEscapedAll(.stdout(), io, label);
    try std.Io.File.writeStreamingAll(.stdout(), io, "\n");
    var buffer: [1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &buffer);
    const output = &stdout_writer.interface;
    const result = stats.result();
    const ratios = stats.ratios();
    // The fields and their order are the JSON document's.
    inline for (@typeInfo(zfastq.StatsResult).@"struct".fields) |field| {
        const value = @field(result, field.name);
        switch (field.type) {
            u64 => try writeUnsignedField(output, field.name, value),
            ?u64 => try writeOptionalUnsignedField(output, field.name, value),
            ?f64 => try writeRatioField(output, field.name, @field(ratios, field.name)),
            else => @compileError("no text form for " ++ field.name),
        }
    }
    try output.flush();
}

fn writeUnsignedField(output: *std.Io.Writer, name: []const u8, value: u64) !void {
    try output.print("{s}: {d}\n", .{ name, value });
}

fn writeOptionalUnsignedField(output: *std.Io.Writer, name: []const u8, value: ?u64) !void {
    if (value) |number| return writeUnsignedField(output, name, number);
    try output.print("{s}: -\n", .{name});
}

/// Rounds the exact ratio to six decimals, so text output never depends on floating point.
fn writeRatioField(
    output: *std.Io.Writer,
    name: []const u8,
    ratio: statistics.Ratio,
) !void {
    const numerator = ratio.numerator;
    const denominator = ratio.denominator;
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

const TestStream = struct {
    buffer: []u8,
    reader: io_layer.ByteSourceReader,

    fn init(source: zfastq.io.ByteSource) !TestStream {
        const buffer = try std.testing.allocator.alloc(
            u8,
            zfastq.limits.DEFAULT_READER_BUFFER_BYTES,
        );
        return .{ .buffer = buffer, .reader = .init(source, buffer) };
    }

    fn deinit(self: *TestStream) void {
        std.testing.allocator.free(self.buffer);
    }
};

/// A `std.testing.tmpDir` and a file in it, relative to the repository root.
const TEST_TEMP_DIR = ".zig-cache/tmp/{s}";
const TEST_TEMP_PATH = TEST_TEMP_DIR ++ "/{s}";

fn testTempPath(allocator: std.mem.Allocator, sub_path: []const u8, name: []const u8) ![]u8 {
    return std.fmt.allocPrint(allocator, TEST_TEMP_PATH, .{ sub_path, name });
}

/// Returns the per-input command failure; a structural pair failure is a test error here.
fn commandFailure(failure: PairCommandFailure) error{ExpectedCommandFailure}!@FieldType(PairCommandFailure, "command") {
    return switch (failure) {
        .command => |command| command,
        .pair => error.ExpectedCommandFailure,
    };
}

/// The rejection exact sampling reports when an input's snapshot changed before or during sampling.
fn expectInputChanged(failure: CommandFailure) !void {
    try std.testing.expectEqualStrings("input_changed", failure.code);
    try std.testing.expectEqualStrings("input changed during exact sampling", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
}

test "[unit] - [interleaved staging]: releases oversized slack before a smaller mate" {
    var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    const allocator = tracking.allocator();
    var staging: std.ArrayList(u8) = .empty;
    defer staging.deinit(allocator);
    const old_capacity = zfastq.limits.DEFAULT_READER_BUFFER_BYTES * 2;
    try staging.ensureTotalCapacityPrecise(allocator, old_capacity);

    try stageCanonicalRecord(
        allocator,
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
    try std.testing.expectEqual(old_capacity, tracking.freed_bytes);
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

test "[edge] - [stats command]: streamed totals match fixed totals at each buffer size" {
    // Two reads with 4 bases; quality bytes 0, 20, 40, 30 sum to 90.
    const valid = "@r\r\nAc\r\n+\r\n!5\r\n@s\nNG\n+\nI?\n";
    // The space is below the quality range; it is the second record's first quality byte.
    const invalid = "@r\nAc\n+\n!5\n@s\nNG\n+\n I\n";
    for ([_]usize{ 1, 4, 7, 64 }) |buffer_len| {
        var storage: [64]u8 = undefined;
        var source = zfastq.io.plain.SliceSource.init(valid);
        var stream = io_layer.ByteSourceReader.init(source.byteSource(), storage[0..buffer_len]);
        var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
            .fail_index = 0,
        });

        const outcome = collectStats(failing.allocator(), &stream.interface, .{});
        const result = switch (outcome) {
            .success => |stats| stats.result(),
            .failure => return error.UnexpectedFailure,
        };
        try std.testing.expect(!failing.has_induced_failure);
        try std.testing.expectEqual(@as(u64, 2), result.reads);
        try std.testing.expectEqual(@as(u64, 4), result.bases);
        try std.testing.expectEqual(@as(u64, 1), result.a);
        try std.testing.expectEqual(@as(u64, 1), result.c);
        try std.testing.expectEqual(@as(u64, 1), result.g);
        try std.testing.expectEqual(@as(u64, 1), result.n);
        try std.testing.expectEqual(@as(u64, 90), result.quality_sum);
        try std.testing.expectEqual(@as(u64, 3), result.q20_bases);
        try std.testing.expectEqual(@as(u64, 2), result.q30_bases);

        var invalid_source = zfastq.io.plain.SliceSource.init(invalid);
        var invalid_stream = io_layer.ByteSourceReader.init(invalid_source.byteSource(), storage[0..buffer_len]);
        const failure = switch (collectStats(std.testing.allocator, &invalid_stream.interface, .{})) {
            .success => return error.ExpectedFailure,
            .failure => |failure| failure,
        };
        try std.testing.expectEqualStrings("S006", failure.code);
        try std.testing.expectEqualStrings("quality byte must be ASCII 33 through 126", failure.message);
        try std.testing.expectEqual(@as(?u64, 1), failure.record_index);
        try std.testing.expectEqual(@as(?u64, 19), failure.byte_offset);
        try std.testing.expectEqual(@as(?u3, 4), failure.line_in_record);
    }
}

/// A standard I/O runtime over the testing runtime, which serves every operation not overridden
/// here. Paths open relative to `dir`. It counts
/// opened and closed files, positional writes and streaming reads. Standard output is captured
/// when `stdout` is set; standard error is captured into `stderr` or discarded. The
/// `fail_write`th positional write (1-based; 0 never fails) returns `error.NoSpaceLeft`.
const FakeIo = struct {
    dir: std.Io.Dir,
    stdout: ?*std.Io.Writer = null,
    stderr: ?*std.Io.Writer = null,
    /// Replaces the first opened path as soon as it opens.
    replace_first_open: ?FirstOpenReplacement = null,
    fail_write: usize = 0,
    /// The nth streaming read (1-based; 0 never fails) returns `error.InputOutput`.
    fail_read: usize = 0,
    opened: usize = 0,
    closed: usize = 0,
    writes: usize = 0,
    reads: usize = 0,
    vtable: @TypeOf(std.Io.failing.vtable.*) = undefined,

    const FirstOpenReplacement = union(enum) {
        /// Renames the path aside and writes these bytes under its name.
        rewrite: []const u8,
        /// Renames the path aside and hard-links this source under its name.
        link: []const u8,
    };

    fn io(self: *FakeIo) std.Io {
        self.vtable = std.testing.io.vtable.*;
        self.vtable.dirOpenFile = openFile;
        self.vtable.dirCreateFile = createFile;
        self.vtable.fileClose = closeFiles;
        self.vtable.fileReadPositional = readPositional;
        self.vtable.fileWritePositional = writePositional;
        self.vtable.operate = operate;
        return .{ .userdata = self, .vtable = &self.vtable };
    }

    fn openFile(
        ctx: ?*anyopaque,
        _: std.Io.Dir,
        path: []const u8,
        options: std.Io.Dir.OpenFileOptions,
    ) std.Io.File.OpenError!std.Io.File {
        const self: *FakeIo = @ptrCast(@alignCast(ctx.?));
        const file = try self.dir.openFile(std.testing.io, path, options);
        errdefer file.close(std.testing.io);
        self.opened += 1;
        if (self.opened == 1) {
            if (self.replace_first_open) |replacement| {
                self.replaceFirstInput(path, replacement) catch return error.Unexpected;
            }
        }
        return file;
    }

    fn replaceFirstInput(self: *FakeIo, path: []const u8, replacement: FirstOpenReplacement) !void {
        const io_runtime = std.testing.io;
        var saved_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
        const saved = try std.fmt.bufPrint(&saved_buffer, "{s}.saved", .{path});
        try self.dir.rename(path, self.dir, saved, io_runtime);
        switch (replacement) {
            .rewrite => |bytes| try self.dir.writeFile(io_runtime, .{ .sub_path = path, .data = bytes }),
            .link => |source| try self.dir.hardLink(source, self.dir, path, io_runtime, .{}),
        }
    }

    fn createFile(
        ctx: ?*anyopaque,
        _: std.Io.Dir,
        path: []const u8,
        options: std.Io.Dir.CreateFileOptions,
    ) std.Io.File.OpenError!std.Io.File {
        const self: *FakeIo = @ptrCast(@alignCast(ctx.?));
        return self.dir.createFile(std.testing.io, path, options);
    }

    fn closeFiles(ctx: ?*anyopaque, files: []const std.Io.File) void {
        const self: *FakeIo = @ptrCast(@alignCast(ctx.?));
        for (files) |file| file.close(std.testing.io);
        self.closed += files.len;
    }

    fn readPositional(
        _: ?*anyopaque,
        file: std.Io.File,
        data: []const []u8,
        offset: u64,
    ) std.Io.File.ReadPositionalError!usize {
        const io_runtime = std.testing.io;
        return io_runtime.vtable.fileReadPositional(io_runtime.userdata, file, data, offset);
    }

    fn writePositional(
        ctx: ?*anyopaque,
        file: std.Io.File,
        header: []const u8,
        data: []const []const u8,
        splat: usize,
        offset: u64,
    ) std.Io.File.WritePositionalError!usize {
        const self: *FakeIo = @ptrCast(@alignCast(ctx.?));
        self.writes += 1;
        if (self.writes == self.fail_write) return error.NoSpaceLeft;
        const io_runtime = std.testing.io;
        return io_runtime.vtable.fileWritePositional(io_runtime.userdata, file, header, data, splat, offset);
    }

    fn operate(ctx: ?*anyopaque, operation: std.Io.Operation) std.Io.Cancelable!std.Io.Operation.Result {
        const self: *FakeIo = @ptrCast(@alignCast(ctx.?));
        switch (operation) {
            .file_read_streaming => {
                self.reads += 1;
                if (self.reads == self.fail_read) return .{ .file_read_streaming = error.InputOutput };
            },
            .file_write_streaming => |write| {
                if (write.file.handle == std.Io.File.stdout().handle) {
                    if (self.stdout) |stdout| return .{ .file_write_streaming = stdout.writeSplatHeader(
                        write.header,
                        write.data,
                        write.splat,
                    ) catch error.NoSpaceLeft };
                } else if (write.file.handle == std.Io.File.stderr().handle) {
                    if (self.stderr) |stderr| return .{ .file_write_streaming = stderr.writeSplatHeader(
                        write.header,
                        write.data,
                        write.splat,
                    ) catch error.NoSpaceLeft };
                    var discarded = write.header.len;
                    for (write.data) |data| discarded += data.len * write.splat;
                    return .{ .file_write_streaming = discarded };
                }
            },
            else => {},
        }
        return std.testing.io.operate(operation);
    }
};

test "[failure] - [stats command]: a closed standard output keeps the partial report and the limit status" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const window = io_layer.DEFAULT_READER_BUFFER_BYTES;
    const long_line = try std.testing.allocator.alloc(u8, window + 2);
    defer std.testing.allocator.free(long_line);
    try tmp.dir.writeFile(io, .{ .sub_path = "small.fastq", .data = "@r\nA\n+\n!\n" });
    @memset(long_line, 'x');
    long_line[0] = '@';
    try tmp.dir.writeFile(io, .{ .sub_path = "limit.fastq", .data = long_line });

    // Stats allocates nothing; the bounded allocator proves it stays untouched.
    const storage = try std.testing.allocator.alloc(u8, 2 * window);
    defer std.testing.allocator.free(storage);
    const block = "input: small.fastq\n" ++
        "reads: 1\nbases: 1\nmin_length: 1\nmax_length: 1\nmean_length: 1.000000\n" ++
        "a: 1\nc: 0\ng: 0\nt: 0\nn: 0\nother_bases: 0\ngc_fraction: 0.000000\n" ++
        "quality_sum: 0\nmean_quality: 0.000000\nq20_bases: 0\nq20_fraction: 0.000000\n" ++
        "q30_bases: 0\nq30_fraction: 0.000000\n";
    // Standard output holds the first report and a partial second one, then fails.
    const partial = block ++ "\ninput: small.fastq\n";
    var stdout_buffer: [partial.len]u8 = undefined;
    var stdout: std.Io.Writer = .fixed(&stdout_buffer);
    var stderr_buffer: [256]u8 = undefined;
    var stderr: std.Io.Writer = .fixed(&stderr_buffer);
    var fake: FakeIo = .{ .dir = tmp.dir, .stdout = &stdout, .stderr = &stderr };
    var bounded = std.heap.FixedBufferAllocator.init(storage);

    const inputs = [_][]const u8{ "small.fastq", "limit.fastq", "small.fastq", "small.fastq" };
    const status = runStats(
        fake.io(),
        bounded.allocator(),
        &inputs,
        .{ .max_line_bytes = window + 1 },
        false,
    );

    try std.testing.expectEqual(@as(u8, 4), status);
    try std.testing.expectEqualStrings(partial, stdout.buffered());
    try std.testing.expectEqualStrings(
        "error: limit.fastq: line length limit exceeded\n" ++
            "error: standard output: I/O error\n",
        stderr.buffered(),
    );
    try std.testing.expectEqual(@as(usize, 3), fake.opened);
    try std.testing.expectEqual(fake.opened, fake.closed);
    try std.testing.expectEqual(@as(usize, 0), bounded.end_index);
}

fn snapshotTestFile(io: std.Io, path: []const u8) !FileSnapshot {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);
    return fileSnapshot(file, io);
}

test "[failure] - [exact sample]: selection allocation failure precedes later block format error" {
    var record_inputs: [2]RecordInput = undefined;
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
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "records.fastq" });
    var selector = sampling.ExactSelector.init(1, 11);
    defer selector.deinit(std.testing.allocator);
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 0,
    });

    const result = sampleExactFirstPass(
        io,
        failing.allocator(),
        &record_inputs[0],
        path,
        &selector,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
        },
    );
    const failure = switch (result) {
        .success => return error.ExpectedFailure,
        .failure => |failure| failure,
    };

    try std.testing.expect(failing.has_induced_failure);
    try std.testing.expectEqualStrings("out_of_memory", failure.code);
    try std.testing.expectEqualStrings("out of memory", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expectEqual(@as(u64, 2), selector.record_count);
}

test "[failure] - [exact sample]: final record-count change keeps a valid prefix" {
    var record_inputs: [2]RecordInput = undefined;
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
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, name });
    const snapshot = try snapshotTestFile(io, path);

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactSecondPass(
        false,
        io,
        std.testing.allocator,
        &record_inputs[0],
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

    try expectInputChanged(failure);
    try std.testing.expectEqualStrings("@one\nA\n+\n!\n", sink.written());
}

test "[failure] - [exact sample]: changed metadata stops the second pass before output" {
    var record_inputs: [2]RecordInput = undefined;
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
        &record_inputs[0],
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

    try expectInputChanged(failure);
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
    // A connected FIFO keeps a blocking open of the path from hanging this in-process test.
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
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, name });
    var input: RecordInput = undefined;
    const failure = switch (initExactInput(&input, io, path, snapshot, null)) {
        .failure => |failure| failure,
        .success => {
            input.deinit(io);
            return error.ExpectedFailure;
        },
    };
    try expectInputChanged(failure);
}

test "[failure] - [exact sample]: a failed read in the first pass reports io error before any selection" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "records.fastq", .data = "@one\nA\n+\n!\n" });
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "records.fastq" });
    var fake: FakeIo = .{ .dir = std.Io.Dir.cwd(), .fail_read = 1 };
    var record_input: RecordInput = undefined;
    var selector = sampling.ExactSelector.init(1, 11);
    defer selector.deinit(std.testing.allocator);

    const failure = switch (sampleExactFirstPass(
        fake.io(),
        std.testing.allocator,
        &record_input,
        path,
        &selector,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
        },
    )) {
        .success => return error.ExpectedFailure,
        .failure => |failure| failure,
    };
    try std.testing.expectEqualStrings("io_error", failure.code);
    try std.testing.expectEqualStrings("I/O error", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expectEqual(@as(u64, 0), selector.record_count);
}

/// The gzip magic, deflate method, no flags, a zero modification time and extra flags, unknown
/// operating system, then the final stored block's header byte (BFINAL set, BTYPE 00).
const gzip_stored_member_header = [_]u8{ 0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255, 1 };

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
        try output.interface.writeAll(&gzip_stored_member_header);
        // A stored block stores LEN and its one's complement NLEN before its bytes.
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
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try testTempPath(allocator, &tmp.sub_path, "input");
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
                const first = sampleExactFirstPass(
                    io,
                    allocator,
                    &record_inputs[0],
                    path,
                    &selector,
                    options,
                );
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
                    &record_inputs[0],
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
                    try expectInputChanged(failure.?);
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
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path1 = try testTempPath(allocator, &tmp.sub_path, "r1");
    defer allocator.free(path1);
    const path2 = try testTempPath(allocator, &tmp.sub_path, "r2");
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
                    const first = sampleExactPairFirstPass(
                        io,
                        allocator,
                        &record_inputs,
                        inputs,
                        &selector,
                        options,
                    );
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
                        &record_inputs,
                        inputs,
                        &writer,
                        if (selection == .all) .empty else .{
                            .low_words = &.{if (selection == .selected) 2 else 1},
                            .middle_bytes = &.{0},
                        },
                        snapshots,
                        selector.record_count,
                        recordStagingLimit(options.max_line_bytes),
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
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try testTempPath(allocator, &tmp.sub_path, "pairs");
    defer allocator.free(path);
    const small_prefix = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n";
    // The first pair starts this many bytes before the reader buffer ends, so its mate one crosses the refill.
    const refill_margin = 300;
    const options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 2,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    const staging_limit = recordStagingLimit(options.max_line_bytes);
    // The name change lands inside the 192-byte identifier that both mates share.
    const header_len = 192;
    const name_change_offset = 180;
    const Change = enum { valid, alphabet2, name, structure, write };
    for ([_]InterleavedFirstRecordStorage{ .reader, .retained }) |storage| {
        var input: std.ArrayList(u8) = .empty;
        defer input.deinit(allocator);
        if (storage == .reader) {
            // Put the selected pair across the transport refill while both mates stay small.
            try input.appendSlice(allocator, "@a/1 ");
            try input.appendNTimes(allocator, 'x', zfastq.limits.DEFAULT_READER_BUFFER_BYTES - refill_margin - small_prefix.len - 1);
            try input.appendSlice(allocator, small_prefix[4..]);
        } else try input.appendSlice(allocator, small_prefix);
        const prefix_len = input.items.len;
        var header_starts: [2]usize = undefined;
        var sequence_starts: [2]usize = undefined;
        var plus_starts: [2]usize = undefined;
        for (0..2) |mate| {
            const field_len: usize = if (storage == .retained)
                zfastq.limits.DEFAULT_READER_BUFFER_BYTES
            else
                8;
            try input.append(allocator, '@');
            header_starts[mate] = input.items.len;
            try input.appendNTimes(allocator, 'b', header_len);
            try input.appendSlice(allocator, if (mate == 0) "/1\n" else "/2\n");
            sequence_starts[mate] = input.items.len;
            try input.appendNTimes(allocator, if (mate == 0) 'A' else 'T', field_len);
            try input.append(allocator, '\n');
            plus_starts[mate] = input.items.len;
            try input.appendSlice(allocator, "+\n");
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
            if (try fastq.nextValidatedMatePreservingRecord(&reader, &first_record, &validator)) |_| {
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
        const first = sampleExactPairFirstPass(
            io,
            allocator,
            &record_inputs,
            &.{path},
            &selector,
            options,
        );
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
                    .alphabet2 => changed[sequence_starts[1]] = '.',
                    .name => changed[header_starts[1] + name_change_offset] = 'c',
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
                    &record_inputs[0],
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
                    if (change == .structure) {
                        try std.testing.expectEqual(@as(?u64, 3), details.record_index);
                        try std.testing.expectEqual(@as(?u3, 3), details.line_in_record);
                        try std.testing.expectEqual(@as(?u64, plus_starts[1]), details.byte_offset);
                    }
                }
                try std.testing.expectEqualStrings(expected_prefix, sink.written());
            }
        }
    }
}

test "[integration] - [interleaved exact sample]: preserves the first mate parse failure" {
    // Every line of this input fits; only the changed plus line fails the second pass.
    const small_line_limit = 64;
    // The byte offset the diagnostic reports for the malformed plus line in `changed`.
    const plus_line_offset: u64 = 40;
    var record_inputs: [2]RecordInput = undefined;
    const SnapshotChangingSink = struct {
        inner: io_layer.SliceSink,
        path: []const u8,
        mtime: ?i96,
        changed: bool = false,

        fn write(ctx: *anyopaque, bytes: []const u8) error{WriteFailed}!void {
            const self: *@This() = @ptrCast(@alignCast(ctx));
            if (self.mtime) |mtime| {
                const file = std.Io.Dir.cwd().openFile(std.testing.io, self.path, .{}) catch
                    return error.WriteFailed;
                defer file.close(std.testing.io);
                file.setTimestamps(std.testing.io, .{
                    .modify_timestamp = .{ .new = .{ .nanoseconds = mtime } },
                }) catch return error.WriteFailed;
                self.mtime = null;
                self.changed = true;
            }
            try self.inner.byteSink().write(bytes);
        }
    };
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try testTempPath(allocator, &tmp.sub_path, "pairs");
    defer allocator.free(path);
    const original = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n@b/1\nA\n+\n!\n@b/2\nT\n+\n#\n";
    const changed = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n@b/1\nA\n+\n!\n@b/2\nT\nx\n#\n";
    const options = SampleOptions{
        .max_line_bytes = small_line_limit,
        .alphabet = .iupac,
        .fraction = null,
        .count = 2,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    const Case = enum { all, second_only, metadata_changed };
    for ([_]bool{ false, true }) |gzip| {
        for (std.enums.values(Case)) |case| {
            errdefer std.debug.print("gzip={} case={s}\n", .{ gzip, @tagName(case) });
            try writeExactTestInput(io, path, original, gzip, null);
            const snapshot = try snapshotTestFile(io, path);
            try writeExactTestInput(io, path, changed, gzip, snapshot);
            var output: [original.len]u8 = undefined;
            var sink = SnapshotChangingSink{
                .inner = io_layer.SliceSink.init(&output),
                .path = path,
                .mtime = if (case == .metadata_changed) snapshot.mtime_nanoseconds + std.time.ns_per_s else null,
            };
            var writer = zfastq.Writer.init(.{
                .ctx = &sink,
                .vtable = &.{ .write = SnapshotChangingSink.write },
            });
            const select_all = case != .second_only;
            const failure = (try sampleExactInterleavedSecondPass(
                select_all,
                io,
                allocator,
                &record_inputs[0],
                path,
                &writer,
                if (select_all) .empty else .{ .low_words = &.{2}, .middle_bytes = &.{0} },
                snapshot,
                2,
                recordStagingLimit(options.max_line_bytes),
                options,
            )).?;
            try std.testing.expect(failure == .command);
            try std.testing.expectEqual(@as(u1, 0), failure.command.input_index);
            const details = failure.command.details;
            if (case == .metadata_changed) {
                try std.testing.expect(sink.changed);
                try expectInputChanged(details);
                try std.testing.expect(details.record_index == null);
                try std.testing.expect(details.line_in_record == null);
                try std.testing.expect(details.byte_offset == null);
            } else {
                try std.testing.expectEqualStrings("S001", details.code);
                try std.testing.expectEqualStrings("plus line must start with '+'", details.message);
                try std.testing.expectEqual(@as(u8, 1), details.exit_code);
                try std.testing.expectEqual(@as(?u64, 3), details.record_index);
                try std.testing.expectEqual(@as(?u3, 3), details.line_in_record);
                try std.testing.expectEqual(@as(?u64, plus_line_offset), details.byte_offset);
            }
            try std.testing.expectEqualStrings(if (select_all) original[0..22] else "", sink.inner.written());
        }
    }
}

test "[failure] - [paired exact sample]: each input snapshot is checked independently" {
    var record_inputs: [2]RecordInput = undefined;
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
    const path1 = try std.fmt.bufPrint(&path1_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r1.fastq" });
    var path2_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path2 = try std.fmt.bufPrint(&path2_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r2.fastq" });
    var pairs_path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const pairs_path = try std.fmt.bufPrint(&pairs_path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "pairs.fastq" });
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
            &record_inputs,
            &inputs,
            &writer,
            .empty,
            snapshots,
            1,
            options,
        )).?;
        const command = try commandFailure(failure);
        try std.testing.expectEqual(@as(u1, @intCast(changed_input)), command.input_index);
        try expectInputChanged(command.details);
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
        &record_inputs[0],
        pairs_path,
        &writer,
        .empty,
        changed_pair_snapshot,
        1,
        recordStagingLimit(zfastq.limits.DEFAULT_MAX_LINE_BYTES),
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
            .pair_mode = .interleaved,
        },
    )).?;
    const command = try commandFailure(failure);
    try expectInputChanged(command.details);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}

test "[integration] - [paired exact sample]: the output pass checks structure and count" {
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var path1_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path1 = try std.fmt.bufPrint(&path1_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r1.fastq" });
    var path2_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path2 = try std.fmt.bufPrint(&path2_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r2.fastq" });
    var pairs_path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const pairs_path = try std.fmt.bufPrint(&pairs_path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "pairs.fastq" });
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
        &record_inputs,
        &inputs,
        &paired_writer,
        .empty,
        .{ snapshot1, snapshot2 },
        2,
        paired_options,
    )).?;
    const count_command = try commandFailure(count_failure);
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
        &record_inputs[0],
        pairs_path,
        &interleaved_writer,
        .empty,
        pair_snapshot,
        1,
        recordStagingLimit(zfastq.limits.DEFAULT_MAX_LINE_BYTES),
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
            .pair_mode = .interleaved,
        },
    )).?;
    const odd_command = try commandFailure(odd_failure);
    try expectInputChanged(odd_command.details);
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
        &record_inputs,
        &inputs,
        &structural_writer,
        .empty,
        .{ structural_snapshot1, structural_snapshot2 },
        1,
        paired_options,
    )).?;
    const structural_command = try commandFailure(structural_failure);
    try std.testing.expectEqual(@as(u1, 1), structural_command.input_index);
    try std.testing.expectEqualStrings("S001", structural_command.details.code);
    try std.testing.expectEqual(@as(u8, 1), structural_command.details.exit_code);
    try std.testing.expectEqual(@as(usize, 0), structural_sink.written().len);
}

test "[unit] - [paired identifiers]: storage reuses its largest allocation" {
    // One byte past the reader buffer, so the first long store must grow the storage.
    const name = try std.testing.allocator.alloc(u8, zfastq.limits.DEFAULT_READER_BUFFER_BYTES + 1);
    defer std.testing.allocator.free(name);
    @memset(name, 'r');
    name[name.len - 1] = 'z';
    var tracking = std.testing.FailingAllocator.init(std.testing.allocator, .{});
    {
        var storage: std.ArrayList(u8) = .empty;
        defer storage.deinit(tracking.allocator());
        try storeNormalizedId(tracking.allocator(), &storage, name[0..3]);
        try storeNormalizedId(tracking.allocator(), &storage, name);
        const allocated = tracking.allocated_bytes;
        for ([_]usize{ 7, name.len, 0, name.len }) |len| {
            try storeNormalizedId(tracking.allocator(), &storage, name[0..len]);
            try std.testing.expectEqualStrings(name[0..len], storage.items);
            try std.testing.expectEqual(allocated, tracking.allocated_bytes);
        }
    }
    try std.testing.expectEqual(tracking.allocated_bytes, tracking.freed_bytes);

    // The second allocation grows the storage for the long name; failing it reports out of memory.
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 1 });
    var storage: std.ArrayList(u8) = .empty;
    defer storage.deinit(failing.allocator());
    try storeNormalizedId(failing.allocator(), &storage, name[0..3]);
    try std.testing.expectError(error.OutOfMemory, storeNormalizedId(failing.allocator(), &storage, name));
}

fn exerciseExactAllocations(allocator: std.mem.Allocator) !void {
    const r1 = "@a/1\nA\n+\n!\n@b/1\nC\n+\n#\n@c/1\nG\n+\n$\n";
    const r2 = "@a/2\nT\n+\n!\n@b/2\nG\n+\n#\n@c/2\nC\n+\n$\n";
    const interleaved =
        "@a/1\nA\n+\n!\n@a/2\nT\n+\n!\n" ++
        "@b/1\nC\n+\n#\n@b/2\nG\n+\n#\n" ++
        "@c/1\nG\n+\n$\n@c/2\nC\n+\n$\n";
    // Outputs for exact counts 0 to 4 with seed 11, pinned from the three inputs above.
    const record_a1 = "@a/1\nA\n+\n!\n";
    const record_b1 = "@b/1\nC\n+\n#\n";
    const record_c1 = "@c/1\nG\n+\n$\n";
    const pair_a = "@a/1\nA\n+\n!\n@a/2\nT\n+\n!\n";
    const pair_b = "@b/1\nC\n+\n#\n@b/2\nG\n+\n#\n";
    const pair_c = "@c/1\nG\n+\n$\n@c/2\nC\n+\n$\n";
    const none_outputs = [_][]const u8{
        "",
        record_a1,
        record_a1 ++ record_c1,
        record_a1 ++ record_b1 ++ record_c1,
        record_a1 ++ record_b1 ++ record_c1,
    };
    const pair_outputs = [_][]const u8{
        "",
        pair_a,
        pair_a ++ pair_c,
        pair_a ++ pair_b ++ pair_c,
        pair_a ++ pair_b ++ pair_c,
    };
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const names = [_][]const u8{ "r1.fastq", "r2.fastq", "interleaved.fastq" };
    var paths: [3][]const u8 = undefined;
    var path_count: usize = 0;
    defer for (paths[0..path_count]) |path| std.testing.allocator.free(path);
    for (names, [_][]const u8{ r1, r2, interleaved }, 0..) |name, bytes, index| {
        try tmp.dir.writeFile(io, .{ .sub_path = name, .data = bytes });
        paths[index] = try testTempPath(std.testing.allocator, &tmp.sub_path, name);
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
                if (try sampleExactPairs(io, allocator, inputs, &writer, count, recordStagingLimit(options.max_line_bytes), options)) |failure| {
                    try std.testing.expectEqual(descriptors, try openDescriptorCount());
                    return pairAllocationFailure(failure);
                }
            }
            try std.testing.expectEqual(descriptors, try openDescriptorCount());
            const expected = switch (mode) {
                .none => none_outputs[@intCast(count)],
                .paired, .interleaved => pair_outputs[@intCast(count)],
            };
            try std.testing.expectEqualStrings(expected, sink.written());
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
        const r1_path = try testTempPath(std.testing.allocator, &tmp.sub_path, "r1");
        defer std.testing.allocator.free(r1_path);
        const r2_path = try testTempPath(std.testing.allocator, &tmp.sub_path, "r2");
        defer std.testing.allocator.free(r2_path);
        // The first open replaces the path it opened; the aliased pair must be refused before any read.
        var fake: FakeIo = .{
            .dir = std.Io.Dir.cwd(),
            .replace_first_open = if (aliases) .{ .rewrite = "replacement" } else .{ .link = r2_path },
        };
        var output: [128]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        const descriptors = try openDescriptorCount();

        const failure = try interleaveInputs(
            fake.io(),
            std.testing.allocator,
            &.{ r1_path, r2_path },
            &writer,
            .{ .max_line_bytes = 8192, .alphabet = .iupac, .pair_name_policy = .illumina },
        );
        try writer.flush();

        if (aliases) {
            try std.testing.expectEqual(@as(u8, 2), failure.?.exitCode());
            try std.testing.expectEqualStrings("same_input", failure.?.command.details.code);
            try std.testing.expectEqual(@as(usize, 0), fake.reads);
        } else {
            try std.testing.expect(failure == null);
        }
        try std.testing.expectEqualStrings(if (aliases) "" else r1 ++ r2, sink.written());
        try std.testing.expectEqual(@as(usize, 1), fake.opened);
        try std.testing.expectEqual(@as(usize, 2), fake.closed);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
    }
}

test "[integration] - [output aliases]: rejection closes inputs on both exact passes" {
    var record_inputs: [2]RecordInput = undefined;
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
        paths[index] = try testTempPath(allocator, &tmp.sub_path, name);
        const file = try tmp.dir.openFile(io, name, .{});
        defer file.close(io);
        snapshots[index] = try fileSnapshot(file, io);
        identities[index] = try fileIdentity(io, file);
    }
    const descriptors = try openDescriptorCount();
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
    const single_failure = (try sampleExactSecondPass(
        true,
        io,
        allocator,
        &record_inputs[0],
        paths[0],
        &writer,
        .empty,
        snapshots[0],
        1,
        options,
    )).?;
    try expectInputChanged(single_failure);
    const interleaved_failure = (try sampleExactInterleavedSecondPass(
        true,
        io,
        allocator,
        &record_inputs[0],
        paths[0],
        &writer,
        .empty,
        snapshots[0],
        1,
        recordStagingLimit(options.max_line_bytes),
        options,
    )).?;
    try std.testing.expectEqualStrings("input_changed", interleaved_failure.command.details.code);
    try std.testing.expectEqual(@as(u8, 3), interleaved_failure.exitCode());
    try std.testing.expectEqual(descriptors, try openDescriptorCount());
    for (identities, 0..) |identity, side| {
        options.output_identity = identity;
        const paired_stream = initPairedRecordInputs(&input1, &input2, io, &paths, identity).?;
        try std.testing.expectEqualStrings("same_output", paired_stream.command.details.code);
        try std.testing.expectEqual(@as(u1, @intCast(side)), paired_stream.command.input_index);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
        const paired_failure = (try sampleExactPairedSecondPass(
            true,
            io,
            allocator,
            &record_inputs,
            &paths,
            &writer,
            .empty,
            snapshots,
            1,
            options,
        )).?;
        try std.testing.expectEqualStrings("input_changed", paired_failure.command.details.code);
        try std.testing.expectEqual(@as(u8, 3), paired_failure.exitCode());
        try std.testing.expectEqual(@as(u1, @intCast(side)), paired_failure.command.input_index);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
    }
    try writer.flush();
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}

test "[failure] - [input resources]: failed opens close owned files and keep borrowed files open" {
    const io = std.testing.io;
    // Standard error is discarded, and the outputs open relative to the working directory like the other paths.
    var quiet: FakeIo = .{ .dir = std.Io.Dir.cwd() };
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "r1", .data = "@a/1\nA\n+\n!\n" });
    try tmp.dir.createDir(io, "directory", .default_dir);
    const base = try std.fmt.allocPrint(std.testing.allocator, TEST_TEMP_DIR, .{tmp.sub_path});
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
    var input: RecordInput = undefined;
    try std.testing.expect(initRecordInput(&input, io, path, null) == null);
    input.deinit(io);
    const directory_failure = initRecordInput(&input, io, directory, null).?;
    try std.testing.expectEqualStrings("io_error", directory_failure.code);
    try std.testing.expectEqualStrings("input is a directory", directory_failure.message);
    const usage_failure = initExactInput(&input, io, directory, null, null).failure;
    try std.testing.expectEqualStrings("usage", usage_failure.code);
    try std.testing.expectEqualStrings("exact-count sampling requires a regular file", usage_failure.message);
    const changed_failure = initExactInput(&input, io, path, changed, null).failure;
    try expectInputChanged(changed_failure);
    const missing_failure = initExactInput(&input, io, missing, changed, null).failure;
    try expectInputChanged(missing_failure);
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
        // Exact-count sampling rejects a directory as an invalid request.
        const exact_expected = if (index >= 3) "usage" else expected;
        const exact_failure = initExactPairedInputs(&input1, &input2, io, &paths, null, null).failure;
        try std.testing.expectEqualStrings(exact_expected, exact_failure.command.details.code);
        try std.testing.expectEqual(descriptors, try openDescriptorCount());
    }

    {
        const borrowed = try tmp.dir.openFile(io, "r1", .{});
        defer borrowed.close(io);
        try input.init(io, borrowed, false);
        input.deinit(io);
        try std.testing.expectEqual(@as(u64, "@a/1\nA\n+\n!\n".len), (try borrowed.stat(io)).size);
    }
    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const inputs = [_][]const u8{ path, missing };
    const failure = (try interleaveInputs(io, std.testing.allocator, &inputs, &writer, .{
        .max_line_bytes = 8192,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    })).?;
    try std.testing.expectEqual(@as(u1, 1), failure.command.input_index);
    try std.testing.expectEqualStrings("io_error", failure.command.details.code);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
    try std.testing.expectEqual(@as(u8, 3), runDeinterleave(quiet.io(), std.testing.allocator, &.{path}, output_path, path, .{
        .max_line_bytes = 8192,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    }));
    try std.testing.expectEqual(descriptors, try openDescriptorCount());
    // Output 2 already existed, so the output 1 it created is gone again.
    try std.testing.expectError(error.FileNotFound, tmp.dir.statFile(io, "out1", .{}));
}

test "[failure] - [paired inputs]: an unusable file handle is reported at its input index" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "r1", .data = "@a/1\nA\n+\n!\n" });
    const borrowed = try tmp.dir.openFile(io, "r1", .{});
    defer borrowed.close(io);
    const invalid: std.Io.File = .{ .handle = -1, .flags = .{ .nonblocking = false } };
    for ([_][2]std.Io.File{ .{ invalid, borrowed }, .{ borrowed, invalid } }, 0..) |files, side| {
        const failure = pairedFilesFailure(io, files[0], files[1], null).?;
        try std.testing.expectEqual(@as(u8, 3), failure.exitCode());
        try std.testing.expectEqual(@as(u1, @intCast(side)), failure.command.input_index);
        try std.testing.expectEqualStrings("failed to inspect file", failure.command.details.message);
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

test "[failure] - [exact sample]: an allocation failure in exact sampling releases memory and descriptors" {
    try exerciseExactAllocations(std.testing.allocator);

    // Every allocation in this body is first-pass reservoir growth, and the second pass allocates
    // nothing for these inputs, so index 0 (the first reservoir growth) is the first induced failure.
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    try std.testing.expectError(error.OutOfMemory, exerciseExactAllocations(failing.allocator()));
    try std.testing.expect(failing.has_induced_failure);
}

test "[edge] - [paired exact sample]: target above pair count keeps indexes implicit" {
    var source1 = io_layer.SliceSource.init("@a/1\nA\n+\n!\n@b/1\nC\n+\n#\n");
    var source2 = io_layer.SliceSource.init("@a/2\nT\n+\n!\n@b/2\nG\n+\n#\n");
    var selector = sampling.ExactSelector.init(3, 11);
    defer selector.deinit(std.testing.allocator);

    var stream1 = try TestStream.init(source1.byteSource());
    defer stream1.deinit();
    var stream2 = try TestStream.init(source2.byteSource());
    defer stream2.deinit();
    const failure = checkPairedSources(
        std.testing.allocator,
        &stream1.reader.interface,
        &stream2.reader.interface,
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
}

test "[failure] - [paired exact sample]: a read failure follows complete earlier pairs" {
    const first1 = "@a/1\nA\n+\n!\n";
    const first2 = "@a/2\nT\n+\n!\n";
    // The first input fails right after its first complete pair.
    var source1: ScriptedSource = .{
        .data = first1 ++ "@b/1\nC\n+\n#\n",
        .fail_at = first1.len,
    };
    var source2 = io_layer.SliceSource.init(first2 ++ "@b/2\nG\n+\n#\n");
    var selector = sampling.ExactSelector.init(1, 11);
    defer selector.deinit(std.testing.allocator);

    var stream1 = try TestStream.init(source1.byteSource());
    defer stream1.deinit();
    var stream2 = try TestStream.init(source2.byteSource());
    defer stream2.deinit();
    const failure = checkPairedSources(
        std.testing.allocator,
        &stream1.reader.interface,
        &stream2.reader.interface,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .pair_mode = .paired,
            .pair_name_policy = .illumina,
        },
        &selector,
    ).?;

    const command = try commandFailure(failure);
    try std.testing.expectEqual(@as(u1, 0), command.input_index);
    try std.testing.expectEqualStrings("io_error", command.details.code);
    try std.testing.expectEqual(@as(u8, 3), command.details.exit_code);
    try std.testing.expectEqual(@as(u64, 1), selector.record_count);
}

test "[property] - [interleaved fraction sample]: source chunks preserve selected pairs" {
    const data = "@a/1 left\r\nAC\r\n+left\r\n!#\r\n@a/2 right\nGT\n+right\n$%\n" ++
        "@b/1\nG\n+\n&\n@b/2\nC\n+\n'\n" ++
        "@c/1\nTTA\n+one\n()*\n@c/2\nAAT\n+two\n+,-";
    // Seed 11 at one half selects pair indexes 0 and 2.
    const expected = "@a/1 left\nAC\n+left\n!#\n@a/2 right\nGT\n+right\n$%\n" ++
        "@c/1\nTTA\n+one\n()*\n@c/2\nAAT\n+two\n+,-\n";
    // Reads of one byte, two bytes, and the whole input.
    for ([_]usize{ 1, 2, data.len }) |chunk_len| {
        var source: ScriptedSource = .{ .data = data, .first_chunk = chunk_len, .chunk_len = chunk_len };
        var output: [expected.len]u8 = undefined;
        var sink = io_layer.SliceSink.init(&output);
        var writer = zfastq.Writer.init(sink.byteSink());
        var selector = sampling.Selector.init(.{ .probability = 0.5 }, 11);
        var selection: PairOutputSelector = .{ .fraction = &selector };

        var stream = try TestStream.init(source.byteSource());
        defer stream.deinit();
        const failure = try sampleInterleavedSource(
            std.testing.allocator,
            &stream.reader.interface,
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

test "[edge] - [single fraction sample]: fraction zero reads without a record reader" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(std.testing.allocator);
    try input.appendSlice(std.testing.allocator, "@large ");
    try input.appendNTimes(std.testing.allocator, 'x', io_layer.DEFAULT_READER_BUFFER_BYTES + 1);
    const header_len = input.items.len;
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "input.fastq" });

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

test "[edge] - [paired fraction sample]: fraction zero reads buffered mate one without allocating" {
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var input_bytes: std.ArrayList(u8) = .empty;
    defer input_bytes.deinit(std.testing.allocator);
    // Mate one is a 128 KiB header that fits in the reader buffer, so it is read in place.
    const long_header_len = 128 * 1024;
    try input_bytes.appendSlice(std.testing.allocator, "@large/1 ");
    try input_bytes.appendNTimes(std.testing.allocator, 'x', long_header_len);
    try input_bytes.appendSlice(
        std.testing.allocator,
        "\nA\n+\n!\n@large/2\nT\n+\n#\n",
    );
    try tmp.dir.writeFile(io, .{
        .sub_path = "pairs.fastq",
        .data = input_bytes.items,
    });
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "pairs.fastq" });
    var output: [1]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    var selector = sampling.Selector.init(.none, 11);
    var output_selector: PairOutputSelector = .{ .fraction = &selector };
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 0,
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

/// A byte source over `data`. Its first read returns at most `first_chunk` bytes, later reads
/// at most `chunk_len`, and a read at `fail_at` returns `error.ReadFailed` instead of data.
const ScriptedSource = struct {
    data: []const u8,
    first_chunk: usize = std.math.maxInt(usize),
    chunk_len: usize = std.math.maxInt(usize),
    fail_at: ?usize = null,
    position: usize = 0,

    fn byteSource(self: *ScriptedSource) zfastq.io.ByteSource {
        return .{ .ctx = self, .vtable = &vtable };
    }

    const vtable = zfastq.io.ByteSource.VTable{ .read = read };

    fn read(ctx: *anyopaque, dest: []u8) error{ReadFailed}!usize {
        const self: *ScriptedSource = @ptrCast(@alignCast(ctx));
        const end = self.fail_at orelse self.data.len;
        if (self.position == self.fail_at) return error.ReadFailed;
        const limit = if (self.position == 0) self.first_chunk else self.chunk_len;
        const amount = @min(dest.len, limit, end - self.position);
        @memcpy(dest[0..amount], self.data[self.position..][0..amount]);
        self.position += amount;
        return amount;
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
        // The second read point lands 7 bytes into the second mate's header.
        const into_second_header = 7;
        for ([_]usize{ mates[0].len, mates[0].len + into_second_header }) |first_chunk| {
            for (0..4) |command| {
                var source: ScriptedSource = .{ .data = input, .first_chunk = first_chunk, .chunk_len = 3 };
                var stream = try TestStream.init(source.byteSource());
                defer stream.deinit();
                const bytes = &stream.reader.interface;
                var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
                    .fail_index = 1,
                });
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
                    1, 2 => try sampleInterleavedSource(failing.allocator(), bytes, &writer1, &selection, recordStagingLimit(64), .{
                        .max_line_bytes = 64,
                        .alphabet = .iupac,
                        .fraction = if (command == 1) .none else .all,
                        .count = null,
                        .seed = 11,
                        .pair_mode = .interleaved,
                        .pair_name_policy = .illumina,
                    }),
                    3 => try deinterleaveSource(failing.allocator(), bytes, &writer1, &writer2, recordStagingLimit(64), .{
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

const DeinterleaveTestSink = struct {
    buffer: [256]u8 = undefined,
    capacity: usize = 256,
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
        if (end > self.capacity) return error.WriteFailed;
        @memcpy(self.buffer[self.length..end], bytes);
        self.length = end;
    }

    fn flush(ctx: *anyopaque) error{WriteFailed}!void {
        const self: *DeinterleaveTestSink = @ptrCast(@alignCast(ctx));
        self.flush_count += 1;
        if (self.fail_flush) return error.WriteFailed;
    }
};

test "[failure] - [interleave]: a write failure after the first mate keeps only that mate" {
    const mate1 = "@pair/1\nA\n+left\n!\n";
    const mate2 = "@pair/2\nT\n+right\n#\n";
    var source1 = io_layer.SliceSource.init(mate1);
    var source2 = io_layer.SliceSource.init(mate2);
    var stream1 = try TestStream.init(source1.byteSource());
    defer stream1.deinit();
    var stream2 = try TestStream.init(source2.byteSource());
    defer stream2.deinit();
    // The sink holds exactly the first mate, so the second mate's first write is refused.
    var sink = DeinterleaveTestSink{ .capacity = mate1.len };
    var writer = zfastq.Writer.init(sink.byteSink());
    var selection: PairOutputSelector = .all;

    try std.testing.expectError(error.WriteFailed, interleaveSources(
        std.testing.allocator,
        &stream1.reader.interface,
        &stream2.reader.interface,
        &writer,
        &selection,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .pair_name_policy = .illumina,
        },
    ));
    try std.testing.expectEqualStrings(mate1, sink.buffer[0..sink.length]);
}

test "[failure] - [deinterleave]: names the output whose write failed and writes nothing partial" {
    const lf_input = "@pair/1\nA\n+left\n!\n@pair/2\nT\n+right\n#\n";
    const crlf_first_mate = "@pair/1\r\nR\r\n+left\r\n!\r\n";
    const crlf_input = crlf_first_mate ++ "@pair/2\r\nT\r\n+right\r\n#\r\n";
    const first_mate = "@pair/1\nA\n+left\n!\n";
    const cases = [_]struct {
        input: []const u8,
        fail_output1: bool = false,
        fail_output2: bool = false,
        output1_capacity: usize = 256,
        max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        first_chunk: usize = std.math.maxInt(usize),
        chunk_len: usize = std.math.maxInt(usize),
        expected_error: anyerror,
        expected_output1: []const u8,
    }{
        .{
            .input = lf_input,
            .fail_output1 = true,
            .expected_error = error.Output1WriteFailed,
            .expected_output1 = "",
        },
        .{
            .input = lf_input,
            .fail_output2 = true,
            .expected_error = error.Output2WriteFailed,
            .expected_output1 = first_mate,
        },
        // Canonical staging submits the complete mate as one write, so a sink one byte short accepts none of it.
        .{
            .input = crlf_input,
            .output1_capacity = first_mate.len - 1,
            .max_line_bytes = 64,
            .first_chunk = crlf_first_mate.len,
            .chunk_len = 3,
            .expected_error = error.Output1WriteFailed,
            .expected_output1 = "",
        },
    };
    for (cases) |case| {
        const options = DeinterleaveOptions{
            .max_line_bytes = case.max_line_bytes,
            .alphabet = .iupac,
            .pair_name_policy = .illumina,
        };
        var source: ScriptedSource = .{
            .data = case.input,
            .first_chunk = case.first_chunk,
            .chunk_len = case.chunk_len,
        };
        var stream = try TestStream.init(source.byteSource());
        defer stream.deinit();
        var sink1 = DeinterleaveTestSink{ .fail_write = case.fail_output1, .capacity = case.output1_capacity };
        var sink2 = DeinterleaveTestSink{ .fail_write = case.fail_output2 };
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());

        try std.testing.expectError(case.expected_error, deinterleaveSource(
            std.testing.allocator,
            &stream.reader.interface,
            &writer1,
            &writer2,
            recordStagingLimit(options.max_line_bytes),
            options,
        ));
        try std.testing.expectEqualStrings(case.expected_output1, sink1.buffer[0..sink1.length]);
        try std.testing.expectEqual(@as(usize, 0), sink2.length);
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
            var stderr: std.Io.Writer = .fixed(&stderr_buffer);
            var fake: FakeIo = .{
                .dir = tmp.dir,
                .stderr = &stderr,
                .fail_write = fail_write,
            };

            const status = runDeinterleave(
                fake.io(),
                std.testing.allocator,
                &.{"input"},
                "out1",
                "out2",
                .{ .max_line_bytes = 8, .alphabet = .iupac, .pair_name_policy = .illumina },
            );

            const expected_status = if (fail_write == 0) case.exit_code else @max(case.exit_code, IO_FAILURE.exit_code);
            const expected_stderr = if (fail_write == 0) case.stderr else try std.fmt.allocPrint(
                allocator,
                "{s}error: out{d}: I/O error\n",
                .{ case.stderr, fail_write },
            );
            try std.testing.expectEqual(expected_status, status);
            try std.testing.expectEqualStrings(expected_stderr, stderr.buffered());
            try std.testing.expectEqual(@as(usize, if (fail_write == 1) 1 else 2), fake.writes);
            try std.testing.expectEqual(@as(usize, 3), fake.closed);
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
    var stream = try TestStream.init(source.byteSource());
    defer stream.deinit();
    var sink1 = DeinterleaveTestSink{};
    var sink2 = DeinterleaveTestSink{};
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 1,
    });
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };

    const failure = (try deinterleaveSource(
        failing.allocator(),
        &stream.reader.interface,
        &writer1,
        &writer2,
        recordStagingLimit(options.max_line_bytes),
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
}

/// Mate one's record is 2 * first_len + 12 bytes, two bytes short of two reader buffers.
const fallback_first_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES - 7;
const fallback_second_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES;
/// Header, both mate lines and the plus line of mate one, with the two sequence and quality fields.
const fallback_first_record_len = "@pair/1\n".len + fallback_first_len + "\n+\n".len + fallback_first_len + "\n".len;

fn fallbackPairInput(allocator: std.mem.Allocator) !std.ArrayList(u8) {
    var input: std.ArrayList(u8) = .empty;
    errdefer input.deinit(allocator);
    try input.appendSlice(allocator, "@pair/1\n");
    try input.appendNTimes(allocator, 'A', fallback_first_len);
    try input.appendSlice(allocator, "\n+\n");
    try input.appendNTimes(allocator, '!', fallback_first_len);
    try input.appendSlice(allocator, "\n@pair/2\n");
    try input.appendNTimes(allocator, 'T', fallback_second_len);
    try input.appendSlice(allocator, "\n+\n");
    try input.appendNTimes(allocator, '#', fallback_second_len);
    try input.append(allocator, '\n');
    return input;
}

test "[edge] - [deinterleave]: mate one is retained across the reader refill" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);

    var source = io_layer.SliceSource.init(input.items);
    var reader = try zfastq.Reader.init(std.testing.allocator, source.byteSource(), .{});
    defer reader.deinit();
    var retained: fastq.RetainedRecordStorage = .{};
    defer retained.deinit(std.testing.allocator);
    var staging: std.ArrayList(u8) = .empty;
    defer staging.deinit(std.testing.allocator);
    var span: ?[]const u8 = null;
    var validator = fastq.AdaptiveRecordValidator.init(.{});
    const record1 = (try fastq.nextWithSpan(&reader, &span)).?;
    try std.testing.expect((try fastq.nextBufferedWithSpan(&reader, &span)) == null);
    const next = try nextAfterPreservingInterleavedMate1(
        std.testing.allocator,
        &reader,
        &retained,
        &staging,
        record1,
        null,
        recordStagingLimit(reader.options.max_line_bytes),
        &validator,
    );
    try std.testing.expect(next.first_storage == .retained);
    try std.testing.expect(std.mem.allEqual(u8, record1.sequence, 'A'));
    try std.testing.expectEqual(fallback_second_len, next.record.?.record.sequence.len);
    try std.testing.expect(std.mem.allEqual(u8, next.record.?.record.sequence, 'T'));
    try std.testing.expect(next.record.?.semantic_error == null);
}

test "[failure] - [deinterleave]: retained pair storage reports out of memory on a later allocation" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);
    try exerciseRetainedPairAllocations(std.testing.allocator, input.items);

    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{
        .fail_index = 1,
    });
    try std.testing.expectError(
        error.OutOfMemory,
        exerciseRetainedPairAllocations(failing.allocator(), input.items),
    );
    try std.testing.expect(failing.has_induced_failure);
}

test "[failure] - [deinterleave]: a failed output write or a one-byte staging limit emits no partial pair" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);
    const output = try std.testing.allocator.alloc(u8, input.items.len);
    defer std.testing.allocator.free(output);
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };

    for (0..3) |failure_mode| {
        var source = io_layer.SliceSource.init(input.items);
        var stream = try TestStream.init(source.byteSource());
        defer stream.deinit();
        var sink1 = io_layer.SliceSink.init(if (failure_mode == 0) output[0..0] else output);
        var sink2 = io_layer.SliceSink.init(output[0..0]);
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());
        const result = deinterleaveSource(
            std.testing.allocator,
            &stream.reader.interface,
            &writer1,
            &writer2,
            if (failure_mode == 2) 1 else recordStagingLimit(options.max_line_bytes),
            options,
        );
        if (failure_mode < 2) {
            try std.testing.expectError(
                if (failure_mode == 0) error.Output1WriteFailed else error.Output2WriteFailed,
                result,
            );
            if (failure_mode == 0) try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
            if (failure_mode == 1) {
                try std.testing.expectEqualStrings(input.items[0..fallback_first_record_len], sink1.written());
            }
        } else {
            const failure = (try result).?;
            try std.testing.expectEqualStrings("arithmetic_limit", failure.command.details.code);
            try std.testing.expectEqual(@as(usize, 0), sink1.written().len);
            try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
        }
    }
}

test "[failure] - [deinterleave]: malformed pairs fail with their diagnostic before any output" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);
    const options = DeinterleaveOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .pair_name_policy = .illumina,
    };
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
        var stream = try TestStream.init(source.byteSource());
        defer stream.deinit();
        var sink1 = DeinterleaveTestSink{};
        var sink2 = DeinterleaveTestSink{};
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());
        const failure = (try deinterleaveSource(
            std.testing.allocator,
            &stream.reader.interface,
            &writer1,
            &writer2,
            recordStagingLimit(options.max_line_bytes),
            options,
        )).?;
        switch (failure_mode) {
            0, 4 => {
                try std.testing.expectEqualStrings("S005", failure.command.details.code);
                try std.testing.expectEqualStrings("sequence and quality lengths differ", failure.command.details.message);
            },
            1 => {
                try std.testing.expectEqualStrings("S006", failure.command.details.code);
                try std.testing.expectEqualStrings("quality byte must be ASCII 33 through 126", failure.command.details.message);
            },
            2 => {
                try std.testing.expect(failure.pair == .name_mismatch);
                try std.testing.expectEqualStrings("paired identifiers or mate markers do not match", failure.pair.message());
            },
            3 => {
                try std.testing.expect(failure.pair == .count_mismatch);
                try std.testing.expectEqualStrings("paired input is missing a mate", failure.pair.message());
            },
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
    var stream = try TestStream.init(source.byteSource());
    defer stream.deinit();
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
        &stream.reader.interface,
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
