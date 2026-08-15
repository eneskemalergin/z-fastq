//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("root.zig");
const fastq = @import("fastq.zig");
const io_layer = @import("io.zig");
const sampling = @import("sample.zig");

const INVALID_QUALITY_MESSAGE = "quality byte must be ASCII 33 through 126";
const PAIR_DIAGNOSTIC_PREFIX_BYTES = 128;

const USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count    Count records in plain or gzip FASTQ inputs
    \\  stats    Report aggregate FASTQ statistics
    \\  check    Validate FASTQ structure, sequence alphabet, and quality range
    \\  sample   Select records by deterministic probability or exact count
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
    \\Check options:
    \\  --paired             Validate two inputs as paired reads
    \\  --interleaved        Validate consecutive records as paired reads
    \\  --pair-names POLICY  Select illumina (default) or exact pair names
    \\
    \\Sample options:
    \\  --fraction P         Use 0, 1, 0.DIGITS, or 1.ZEROES
    \\  --count K            Select exactly min(K, records) from a file
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
    \\
;

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const gpa = init.gpa;

    var args = std.process.Args.Iterator.initAllocator(init.minimal.args, gpa) catch {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
        std.process.exit(3);
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
    else {
        std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown command: ") catch {};
        writeEscaped(.stderr(), io, cmd);
        std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
        printUsageAndExit(io);
    };

    var max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES;
    var alphabet: zfastq.Alphabet = .iupac;
    var pair_mode: PairMode = .none;
    var pair_name_policy: PairNamePolicy = .illumina;
    var pair_names_set = false;
    var json_output = false;
    var fraction: ?sampling.Fraction = null;
    var sample_count: ?u64 = null;
    var sample_seed: u64 = 11;
    var positional = std.ArrayList([]const u8).empty;
    defer positional.deinit(gpa);

    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--")) {
            while (args.next()) |path| {
                positional.append(gpa, path) catch {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: out of memory\n",
                    ) catch {};
                    std.process.exit(3);
                };
            }
            break;
        }
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            printHelpAndExit(io);
        }
        if (std.mem.eql(u8, arg, "--max-line-bytes")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --max-line-bytes requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            max_line_bytes = std.fmt.parseInt(usize, value, 10) catch |err| switch (err) {
                error.Overflow => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: --max-line-bytes exceeds supported limit\n",
                    ) catch {};
                    std.process.exit(4);
                },
                error.InvalidCharacter => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: invalid --max-line-bytes value\n",
                    ) catch {};
                    std.process.exit(2);
                },
            };
            continue;
        }
        if ((command == .stats or command == .check) and std.mem.eql(u8, arg, "--json")) {
            json_output = true;
            continue;
        }
        if ((command == .check or command == .sample) and
            std.mem.eql(u8, arg, "--alphabet"))
        {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --alphabet requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            alphabet = if (std.mem.eql(u8, value, "iupac"))
                .iupac
            else if (std.mem.eql(u8, value, "acgtn"))
                .acgtn
            else {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --alphabet must be iupac or acgtn\n",
                ) catch {};
                std.process.exit(2);
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--fraction")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --fraction requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            fraction = sampling.Fraction.parse(value) catch {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: invalid --fraction value\n",
                ) catch {};
                std.process.exit(2);
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--count")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --count requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            sample_count = sampling.parseCount(value) catch |err| switch (err) {
                error.InvalidCount => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: invalid --count value\n",
                    ) catch {};
                    std.process.exit(2);
                },
                error.Overflow => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: --count exceeds supported limit\n",
                    ) catch {};
                    std.process.exit(4);
                },
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--seed")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --seed requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            sample_seed = sampling.parseSeed(value) catch |err| switch (err) {
                error.InvalidSeed => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: invalid --seed value\n",
                    ) catch {};
                    std.process.exit(2);
                },
                error.Overflow => {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: --seed exceeds supported limit\n",
                    ) catch {};
                    std.process.exit(4);
                },
            };
            continue;
        }
        if (command == .check and std.mem.eql(u8, arg, "--paired")) {
            if (pair_mode == .interleaved) {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --paired and --interleaved are mutually exclusive\n",
                ) catch {};
                std.process.exit(2);
            }
            pair_mode = .paired;
            continue;
        }
        if (command == .check and std.mem.eql(u8, arg, "--interleaved")) {
            if (pair_mode == .paired) {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --paired and --interleaved are mutually exclusive\n",
                ) catch {};
                std.process.exit(2);
            }
            pair_mode = .interleaved;
            continue;
        }
        if (command == .check and std.mem.eql(u8, arg, "--pair-names")) {
            const value = args.next() orelse {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --pair-names requires a value\n",
                ) catch {};
                std.process.exit(2);
            };
            pair_name_policy = if (std.mem.eql(u8, value, "illumina"))
                .illumina
            else if (std.mem.eql(u8, value, "exact"))
                .exact
            else {
                std.Io.File.writeStreamingAll(
                    .stderr(),
                    io,
                    "error: --pair-names must be illumina or exact\n",
                ) catch {};
                std.process.exit(2);
            };
            pair_names_set = true;
            continue;
        }
        if (arg.len > 1 and std.mem.startsWith(u8, arg, "-")) {
            std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown ") catch {};
            std.Io.File.writeStreamingAll(.stderr(), io, @tagName(command)) catch {};
            std.Io.File.writeStreamingAll(.stderr(), io, " option: ") catch {};
            writeEscaped(.stderr(), io, arg);
            std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
            std.process.exit(2);
        }
        positional.append(gpa, arg) catch {
            std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
            std.process.exit(3);
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
            },
        ),
    };
    std.process.exit(code);
}

const Command = enum { count, stats, check, sample };

fn printUsageAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, USAGE) catch {};
    std.process.exit(2);
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
    if (inputs.len == 0) {
        std.Io.File.writeStreamingAll(
            .stderr(),
            io,
            "error: count requires at least one input\n",
        ) catch {};
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

    var exit_code: u8 = 0;
    for (inputs) |input| {
        const count_result = if (std.mem.eql(u8, input, "-"))
            countStdin(io, options)
        else
            countFile(io, input, options);
        const count = count_result catch |err| switch (err) {
            error.Io => {
                exit_code = @max(exit_code, 3);
                continue;
            },
            error.Format => {
                exit_code = @max(exit_code, 1);
                continue;
            },
            error.Limit => {
                printPathError(io, input, "line length limit exceeded");
                exit_code = @max(exit_code, 4);
                continue;
            },
            error.OutOfMemory => {
                std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                return 3;
            },
        };

        printCount(io, count) catch return @max(exit_code, 3);
    }

    return exit_code;
}

const CountError = error{
    Io,
    Format,
    Limit,
    OutOfMemory,
};

fn countFile(
    io: std.Io,
    path: []const u8,
    options: InputOptions,
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

    return countInput(io, path, file, options);
}

fn countStdin(io: std.Io, options: InputOptions) CountError!u64 {
    return countInput(io, "-", std.Io.File.stdin(), options);
}

fn countInput(
    io: std.Io,
    label: []const u8,
    file: std.Io.File,
    options: InputOptions,
) CountError!u64 {
    var sniff_buffer: [2]u8 = undefined;
    var sniff_reader = file.readerStreaming(io, &sniff_buffer);
    const prefix = sniff_reader.interface.peek(2) catch |err| switch (err) {
        error.EndOfStream => null,
        error.ReadFailed => {
            printPathError(io, label, "I/O error");
            return error.Io;
        },
    };

    if (prefix) |bytes| {
        if (std.mem.eql(u8, bytes, &.{ 0x1f, 0x8b })) {
            var gzip_buffer: [64 * 1024]u8 = undefined;
            // The gzip reader must replay the magic consumed by the sniff reader.
            @memcpy(gzip_buffer[0..2], bytes);
            var gzip_reader = file.readerStreaming(io, &gzip_buffer);
            gzip_reader.pos = 2;
            gzip_reader.interface.end = 2;
            var gzip_source = io_layer.GzipSource.init(&gzip_reader.interface);
            return countGzipSource(io, label, &gzip_source, options);
        }
    }

    var plain_source = zfastq.io.plain.ReaderSource.init(&sniff_reader.interface);
    return countSource(io, label, plain_source.byteSource(), options);
}

fn countGzipSource(
    io: std.Io,
    label: []const u8,
    source: *io_layer.GzipSource,
    options: InputOptions,
) CountError!u64 {
    const scan_options = zfastq.count_scan.Options{
        .max_line_bytes = options.max_line_bytes,
    };
    var scanner = zfastq.count_scan.Scanner.init(scan_options);
    var decompressor_buffer: [io_layer.COUNT_DECOMPRESS_BUFFER_BYTES]u8 = undefined;
    while (io_layer.readGzipChunk(source, &decompressor_buffer) catch {
        printPathError(io, label, "I/O error");
        return error.Io;
    }) |decoded| {
        _ = scanner.feed(decoded) catch |err| {
            return mapScanError(io, label, &scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapScanError(io, label, &scanner, err);
    };
    return scanner.record_index;
}

fn countSource(
    io: std.Io,
    label: []const u8,
    source: zfastq.io.ByteSource,
    options: InputOptions,
) CountError!u64 {
    const scan_options = zfastq.count_scan.Options{
        .max_line_bytes = options.max_line_bytes,
    };
    var scanner = zfastq.count_scan.Scanner.init(scan_options);
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const n = source.read(&buf) catch {
            printPathError(io, label, "I/O error");
            return error.Io;
        };
        if (n == 0) break;
        _ = scanner.feed(buf[0..n]) catch |err| {
            return mapScanError(io, label, &scanner, err);
        };
    }
    scanner.finishEof() catch |err| {
        return mapScanError(io, label, &scanner, err);
    };
    return scanner.record_index;
}

fn mapScanError(
    io: std.Io,
    path: []const u8,
    scanner: *zfastq.count_scan.Scanner,
    err: zfastq.ReaderError,
) CountError {
    return switch (err) {
        error.S001InvalidPlusLine,
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

const CommandFailure = struct {
    code: []const u8,
    message: []const u8,
    exit_code: u8,
    record_index: ?u64 = null,
    byte_offset: ?u64 = null,
    line_in_record: ?u3 = null,

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

fn validateCommandRecord(
    record: zfastq.Record,
    offsets: zfastq.RecordOffsets,
    record_index: u64,
    alphabet: zfastq.Alphabet,
) ?CommandFailure {
    const semantic_error = zfastq.validateRecord(
        record,
        .{ .alphabet = alphabet },
    ) orelse return null;
    const field_offset = switch (semantic_error.field) {
        .sequence => offsets.sequence,
        .quality => offsets.quality,
    };
    const relative_offset = std.math.cast(u64, semantic_error.byte_index) orelse {
        return CommandFailure.plain(
            "arithmetic_limit",
            "input location exceeds supported limit",
            4,
        );
    };
    const byte_offset = std.math.add(u64, field_offset, relative_offset) catch {
        return CommandFailure.plain(
            "arithmetic_limit",
            "input location exceeds supported limit",
            4,
        );
    };
    return CommandFailure.lint(.{
        .code = semantic_error.code,
        .message = semantic_error.message,
        .record_index = record_index,
        .byte_offset = byte_offset,
        .line_in_record = switch (semantic_error.field) {
            .sequence => 2,
            .quality => 4,
        },
    });
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

const PairNamePolicy = enum {
    illumina,
    exact,
};

const CheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
};

const PairedCheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
    pair_mode: PairMode = .none,
    pair_name_policy: PairNamePolicy = .illumina,
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

const PairName = struct {
    first_token: []const u8,
    normalized_id: []const u8,
    mate_markers: u2,
    first_mate_marker: ?u2,
};

const PairRecordDiagnostic = struct {
    record_index: u64,
    byte_offset: u64,
    first_token: BoundedBytes,
    normalized_id: BoundedBytes,
    mate_markers: u2,

    fn init(name: PairName, record_index: u64, byte_offset: u64) PairRecordDiagnostic {
        return .{
            .record_index = record_index,
            .byte_offset = byte_offset,
            .first_token = .init(name.first_token),
            .normalized_id = .init(name.normalized_id),
            .mate_markers = name.mate_markers,
        };
    }

    fn initStored(
        normalized_id: []const u8,
        first_token_len: usize,
        first_mate_marker: ?u2,
        mate_markers: u2,
        record_index: u64,
        byte_offset: u64,
    ) PairRecordDiagnostic {
        var first_token = BoundedBytes.init(normalized_id);
        first_token.full_len = first_token_len;
        if (first_mate_marker) |marker| {
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
            .normalized_id = .init(normalized_id),
            .mate_markers = mate_markers,
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

const CheckFailure = union(enum) {
    command: struct {
        input_index: u1,
        details: CommandFailure,
    },
    pair: PairFailure,

    fn exitCode(self: CheckFailure) u8 {
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

// Keep check code generation out of the measured count dispatch path.
noinline fn runCheck(
    io: std.Io,
    inputs: []const []const u8,
    options: CheckOptions,
    json_output: bool,
) u8 {
    if (validateRecordCommandInputs(io, .check, inputs)) |exit_code| return exit_code;
    if (json_output) return runCheckJson(io, inputs, options);

    var exit_code: u8 = 0;
    for (inputs) |input| {
        const failure = if (std.mem.eql(u8, input, "-"))
            checkStdin(io, options)
        else
            checkFile(io, input, options);
        if (failure) |details| {
            printCommandFailure(io, input, details);
            exit_code = @max(exit_code, details.exit_code);
        }
    }
    return exit_code;
}

// Keep paired check code generation out of the measured single-input path.
noinline fn runPairedCheckCommand(
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
    if (failure) |details| printCheckFailure(io, inputs, options.pair_mode, details);
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
        const failure = if (std.mem.eql(u8, input, "-"))
            checkStdin(io, options)
        else
            checkFile(io, input, options);
        writeCheckJsonResult(&json, input, failure) catch return 3;
        if (failure) |details| exit_code = @max(exit_code, details.exit_code);
    }
    finishJsonDocument(&json) catch return 3;
    stdout_writer.interface.flush() catch return 3;
    return exit_code;
}

fn runPairedCheckJson(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: PairedCheckOptions,
) u8 {
    var stdout_buffer: [16 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var json: std.json.Stringify = .{ .writer = &stdout_writer.interface };

    beginJsonDocument(&json, "z-fastq/check-v1") catch return 3;
    const failure = checkPairMode(io, allocator, inputs, options);
    writePairedCheckJsonResult(&json, inputs, options.pair_mode, failure) catch return 3;
    finishJsonDocument(&json) catch return 3;
    stdout_writer.interface.flush() catch return 3;
    return if (failure) |details| details.exitCode() else 0;
}

fn checkPairMode(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: PairedCheckOptions,
) ?CheckFailure {
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
) ?CheckFailure {
    var input1: RecordInput = undefined;
    if (initCheckInput(&input1, io, inputs[0])) |failure| {
        return .{ .command = .{ .input_index = 0, .details = failure } };
    }
    defer input1.deinit(io);

    var input2: RecordInput = undefined;
    if (initCheckInput(&input2, io, inputs[1])) |failure| {
        return .{ .command = .{ .input_index = 1, .details = failure } };
    }
    defer input2.deinit(io);

    var reader1 = zfastq.Reader.init(
        allocator,
        input1.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return commandCheckFailure(0, "out_of_memory", "out of memory", 3);
    defer reader1.deinit();
    var reader2 = zfastq.Reader.init(
        allocator,
        input2.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return commandCheckFailure(1, "out_of_memory", "out of memory", 3);
    defer reader2.deinit();

    while (true) {
        const record1 = reader1.next() catch |err| {
            return .{ .command = .{
                .input_index = 0,
                .details = mapReaderFailure(&reader1, err),
            } };
        };
        const record2 = reader2.next() catch |err| {
            return .{ .command = .{
                .input_index = 1,
                .details = mapReaderFailure(&reader2, err),
            } };
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
        if (validatePairRecord(record1.?, offsets1, record_index1, options)) |failure| {
            return .{ .command = .{ .input_index = 0, .details = failure } };
        }
        if (validatePairRecord(record2.?, offsets2, record_index2, options)) |failure| {
            return .{ .command = .{ .input_index = 1, .details = failure } };
        }

        const name1 = parsePairName(record1.?.header, options.pair_name_policy);
        const name2 = parsePairName(record2.?.header, options.pair_name_policy);
        if (!pairNamesMatch(name1, name2)) {
            return .{ .pair = .{ .name_mismatch = .{
                .pair_index = record_index1,
                .records = .{
                    .init(name1, record_index1, offsets1.header),
                    .init(name2, record_index2, offsets2.header),
                },
            } } };
        }
    }
}

fn checkInterleaved(
    io: std.Io,
    allocator: std.mem.Allocator,
    input_label: []const u8,
    options: PairedCheckOptions,
) ?CheckFailure {
    var input: RecordInput = undefined;
    if (initCheckInput(&input, io, input_label)) |failure| {
        return .{ .command = .{ .input_index = 0, .details = failure } };
    }
    defer input.deinit(io);

    var reader = zfastq.Reader.init(
        allocator,
        input.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return commandCheckFailure(0, "out_of_memory", "out of memory", 3);
    defer reader.deinit();
    var normalized_id: std.ArrayList(u8) = .empty;
    defer normalized_id.deinit(allocator);

    while (true) {
        const record1 = reader.next() catch |err| {
            return .{ .command = .{
                .input_index = 0,
                .details = mapReaderFailure(&reader, err),
            } };
        } orelse return null;
        const record_index1 = reader.recordIndex() - 1;
        const offsets1 = reader.currentRecordOffsets().?;
        const semantic1 = validatePairRecord(record1, offsets1, record_index1, options);

        var first_token_len: usize = 0;
        var first_mate_marker: ?u2 = null;
        var mate1_markers: u2 = 0;
        if (semantic1 == null) {
            const name1 = parsePairName(record1.header, options.pair_name_policy);
            normalized_id.clearRetainingCapacity();
            normalized_id.ensureTotalCapacityPrecise(allocator, name1.normalized_id.len) catch {
                return commandCheckFailure(0, "out_of_memory", "out of memory", 3);
            };
            normalized_id.appendSliceAssumeCapacity(name1.normalized_id);
            first_token_len = name1.first_token.len;
            first_mate_marker = name1.first_mate_marker;
            mate1_markers = name1.mate_markers;
        }

        const record2 = reader.next() catch |err| {
            return .{ .command = .{
                .input_index = 0,
                .details = mapReaderFailure(&reader, err),
            } };
        } orelse return .{ .pair = .{ .count_mismatch = .{
            .pair_index = record_index1 / 2,
            .remaining_side = 0,
            .record_indexes = .{
                record_index1,
                if (record_index1 == 0) null else record_index1 - 1,
            },
        } } };

        if (semantic1) |failure| {
            return .{ .command = .{ .input_index = 0, .details = failure } };
        }

        const record_index2 = reader.recordIndex() - 1;
        const offsets2 = reader.currentRecordOffsets().?;
        if (validatePairRecord(record2, offsets2, record_index2, options)) |failure| {
            return .{ .command = .{ .input_index = 0, .details = failure } };
        }

        const name1: PairName = .{
            .first_token = normalized_id.items,
            .normalized_id = normalized_id.items,
            .mate_markers = mate1_markers,
            .first_mate_marker = first_mate_marker,
        };
        const name2 = parsePairName(record2.header, options.pair_name_policy);
        if (!pairNamesMatch(name1, name2)) {
            return .{ .pair = .{ .name_mismatch = .{
                .pair_index = record_index1 / 2,
                .records = .{
                    .initStored(
                        normalized_id.items,
                        first_token_len,
                        first_mate_marker,
                        mate1_markers,
                        record_index1,
                        offsets1.header,
                    ),
                    .init(name2, record_index2, offsets2.header),
                },
            } } };
        }
    }
}

fn initCheckInput(
    input: *RecordInput,
    io: std.Io,
    label: []const u8,
) ?CommandFailure {
    if (std.mem.eql(u8, label, "-")) {
        input.init(io, .stdin(), false) catch {
            return CommandFailure.plain("io_error", "I/O error", 3);
        };
        return null;
    }

    const file = std.Io.Dir.cwd().openFile(io, label, .{}) catch |err| switch (err) {
        error.FileNotFound => return CommandFailure.plain("io_error", "file not found", 3),
        else => return CommandFailure.plain("io_error", "failed to open file", 3),
    };
    input.init(io, file, true) catch {
        file.close(io);
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    return null;
}

fn commandCheckFailure(
    input_index: u1,
    code: []const u8,
    message: []const u8,
    exit_code: u8,
) CheckFailure {
    return .{ .command = .{
        .input_index = input_index,
        .details = CommandFailure.plain(code, message, exit_code),
    } };
}

fn validatePairRecord(
    record: zfastq.Record,
    offsets: zfastq.RecordOffsets,
    record_index: u64,
    options: PairedCheckOptions,
) ?CommandFailure {
    return validateCommandRecord(record, offsets, record_index, options.alphabet);
}

fn lastRecordIndex(reader: *const zfastq.Reader) ?u64 {
    return if (reader.recordIndex() == 0) null else reader.recordIndex() - 1;
}

fn parsePairName(header: []const u8, policy: PairNamePolicy) PairName {
    var first_end: usize = 0;
    while (first_end < header.len and header[first_end] != ' ' and header[first_end] != '\t') {
        first_end += 1;
    }
    const first_token = header[0..first_end];
    if (policy == .exact) {
        return .{
            .first_token = first_token,
            .normalized_id = first_token,
            .mate_markers = 0,
            .first_mate_marker = null,
        };
    }

    var normalized_id = first_token;
    var mate_markers: u2 = 0;
    const first_mate_marker = terminalMateMarker(first_token);
    if (first_mate_marker) |mate| {
        mate_markers |= mateMask(mate);
        normalized_id = first_token[0 .. first_token.len - 2];
    }

    var second_start = first_end;
    while (second_start < header.len and
        (header[second_start] == ' ' or header[second_start] == '\t'))
    {
        second_start += 1;
    }
    var second_end = second_start;
    while (second_end < header.len and header[second_end] != ' ' and
        header[second_end] != '\t')
    {
        second_end += 1;
    }
    const second_token = header[second_start..second_end];
    if (second_token.len >= 2 and second_token[1] == ':' and
        (second_token[0] == '1' or second_token[0] == '2'))
    {
        mate_markers |= mateMask(@intCast(second_token[0] - '0'));
    }
    if (terminalMateMarker(second_token)) |mate| {
        mate_markers |= mateMask(mate);
    }
    return .{
        .first_token = first_token,
        .normalized_id = normalized_id,
        .mate_markers = mate_markers,
        .first_mate_marker = first_mate_marker,
    };
}

fn terminalMateMarker(token: []const u8) ?u2 {
    if (token.len < 2 or token[token.len - 2] != '/') return null;
    return switch (token[token.len - 1]) {
        '1' => 1,
        '2' => 2,
        else => null,
    };
}

fn mateMask(mate: u2) u2 {
    return if (mate == 1) 0b01 else 0b10;
}

fn pairNamesMatch(name1: PairName, name2: PairName) bool {
    if (!std.mem.eql(u8, name1.normalized_id, name2.normalized_id)) return false;
    if (name1.mate_markers == 0b11 or name2.mate_markers == 0b11) return false;
    if ((name1.mate_markers == 0) != (name2.mate_markers == 0)) return false;
    if (name1.mate_markers == 0) return true;
    return name1.mate_markers == 0b01 and name2.mate_markers == 0b10;
}

fn checkFile(
    io: std.Io,
    path: []const u8,
    options: CheckOptions,
) ?CommandFailure {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return CommandFailure.plain("io_error", "file not found", 3),
        else => return CommandFailure.plain("io_error", "failed to open file", 3),
    };

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);
    return checkSource(input.byteSource(), options);
}

fn checkStdin(
    io: std.Io,
    options: CheckOptions,
) ?CommandFailure {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);
    return checkSource(input.byteSource(), options);
}

fn checkSource(
    source: zfastq.io.ByteSource,
    options: CheckOptions,
) ?CommandFailure {
    var scanner = fastq.CheckScanner.init(
        .{ .max_line_bytes = options.max_line_bytes },
        .{ .alphabet = options.alphabet },
    );
    var buf: [zfastq.limits.COUNT_READ_BUFFER_BYTES]u8 = undefined;
    while (true) {
        const n = source.read(&buf) catch {
            return CommandFailure.plain("io_error", "I/O error", 3);
        };
        if (n == 0) break;
        _ = scanner.feed(buf[0..n]) catch |err| {
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
        error.Format => if (scanner.takeLastError()) |details|
            CommandFailure.lint(details)
        else
            CommandFailure.plain("io_error", "validation failed without details", 3),
        error.LineTooLong => CommandFailure.plain(
            "line_limit",
            "line length limit exceeded",
            4,
        ),
        error.ArithmeticLimit => CommandFailure.plain(
            "arithmetic_limit",
            "input location exceeds supported limit",
            4,
        ),
    };
}

// --- Stats command ---

// Keep stats code generation out of the measured count dispatch path.
noinline fn runStats(
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
        const outcome = if (std.mem.eql(u8, input, "-"))
            statsStdin(io, allocator, options)
        else
            statsFile(io, allocator, input, options);
        const stats = switch (outcome) {
            .success => |stats| stats,
            .failure => |failure| {
                if (std.mem.eql(u8, failure.code, "out_of_memory")) {
                    std.Io.File.writeStreamingAll(
                        .stderr(),
                        io,
                        "error: out of memory\n",
                    ) catch {};
                    return 3;
                }
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
        const outcome = if (std.mem.eql(u8, input, "-"))
            statsStdin(io, allocator, options)
        else
            statsFile(io, allocator, input, options);
        writeStatsJsonResult(&json, input, outcome) catch return 3;
        switch (outcome) {
            .success => {},
            .failure => |failure| exit_code = @max(exit_code, failure.exit_code),
        }
    }
    finishJsonDocument(&json) catch return 3;
    stdout_writer.interface.flush() catch return 3;
    return exit_code;
}

fn statsFile(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    options: InputOptions,
) StatsOutcome {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .{ .failure = CommandFailure.plain(
            "io_error",
            "file not found",
            3,
        ) },
        else => return .{ .failure = CommandFailure.plain(
            "io_error",
            "failed to open file",
            3,
        ) },
    };

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return .{ .failure = CommandFailure.plain("io_error", "I/O error", 3) };
    };
    defer input.deinit(io);
    return collectStats(allocator, input.byteSource(), options);
}

fn statsStdin(
    io: std.Io,
    allocator: std.mem.Allocator,
    options: InputOptions,
) StatsOutcome {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        return .{ .failure = CommandFailure.plain("io_error", "I/O error", 3) };
    };
    defer input.deinit(io);
    return collectStats(allocator, input.byteSource(), options);
}

const RecordInput = struct {
    file: std.Io.File,
    owns_file: bool,
    read_buffer: [64 * 1024]u8,
    file_reader: std.Io.File.Reader,
    source: union(enum) {
        plain: io_layer.ReaderSource,
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
        self.source = .{ .plain = io_layer.ReaderSource.init(&self.file_reader.interface) };
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
};

fn collectStats(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    options: InputOptions,
) StatsOutcome {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return .{ .failure = CommandFailure.plain(
        "out_of_memory",
        "out of memory",
        3,
    ) };
    defer reader.deinit();

    var stats: zfastq.Stats = .{};
    while (reader.next() catch |err| {
        return .{ .failure = mapReaderFailure(&reader, err) };
    }) |record| {
        stats.addRecord(record) catch |err| switch (err) {
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
                const relative_offset = std.math.cast(u64, quality_error.byte_index) orelse
                    return .{ .failure = CommandFailure.plain(
                        "arithmetic_limit",
                        "statistics arithmetic limit exceeded",
                        4,
                    ) };
                const byte_offset = std.math.add(
                    u64,
                    offsets.quality,
                    relative_offset,
                ) catch return .{ .failure = CommandFailure.plain(
                    "arithmetic_limit",
                    "statistics arithmetic limit exceeded",
                    4,
                ) };
                return .{ .failure = CommandFailure.lint(.{
                    .code = .s006_invalid_quality_range,
                    .message = INVALID_QUALITY_MESSAGE,
                    .record_index = reader.recordIndex() - 1,
                    .byte_offset = byte_offset,
                    .line_in_record = 4,
                }) };
            },
            error.S005LengthMismatch => {
                const offsets = reader.currentRecordOffsets() orelse {
                    return .{ .failure = CommandFailure.plain(
                        "io_error",
                        "record location is unavailable",
                        3,
                    ) };
                };
                return .{ .failure = CommandFailure.lint(.{
                    .code = .s005_length_mismatch,
                    .message = "sequence and quality lengths differ",
                    .record_index = reader.recordIndex() - 1,
                    .byte_offset = offsets.quality,
                    .line_in_record = 4,
                }) };
            },
            error.Overflow => return .{ .failure = CommandFailure.plain(
                "arithmetic_limit",
                "statistics arithmetic limit exceeded",
                4,
            ) },
        };
    }
    return .{ .success = stats };
}

fn mapReaderFailure(
    reader: *zfastq.Reader,
    err: zfastq.ReaderError,
) CommandFailure {
    return switch (err) {
        error.S001InvalidPlusLine,
        error.S003InvalidHeader,
        error.S004TruncatedRecord,
        error.S005LengthMismatch,
        => if (reader.takeLastError()) |details|
            CommandFailure.lint(details)
        else
            CommandFailure.plain("io_error", "validation failed without details", 3),
        error.LineTooLong => CommandFailure.plain(
            "line_limit",
            "line length limit exceeded",
            4,
        ),
        error.OutOfMemory => CommandFailure.plain("out_of_memory", "out of memory", 3),
        error.Io => CommandFailure.plain("io_error", "I/O error", 3),
    };
}

// --- Sample command ---

const SampleOptions = struct {
    max_line_bytes: usize,
    alphabet: zfastq.Alphabet,
    fraction: ?sampling.Fraction,
    count: ?u64,
    seed: u64,
};

const SampleMode = union(enum) {
    fraction: sampling.Fraction,
    count: u64,
};

noinline fn runSample(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: SampleOptions,
) u8 {
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

    var stdout_buffer: [64 * 1024]u8 = undefined;
    var stdout_writer = std.Io.File.stdout().writerStreaming(io, &stdout_buffer);
    var sink_adapter = io_layer.WriterSink.init(&stdout_writer.interface);
    var writer = zfastq.Writer.init(sink_adapter.byteSink());
    const failure = (switch (mode) {
        .fraction => |fraction| blk: {
            var selector = sampling.Selector.init(fraction, options.seed);
            break :blk if (std.mem.eql(u8, input, "-"))
                sampleFractionStdin(io, allocator, &writer, &selector, options)
            else
                sampleFractionFile(io, allocator, input, &writer, &selector, options);
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

fn sampleFractionFile(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    writer: *zfastq.Writer,
    selector: *sampling.Selector,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return CommandFailure.plain("io_error", "file not found", 3),
        else => return CommandFailure.plain("io_error", "failed to open file", 3),
    };

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);
    return sampleFractionSource(allocator, input.byteSource(), writer, selector, options);
}

fn sampleFractionStdin(
    io: std.Io,
    allocator: std.mem.Allocator,
    writer: *zfastq.Writer,
    selector: *sampling.Selector,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    var input: RecordInput = undefined;
    input.init(io, .stdin(), false) catch {
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);
    return sampleFractionSource(allocator, input.byteSource(), writer, selector, options);
}

fn sampleFractionSource(
    allocator: std.mem.Allocator,
    source: zfastq.io.ByteSource,
    writer: *zfastq.Writer,
    selector: *sampling.Selector,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    var reader = zfastq.Reader.init(
        allocator,
        source,
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return CommandFailure.plain("out_of_memory", "out of memory", 3);
    defer reader.deinit();

    while (true) switch (nextValidatedSampleRecord(&reader, options.alphabet)) {
        .done => return null,
        .failure => |failure| return failure,
        .record => |validated| {
            if (!selector.selectRecord()) continue;
            if (validated.canonical_span) |span| {
                fastq.writeCanonicalRecordSpan(writer, span) catch return error.WriteFailed;
            } else {
                fastq.writeValidatedRecord(writer, validated.record) catch
                    return error.WriteFailed;
            }
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
    alphabet: zfastq.Alphabet,
) SampleRecordOutcome {
    var canonical_span: ?[]const u8 = null;
    const record = fastq.nextWithoutId(reader, &canonical_span) catch |err| {
        return .{ .failure = mapReaderFailure(reader, err) };
    } orelse return .done;
    const offsets = reader.currentRecordOffsets() orelse return .{ .failure = CommandFailure.plain("io_error", "record location is unavailable", 3) };
    if (validateCommandRecord(
        record,
        offsets,
        reader.recordIndex() - 1,
        alphabet,
    )) |failure| {
        return .{ .failure = failure };
    }
    return .{ .record = .{
        .record = record,
        .canonical_span = canonical_span,
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
        record_count: u64,
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
    selector.finish(allocator) catch {
        return CommandFailure.plain("out_of_memory", "out of memory", 3);
    };
    if (selector.indexes.items.len == 0) return null;

    return sampleExactSecondPass(
        io,
        allocator,
        path,
        writer,
        selector.indexes.items,
        completed.snapshot,
        completed.record_count,
        options,
    );
}

fn sampleExactFirstPass(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    selector: *sampling.ExactSelector,
    options: SampleOptions,
) ExactFirstPass {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return .{ .failure = CommandFailure.plain(
            "io_error",
            "file not found",
            3,
        ) },
        else => return .{ .failure = CommandFailure.plain(
            "io_error",
            "failed to open file",
            3,
        ) },
    };
    const initial = fileSnapshot(file, io) catch |err| {
        file.close(io);
        return .{ .failure = switch (err) {
            error.NotRegularFile => CommandFailure.plain(
                "io_error",
                "exact-count sampling requires a regular file",
                3,
            ),
            else => CommandFailure.plain("io_error", "failed to inspect file", 3),
        } };
    };

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return .{ .failure = CommandFailure.plain("io_error", "I/O error", 3) };
    };
    defer input.deinit(io);

    var reader = zfastq.Reader.init(
        allocator,
        input.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return .{ .failure = CommandFailure.plain(
        "out_of_memory",
        "out of memory",
        3,
    ) };
    defer reader.deinit();

    var failure: ?CommandFailure = null;
    while (true) switch (nextValidatedSampleRecord(&reader, options.alphabet)) {
        .done => break,
        .failure => |details| {
            failure = details;
            break;
        },
        .record => selector.considerRecord(allocator, reader.recordIndex()) catch |err| {
            failure = switch (err) {
                error.OutOfMemory => CommandFailure.plain("out_of_memory", "out of memory", 3),
                error.Overflow => CommandFailure.plain(
                    "arithmetic_limit",
                    "sample index storage exceeds supported limit",
                    4,
                ),
            };
            break;
        },
    };

    const final = fileSnapshot(input.file, io) catch return .{ .failure = CommandFailure.plain("io_error", "failed to inspect file", 3) };
    if (!sameFileSnapshot(initial, final)) return .{ .failure = inputChangedFailure() };
    if (failure) |details| return .{ .failure = details };
    return .{ .success = .{
        .snapshot = initial,
        .record_count = reader.recordIndex(),
    } };
}

fn sampleExactSecondPass(
    io: std.Io,
    allocator: std.mem.Allocator,
    path: []const u8,
    writer: *zfastq.Writer,
    selected: []const u64,
    expected_snapshot: FileSnapshot,
    expected_count: u64,
    options: SampleOptions,
) error{WriteFailed}!?CommandFailure {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| switch (err) {
        error.FileNotFound => return inputChangedFailure(),
        else => return CommandFailure.plain("io_error", "failed to reopen file", 3),
    };
    const initial = fileSnapshot(file, io) catch |err| {
        file.close(io);
        return switch (err) {
            error.NotRegularFile => inputChangedFailure(),
            else => CommandFailure.plain("io_error", "failed to inspect file", 3),
        };
    };
    if (!sameFileSnapshot(expected_snapshot, initial)) {
        file.close(io);
        return inputChangedFailure();
    }

    var input: RecordInput = undefined;
    input.init(io, file, true) catch {
        file.close(io);
        return CommandFailure.plain("io_error", "I/O error", 3);
    };
    defer input.deinit(io);

    var reader = zfastq.Reader.init(
        allocator,
        input.byteSource(),
        .{ .max_line_bytes = options.max_line_bytes },
    ) catch return CommandFailure.plain("out_of_memory", "out of memory", 3);
    defer reader.deinit();

    var selected_cursor: usize = 0;
    var failure: ?CommandFailure = null;
    while (true) switch (nextValidatedSampleRecord(&reader, options.alphabet)) {
        .done => break,
        .failure => |details| {
            failure = details;
            break;
        },
        .record => |validated| {
            if (selected_cursor == selected.len or
                reader.recordIndex() != selected[selected_cursor]) continue;
            if (validated.canonical_span) |span| {
                fastq.writeCanonicalRecordSpan(writer, span) catch return error.WriteFailed;
            } else {
                fastq.writeValidatedRecord(writer, validated.record) catch
                    return error.WriteFailed;
            }
            selected_cursor += 1;
        },
    };

    const final = fileSnapshot(input.file, io) catch
        return CommandFailure.plain("io_error", "failed to inspect file", 3);
    if (!sameFileSnapshot(expected_snapshot, final)) return inputChangedFailure();
    if (failure) |details| return details;
    if (reader.recordIndex() != expected_count or selected_cursor != selected.len) {
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
    failure: ?CheckFailure,
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
    var run_start: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte >= 0x20 and byte <= 0x7e and byte != '\\') continue;

        try std.json.Stringify.encodeJsonStringChars(
            bytes[run_start..index],
            .{},
            json.writer,
        );
        if (byte == '\\') {
            try std.json.Stringify.encodeJsonStringChars("\\\\", .{}, json.writer);
        } else {
            const escaped = [4]u8{ '\\', 'x', HEX[byte >> 4], HEX[byte & 0x0f] };
            try std.json.Stringify.encodeJsonStringChars(&escaped, .{}, json.writer);
        }
        run_start = index + 1;
    }
    try std.json.Stringify.encodeJsonStringChars(bytes[run_start..], .{}, json.writer);
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
    try writeUnsignedField(io, "reads", result.reads);
    try writeUnsignedField(io, "bases", result.bases);
    try writeOptionalUnsignedField(io, "min_length", result.min_length);
    try writeOptionalUnsignedField(io, "max_length", result.max_length);
    try writeRatioField(io, "mean_length", result.bases, result.reads);
    try writeUnsignedField(io, "a", result.a);
    try writeUnsignedField(io, "c", result.c);
    try writeUnsignedField(io, "g", result.g);
    try writeUnsignedField(io, "t", result.t);
    try writeUnsignedField(io, "n", result.n);
    try writeUnsignedField(io, "other_bases", result.other_bases);
    try writeRatioField(
        io,
        "gc_fraction",
        @as(u128, result.g) + result.c,
        @as(u128, result.a) + result.c + result.g + result.t,
    );
    try writeUnsignedField(io, "quality_sum", result.quality_sum);
    try writeRatioField(io, "mean_quality", result.quality_sum, result.bases);
    try writeUnsignedField(io, "q20_bases", result.q20_bases);
    try writeRatioField(io, "q20_fraction", result.q20_bases, result.bases);
    try writeUnsignedField(io, "q30_bases", result.q30_bases);
    try writeRatioField(io, "q30_fraction", result.q30_bases, result.bases);
}

fn writeUnsignedField(io: std.Io, name: []const u8, value: u64) !void {
    var buf: [96]u8 = undefined;
    const line = try std.fmt.bufPrint(&buf, "{s}: {d}\n", .{ name, value });
    try std.Io.File.writeStreamingAll(.stdout(), io, line);
}

fn writeOptionalUnsignedField(io: std.Io, name: []const u8, value: ?u64) !void {
    if (value) |number| return writeUnsignedField(io, name, number);
    var buf: [64]u8 = undefined;
    const line = try std.fmt.bufPrint(&buf, "{s}: -\n", .{name});
    try std.Io.File.writeStreamingAll(.stdout(), io, line);
}

fn writeRatioField(
    io: std.Io,
    name: []const u8,
    numerator: u128,
    denominator: u128,
) !void {
    var buf: [96]u8 = undefined;
    if (denominator == 0) {
        const line = try std.fmt.bufPrint(&buf, "{s}: -\n", .{name});
        return std.Io.File.writeStreamingAll(.stdout(), io, line);
    }
    const scale = 1_000_000;
    const rounded = (numerator * scale + denominator / 2) / denominator;
    const line = try std.fmt.bufPrint(
        &buf,
        "{s}: {d}.{d:0>6}\n",
        .{ name, rounded / scale, rounded % scale },
    );
    try std.Io.File.writeStreamingAll(.stdout(), io, line);
}

fn printPathError(io: std.Io, path: []const u8, message: []const u8) void {
    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
}

fn printCommandFailure(io: std.Io, path: []const u8, failure: CommandFailure) void {
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

fn printCheckFailure(
    io: std.Io,
    inputs: []const []const u8,
    pair_mode: PairMode,
    failure: CheckFailure,
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

fn printParseError(io: std.Io, path: []const u8, details: zfastq.ParseError) void {
    std.Io.File.writeStreamingAll(.stderr(), io, "error: ") catch {};
    writeEscaped(.stderr(), io, path);
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, zfastq.codeTag(details.code)) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, ": ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, details.message) catch {};

    var buf: [96]u8 = undefined;
    const suffix = std.fmt.bufPrint(
        &buf,
        " (record {d}, line {d}, offset {d})\n",
        .{
            details.record_index,
            details.line_in_record,
            details.byte_offset,
        },
    ) catch return;
    std.Io.File.writeStreamingAll(.stderr(), io, suffix) catch {};
}

const HEX = "0123456789ABCDEF";

fn writeEscaped(file: std.Io.File, io: std.Io, bytes: []const u8) void {
    writeEscapedAll(file, io, bytes) catch {};
}

fn writeEscapedAll(file: std.Io.File, io: std.Io, bytes: []const u8) !void {
    var run_start: usize = 0;
    for (bytes, 0..) |byte, index| {
        if (byte >= 0x20 and byte <= 0x7e and byte != '\\') continue;

        try std.Io.File.writeStreamingAll(file, io, bytes[run_start..index]);
        if (byte == '\\') {
            try std.Io.File.writeStreamingAll(file, io, "\\\\");
        } else {
            const escaped = [4]u8{ '\\', 'x', HEX[byte >> 4], HEX[byte & 0x0f] };
            try std.Io.File.writeStreamingAll(file, io, &escaped);
        }
        run_start = index + 1;
    }
    try std.Io.File.writeStreamingAll(file, io, bytes[run_start..]);
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
    var path_buffer: [std.fs.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(
        &path_buffer,
        ".zig-cache/tmp/{s}/{s}",
        .{ tmp.sub_path, name },
    );
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    const snapshot = try fileSnapshot(file, io);
    file.close(io);

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactSecondPass(
        io,
        std.testing.allocator,
        path,
        &writer,
        &.{1},
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
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    var snapshot = try fileSnapshot(file, io);
    file.close(io);
    snapshot.size += 1;

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactSecondPass(
        io,
        std.testing.allocator,
        path,
        &writer,
        &.{1},
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

test "[failure] - [exact sample]: the second pass repeats semantic validation" {
    const io = std.testing.io;
    const path = "tests/data/synthetic/bad_alphabet.fastq";
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    const snapshot = try fileSnapshot(file, io);
    file.close(io);

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactSecondPass(
        io,
        std.testing.allocator,
        path,
        &writer,
        &.{1},
        snapshot,
        1,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .fraction = null,
            .count = 1,
            .seed = 11,
        },
    )).?;

    try std.testing.expectEqualStrings("S002", failure.code);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}
