//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("root.zig");
const fastq = @import("fastq.zig");
const io_layer = @import("io.zig");
const pairing = @import("pair.zig");
const sampling = @import("sample.zig");
const statistics = @import("stats.zig");
const qc = @import("qc.zig");
const preprocess = @import("preprocess.zig");

pub const std_options: std.Options = .{ .networking = false };

const Command = enum { count, stats, check, sample, interleave, deinterleave, process };

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

fn argumentError(io: std.Io, message: []const u8, exit_code: u8) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
    std.process.exit(exit_code);
}

fn printUsageAndExit(io: std.Io) noreturn {
    argumentError(io, USAGE, 2);
}

fn writeStderr(io: std.Io, message: []const u8) void {
    std.Io.File.writeStreamingAll(.stderr(), io, message) catch {};
}

fn usageError(io: std.Io, message: []const u8) u8 {
    writeStderr(io, message);
    return 2;
}

fn printHelpAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, USAGE) catch std.process.exit(stdoutFailed(io, 0));
    std.process.exit(0);
}

const GivenOptions = struct {
    max_line_bytes: bool = false,
    alphabet: bool = false,
    fraction: bool = false,
    count: bool = false,
    seed: bool = false,
    pair_names: bool = false,
    out1: bool = false,
    out2: bool = false,
    report: bool = false,
    contaminants: bool = false,
    max_pairs: bool = false,
    adapter1: bool = false,
    adapter2: bool = false,
    trim_poly_g: bool = false,
    trim_tail_window: bool = false,
    trim_tail_quality: bool = false,
    min_length: bool = false,
    max_n: bool = false,
    low_quality: bool = false,
    max_low_quality_percent: bool = false,
};

fn optionValue(
    io: std.Io,
    args: *std.process.Args.Iterator,
    given: *bool,
    comptime option: []const u8,
) []const u8 {
    const value = args.next() orelse argumentError(io, "error: " ++ option ++ " requires a value\n", 2);
    if (given.*) argumentError(io, "error: " ++ option ++ " may appear only once\n", 2);
    given.* = true;
    return value;
}

fn adapterOption(
    io: std.Io,
    args: *std.process.Args.Iterator,
    given: *bool,
    comptime option: []const u8,
) []const u8 {
    const message = std.fmt.comptimePrint("error: {s} must be {d} to {d} bases of A, C, G, and T\n", .{
        option,
        preprocess.MIN_ADAPTER_BASES,
        preprocess.MAX_ADAPTER_BASES,
    });
    const value = optionValue(io, args, given, option);
    if (value.len < preprocess.MIN_ADAPTER_BASES or value.len > preprocess.MAX_ADAPTER_BASES) {
        argumentError(io, message, 2);
    }
    for (value) |base| switch (base) {
        'A', 'C', 'G', 'T' => {},
        else => argumentError(io, message, 2),
    };
    return value;
}

fn rangeOption(
    io: std.Io,
    args: *std.process.Args.Iterator,
    given: *bool,
    comptime T: type,
    comptime option: []const u8,
    comptime min: T,
    comptime max: T,
) T {
    const message = std.fmt.comptimePrint("error: {s} must be {d} to {d}\n", .{ option, min, max });
    const value = sampling.parseDecimal(T, optionValue(io, args, given, option)) catch argumentError(io, message, 2);
    if (value < min or value > max) argumentError(io, message, 2);
    return value;
}

const ReadOptionValues = struct {
    adapters: [2]?[]const u8 = .{ null, null },
    poly_g: ?u32 = null,
    tail_window: ?u16 = null,
    tail_quality: ?u8 = null,
    min_length: ?u32 = null,
    max_n: ?u8 = null,
    low_quality: ?u8 = null,
    max_low_quality_percent: ?u8 = null,
};

fn processReadOptions(io: std.Io, values: ReadOptionValues) preprocess.Options {
    if ((values.adapters[0] == null) != (values.adapters[1] == null)) {
        argumentError(io, "error: --adapter1 and --adapter2 must be given together\n", 2);
    }
    if ((values.tail_window == null) != (values.tail_quality == null)) {
        argumentError(io, "error: --trim-tail-window and --trim-tail-quality must be given together\n", 2);
    }
    if ((values.low_quality == null) != (values.max_low_quality_percent == null)) {
        argumentError(io, "error: --low-quality and --max-low-quality-percent must be given together\n", 2);
    }
    return .{
        .adapters = if (values.adapters[0]) |adapter1| .{ adapter1, values.adapters[1].? } else null,
        .poly_g = values.poly_g,
        .tail = if (values.tail_window) |window| .{ .window = window, .quality = values.tail_quality.? } else null,
        .min_length = values.min_length,
        .max_n = values.max_n,
        .low_quality = if (values.low_quality) |quality| .{ .quality = quality, .percent = values.max_low_quality_percent.? } else null,
    };
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
    else if (std.mem.eql(u8, cmd, "process"))
        .process
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
    var report_path: ?[]const u8 = null;
    var contaminants_path: ?[]const u8 = null;
    var pair_duplication = false;
    var gzip_output = false;
    var max_pairs: ?u64 = null;
    var read_options: ReadOptionValues = .{};
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
            const value = optionValue(io, &args, &given.max_line_bytes, "--max-line-bytes");
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
            command == .deinterleave or command == .process) and
            std.mem.eql(u8, arg, "--alphabet"))
        {
            const value = optionValue(io, &args, &given.alphabet, "--alphabet");
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
            const value = optionValue(io, &args, &given.fraction, "--fraction");
            fraction = sampling.Fraction.parse(value) catch {
                argumentError(io, "error: invalid --fraction value\n", 2);
            };
            continue;
        }
        if (command == .sample and std.mem.eql(u8, arg, "--count")) {
            const value = optionValue(io, &args, &given.count, "--count");
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
            const value = optionValue(io, &args, &given.seed, "--seed");
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
            command == .deinterleave or command == .process) and
            std.mem.eql(u8, arg, "--pair-names"))
        {
            const value = optionValue(io, &args, &given.pair_names, "--pair-names");
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
        if ((command == .deinterleave or command == .process) and
            (std.mem.eql(u8, arg, "--out1") or std.mem.eql(u8, arg, "--out2")))
        {
            if (std.mem.eql(u8, arg, "--out1")) {
                output1 = optionValue(io, &args, &given.out1, "--out1");
            } else {
                output2 = optionValue(io, &args, &given.out2, "--out2");
            }
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--gzip")) {
            gzip_output = true;
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--report")) {
            report_path = optionValue(io, &args, &given.report, "--report");
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--contaminants")) {
            contaminants_path = optionValue(io, &args, &given.contaminants, "--contaminants");
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--pair-duplication")) {
            pair_duplication = true;
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--max-pairs")) {
            const value = optionValue(io, &args, &given.max_pairs, "--max-pairs");
            max_pairs = sampling.parseDecimal(u64, value) catch |err| switch (err) {
                error.Overflow => {
                    argumentError(io, "error: --max-pairs exceeds supported limit\n", 4);
                },
                error.InvalidCharacter => {
                    argumentError(io, "error: invalid --max-pairs value\n", 2);
                },
            };
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--adapter1")) {
            read_options.adapters[0] = adapterOption(io, &args, &given.adapter1, "--adapter1");
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--adapter2")) {
            read_options.adapters[1] = adapterOption(io, &args, &given.adapter2, "--adapter2");
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--trim-poly-g")) {
            read_options.poly_g = rangeOption(io, &args, &given.trim_poly_g, u32, "--trim-poly-g", 1, std.math.maxInt(u32));
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--trim-tail-window")) {
            read_options.tail_window = rangeOption(io, &args, &given.trim_tail_window, u16, "--trim-tail-window", 1, preprocess.MAX_TAIL_WINDOW);
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--trim-tail-quality")) {
            read_options.tail_quality = rangeOption(io, &args, &given.trim_tail_quality, u8, "--trim-tail-quality", 1, preprocess.MAX_TAIL_QUALITY);
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--min-length")) {
            read_options.min_length = rangeOption(io, &args, &given.min_length, u32, "--min-length", 0, std.math.maxInt(u32));
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--max-n")) {
            read_options.max_n = rangeOption(io, &args, &given.max_n, u8, "--max-n", 0, preprocess.MAX_N_LIMIT);
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--low-quality")) {
            read_options.low_quality = rangeOption(io, &args, &given.low_quality, u8, "--low-quality", 0, preprocess.MAX_LOW_QUALITY);
            continue;
        }
        if (command == .process and std.mem.eql(u8, arg, "--max-low-quality-percent")) {
            read_options.max_low_quality_percent = rangeOption(io, &args, &given.max_low_quality_percent, u8, "--max-low-quality-percent", 0, preprocess.MAX_LOW_QUALITY_PERCENT);
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
    if (command == .process and report_path == null) {
        if (contaminants_path != null) argumentError(io, "error: --contaminants requires --report\n", 2);
        if (pair_duplication) argumentError(io, "error: --pair-duplication requires --report\n", 2);
    }
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
        .process => runProcess(io, gpa, positional.items, .{
            .max_line_bytes = max_line_bytes,
            .alphabet = alphabet,
            .pair_name_policy = pair_name_policy,
            .output1 = output1,
            .output2 = output2,
            .report = report_path,
            .contaminants = contaminants_path,
            .pair_duplication = pair_duplication,
            .gzip = gzip_output,
            .max_pairs = max_pairs,
            .trim = processReadOptions(io, read_options),
        }),
    };
    std.process.exit(code);
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

const STATISTICS_LIMIT = CommandFailure.plain("arithmetic_limit", "statistics arithmetic limit exceeded", 4);

const INSPECT_FAILURE = CommandFailure.plain("io_error", "failed to inspect file", 3);

const RECORD_LOCATION_UNAVAILABLE = CommandFailure.plain("io_error", "record location is unavailable", 3);

const OUT_OF_MEMORY_LINE = "error: " ++ OUT_OF_MEMORY.message ++ "\n";

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

const InputOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
};

const FileIdentity = struct {
    device_major: u32,
    device_minor: u32,
    inode: u64,

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
        reader.buffer = &self.read_buffer;
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

fn inputTypeFailure(io: std.Io, file: std.Io.File) ?CommandFailure {
    const stat = statRecordFile(io, file) catch return null;
    if (stat.mode & std.os.linux.S.IFMT == std.os.linux.S.IFDIR) {
        return CommandFailure.plain("io_error", "input is a directory", 3);
    }
    return null;
}

fn namesStdin(label: []const u8) bool {
    return std.mem.eql(u8, label, "-");
}

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
        error.BadFileDescriptor => return null,
        else => |other| return other,
    };
    if (stat.mode & std.os.linux.S.IFMT != std.os.linux.S.IFREG) return null;
    return .of(stat);
}

fn checkedStdoutIdentity(io: std.Io) error{Reported}!?FileIdentity {
    return stdoutFileIdentity(io) catch {
        if (canReportInputFailure(io, null)) printOutputFailure(io);
        return error.Reported;
    };
}

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
                return usageError(io, "error: standard input may appear at most once\n");
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
        error.S005LengthMismatch => @panic("Reader returned unequal sequence and quality lengths"),
        error.Overflow => return STATISTICS_LIMIT,
    }
}

const CheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
};

const PairedCheckOptions = struct {
    max_line_bytes: usize = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    alphabet: zfastq.Alphabet = .iupac,
    pair_mode: PairMode = .none,
    pair_name_policy: pairing.NamePolicy = .illumina,
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
        return usageError(io, "error: --pair-names requires --paired or --interleaved\n");
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
        return usageError(io, message);
    }
    if (pair_mode == .paired and stdinRepeats(inputs)) {
        return usageError(io, "error: paired inputs may contain standard input at most once\n");
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
    end,
    pair,
    failure: PairCommandFailure,
};

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

const ExactUnit = enum { done, end, failed };

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
            return usageError(io, "error: --fraction and --count are mutually exclusive\n");
        }
        break :blk .{ .fraction = fraction };
    } else if (options.count) |count|
        .{ .count = count }
    else {
        return usageError(io, "error: sample requires --fraction P or --count K\n");
    };
    if (options.pair_mode != .none) {
        return runPairSampleCommand(io, allocator, inputs, mode, options);
    }
    if (options.pair_names_set) {
        return usageError(io, "error: --pair-names requires --paired or --interleaved\n");
    }
    if (inputs.len != 1) {
        return usageError(io, "error: sample requires exactly one input\n");
    }
    const input = inputs[0];
    if (mode == .count and namesStdin(input)) {
        return usageError(io, "error: exact-count sampling requires a file path\n");
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
        return usageError(io, message);
    }
    switch (mode) {
        .count => {
            for (inputs) |input| {
                if (!namesStdin(input)) continue;
                return usageError(io, "error: paired exact-count sampling requires file paths\n");
            }
        },
        .fraction => if (options.pair_mode == .paired and stdinRepeats(inputs)) {
            return usageError(io, "error: paired sample inputs may contain standard input at most once\n");
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
    if (inputs.len != 2) return usageError(io, "error: interleave requires exactly two inputs\n");
    if (stdinRepeats(inputs)) {
        return usageError(io, "error: interleave inputs may contain standard input at most once\n");
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

const DeinterleaveOptions = struct {
    max_line_bytes: usize,
    alphabet: zfastq.Alphabet,
    pair_name_policy: pairing.NamePolicy,
};

const DeinterleaveWriteError = error{
    Output1WriteFailed,
    Output2WriteFailed,
};

const RecordOutput = struct {
    path: []const u8,
    file: ?std.Io.File = null,
    write_buffer: [RECORD_OUTPUT_BUFFER_BYTES]u8 = undefined,
    sink: io_layer.FileSink = undefined,
    gzip: ?*io_layer.GzipSink = null,
    writer: zfastq.Writer = undefined,

    fn init(path: []const u8) RecordOutput {
        return .{ .path = path };
    }

    fn create(self: *RecordOutput, io: std.Io) ?CommandFailure {
        if (createExclusiveOutput(io, self.path, &self.file)) |failure| return failure;
        const file = self.file.?;
        self.sink = io_layer.FileSink.init(io, file, &self.write_buffer);
        if (self.gzip) |gzip| {
            gzip.init(&self.sink.writer) catch return IO_FAILURE;
            self.writer = zfastq.Writer.init(gzip.byteSink());
        } else {
            self.writer = zfastq.Writer.init(self.sink.byteSink());
        }
        return null;
    }

    fn finish(self: *RecordOutput) error{WriteFailed}!void {
        if (self.gzip) |gzip| try gzip.finish();
        self.sink.writer.flush() catch return error.WriteFailed;
    }

    fn close(self: *RecordOutput, io: std.Io) void {
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
    if (inputTypeFailure(io, input_file)) |failure| {
        printCommandFailure(io, input_label, failure);
        return failure.exit_code;
    }

    var output1 = RecordOutput.init(paths[0]);
    defer output1.close(io);
    var output2 = RecordOutput.init(paths[1]);
    defer output2.close(io);
    if (output1.create(io)) |failure| {
        printCommandFailure(io, output1.path, failure);
        return failure.exit_code;
    }
    if (output2.create(io)) |failure| {
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
        writeStderr(io, "error: deinterleave requires exactly one input\n");
        return null;
    }
    const path1 = output1_path orelse {
        writeStderr(io, "error: deinterleave requires --out1\n");
        return null;
    };
    const path2 = output2_path orelse {
        writeStderr(io, "error: deinterleave requires --out2\n");
        return null;
    };
    if (namesStdin(path1) or namesStdin(path2)) {
        writeStderr(io, "error: deinterleave output paths cannot be standard output\n");
        return null;
    }
    if (std.mem.eql(u8, path1, path2)) {
        writeStderr(io, "error: deinterleave output paths must differ\n");
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

const ProcessOptions = struct {
    max_line_bytes: usize,
    alphabet: zfastq.Alphabet,
    pair_name_policy: pairing.NamePolicy,
    output1: ?[]const u8,
    output2: ?[]const u8,
    report: ?[]const u8,
    contaminants: ?[]const u8 = null,
    pair_duplication: bool = false,
    gzip: bool,
    max_pairs: ?u64,
    trim: preprocess.Options,
};

const ProcessCounts = struct {
    read: u64 = 0,
    written: u64 = 0,
    dropped: u64 = 0,
    input_complete: bool = false,
};

fn runProcess(
    io: std.Io,
    allocator: std.mem.Allocator,
    inputs: []const []const u8,
    options: ProcessOptions,
) u8 {
    const paths = validateProcessArguments(io, inputs, options) orelse return 2;
    const report_to_stdout = if (options.report) |path| namesStdin(path) else false;
    const output_identity: ?FileIdentity = if (report_to_stdout)
        checkedStdoutIdentity(io) catch return 3
    else
        null;

    var input1: RecordInput = undefined;
    var input2: RecordInput = undefined;
    if (initPairedRecordInputs(&input1, &input2, io, inputs, output_identity)) |failure| {
        printPairCommandFailure(io, inputs, .paired, failure);
        return failure.exitCode();
    }
    defer input1.deinit(io);
    defer input2.deinit(io);
    const snapshots = [2]?FileSnapshot{
        fileSnapshot(input1.file, io) catch null,
        fileSnapshot(input2.file, io) catch null,
    };

    var contaminants: ?qc.ContaminantList = null;
    defer if (contaminants) |*list| list.deinit();
    if (options.contaminants) |path| {
        if (loadContaminants(io, allocator, path, &contaminants)) |code| return code;
    }
    var pair_duplication: preprocess.PairDuplication = undefined;
    std.debug.assert(!options.pair_duplication or options.report != null);
    if (options.pair_duplication) pair_duplication.init(allocator) catch {
        writeStderr(io, OUT_OF_MEMORY_LINE);
        return OUT_OF_MEMORY.exit_code;
    };
    defer if (options.pair_duplication) pair_duplication.deinit(allocator);

    var gzip_sinks: [2]?*io_layer.GzipSink = .{ null, null };
    defer for (gzip_sinks) |slot| if (slot) |sink| allocator.destroy(sink);
    if (options.gzip) for (&gzip_sinks) |*slot| {
        slot.* = allocator.create(io_layer.GzipSink) catch {
            writeStderr(io, OUT_OF_MEMORY_LINE);
            return OUT_OF_MEMORY.exit_code;
        };
    };

    var outputs = [2]RecordOutput{ .init(paths[0]), .init(paths[1]) };
    defer for (&outputs) |*output| output.close(io);
    for (&outputs, gzip_sinks, 0..) |*output, gzip_sink, index| {
        output.gzip = gzip_sink;
        if (output.create(io)) |failure| {
            removeCreatedOutputs(io, outputs[0..index]);
            printCommandFailure(io, output.path, failure);
            return failure.exit_code;
        }
    }
    var report_file: ?std.Io.File = null;
    defer if (report_file) |file| file.close(io);
    if (options.report) |path| if (!report_to_stdout) {
        if (createExclusiveOutput(io, path, &report_file)) |failure| {
            removeCreatedOutputs(io, &outputs);
            printCommandFailure(io, path, failure);
            return failure.exit_code;
        }
    };

    var collectors: ProcessCollectors = .{ .pair_duplication = if (options.pair_duplication) &pair_duplication else null };
    defer collectors.deinit(allocator);
    var counts: ProcessCounts = .{};
    var failed_output: ?u1 = null;
    const processed = processSources(
        allocator,
        .{ input1.stream(), input2.stream() },
        .{ &outputs[0].writer, &outputs[1].writer },
        &collectors,
        &counts,
        options,
    ) catch |err| failure: {
        failed_output = if (err == error.Output1WriteFailed) 0 else 1;
        break :failure null;
    };

    var input_failure = processed;
    if (collectors.pair_duplication) |pairs| pairs.flush(allocator) catch |err| {
        if (input_failure == null) input_failure = pairCommandFailure(0, switch (err) {
            error.OutOfMemory => OUT_OF_MEMORY,
            error.Overflow => STATISTICS_LIMIT,
        });
    };
    const changes_reads = options.trim.changesReads();
    var stats: [2][2]zfastq.Stats = @splat(.{ .{}, .{} });
    if (options.report != null) finishCollectors(allocator, &collectors.qc, changes_reads, &stats, &input_failure);

    var exit_code: u8 = 0;
    if (input_failure) |failure| {
        printPairCommandFailure(io, inputs, .paired, failure);
        exit_code = failure.exitCode();
    }
    if (failed_output) |index| {
        printPathError(io, outputs[index].path, IO_FAILURE.message);
        exit_code = @max(exit_code, IO_FAILURE.exit_code);
    }
    for (&outputs, 0..) |*output, index| {
        if (failed_output != null and failed_output.? == index) continue;
        output.finish() catch {
            printPathError(io, output.path, IO_FAILURE.message);
            exit_code = @max(exit_code, IO_FAILURE.exit_code);
            if (failed_output == null) failed_output = @intCast(index);
        };
    }

    const report_path = options.report orelse return exit_code;
    const report: ProcessReport = .{
        .inputs = inputs,
        .snapshots = snapshots,
        .outputs = paths,
        .options = options,
        .status_ok = exit_code == 0,
        .counts = counts,
        .stats = &stats,
        .qc = &collectors.qc,
        .changes_reads = changes_reads,
        .trimming = &collectors.trimming,
        .filtering = &collectors.filtering,
        .inserts = &collectors.inserts,
        .contaminants = if (contaminants) |*list| list else null,
        .pair_duplication = collectors.pair_duplication,
        .input_failure = input_failure,
        .failed_output = if (failed_output) |index| paths[index] else null,
    };
    var report_buffer: [JSON_OUTPUT_BUFFER_BYTES]u8 = undefined;
    var report_writer = if (report_file) |file|
        file.writerStreaming(io, &report_buffer)
    else
        std.Io.File.stdout().writerStreaming(io, &report_buffer);
    var json: std.json.Stringify = .{ .writer = &report_writer.interface };
    writeProcessReport(&json, report) catch return reportFailed(io, report_path, report_to_stdout, exit_code);
    report_writer.interface.flush() catch return reportFailed(io, report_path, report_to_stdout, exit_code);
    return exit_code;
}

const ProcessCollectors = struct {
    qc: [2]qc.MatePopulations = .{ .{}, .{} },
    trimming: [2]preprocess.MateCounts = .{ .{}, .{} },
    filtering: preprocess.FilterCounts = .{},
    inserts: preprocess.InsertCounts = .{},
    pair_duplication: ?*preprocess.PairDuplication = null,

    fn deinit(self: *ProcessCollectors, allocator: std.mem.Allocator) void {
        for (&self.qc) |*populations| populations.deinit(allocator);
    }
};

fn finishCollectors(
    allocator: std.mem.Allocator,
    collectors: *[2]qc.MatePopulations,
    changes_reads: bool,
    stats: *[2][2]zfastq.Stats,
    input_failure: *?PairCommandFailure,
) void {
    for (collectors, 0..) |*populations, index| {
        populations.finish(allocator) catch {
            if (input_failure.* == null) input_failure.* = pairCommandFailure(@intCast(index), OUT_OF_MEMORY);
        };
        const retained = if (changes_reads) &populations.retained else &populations.raw;
        for ([_]*const qc.MateQc{ &populations.raw, retained }, stats) |collector, *population_stats| {
            population_stats[index] = collector.stats() catch |err| stats: {
                if (input_failure.* == null) input_failure.* = pairCommandFailure(@intCast(index), qcFailure(err));
                break :stats .{};
            };
        }
    }
}

fn loadContaminants(io: std.Io, allocator: std.mem.Allocator, path: []const u8, list: *?qc.ContaminantList) ?u8 {
    const text = std.Io.Dir.cwd().readFileAlloc(io, path, allocator, .limited(qc.CONTAMINANT_FILE_LIMIT)) catch |err| {
        const failure: CommandFailure = switch (err) {
            error.OutOfMemory => OUT_OF_MEMORY,
            error.StreamTooLong => CommandFailure.plain("file_limit", "contaminant list exceeds 1 MiB", 4),
            error.FileNotFound => CommandFailure.plain("io_error", "file not found", 3),
            else => CommandFailure.plain("io_error", "failed to read file", 3),
        };
        printPathError(io, path, failure.message);
        return failure.exit_code;
    };
    defer allocator.free(text);
    var failed_line: usize = 0;
    list.* = qc.ContaminantList.parse(allocator, text, &failed_line) catch |err| switch (err) {
        error.OutOfMemory => {
            writeStderr(io, OUT_OF_MEMORY_LINE);
            return OUT_OF_MEMORY.exit_code;
        },
        error.InvalidBase => {
            var buffer: [96]u8 = undefined;
            printPathError(io, path, std.fmt.bufPrint(&buffer, "line {d}: contaminant sequence has a base other than A, C, G or T", .{failed_line}) catch unreachable);
            return 2;
        },
    };
    return null;
}

fn qcFailure(err: (error{OutOfMemory} || statistics.StatsError)) CommandFailure {
    return switch (err) {
        error.OutOfMemory => OUT_OF_MEMORY,
        error.Overflow => STATISTICS_LIMIT,
        error.S005LengthMismatch, error.S006InvalidQuality => CommandFailure.plain("io_error", "quality validation failed without details", 3),
    };
}

fn reportFailed(io: std.Io, path: []const u8, to_stdout: bool, exit_code: u8) u8 {
    if (to_stdout) return stdoutFailed(io, exit_code);
    printPathError(io, path, IO_FAILURE.message);
    return @max(exit_code, IO_FAILURE.exit_code);
}

fn validateProcessArguments(
    io: std.Io,
    inputs: []const []const u8,
    options: ProcessOptions,
) ?[2][]const u8 {
    if (inputs.len != 2) {
        writeStderr(io, "error: process requires exactly two inputs\n");
        return null;
    }
    if (stdinRepeats(inputs)) {
        writeStderr(io, "error: process inputs may contain standard input at most once\n");
        return null;
    }
    const path1 = options.output1 orelse {
        writeStderr(io, "error: process requires --out1\n");
        return null;
    };
    const path2 = options.output2 orelse {
        writeStderr(io, "error: process requires --out2\n");
        return null;
    };
    if (namesStdin(path1) or namesStdin(path2)) {
        writeStderr(io, "error: process output paths cannot be standard output\n");
        return null;
    }
    if (!options.gzip) {
        if (std.mem.endsWith(u8, path1, ".gz")) {
            writeStderr(io, "error: --out1 names a .gz file but --gzip is not given\n");
            return null;
        }
        if (std.mem.endsWith(u8, path2, ".gz")) {
            writeStderr(io, "error: --out2 names a .gz file but --gzip is not given\n");
            return null;
        }
    }
    const report = options.report orelse "-";
    if (std.mem.eql(u8, path1, path2) or
        std.mem.eql(u8, report, path1) or
        std.mem.eql(u8, report, path2))
    {
        writeStderr(io, "error: process output paths must differ\n");
        return null;
    }
    return .{ path1, path2 };
}

fn createExclusiveOutput(io: std.Io, path: []const u8, created: *?std.Io.File) ?CommandFailure {
    created.* = std.Io.Dir.cwd().createFile(io, path, .{ .exclusive = true }) catch |err| {
        return if (err == error.PathAlreadyExists)
            CommandFailure.plain("io_error", "output path already exists", 3)
        else
            CommandFailure.plain("io_error", "failed to create output", 3);
    };
    return null;
}

fn removeCreatedOutputs(io: std.Io, outputs: []RecordOutput) void {
    for (outputs) |*output| {
        output.close(io);
        std.Io.Dir.cwd().deleteFile(io, output.path) catch {};
    }
}

const ProcessPairs = struct {
    allocator: std.mem.Allocator,
    reader1: *zfastq.Reader,
    reader2: *zfastq.Reader,
    validator1: *fastq.AdaptiveRecordValidator,
    validator2: *fastq.AdaptiveRecordValidator,
    writers: [2]*zfastq.Writer,
    collectors: *ProcessCollectors,
    counts: *ProcessCounts,
    max_pairs: ?u64,
    policy: pairing.NamePolicy,
    trimmer: preprocess.Trimmer,
    collect_qc: bool,
    trim_failure: ?PairCommandFailure = null,

    fn atCap(self: *const ProcessPairs) bool {
        const cap = self.max_pairs orelse return false;
        return self.counts.read == cap;
    }

    inline fn emit(
        self: *ProcessPairs,
        records: [2]zfastq.Record,
        spans: [2]?[]const u8,
    ) DeinterleaveWriteError!?PairCommandFailure {
        const pairs = self.collectors.pair_duplication orelse return self.emitPair(records, spans);
        const hash = pairs.hash(records[0].sequence, records[1].sequence);
        if (try self.emitPair(records, spans)) |failure| return failure;
        pairs.push(self.allocator, hash) catch |err| return pairCommandFailure(0, switch (err) {
            error.OutOfMemory => OUT_OF_MEMORY,
            error.Overflow => STATISTICS_LIMIT,
        });
        return null;
    }

    inline fn emitPair(
        self: *ProcessPairs,
        records: [2]zfastq.Record,
        spans: [2]?[]const u8,
    ) DeinterleaveWriteError!?PairCommandFailure {
        self.counts.read += 1;
        if (self.collect_qc) for (&self.collectors.qc, records, 0..) |*collector, record, index| {
            collector.addRead(self.allocator, record.sequence, record.quality) catch |err| {
                return pairCommandFailure(@intCast(index), qcFailure(err));
            };
        };
        if (self.trimmer.options.changesReads()) {
            self.emitTransformed(records, spans) catch |err| switch (err) {
                error.TrimFailed => return self.trim_failure,
                error.Output1WriteFailed, error.Output2WriteFailed => |write_error| return write_error,
            };
            return null;
        }
        writeCheckedRecord(self.writers[0], records[0], spans[0]) catch return error.Output1WriteFailed;
        writeCheckedRecord(self.writers[1], records[1], spans[1]) catch return error.Output2WriteFailed;
        self.counts.written += 1;
        return null;
    }

    fn emitTransformed(
        self: *ProcessPairs,
        records: [2]zfastq.Record,
        spans: [2]?[]const u8,
    ) (DeinterleaveWriteError || error{TrimFailed})!void {
        const result = self.trimmer.trimPair(
            .{ records[0].sequence, records[1].sequence },
            .{ records[0].quality, records[1].quality },
        );
        self.collectors.inserts.add(&result) catch {
            self.trim_failure = pairCommandFailure(0, STATISTICS_LIMIT);
            return error.TrimFailed;
        };
        for (&self.collectors.trimming, result.removed, 0..) |*counts, removed, index| {
            counts.add(removed) catch {
                self.trim_failure = pairCommandFailure(@intCast(index), STATISTICS_LIMIT);
                return error.TrimFailed;
            };
        }
        const decision = preprocess.decide(
            self.trimmer.options,
            &result,
            .{ records[0].sequence[0..result.kept[0]], records[1].sequence[0..result.kept[1]] },
            .{ records[0].quality[0..result.kept[0]], records[1].quality[0..result.kept[1]] },
        );
        if (decision.pair != null) {
            self.collectors.filtering.add(decision) catch {
                self.trim_failure = pairCommandFailure(0, STATISTICS_LIMIT);
                return error.TrimFailed;
            };
            if (self.collect_qc) for (&self.collectors.qc, records, 0..) |*collector, record, index| {
                collector.dropRetained(self.allocator, record.sequence, record.quality) catch |err| {
                    self.trim_failure = pairCommandFailure(@intCast(index), qcFailure(err));
                    return error.TrimFailed;
                };
            };
            self.counts.dropped += 1;
            return;
        }
        var trimmed = records;
        var trimmed_spans = spans;
        for (&trimmed, &trimmed_spans, result.kept, 0..) |*record, *span, kept, index| {
            if (self.collect_qc) self.collectors.qc[index].addRetained(self.allocator, record.sequence, record.quality, kept) catch |err| {
                self.trim_failure = pairCommandFailure(@intCast(index), qcFailure(err));
                return error.TrimFailed;
            };
            if (kept != record.sequence.len) span.* = null;
            record.sequence = record.sequence[0..kept];
            record.quality = record.quality[0..kept];
        }
        writeCheckedRecord(self.writers[0], trimmed[0], trimmed_spans[0]) catch return error.Output1WriteFailed;
        writeCheckedRecord(self.writers[1], trimmed[1], trimmed_spans[1]) catch return error.Output2WriteFailed;
        self.counts.written += 1;
    }
};

noinline fn processSpanPairs(
    pairs: *ProcessPairs,
    pending1: *?fastq.RecordSpan,
) DeinterleaveWriteError!?PairCommandFailure {
    while (!pairs.atCap()) {
        const mate1 = fastq.nextRecordSpan(pairs.reader1, pairs.validator1) orelse return null;
        const mate2 = fastq.nextRecordSpan(pairs.reader2, pairs.validator2) orelse {
            pending1.* = mate1;
            return null;
        };
        if (!pairing.headersMatch(mate1.header(), mate2.header(), pairs.policy)) {
            const record_index1 = pairs.reader1.recordIndex() - 1;
            return headerNameMismatch(
                mate1.header(),
                mate2.header(),
                pairs.policy,
                record_index1,
                .{ record_index1, pairs.reader2.recordIndex() - 1 },
                .{ pairs.reader1.currentRecordOffsets().?.header, pairs.reader2.currentRecordOffsets().?.header },
            );
        }
        if (!spanIsWritable(mate1)) return pairCommandFailure(0, unwritableRecordFailure());
        if (!spanIsWritable(mate2)) return pairCommandFailure(1, unwritableRecordFailure());
        if (try pairs.emit(
            .{ mate1.record(), mate2.record() },
            .{ if (mate1.canonical) mate1.bytes else null, if (mate2.canonical) mate2.bytes else null },
        )) |failure| return failure;
    }
    return null;
}

fn processSources(
    allocator: std.mem.Allocator,
    sources: [2]*std.Io.Reader,
    writers: [2]*zfastq.Writer,
    collectors: *ProcessCollectors,
    counts: *ProcessCounts,
    options: ProcessOptions,
) DeinterleaveWriteError!?PairCommandFailure {
    if (options.report != null and options.trim.changesReads()) {
        for (&collectors.qc) |*populations| populations.enableRetained();
    }
    const reader_options: fastq.ReaderOptions = .{ .max_line_bytes = options.max_line_bytes };
    var reader1 = fastq.initStreamReader(allocator, sources[0], reader_options);
    defer reader1.deinit();
    var reader2 = fastq.initStreamReader(allocator, sources[1], reader_options);
    defer reader2.deinit();
    var validator1 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var validator2 = fastq.AdaptiveRecordValidator.init(.{ .alphabet = options.alphabet });
    var pairs: ProcessPairs = .{
        .allocator = allocator,
        .reader1 = &reader1,
        .reader2 = &reader2,
        .validator1 = &validator1,
        .validator2 = &validator2,
        .writers = writers,
        .collectors = collectors,
        .counts = counts,
        .max_pairs = options.max_pairs,
        .policy = options.pair_name_policy,
        .trimmer = .init(options.trim),
        .collect_qc = options.report != null,
    };

    while (true) {
        var span1: ?fastq.RecordSpan = null;
        if (try processSpanPairs(&pairs, &span1)) |failure| return failure;
        if (pairs.atCap()) return null;
        var next1 = if (span1) |mate1| fastq.validatedSpan(mate1) else fastq.nextValidatedRecord(&reader1, &validator1) catch |err| {
            return pairCommandFailure(0, mapReaderFailure(&reader1, err));
        };
        var next2 = fastq.nextValidatedRecord(&reader2, &validator2) catch |err| {
            return pairCommandFailure(1, mapReaderFailure(&reader2, err));
        };
        const maybe1: ?*const fastq.ValidatedRecord = if (next1) |*value| value else null;
        const maybe2: ?*const fastq.ValidatedRecord = if (next2) |*value| value else null;
        switch (verifyPair(
            &reader1,
            &reader2,
            if (maybe1) |mate1| mate1.validatedHeader() else null,
            if (maybe2) |mate2| mate2.validatedHeader() else null,
            options.pair_name_policy,
        )) {
            .end => {
                counts.input_complete = true;
                return null;
            },
            .failure => |failure| return failure,
            .pair => {},
        }
        const mate1 = maybe1.?;
        const mate2 = maybe2.?;
        if (recordHasUnwritableEnding(mate1.record, mate1.canonical_span)) {
            return pairCommandFailure(0, unwritableRecordFailure());
        }
        if (recordHasUnwritableEnding(mate2.record, mate2.canonical_span)) {
            return pairCommandFailure(1, unwritableRecordFailure());
        }
        if (try pairs.emit(
            .{ mate1.record, mate2.record },
            .{ mate1.canonical_span, mate2.canonical_span },
        )) |failure| return failure;
    }
}

const ProcessReport = struct {
    inputs: []const []const u8,
    snapshots: [2]?FileSnapshot,
    outputs: [2][]const u8,
    options: ProcessOptions,
    status_ok: bool,
    counts: ProcessCounts,
    stats: *const [2][2]zfastq.Stats,
    qc: *const [2]qc.MatePopulations,
    changes_reads: bool,
    trimming: *const [2]preprocess.MateCounts,
    filtering: *const preprocess.FilterCounts,
    inserts: *const preprocess.InsertCounts,
    contaminants: ?*const qc.ContaminantList = null,
    pair_duplication: ?*const preprocess.PairDuplication = null,
    input_failure: ?PairCommandFailure,
    failed_output: ?[]const u8,
};

fn writeProcessReport(json: *std.json.Stringify, report: ProcessReport) !void {
    try writeJsonHeader(json, "z-fastq/process-v1");
    try json.objectField("options");
    try json.beginObject();
    try json.objectField("gzip");
    try json.write(report.options.gzip);
    try json.objectField("max_pairs");
    try json.write(report.options.max_pairs);
    try json.objectField("pair_names");
    try json.write(@tagName(report.options.pair_name_policy));
    try json.objectField("alphabet");
    try json.write(@tagName(report.options.alphabet));
    try json.objectField("max_line_bytes");
    try json.write(report.options.max_line_bytes);
    const trim = report.options.trim;
    try json.objectField("adapter1");
    try json.write(if (trim.adapters) |adapters| adapters[0] else null);
    try json.objectField("adapter2");
    try json.write(if (trim.adapters) |adapters| adapters[1] else null);
    try json.objectField("trim_poly_g");
    try json.write(trim.poly_g);
    try json.objectField("trim_tail_window");
    try json.write(if (trim.tail) |tail| tail.window else null);
    try json.objectField("trim_tail_quality");
    try json.write(if (trim.tail) |tail| tail.quality else null);
    try json.objectField("min_length");
    try json.write(trim.min_length);
    try json.objectField("max_n");
    try json.write(trim.max_n);
    try json.objectField("low_quality");
    try json.write(if (trim.low_quality) |low_quality| low_quality.quality else null);
    try json.objectField("max_low_quality_percent");
    try json.write(if (trim.low_quality) |low_quality| low_quality.percent else null);
    try json.objectField("contaminants");
    try json.write(report.options.contaminants);
    try json.objectField("pair_duplication");
    try json.write(report.options.pair_duplication);
    try json.endObject();

    try json.objectField("inputs");
    try json.beginArray();
    for (report.inputs, report.snapshots) |path, snapshot| {
        try json.beginObject();
        try json.objectField("path");
        try writeEscapedJsonString(json, path);
        try json.objectField("size");
        try json.write(if (snapshot) |known| known.size else null);
        try json.objectField("mtime_ns");
        try json.write(if (snapshot) |known| known.mtime_nanoseconds else null);
        try json.endObject();
    }
    try json.endArray();
    try json.objectField("outputs");
    try json.beginArray();
    for (report.outputs) |path| {
        try json.beginObject();
        try json.objectField("path");
        try writeEscapedJsonString(json, path);
        try json.endObject();
    }
    try json.endArray();

    try json.objectField("status");
    try json.write(if (report.status_ok) "ok" else "error");
    try json.objectField("input_complete");
    try json.write(report.counts.input_complete);
    try json.objectField("pairs");
    try json.beginObject();
    try json.objectField("read");
    try json.write(report.counts.read);
    try json.objectField("written");
    try json.write(report.counts.written);
    try json.objectField("dropped");
    try json.write(report.counts.dropped);
    try json.endObject();
    try json.objectField("mates");
    try json.beginObject();
    for ([_][]const u8{ "raw", "retained" }, report.stats) |population, *population_stats| {
        try json.objectField(population);
        try json.beginArray();
        for (population_stats) |*mate_stats| {
            try json.beginObject();
            try writeStatsResultFields(json, mate_stats);
            try json.endObject();
        }
        try json.endArray();
    }
    try json.endObject();
    try json.objectField("qc");
    try json.beginObject();
    for ([_][]const u8{ "raw", "retained" }, [_]bool{ false, report.changes_reads }) |population, retained| {
        try json.objectField(population);
        try json.beginArray();
        for (report.qc) |*populations| {
            try populations.writeJson(json, retained, report.contaminants);
        }
        try json.endArray();
    }
    try json.endObject();
    try json.objectField("trimming");
    try json.beginArray();
    for (report.trimming) |*mate_counts| try writeTrimCounts(json, trim, mate_counts);
    try json.endArray();
    try json.objectField("filtering");
    if (report.changes_reads) try writeFilterCounts(json, report.filtering) else try json.write(null);
    try json.objectField("inserts");
    if (trim.adapters != null) try writeInsertCounts(json, report.inserts) else try json.write(null);
    try json.objectField("pair_duplication");
    if (report.pair_duplication) |pairs| try pairs.writeJson(json) else try json.write(null);

    const failed_input: ?[]const u8 = if (report.input_failure) |failure| switch (failure) {
        .command => |details| report.inputs[details.input_index],
        .pair => null,
    } else null;
    try json.objectField("failed_input");
    if (failed_input) |path| try writeEscapedJsonString(json, path) else try json.write(null);
    try json.objectField("failed_output");
    if (report.failed_output) |path| try writeEscapedJsonString(json, path) else try json.write(null);
    if (report.input_failure) |failure| switch (failure) {
        .command => |details| try writeJsonFailure(json, details.details),
        .pair => |details| try writePairJsonFailure(json, details),
    } else if (report.failed_output != null) {
        try writeJsonFailure(json, IO_FAILURE);
    } else {
        try json.objectField("error");
        try json.write(null);
    }
    try json.endObject();
    try json.writer.writeByte('\n');
}

fn writeTrimCounts(
    json: *std.json.Stringify,
    trim: preprocess.Options,
    counts: *const preprocess.MateCounts,
) !void {
    try json.beginObject();
    for (std.enums.values(preprocess.Step), counts.steps) |step, step_counts| {
        const enabled = switch (step) {
            .quality_tail => trim.tail != null,
            .poly_g => trim.poly_g != null,
            .overlap, .adapter => trim.adapters != null,
        };
        try json.objectField(@tagName(step));
        if (!enabled) {
            try json.write(null);
            continue;
        }
        try json.beginObject();
        try json.objectField("reads");
        try json.write(step_counts.reads);
        try json.objectField("bases");
        try json.write(step_counts.bases);
        try json.endObject();
    }
    try json.endObject();
}

fn writeInsertCounts(json: *std.json.Stringify, counts: *const preprocess.InsertCounts) !void {
    try json.beginObject();
    try json.objectField("pairs");
    try json.write(counts.pairs);
    try json.objectField("overlapping");
    try json.write(counts.overlapping);
    try json.objectField("unknown");
    try json.write(counts.unknown());
    try json.objectField("peak");
    try json.write(counts.peak());
    try json.objectField("median");
    try json.write(counts.median());
    try json.objectField("histogram");
    try json.beginArray();
    for (counts.histogram[0..preprocess.INSERT_BINS]) |count| try json.write(count);
    try json.endArray();
    try json.endObject();
}

fn writeFilterCounts(json: *std.json.Stringify, counts: *const preprocess.FilterCounts) !void {
    try json.beginObject();
    try json.objectField("dropped_pairs");
    try writeReasonCounts(json, &counts.dropped_pairs);
    try json.objectField("failed_reads");
    try json.beginArray();
    for (&counts.failed_reads) |*mate| try writeReasonCounts(json, mate);
    try json.endArray();
    try json.endObject();
}

fn writeReasonCounts(json: *std.json.Stringify, counts: *const [preprocess.REASON_COUNT]u64) !void {
    try json.beginObject();
    for (std.enums.values(preprocess.Reason), counts) |reason, count| {
        try json.objectField(@tagName(reason));
        try json.write(count);
    }
    try json.endObject();
}

fn writeJsonHeader(json: *std.json.Stringify, schema: []const u8) !void {
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
}

fn beginJsonDocument(json: *std.json.Stringify, schema: []const u8) !void {
    try writeJsonHeader(json, schema);
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
            try json.objectField("status");
            try json.write("ok");
            try writeStatsResultFields(json, &stats);
        },
    }
    try json.endObject();
}

fn writeStatsResultFields(json: *std.json.Stringify, stats: *const zfastq.Stats) !void {
    const result = stats.result();
    inline for (@typeInfo(zfastq.StatsResult).@"struct".fields) |field| {
        try json.objectField(field.name);
        try json.write(@field(result, field.name));
    }
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

const TEST_TEMP_DIR = ".zig-cache/tmp/{s}";
const TEST_TEMP_PATH = TEST_TEMP_DIR ++ "/{s}";

fn testTempPath(allocator: std.mem.Allocator, sub_path: []const u8, name: []const u8) ![]u8 {
    return std.fmt.allocPrint(allocator, TEST_TEMP_PATH, .{ sub_path, name });
}

fn commandFailure(failure: PairCommandFailure) error{ExpectedCommandFailure}!@FieldType(PairCommandFailure, "command") {
    return switch (failure) {
        .command => |command| command,
        .pair => error.ExpectedCommandFailure,
    };
}

fn expectInputChanged(failure: CommandFailure) !void {
    try std.testing.expectEqualStrings("input_changed", failure.code);
    try std.testing.expectEqualStrings("input changed during exact sampling", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
}

const TEST_SAMPLE_OPTIONS: SampleOptions = .{
    .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    .alphabet = .iupac,
    .fraction = null,
    .count = 1,
    .seed = 11,
};

const TEST_DEINTERLEAVE_OPTIONS: DeinterleaveOptions = .{
    .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
    .alphabet = .iupac,
    .pair_name_policy = .illumina,
};

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
    const valid = "@r\r\nAc\r\n+\r\n!5\r\n@s\nNG\n+\nI?\n";
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

const FakeIo = struct {
    dir: std.Io.Dir,
    stdout: ?*std.Io.Writer = null,
    stderr: ?*std.Io.Writer = null,
    replace_first_open: ?FirstOpenReplacement = null,
    fail_write: usize = 0,
    fail_read: usize = 0,
    opened: usize = 0,
    closed: usize = 0,
    writes: usize = 0,
    reads: usize = 0,
    vtable: @TypeOf(std.Io.failing.vtable.*) = undefined,

    const FirstOpenReplacement = union(enum) {
        rewrite: []const u8,
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

    const storage = try std.testing.allocator.alloc(u8, 2 * window);
    defer std.testing.allocator.free(storage);
    const block = "input: small.fastq\n" ++
        "reads: 1\nbases: 1\nmin_length: 1\nmax_length: 1\nmean_length: 1.000000\n" ++
        "a: 1\nc: 0\ng: 0\nt: 0\nn: 0\nother_bases: 0\ngc_fraction: 0.000000\n" ++
        "quality_sum: 0\nmean_quality: 0.000000\nq20_bases: 0\nq20_fraction: 0.000000\n" ++
        "q30_bases: 0\nq30_fraction: 0.000000\n";

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
        TEST_SAMPLE_OPTIONS,
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
        TEST_SAMPLE_OPTIONS,
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
        TEST_SAMPLE_OPTIONS,
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
        TEST_SAMPLE_OPTIONS,
    )) {
        .success => return error.ExpectedFailure,
        .failure => |failure| failure,
    };
    try std.testing.expectEqualStrings("io_error", failure.code);
    try std.testing.expectEqualStrings("I/O error", failure.message);
    try std.testing.expectEqual(@as(u8, 3), failure.exit_code);
    try std.testing.expectEqual(@as(u64, 0), selector.record_count);
}

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

const InterleavedStorageFixture = struct {
    bytes: std.ArrayList(u8),
    prefix_len: usize,
    header_starts: [2]usize,
    sequence_starts: [2]usize,
    plus_starts: [2]usize,

    fn deinit(self: *InterleavedStorageFixture, allocator: std.mem.Allocator) void {
        self.bytes.deinit(allocator);
    }
};

fn buildInterleavedStorageFixture(
    allocator: std.mem.Allocator,
    storage: InterleavedFirstRecordStorage,
) !InterleavedStorageFixture {
    const small_prefix = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n";
    const refill_margin = 300;
    const header_len = 192;
    const staging_limit = recordStagingLimit(zfastq.limits.DEFAULT_MAX_LINE_BYTES);
    var fixture: InterleavedStorageFixture = .{
        .bytes = .empty,
        .prefix_len = undefined,
        .header_starts = undefined,
        .sequence_starts = undefined,
        .plus_starts = undefined,
    };
    errdefer fixture.deinit(allocator);
    if (storage == .reader) {
        try fixture.bytes.appendSlice(allocator, "@a/1 ");
        try fixture.bytes.appendNTimes(
            allocator,
            'x',
            zfastq.limits.DEFAULT_READER_BUFFER_BYTES - refill_margin - small_prefix.len - 1,
        );
        try fixture.bytes.appendSlice(allocator, small_prefix[4..]);
    } else try fixture.bytes.appendSlice(allocator, small_prefix);
    fixture.prefix_len = fixture.bytes.items.len;
    for (0..2) |mate| {
        const field_len: usize = if (storage == .retained)
            zfastq.limits.DEFAULT_READER_BUFFER_BYTES / 2
        else
            8;
        try fixture.bytes.append(allocator, '@');
        fixture.header_starts[mate] = fixture.bytes.items.len;
        try fixture.bytes.appendNTimes(allocator, 'b', header_len);
        try fixture.bytes.appendSlice(allocator, if (mate == 0) "/1\n" else "/2\n");
        fixture.sequence_starts[mate] = fixture.bytes.items.len;
        try fixture.bytes.appendNTimes(allocator, if (mate == 0) 'A' else 'T', field_len);
        try fixture.bytes.append(allocator, '\n');
        fixture.plus_starts[mate] = fixture.bytes.items.len;
        try fixture.bytes.appendSlice(allocator, "+\n");
        try fixture.bytes.appendNTimes(allocator, '!', field_len);
        try fixture.bytes.append(allocator, '\n');
    }

    var source = io_layer.SliceSource.init(fixture.bytes.items);
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
    return fixture;
}

test "[integration] - [interleaved exact sample]: revalidates selected pairs after mate storage changes" {
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try testTempPath(allocator, &tmp.sub_path, "pairs");
    defer allocator.free(path);
    const options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 2,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    const staging_limit = recordStagingLimit(options.max_line_bytes);
    const name_change_offset = 180;
    const Change = enum { valid, alphabet2, name, structure };
    for ([_]InterleavedFirstRecordStorage{ .reader, .retained }) |storage| {
        var fixture = try buildInterleavedStorageFixture(allocator, storage);
        defer fixture.deinit(allocator);
        const input = fixture.bytes.items;
        const prefix = input[0..fixture.prefix_len];
        try writeExactTestInput(io, path, input, false, null);
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
        const changed = try allocator.dupe(u8, input);
        defer allocator.free(changed);
        const output = try allocator.alloc(u8, input.len);
        defer allocator.free(output);
        for ([_]bool{ false, true }) |select_all| {
            for (std.enums.values(Change)) |change| {
                errdefer std.debug.print("storage={s} select_all={} change={s}\n", .{ @tagName(storage), select_all, @tagName(change) });
                @memcpy(changed, input);
                switch (change) {
                    .valid => {},
                    .alphabet2 => changed[fixture.sequence_starts[1]] = '.',
                    .name => changed[fixture.header_starts[1] + name_change_offset] = 'c',
                    .structure => {
                        changed[fixture.sequence_starts[0]] = '.';
                        changed[fixture.plus_starts[1]] = 'x';
                    },
                }
                try writeExactTestInput(io, path, changed, false, snapshot);
                const expected_prefix = if (select_all) prefix else "";
                var sink = io_layer.SliceSink.init(output);
                var writer = zfastq.Writer.init(sink.byteSink());
                const result = try sampleExactInterleavedSecondPass(
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
                if (change == .valid) {
                    try std.testing.expect(result == null);
                    try std.testing.expectEqualStrings(if (select_all) input else input[prefix.len..], sink.written());
                } else {
                    const failure = result orelse return error.ExpectedFailure;
                    try std.testing.expect(failure == .command);
                    const details = failure.command.details;
                    try std.testing.expectEqualStrings(if (change == .structure) "S001" else "input_changed", details.code);
                    try std.testing.expectEqual(@as(u8, if (change == .structure) 1 else 3), details.exit_code);
                    if (change == .structure) {
                        try std.testing.expectEqual(@as(?u64, 3), details.record_index);
                        try std.testing.expectEqual(@as(?u3, 3), details.line_in_record);
                        try std.testing.expectEqual(@as(?u64, fixture.plus_starts[1]), details.byte_offset);
                    }
                    try std.testing.expectEqualStrings(expected_prefix, sink.written());
                }
            }
        }
    }
}

test "[failure] - [interleaved exact sample]: an output write failure keeps the written pairs and reports WriteFailed" {
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try testTempPath(allocator, &tmp.sub_path, "pairs");
    defer allocator.free(path);
    const options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 2,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    const staging_limit = recordStagingLimit(options.max_line_bytes);
    for ([_]InterleavedFirstRecordStorage{ .reader, .retained }) |storage| {
        var fixture = try buildInterleavedStorageFixture(allocator, storage);
        defer fixture.deinit(allocator);
        const input = fixture.bytes.items;
        const prefix = input[0..fixture.prefix_len];
        try writeExactTestInput(io, path, input, false, null);
        const snapshot = try snapshotTestFile(io, path);
        const output = try allocator.alloc(u8, input.len);
        defer allocator.free(output);
        for ([_]bool{ false, true }) |select_all| {
            errdefer std.debug.print("storage={s} select_all={}\n", .{ @tagName(storage), select_all });
            const expected_prefix = if (select_all) prefix else "";
            var sink = io_layer.SliceSink.init(output[0..expected_prefix.len]);
            var writer = zfastq.Writer.init(sink.byteSink());
            try std.testing.expectError(error.WriteFailed, sampleExactInterleavedSecondPass(
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
            ));
            try std.testing.expectEqualStrings(expected_prefix, sink.written());
        }
    }
}

test "[integration] - [interleaved exact sample]: preserves the first mate parse failure" {
    const small_line_limit = 64;
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

test "[failure] - [paired exact sample]: a record count change after output keeps the written prefix" {
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@a/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@a/2\nT\n+\n#\n" });
    var path1_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path1 = try std.fmt.bufPrint(&path1_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r1.fastq" });
    var path2_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path2 = try std.fmt.bufPrint(&path2_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r2.fastq" });
    const inputs = [_][]const u8{ path1, path2 };
    const snapshot1 = try snapshotTestFile(io, path1);
    const snapshot2 = try snapshotTestFile(io, path2);
    var options = TEST_SAMPLE_OPTIONS;
    options.pair_mode = .paired;

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
        .{ snapshot1, snapshot2 },
        2,
        options,
    )).?;

    const command = try commandFailure(failure);
    try std.testing.expectEqualStrings("input_changed", command.details.code);
    try std.testing.expectEqualStrings("@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n", sink.written());
}

test "[failure] - [interleaved exact sample]: an extra mate after the sampled pairs reports input changed" {
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const first_pair = "@a/1\nA\n+\n!\n@a/2\nT\n+\n#\n";
    try tmp.dir.writeFile(io, .{ .sub_path = "pairs.fastq", .data = first_pair ++ "@odd/1\nC\n+\n$\n" });
    var path_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path = try std.fmt.bufPrint(&path_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "pairs.fastq" });
    const snapshot = try snapshotTestFile(io, path);
    var options = TEST_SAMPLE_OPTIONS;
    options.pair_mode = .interleaved;

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const failure = (try sampleExactInterleavedSecondPass(
        true,
        io,
        std.testing.allocator,
        &record_inputs[0],
        path,
        &writer,
        .empty,
        snapshot,
        1,
        recordStagingLimit(options.max_line_bytes),
        options,
    )).?;

    const command = try commandFailure(failure);
    try expectInputChanged(command.details);
    try std.testing.expectEqualStrings(first_pair, sink.written());
}

test "[failure] - [paired exact sample]: a malformed record in the second pass reports its input index" {
    var record_inputs: [2]RecordInput = undefined;
    const io = std.testing.io;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    try tmp.dir.writeFile(io, .{ .sub_path = "r1.fastq", .data = "@bad/1\nA\n+\n!\n" });
    try tmp.dir.writeFile(io, .{ .sub_path = "r2.fastq", .data = "@bad/2\nT\nx\n#\n" });
    var path1_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path1 = try std.fmt.bufPrint(&path1_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r1.fastq" });
    var path2_buffer: [std.Io.Dir.max_path_bytes]u8 = undefined;
    const path2 = try std.fmt.bufPrint(&path2_buffer, TEST_TEMP_PATH, .{ tmp.sub_path, "r2.fastq" });
    const inputs = [_][]const u8{ path1, path2 };
    const snapshot1 = try snapshotTestFile(io, path1);
    const snapshot2 = try snapshotTestFile(io, path2);
    var options = TEST_SAMPLE_OPTIONS;
    options.pair_mode = .paired;

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
        .{ snapshot1, snapshot2 },
        1,
        options,
    )).?;

    const command = try commandFailure(failure);
    try std.testing.expectEqual(@as(u1, 1), command.input_index);
    try std.testing.expectEqualStrings("S001", command.details.code);
    try std.testing.expectEqual(@as(u8, 1), command.details.exit_code);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
}

test "[unit] - [paired identifiers]: storage reuses its largest allocation" {
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
    var failing = std.testing.FailingAllocator.init(std.testing.allocator, .{ .fail_index = 0 });
    try std.testing.expectError(error.OutOfMemory, exerciseExactAllocations(failing.allocator()));
    try std.testing.expect(failing.has_induced_failure);
}

test "[failure] - [interleaved exact sample]: a second-pass allocation for a spanning mate reports out of memory" {
    const second_pass_first_allocation = 1;
    const io = std.testing.io;
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    const path = try testTempPath(allocator, &tmp.sub_path, "pairs");
    defer allocator.free(path);
    var input: std.ArrayList(u8) = .empty;
    defer input.deinit(allocator);
    try input.appendSlice(allocator, "@pair/1\nA\n+");
    try input.appendNTimes(allocator, 'x', zfastq.limits.DEFAULT_READER_BUFFER_BYTES);
    try input.appendSlice(allocator, "\n!\n@pair/2\nT\n+\n#\n");
    try writeExactTestInput(io, path, input.items, false, null);
    var record_inputs: [2]RecordInput = undefined;
    const options = SampleOptions{
        .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
        .alphabet = .iupac,
        .fraction = null,
        .count = 1,
        .seed = 11,
        .pair_mode = .interleaved,
    };
    var selector = sampling.ExactSelector.init(1, options.seed);
    defer selector.deinit(allocator);
    var failing = std.testing.FailingAllocator.init(allocator, .{
        .fail_index = second_pass_first_allocation,
    });

    const first = sampleExactPairFirstPass(
        io,
        failing.allocator(),
        &record_inputs,
        &.{path},
        &selector,
        options,
    );
    try std.testing.expect(first == .success);
    try std.testing.expect(!failing.has_induced_failure);
    try std.testing.expectEqual(@as(usize, second_pass_first_allocation), failing.alloc_index);
    const snapshot = first.success.snapshots.interleaved;

    var output: [64]u8 = undefined;
    var sink = io_layer.SliceSink.init(&output);
    var writer = zfastq.Writer.init(sink.byteSink());
    const descriptors = try openDescriptorCount();
    const failure = (try sampleExactInterleavedSecondPass(
        true,
        io,
        failing.allocator(),
        &record_inputs[0],
        path,
        &writer,
        .empty,
        snapshot,
        1,
        recordStagingLimit(options.max_line_bytes),
        options,
    )).?;

    try std.testing.expect(failing.has_induced_failure);
    const command = try commandFailure(failure);
    try std.testing.expectEqual(@as(u1, 0), command.input_index);
    try std.testing.expectEqualStrings("out_of_memory", command.details.code);
    try std.testing.expectEqual(@as(u8, 3), command.details.exit_code);
    try std.testing.expectEqual(@as(usize, 0), sink.written().len);
    try std.testing.expectEqual(descriptors, try openDescriptorCount());
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
    const expected = "@a/1 left\nAC\n+left\n!#\n@a/2 right\nGT\n+right\n$%\n" ++
        "@c/1\nTTA\n+one\n()*\n@c/2\nAAT\n+two\n+,-\n";
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

test "[failure] - [process]: a mate 2 write failure keeps mate 1 and counts no written pair" {
    const mate1 = "@pair/1\nA\n+left\n!\n";
    const mate2 = "@pair/2\nT\n+right\n#\n";
    const trims = [_]preprocess.Options{ .{}, .{ .adapters = .{ "ACGT", "ACGT" } }, .{ .max_n = 0, .min_length = 1 } };
    for (trims, 0..) |trim, case| {
        errdefer std.debug.print("case {d}\n", .{case});
        var source1 = io_layer.SliceSource.init(mate1);
        var source2 = io_layer.SliceSource.init(mate2);
        var stream1 = try TestStream.init(source1.byteSource());
        defer stream1.deinit();
        var stream2 = try TestStream.init(source2.byteSource());
        defer stream2.deinit();
        var sink1 = DeinterleaveTestSink{};
        var sink2 = DeinterleaveTestSink{ .fail_write = true };
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());
        var collectors: ProcessCollectors = .{};
        defer collectors.deinit(std.testing.allocator);
        var counts: ProcessCounts = .{};

        try std.testing.expectError(error.Output2WriteFailed, processSources(
            std.testing.allocator,
            .{ &stream1.reader.interface, &stream2.reader.interface },
            .{ &writer1, &writer2 },
            &collectors,
            &counts,
            .{
                .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
                .alphabet = .iupac,
                .pair_name_policy = .illumina,
                .output1 = null,
                .output2 = null,
                .report = null,
                .gzip = false,
                .max_pairs = null,
                .trim = trim,
            },
        ));
        try std.testing.expectEqualStrings(mate1, sink1.buffer[0..sink1.length]);
        try std.testing.expectEqual(@as(u64, 1), counts.read);
        try std.testing.expectEqual(@as(u64, 0), counts.written);
        try std.testing.expect(!counts.input_complete);
    }
}

test "[property] - [process]: trimmed, dropped, and long pairs release every allocation" {
    const allocator = std.testing.allocator;
    const Pair = struct { name: []const u8, length: usize, quality: u8 = 'I', kept: bool };
    const pairs = [_]Pair{
        .{ .name = "kept", .length = 60, .kept = true },
        .{ .name = "short", .length = 10, .kept = false },
        .{ .name = "long", .length = 1100, .kept = true },
        .{ .name = "low", .length = 60, .quality = '#', .kept = false },
        .{ .name = "after", .length = 60, .kept = true },
    };
    var inputs: [2]std.ArrayList(u8) = .{ .empty, .empty };
    defer for (&inputs) |*input| input.deinit(allocator);
    var state: u32 = 7;
    for (pairs) |pair| {
        for (&inputs, 1..) |*input, mate| {
            try input.print(allocator, "@{s}/{d}\n", .{ pair.name, mate });
            for (0..pair.length) |_| {
                state = state *% 1103515245 +% 12345;
                try input.append(allocator, "ACGT"[(state >> 16) & 3]);
            }
            try input.appendSlice(allocator, "\n+\n");
            try input.appendNTimes(allocator, pair.quality, pair.length);
            try input.append(allocator, '\n');
        }
    }
    var source1 = io_layer.SliceSource.init(inputs[0].items);
    var source2 = io_layer.SliceSource.init(inputs[1].items);
    var stream1 = try TestStream.init(source1.byteSource());
    defer stream1.deinit();
    var stream2 = try TestStream.init(source2.byteSource());
    defer stream2.deinit();
    var output1: [4096]u8 = undefined;
    var output2: [4096]u8 = undefined;
    var sink1 = io_layer.SliceSink.init(&output1);
    var sink2 = io_layer.SliceSink.init(&output2);
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());
    var collectors: ProcessCollectors = .{};
    defer collectors.deinit(allocator);
    var counts: ProcessCounts = .{};

    try std.testing.expectEqual(@as(?PairCommandFailure, null), try processSources(
        allocator,
        .{ &stream1.reader.interface, &stream2.reader.interface },
        .{ &writer1, &writer2 },
        &collectors,
        &counts,
        .{
            .max_line_bytes = zfastq.limits.DEFAULT_MAX_LINE_BYTES,
            .alphabet = .iupac,
            .pair_name_policy = .illumina,
            .output1 = null,
            .output2 = null,
            .report = "-",
            .gzip = false,
            .max_pairs = null,
            .trim = .{
                .adapters = .{ "CTGTCTCTTATACACATCT", "CTGTCTCTTATACACATCT" },
                .poly_g = 10,
                .tail = .{ .window = 4, .quality = 20 },
                .min_length = 30,
                .max_n = 5,
                .low_quality = .{ .quality = 15, .percent = 40 },
            },
        },
    ));
    try writer1.flush();
    try writer2.flush();
    try std.testing.expect(counts.input_complete);
    try std.testing.expectEqual(ProcessCounts{ .read = 5, .written = 3, .dropped = 2, .input_complete = true }, counts);
    try std.testing.expectEqual(@as(u64, 1), collectors.filtering.dropped_pairs[@intFromEnum(preprocess.Reason.too_short)]);
    try std.testing.expectEqual(@as(u64, 1), collectors.filtering.dropped_pairs[@intFromEnum(preprocess.Reason.low_quality)]);
    var stats: [2][2]zfastq.Stats = @splat(.{ .{}, .{} });
    var failure: ?PairCommandFailure = null;
    finishCollectors(allocator, &collectors.qc, true, &stats, &failure);
    try std.testing.expectEqual(@as(?PairCommandFailure, null), failure);
    for (stats[0], stats[1]) |raw, retained| {
        try std.testing.expectEqual(@as(u64, 5), raw.reads);
        try std.testing.expectEqual(@as(u64, 3), retained.reads);
        try std.testing.expectEqual(@as(u64, 1220), retained.bases);
    }
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
    for (cases, 0..) |case, case_index| {
        errdefer std.debug.print("case={d}\n", .{case_index});
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
    const options = TEST_DEINTERLEAVE_OPTIONS;

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
    const command_failure = (try commandFailure(failure)).details;
    try std.testing.expectEqualStrings("out_of_memory", command_failure.code);
    try std.testing.expectEqual(@as(usize, 0), sink1.length);
    try std.testing.expectEqual(@as(usize, 0), sink2.length);
}

const fallback_first_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES - 7;
const fallback_second_len = zfastq.limits.DEFAULT_READER_BUFFER_BYTES;
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

test "[failure] - [deinterleave]: a failed output write emits no partial pair" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);
    const output = try std.testing.allocator.alloc(u8, input.items.len);
    defer std.testing.allocator.free(output);
    const options = TEST_DEINTERLEAVE_OPTIONS;

    for (0..2) |failed_output| {
        var source = io_layer.SliceSource.init(input.items);
        var stream = try TestStream.init(source.byteSource());
        defer stream.deinit();
        var sink1 = io_layer.SliceSink.init(if (failed_output == 0) output[0..0] else output);
        var sink2 = io_layer.SliceSink.init(output[0..0]);
        var writer1 = zfastq.Writer.init(sink1.byteSink());
        var writer2 = zfastq.Writer.init(sink2.byteSink());

        try std.testing.expectError(
            if (failed_output == 0) error.Output1WriteFailed else error.Output2WriteFailed,
            deinterleaveSource(
                std.testing.allocator,
                &stream.reader.interface,
                &writer1,
                &writer2,
                recordStagingLimit(options.max_line_bytes),
                options,
            ),
        );
        if (failed_output == 0) {
            try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
        } else {
            try std.testing.expectEqualStrings(input.items[0..fallback_first_record_len], sink1.written());
        }
    }
}

test "[failure] - [deinterleave]: a one-byte staging limit reports arithmetic limit before any output" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);
    const output = try std.testing.allocator.alloc(u8, input.items.len);
    defer std.testing.allocator.free(output);
    var source = io_layer.SliceSource.init(input.items);
    var stream = try TestStream.init(source.byteSource());
    defer stream.deinit();
    var sink1 = io_layer.SliceSink.init(output);
    var sink2 = io_layer.SliceSink.init(output[0..0]);
    var writer1 = zfastq.Writer.init(sink1.byteSink());
    var writer2 = zfastq.Writer.init(sink2.byteSink());

    const failure = (try deinterleaveSource(
        std.testing.allocator,
        &stream.reader.interface,
        &writer1,
        &writer2,
        1,
        TEST_DEINTERLEAVE_OPTIONS,
    )).?;

    const command = try commandFailure(failure);
    try std.testing.expectEqualStrings("arithmetic_limit", command.details.code);
    try std.testing.expectEqual(@as(usize, 0), sink1.written().len);
    try std.testing.expectEqual(@as(usize, 0), sink2.written().len);
}

test "[failure] - [deinterleave]: malformed pairs fail with their diagnostic before any output" {
    var input = try fallbackPairInput(std.testing.allocator);
    defer input.deinit(std.testing.allocator);
    const options = TEST_DEINTERLEAVE_OPTIONS;
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
    const options = TEST_DEINTERLEAVE_OPTIONS;
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
