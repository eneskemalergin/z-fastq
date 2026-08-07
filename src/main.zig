//! CLI entry and subcommand dispatcher for z-fastq.

const std = @import("std");
const zfastq = @import("z-fastq");
const cli = @import("cli/mod.zig");

const USAGE =
    \\usage: z-fastq <command> [options] [args...]
    \\
    \\Commands:
    \\  count    Count records in plain FASTQ files
    \\
    \\General options:
    \\  -h, --help       Show this help message
    \\  -V, --version    Print version
    \\  --max-line-bytes N   Override default line length limit
    \\
    \\Count usage:
    \\  z-fastq count <file.fastq> [file2.fastq ...]
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

    var max_line_bytes = zfastq.limits.default_max_line_bytes;
    var positional = std.ArrayList([]const u8).empty;
    defer positional.deinit(gpa);

    if (std.mem.eql(u8, cmd, "count")) {
        while (args.next()) |arg| {
            if (std.mem.eql(u8, arg, "--max-line-bytes")) {
                const value = args.next() orelse {
                    std.Io.File.writeStreamingAll(.stderr(), io, "error: --max-line-bytes requires a value\n") catch {};
                    std.process.exit(2);
                };
                max_line_bytes = std.fmt.parseInt(usize, value, 10) catch {
                    std.Io.File.writeStreamingAll(.stderr(), io, "error: invalid --max-line-bytes value\n") catch {};
                    std.process.exit(2);
                };
                continue;
            }
            positional.append(gpa, arg) catch {
                std.Io.File.writeStreamingAll(.stderr(), io, "error: out of memory\n") catch {};
                std.process.exit(3);
            };
        }

        const code = cli.count.run(io, gpa, positional.items, .{
            .max_line_bytes = max_line_bytes,
        });
        std.process.exit(code);
    }

    std.Io.File.writeStreamingAll(.stderr(), io, "error: unknown command: ") catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, cmd) catch {};
    std.Io.File.writeStreamingAll(.stderr(), io, "\n") catch {};
    printUsageAndExit(io);
}

fn printUsageAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, USAGE) catch {};
    std.process.exit(2);
}

fn printHelpAndExit(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stdout(), io, USAGE) catch {};
    std.process.exit(0);
}

fn printVersionAndExit(io: std.Io) noreturn {
    var buf: [64]u8 = undefined;
    const line = std.fmt.bufPrint(&buf, "z-fastq {s}\n", .{zfastq.VERSION}) catch "z-fastq\n";
    std.Io.File.writeStreamingAll(.stdout(), io, line) catch {};
    std.process.exit(0);
}
