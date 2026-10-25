//! Subset checker: is a sampler's output an ordered subset of its input?
//!
//!     subset [--pairs] INPUT [INPUT2] < OUTPUT
//!
//! OUTPUT is the sampler's four-line FASTQ records. With one INPUT, each
//! output record must be an input record, in input order. With `--pairs`,
//! INPUT is interleaved and kept records must come as whole input pairs. With
//! INPUT and INPUT2 (mate files), OUTPUT is interleaved: odd records come
//! from INPUT, even ones from INPUT2, and both mates of a pair share an index.
//!
//! A record matches by header, sequence, and quality, each without a trailing
//! CR. The walk is a greedy subsequence match, so it decides order correctly
//! even where records repeat. Prints `name value` facts:
//!
//!     kept N           records written (pairs with INPUT2)
//!     order yes|no R   no: output record R is not the next input match
//!     pairs yes|no     mates kept together (with --pairs or INPUT2)
//!     plus N           kept records whose plus line differs from the input's
//!     bare N           kept records whose plus line is a bare `+`
//!     cr N             kept lines whose trailing CR differs from the input's
//!     crlines N        output lines that end in CR
//!     fingerprint H    FNV-1a of the kept record (or pair) indexes
//!
//! The whole output is read even after a record out of order, so the sampler
//! always runs to its own end. Exits 2 on a usage or read error, 0 otherwise.
const std = @import("std");

const line_max = 8 << 20;

const Line = struct {
    bytes: []const u8,
    cr: bool,

    fn of(raw: []const u8) Line {
        const cr = raw.len > 0 and raw[raw.len - 1] == '\r';
        return .{ .bytes = if (cr) raw[0 .. raw.len - 1] else raw, .cr = cr };
    }
};

/// One output record, copied out of the reader's buffer.
const Record = struct {
    lines: [4]std.ArrayList(u8) = .{ .empty, .empty, .empty, .empty },
    cr: [4]bool = .{ false, false, false, false },

    fn read(self: *Record, gpa: std.mem.Allocator, reader: *std.Io.Reader) !bool {
        for (0..4) |i| {
            const raw = try reader.takeDelimiter('\n') orelse {
                if (i == 0) return false;
                return error.TruncatedRecord;
            };
            const line = Line.of(raw);
            self.lines[i].clearRetainingCapacity();
            try self.lines[i].appendSlice(gpa, line.bytes);
            self.cr[i] = line.cr;
        }
        return true;
    }
};

/// The input side: reads records until one matches, counting its index.
const Walk = struct {
    reader: *std.Io.Reader,
    index: u64 = 0,
    plus: u64 = 0,
    cr: u64 = 0,

    /// Advances to the next input record equal to `want`; false at end of input.
    fn find(self: *Walk, want: *const Record) !bool {
        while (true) {
            const h = Line.of(try self.reader.takeDelimiter('\n') orelse return false);
            self.index += 1;
            var same = std.mem.eql(u8, h.bytes, want.lines[0].items);
            var cr_diff: u64 = @intFromBool(h.cr != want.cr[0]);
            var plus_diff = false;
            for (1..4) |i| {
                const l = Line.of(try self.reader.takeDelimiter('\n') orelse return error.TruncatedInput);
                cr_diff += @intFromBool(l.cr != want.cr[i]);
                if (i == 2) {
                    plus_diff = !std.mem.eql(u8, l.bytes, want.lines[2].items);
                } else {
                    same = same and std.mem.eql(u8, l.bytes, want.lines[i].items);
                }
            }
            if (same) {
                self.plus += @intFromBool(plus_diff);
                self.cr += cr_diff;
                return true;
            }
        }
    }
};

fn fnv(hash: *u64, value: u64) void {
    var bytes: [8]u8 = undefined;
    std.mem.writeInt(u64, &bytes, value, .little);
    for (bytes) |b| {
        hash.* ^= b;
        hash.* *%= 0x100000001b3;
    }
}

fn usage(io: std.Io) noreturn {
    std.Io.File.writeStreamingAll(.stderr(), io, "usage: subset [--pairs] INPUT [INPUT2] < OUTPUT\n") catch {};
    std.process.exit(2);
}

fn fail(io: std.Io, what: []const u8, err: anyerror) noreturn {
    var buf: [512]u8 = undefined;
    const text = std.fmt.bufPrint(&buf, "subset: {s}: {s}\n", .{ what, @errorName(err) }) catch "subset: error\n";
    std.Io.File.writeStreamingAll(.stderr(), io, text) catch {};
    std.process.exit(2);
}

pub fn main(init: std.process.Init) void {
    const io = init.io;
    const gpa = init.gpa;
    const args = std.process.Args.toSlice(init.minimal.args, init.arena.allocator()) catch |err| fail(io, "arguments", err);

    var pairs = false;
    var paths: [2][]const u8 = undefined;
    var n_paths: usize = 0;
    for (args[1..]) |arg| {
        if (std.mem.eql(u8, arg, "--pairs")) {
            pairs = true;
        } else {
            if (n_paths == 2) usage(io);
            paths[n_paths] = arg;
            n_paths += 1;
        }
    }
    if (n_paths == 0 or (pairs and n_paths == 2)) usage(io);
    const mates = n_paths == 2;

    var walks: [2]Walk = undefined;
    var readers: [2]std.Io.File.Reader = undefined;
    for (0..n_paths) |i| {
        const file = std.Io.Dir.cwd().openFile(io, paths[i], .{}) catch |err| fail(io, paths[i], err);
        const buf = gpa.alloc(u8, line_max) catch |err| fail(io, "memory", err);
        readers[i] = file.reader(io, buf);
        walks[i] = .{ .reader = &readers[i].interface };
    }
    const out_buf = gpa.alloc(u8, line_max) catch |err| fail(io, "memory", err);
    var out_reader = std.Io.File.stdin().readerStreaming(io, out_buf);

    var record: Record = .{};
    var kept: u64 = 0;
    var ordered = true;
    var bad_at: u64 = 0;
    var paired = true;
    var hash: u64 = 0xcbf29ce484222325;
    var first_index: u64 = 0;
    var bare: u64 = 0;
    var crlines: u64 = 0;
    while (ordered) {
        const more = record.read(gpa, &out_reader.interface) catch |err| switch (err) {
            error.TruncatedRecord => {
                ordered = false;
                bad_at = kept + 1;
                break;
            },
            else => fail(io, "output", err),
        };
        if (!more) break;
        kept += 1;
        if (std.mem.eql(u8, record.lines[2].items, "+")) bare += 1;
        for (record.cr) |c| crlines += @intFromBool(c);
        const walk = &walks[if (mates and kept % 2 == 0) 1 else 0];
        const found = walk.find(&record) catch |err| fail(io, "input", err);
        if (!found) {
            ordered = false;
            bad_at = kept;
            break;
        }
        if (mates or pairs) {
            if (kept % 2 == 1) {
                first_index = walk.index;
                // An interleaved pair starts on an odd input record.
                if (pairs and walk.index % 2 == 0) paired = false;
            } else {
                // The second mate is the first mate's partner.
                const partner = if (mates) first_index else first_index + 1;
                if (walk.index != partner) paired = false;
                fnv(&hash, if (mates) first_index else (first_index + 1) / 2);
            }
        } else {
            fnv(&hash, walk.index);
        }
    }
    if ((mates or pairs) and kept % 2 == 1) paired = false;
    // Read the rest of the output, so the sampler is not stopped by a closed pipe.
    if (!ordered) _ = out_reader.interface.discardRemaining() catch |err| fail(io, "output", err);

    var plus: u64 = 0;
    var cr: u64 = 0;
    for (walks[0..n_paths]) |w| {
        plus += w.plus;
        cr += w.cr;
    }
    var buf: [512]u8 = undefined;
    var w = std.Io.File.stdout().writer(io, &buf);
    const out = &w.interface;
    out.print("kept {d}\n", .{if (mates or pairs) kept / 2 else kept}) catch |err| fail(io, "stdout", err);
    if (ordered) out.writeAll("order yes\n") catch |err| fail(io, "stdout", err) else out.print("order no {d}\n", .{bad_at}) catch |err| fail(io, "stdout", err);
    if (mates or pairs) out.print("pairs {s}\n", .{if (paired) "yes" else "no"}) catch |err| fail(io, "stdout", err);
    out.print("plus {d}\nbare {d}\ncr {d}\ncrlines {d}\nfingerprint {x:0>16}\n", .{ plus, bare, cr, crlines, hash }) catch |err| fail(io, "stdout", err);
    out.flush() catch |err| fail(io, "stdout", err);
}
