//! Benchmark exec helper: give the measured program a fresh stdin, a stdout
//! file, and per-run paths, then replace this process with it.
//!
//!     zexec [-i FILE] [-o FILE] [-r FILE] [-u FILE]... -- PROGRAM [ARG...]
//!
//! `-i` opens FILE as stdin and `-o` creates or truncates FILE as stdout, for
//! each run. "@RUN@" in FILE paths and in ARGs becomes zexec's process id, so
//! every run writes new paths. PROGRAM must contain a "/": there is no PATH
//! search. Arguments without the token pass through untouched. Any setup or
//! exec failure prints one line to stderr and exits 127.
//!
//! Without `-r`, zexec replaces itself with PROGRAM. The runner that started
//! zexec then reports PROGRAM's resources, but a peak RSS below the runner's own
//! footprint at fork is hidden: the high-water mark it inherits survives exec.
//!
//! With `-r`, zexec runs PROGRAM as its child, waits, appends one line to FILE,
//! and ends the way PROGRAM ended (same exit status, same signal). Each `-u`
//! FILE is removed once PROGRAM has ended and its line is written, so a run's
//! output files never pile up across runs. The line is
//! `maxrss_kib user_us system_us minor_faults status wall_us`, tab-separated:
//! the child's own `wait4` usage, then the monotonic time from fork to wait. zexec is static with no libc, so the child starts
//! with a few KiB of inherited pages instead of the runner's, and its peak RSS
//! is its own.
const std = @import("std");
const linux = std.os.linux;

const token = "@RUN@";
const max_args = 4096;

var argv_out: [max_args:null]?[*:0]const u8 = undefined;
var pool: [1 << 18]u8 = undefined;
var pool_used: usize = 0;
var pid_buf: [20]u8 = undefined;
var pid: []const u8 = undefined;
var unlink_paths: [8][*:0]const u8 = undefined;
var unlink_count: usize = 0;

pub fn main(init: std.process.Init.Minimal) noreturn {
    const argv = init.args.vector;
    pid = std.fmt.bufPrint(&pid_buf, "{d}", .{linux.getpid()}) catch unreachable;

    var stdin_path: ?[*:0]const u8 = null;
    var stdout_path: ?[*:0]const u8 = null;
    var report_path: ?[*:0]const u8 = null;
    var i: usize = 1;
    while (true) : (i += 2) {
        if (i >= argv.len) die("usage:", "zexec [-i FILE] [-o FILE] -- PROGRAM [ARG...]", null);
        const arg = std.mem.span(argv[i]);
        if (std.mem.eql(u8, arg, "--")) break;
        if (i + 1 >= argv.len) die("missing value for", arg, null);
        if (std.mem.eql(u8, arg, "-i")) {
            stdin_path = expand(argv[i + 1]);
        } else if (std.mem.eql(u8, arg, "-o")) {
            stdout_path = expand(argv[i + 1]);
        } else if (std.mem.eql(u8, arg, "-r")) {
            report_path = expand(argv[i + 1]);
        } else if (std.mem.eql(u8, arg, "-u")) {
            if (unlink_count == unlink_paths.len) die("too many", "-u files", null);
            unlink_paths[unlink_count] = expand(argv[i + 1]);
            unlink_count += 1;
        } else die("unknown option", arg, null);
    }
    i += 1;
    if (i >= argv.len) die("missing", "PROGRAM", null);
    if (argv.len - i > max_args) die("too many arguments for", "PROGRAM", null);
    const program = std.mem.span(argv[i]);
    if (std.mem.indexOfScalar(u8, program, '/') == null) die("PROGRAM needs a path with '/':", program, null);

    var n: usize = 0;
    while (i < argv.len) : ({
        i += 1;
        n += 1;
    }) argv_out[n] = expand(argv[i]);
    argv_out[n] = null;

    const envp: [*:null]const ?[*:0]const u8 = @ptrCast(init.environ.block.slice.ptr);
    const report = report_path orelse run(program, stdin_path, stdout_path, envp);

    const fd = linux.open(report, .{ .ACCMODE = .WRONLY, .CREAT = true, .APPEND = true, .CLOEXEC = true }, 0o644);
    if (linux.errno(fd) != .SUCCESS) die("open", std.mem.span(report), linux.errno(fd));
    var start: linux.timespec = undefined;
    _ = linux.clock_gettime(.MONOTONIC, &start);
    const forked = linux.fork();
    if (linux.errno(forked) != .SUCCESS) die("fork", program, linux.errno(forked));
    if (forked == 0) {
        // Die with zexec, so a runner timeout cannot leave the program running.
        _ = linux.prctl(@intFromEnum(linux.PR.SET_PDEATHSIG), @intFromEnum(linux.SIG.KILL), 0, 0, 0);
        run(program, stdin_path, stdout_path, envp);
    }

    var status: u32 = 0;
    var usage: linux.rusage = undefined;
    while (true) {
        const waited = linux.wait4(@intCast(forked), &status, 0, &usage);
        if (linux.errno(waited) == .SUCCESS) break;
        if (linux.errno(waited) != .INTR) die("wait", program, linux.errno(waited));
    }
    var end: linux.timespec = undefined;
    _ = linux.clock_gettime(.MONOTONIC, &end);
    const wall_us = @divTrunc((@as(i64, end.sec) - @as(i64, start.sec)) * 1_000_000_000 + (@as(i64, end.nsec) - @as(i64, start.nsec)), 1000);
    var line: [160]u8 = undefined;
    const text = std.fmt.bufPrint(&line, "{d}\t{d}\t{d}\t{d}\t{d}\t{d}\n", .{
        usage.maxrss,
        micros(usage.utime),
        micros(usage.stime),
        usage.minflt,
        status,
        wall_us,
    }) catch unreachable;
    if (linux.write(@intCast(fd), text.ptr, text.len) != text.len) die("write", std.mem.span(report), null);
    // A file the program did not write is not an error here.
    for (unlink_paths[0..unlink_count]) |path| _ = linux.unlinkat(linux.AT.FDCWD, path, 0);

    if (linux.W.IFSIGNALED(status)) {
        const sig = linux.W.TERMSIG(status);
        _ = linux.kill(linux.getpid(), sig);
        linux.exit_group(128 + @as(i32, @intCast(@intFromEnum(sig))));
    }
    linux.exit_group(linux.W.EXITSTATUS(status));
}

/// Applies the redirects and replaces this process with PROGRAM.
fn run(program: []const u8, stdin_path: ?[*:0]const u8, stdout_path: ?[*:0]const u8, envp: [*:null]const ?[*:0]const u8) noreturn {
    if (stdin_path) |path| redirect(path, .{ .ACCMODE = .RDONLY, .CLOEXEC = true }, 0);
    if (stdout_path) |path| redirect(path, .{ .ACCMODE = .WRONLY, .CREAT = true, .TRUNC = true, .CLOEXEC = true }, 1);
    const rc = linux.execve(argv_out[0].?, &argv_out, envp);
    die("exec", program, linux.errno(rc));
}

fn micros(t: linux.timeval) i64 {
    return @as(i64, t.sec) * 1_000_000 + t.usec;
}

/// Returns `arg` itself, or a copy in `pool` with every token replaced by the pid.
fn expand(arg: [*:0]const u8) [*:0]const u8 {
    const text = std.mem.span(arg);
    if (std.mem.indexOf(u8, text, token) == null) return arg;
    const len = std.mem.replacementSize(u8, text, token, pid);
    if (pool.len - pool_used < len + 1) die("arguments too long:", text, null);
    const out = pool[pool_used..][0 .. len + 1];
    _ = std.mem.replace(u8, text, token, pid, out[0..len]);
    out[len] = 0;
    pool_used += len + 1;
    return out[0..len :0].ptr;
}

/// Opens `path` onto descriptor `fd`, leaving no extra descriptor open.
fn redirect(path: [*:0]const u8, flags: linux.O, fd: i32) void {
    const rc = linux.open(path, flags, 0o644);
    if (linux.errno(rc) != .SUCCESS) die("open", std.mem.span(path), linux.errno(rc));
    const opened: i32 = @intCast(rc);
    if (opened == fd) {
        // Already on the target descriptor: only clear close-on-exec.
        _ = linux.fcntl(fd, linux.F.SETFD, 0);
        return;
    }
    const dup = linux.dup2(opened, fd);
    if (linux.errno(dup) != .SUCCESS) die("dup2", std.mem.span(path), linux.errno(dup));
    _ = linux.close(opened);
}

fn die(step: []const u8, detail: []const u8, err: ?linux.E) noreturn {
    var buf: [1024]u8 = undefined;
    const shown = detail[0..@min(detail.len, 800)];
    const msg = if (err) |e|
        std.fmt.bufPrint(&buf, "zexec: {s} {s}: {s}\n", .{ step, shown, @tagName(e) }) catch "zexec: error\n"
    else
        std.fmt.bufPrint(&buf, "zexec: {s} {s}\n", .{ step, shown }) catch "zexec: error\n";
    _ = linux.write(2, msg.ptr, msg.len);
    linux.exit_group(127);
}
