//! Product build: z-fastq exe (ReleaseFast + strip by default) and unit tests.

const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    if (optimize == .ReleaseSmall) {
        std.debug.print(
            "error: ReleaseSmall is unsupported; use -Doptimize=ReleaseFast (strips by default)\n",
            .{},
        );
        std.process.exit(1);
    }
    const strip = b.option(bool, "strip", "Strip debug info") orelse (optimize != .Debug);

    const lib_module = b.createModule(.{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    const exe = b.addExecutable(.{
        .name = "z-fastq",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .strip = strip,
            .imports = &.{
                .{ .name = "z-fastq", .module = lib_module },
            },
        }),
    });

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run z-fastq");
    run_step.dependOn(&run_cmd.step);

    const import_lib = [_]std.Build.Module.Import{
        .{ .name = "z-fastq", .module = lib_module },
    };

    const root_test_module = b.createModule(.{
        .root_source_file = b.path("tests/root_test.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });

    const reader_test_module = b.createModule(.{
        .root_source_file = b.path("tests/reader_test.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });

    const writer_test_module = b.createModule(.{
        .root_source_file = b.path("tests/writer_test.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });

    const count_test_module = b.createModule(.{
        .root_source_file = b.path("tests/count_test.zig"),
        .target = target,
        .optimize = optimize,
    });
    count_test_module.link_libc = true;

    const run_root_test = b.addRunArtifact(b.addTest(.{ .root_module = root_test_module }));
    run_root_test.step.dependOn(b.getInstallStep());

    const run_reader_test = b.addRunArtifact(b.addTest(.{ .root_module = reader_test_module }));
    const run_writer_test = b.addRunArtifact(b.addTest(.{ .root_module = writer_test_module }));

    const run_count_test = b.addRunArtifact(b.addTest(.{ .root_module = count_test_module }));
    run_count_test.step.dependOn(b.getInstallStep());

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_root_test.step);
    test_step.dependOn(&run_reader_test.step);
    test_step.dependOn(&run_writer_test.step);
    test_step.dependOn(&run_count_test.step);
}
