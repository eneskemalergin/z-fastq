//! Build graph for z-fastq.

const std = @import("std");

const KernelBackend = enum { dispatch, portable };

pub fn build(b: *std.Build) void {
    const requested_target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    switch (optimize) {
        .Debug, .ReleaseSafe, .ReleaseFast => {},
        else => {
            std.debug.print(
                "error: z-fastq builds only in Debug, ReleaseSafe, or ReleaseFast\n",
                .{},
            );
            std.process.exit(1);
        },
    }
    if (requested_target.result.cpu.arch != .x86_64 or
        requested_target.result.os.tag != .linux)
    {
        std.debug.print(
            "error: z-fastq currently supports Linux x86-64 builds only\n",
            .{},
        );
        std.process.exit(1);
    }
    const release_build = optimize == .ReleaseSafe or optimize == .ReleaseFast;
    if (b.option(bool, "static", "Confirm a static release build")) |requested_static| {
        if (requested_static != release_build) {
            std.debug.print(
                "error: Debug is native; ReleaseSafe and ReleaseFast are always static\n",
                .{},
            );
            std.process.exit(1);
        }
    }
    const target = if (release_build)
        b.resolveTargetQuery(.{
            .cpu_arch = .x86_64,
            .os_tag = .linux,
            .abi = .musl,
        })
    else
        requested_target;
    const strip = optimize == .ReleaseFast;
    const kernel_backend = b.option(
        KernelBackend,
        "kernel-backend",
        "Zipir CPU kernels: dispatch selects for the running CPU; portable forces the portable kernels",
    ) orelse .dispatch;
    const zipir = b.dependency("zipir", .{
        .target = target,
        .optimize = optimize,
        .@"kernel-backend" = kernel_backend,
    }).module("zipir");
    const package_version = @import("build.zig.zon").version;
    const build_options = b.addOptions();
    build_options.addOption([:0]const u8, "version", package_version);

    const lib_module = b.addModule("z-fastq", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    lib_module.addOptions("build_options", build_options);

    const exe = b.addExecutable(.{
        .name = "z-fastq",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .single_threaded = true,
            .strip = strip,
        }),
    });
    exe.root_module.addOptions("build_options", build_options);

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

    const fastq_test_module = b.createModule(.{
        .root_source_file = b.path("src/fastq.zig"),
        .target = target,
        .optimize = optimize,
    });
    for ([_]*std.Build.Module{ lib_module, exe.root_module, fastq_test_module }) |module| {
        module.addImport("zipir", zipir);
    }

    const reader_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_reader.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });
    const test_options = b.addOptions();
    test_options.addOption([]const u8, "package_version", package_version);

    const writer_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_writer.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });

    const count_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_count.zig"),
        .target = target,
        .optimize = optimize,
    });
    count_test_module.link_libc = true;

    const stats_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_stats.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });
    stats_test_module.link_libc = true;

    const check_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_check.zig"),
        .target = target,
        .optimize = optimize,
    });
    check_test_module.link_libc = true;

    const sample_internal_module = b.createModule(.{
        .root_source_file = b.path("src/sample.zig"),
        .target = target,
        .optimize = optimize,
    });
    const sample_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_sample.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "sample_internal", .module = sample_internal_module },
        },
    });
    sample_test_module.link_libc = true;

    const interleave_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_interleave.zig"),
        .target = target,
        .optimize = optimize,
    });
    interleave_test_module.link_libc = true;

    const deinterleave_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_deinterleave.zig"),
        .target = target,
        .optimize = optimize,
    });
    deinterleave_test_module.link_libc = true;

    for ([_]*std.Build.Module{
        reader_test_module,
        count_test_module,
        stats_test_module,
        check_test_module,
        sample_test_module,
        interleave_test_module,
        deinterleave_test_module,
    }) |module| {
        module.addOptions("test_options", test_options);
    }

    const run_fastq_test = b.addRunArtifact(b.addTest(.{ .root_module = fastq_test_module }));
    const run_main_test = b.addRunArtifact(b.addTest(.{ .root_module = exe.root_module }));
    const run_reader_test = b.addRunArtifact(b.addTest(.{ .root_module = reader_test_module }));
    const run_writer_test = b.addRunArtifact(b.addTest(.{ .root_module = writer_test_module }));

    const run_count_test = b.addRunArtifact(b.addTest(.{ .root_module = count_test_module }));
    run_count_test.step.dependOn(b.getInstallStep());
    const run_stats_test = b.addRunArtifact(b.addTest(.{ .root_module = stats_test_module }));
    run_stats_test.step.dependOn(b.getInstallStep());
    const run_check_test = b.addRunArtifact(b.addTest(.{ .root_module = check_test_module }));
    run_check_test.step.dependOn(b.getInstallStep());
    const run_sample_test = b.addRunArtifact(b.addTest(.{ .root_module = sample_test_module }));
    run_sample_test.step.dependOn(b.getInstallStep());
    const run_interleave_test = b.addRunArtifact(b.addTest(.{
        .root_module = interleave_test_module,
    }));
    run_interleave_test.step.dependOn(b.getInstallStep());
    const run_deinterleave_test = b.addRunArtifact(b.addTest(.{
        .root_module = deinterleave_test_module,
    }));
    run_deinterleave_test.step.dependOn(b.getInstallStep());

    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_fastq_test.step);
    test_step.dependOn(&run_main_test.step);
    test_step.dependOn(&run_reader_test.step);
    test_step.dependOn(&run_writer_test.step);
    test_step.dependOn(&run_count_test.step);
    test_step.dependOn(&run_stats_test.step);
    test_step.dependOn(&run_check_test.step);
    test_step.dependOn(&run_sample_test.step);
    test_step.dependOn(&run_interleave_test.step);
    test_step.dependOn(&run_deinterleave_test.step);
}
