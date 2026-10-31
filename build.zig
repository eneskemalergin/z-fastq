//! Build graph for z-fastq. Every step compiles with LLVM, Debug included: with Zig 0.16's
//! self-hosted x86-64 backend, a Debug build whose target has PCLMUL compiles Zipir's PCLMUL
//! CRC32 kernel inline, and the result corrupts later state in tests.

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
    const arch = requested_target.result.cpu.arch;
    if ((arch != .x86_64 and arch != .aarch64) or requested_target.result.os.tag != .linux) {
        std.debug.print(
            "error: z-fastq currently supports Linux x86-64 and aarch64 builds only\n",
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
            .cpu_arch = arch,
            .os_tag = .linux,
            .abi = .musl,
        })
    else
        requested_target;
    const strip = b.option(bool, "strip", "Strip debug information (default: on for ReleaseFast)") orelse (optimize == .ReleaseFast);
    const kernel_backend = b.option(
        KernelBackend,
        "kernel-backend",
        "SIMD kernels: dispatch picks them for the running CPU; portable forces portable code",
    ) orelse .dispatch;
    const zipir = b.dependency("zipir", .{
        .target = target,
        .optimize = optimize,
        .@"kernel-backend" = kernel_backend,
    }).module("zipir");
    const package_version = @import("build.zig.zon").version;
    const build_options = b.addOptions();
    build_options.addOption([:0]const u8, "version", package_version);
    const link_avx2_object = kernel_backend == .dispatch and arch == .x86_64 and
        !std.Target.x86.featureSetHas(target.result.cpu.features, .avx2);
    build_options.addOption(bool, "avx2_object", link_avx2_object);
    build_options.addOption(bool, "simd_object", false);
    const avx2_object: ?*std.Build.Step.Compile = if (!link_avx2_object) null else blk: {
        var query = target.query;
        query.cpu_model = .baseline;
        query.cpu_features_add = .empty;
        query.cpu_features_sub = .empty;
        query.cpu_features_add.addFeature(@intFromEnum(std.Target.x86.Feature.avx2));
        const object_options = b.addOptions();
        object_options.addOption(bool, "simd_object", true);
        object_options.addOption(bool, "avx2_object", false);
        const module = b.createModule(.{
            .root_source_file = b.path("src/simd.zig"),
            .target = b.resolveTargetQuery(query),
            .optimize = optimize,
            .single_threaded = true,
            .strip = strip,
        });
        module.addOptions("build_options", object_options);
        break :blk b.addObject(.{ .name = "simd_avx2", .root_module = module, .use_llvm = true });
    };

    const lib_module = b.addModule("z-fastq", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });
    lib_module.addOptions("build_options", build_options);

    const exe = b.addExecutable(.{
        .name = "z-fastq",
        .use_llvm = true,
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .single_threaded = true,
            .strip = strip,
            .unwind_tables = if (strip) .none else null,
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

    for ([_]*std.Build.Module{ lib_module, exe.root_module }) |module| {
        module.addImport("zipir", zipir);
        if (avx2_object) |object| module.addObject(object);
    }

    const test_options = b.addOptions();
    test_options.addOption([]const u8, "package_version", package_version);
    const sample_internal_module = b.createModule(.{
        .root_source_file = b.path("src/sample.zig"),
        .target = target,
        .optimize = optimize,
    });
    const Suite = struct {
        path: []const u8,
        imports: []const std.Build.Module.Import = &.{},
        runs_binary: bool = true,
    };
    const suites = [_]Suite{
        .{ .path = "tests/test_reader.zig", .imports = &import_lib, .runs_binary = false },
        .{ .path = "tests/test_writer.zig", .imports = &import_lib, .runs_binary = false },
        .{ .path = "tests/test_count.zig", .imports = &import_lib },
        .{ .path = "tests/test_stats.zig", .imports = &import_lib },
        .{ .path = "tests/test_check.zig" },
        .{ .path = "tests/test_sample.zig", .imports = &.{.{ .name = "sample_internal", .module = sample_internal_module }} },
        .{ .path = "tests/test_interleave.zig" },
        .{ .path = "tests/test_deinterleave.zig" },
        .{ .path = "tests/test_process.zig" },
    };

    const test_step = b.step("test", "Run unit tests");
    const source_tests = b.addRunArtifact(b.addTest(.{ .root_module = exe.root_module, .use_llvm = true }));
    test_step.dependOn(&source_tests.step);
    b.step("test-src", "Run the in-file tests of src/").dependOn(&source_tests.step);
    for (suites) |suite| {
        const module = b.createModule(.{
            .root_source_file = b.path(suite.path),
            .target = target,
            .optimize = optimize,
            .imports = suite.imports,
        });
        module.addOptions("test_options", test_options);
        const run = b.addRunArtifact(b.addTest(.{ .root_module = module, .use_llvm = true }));
        if (suite.runs_binary) {
            module.link_libc = true;
            run.step.dependOn(b.getInstallStep());
        }
        test_step.dependOn(&run.step);
        const name = std.fs.path.stem(suite.path)["test_".len..];
        b.step(b.fmt("test-{s}", .{name}), b.fmt("Run {s}", .{suite.path})).dependOn(&run.step);
    }
}
