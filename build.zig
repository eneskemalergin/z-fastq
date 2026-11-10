//! Build graph for z-fastq. Every step compiles with LLVM, Debug included: with Zig 0.16's
//! self-hosted x86-64 backend, a Debug build whose target has PCLMUL compiles Zipir's PCLMUL
//! CRC32 kernel inline, and the result corrupts later state in tests.

const std = @import("std");

const KernelBackend = enum { dispatch, portable };

pub fn build(b: *std.Build) void {
    const requested_query = b.standardTargetOptionsQueryOnly(.{});
    const requested_target = b.resolveTargetQuery(requested_query);
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
    const os = requested_target.result.os.tag;
    if ((arch != .x86_64 and arch != .aarch64) or (os != .linux and os != .macos)) {
        std.debug.print(
            "error: z-fastq supports x86-64 and aarch64 builds for Linux and macOS only\n",
            .{},
        );
        std.process.exit(1);
    }
    const release_build = optimize == .ReleaseSafe or optimize == .ReleaseFast;
    const static_build = release_build and os == .linux;
    if (b.option(bool, "static", "Confirm a static Linux release build")) |requested_static| {
        if (requested_static != static_build) {
            std.debug.print("error: {s}\n", .{if (os == .macos)
                "macOS builds link libSystem dynamically and are never static"
            else
                "Debug is native; ReleaseSafe and ReleaseFast are always static"});
            std.process.exit(1);
        }
    }
    const target = if (release_build)
        b.resolveTargetQuery(.{
            .cpu_arch = arch,
            .cpu_model = requested_query.cpu_model,
            .cpu_features_add = requested_query.cpu_features_add,
            .cpu_features_sub = requested_query.cpu_features_sub,
            .os_tag = os,
            .abi = if (os == .linux) .musl else null,
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
    const x86_dispatch = kernel_backend == .dispatch and arch == .x86_64;
    const target_features = target.result.cpu.features;
    const link_avx2_object = x86_dispatch and !std.Target.x86.featureSetHas(target_features, .avx2);
    const link_avx512_object = x86_dispatch and
        !std.Target.x86.featureSetHasAll(target_features, .{ .avx512f, .avx512bw, .avx512vl, .avx512dq });
    build_options.addOption(bool, "avx2_object", link_avx2_object);
    build_options.addOption(bool, "avx512_object", link_avx512_object);
    build_options.addOption(bool, "simd_object", false);
    const kernel_objects = [_]?*std.Build.Step.Compile{
        if (link_avx2_object) kernelObject(b, target, optimize, strip, "simd_avx2", .baseline, &.{.avx2}, &.{}) else null,
        if (link_avx512_object) kernelObject(b, target, optimize, strip, "simd_avx512", .{ .explicit = &std.Target.x86.cpu.x86_64_v4 }, &.{}, &.{.prefer_256_bit}) else null,
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
    if (arch == .x86_64 and os == .linux and
        std.Target.x86.featureSetHasAll(target_features, .{ .avx2, .bmi, .bmi2, .fma, .f16c, .lzcnt, .movbe }))
    {
        const require_x86_64_v4 = std.Target.x86.featureSetHasAll(target_features, .{ .avx512f, .avx512bw, .avx512vl, .avx512dq });
        exe.root_module.addObject(x86_64Entry(b, target, optimize, strip, require_x86_64_v4));
        exe.entry = .{ .symbol_name = "zfq_x86_64_entry" };
    }

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
        for (kernel_objects) |kernel_object| if (kernel_object) |object| module.addObject(object);
    }

    const test_exec = b.option(
        []const []const u8,
        "test-exec",
        "Run each test binary, and each CLI process the tests start, through this command (one argument per option)",
    ) orelse &.{};
    const test_exec_cmd: []const ?[]const u8 = blk: {
        const cmd = b.allocator.alloc(?[]const u8, test_exec.len + 1) catch @panic("OOM");
        for (test_exec, 0..) |arg, index| cmd[index] = arg;
        cmd[test_exec.len] = null;
        break :blk cmd;
    };
    const test_options = b.addOptions();
    test_options.addOption([]const u8, "package_version", package_version);
    test_options.addOption([]const []const u8, "cli_exec", test_exec);
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
    const source_test = b.addTest(.{ .root_module = exe.root_module, .use_llvm = true });
    if (test_exec.len > 0) source_test.setExecCmd(test_exec_cmd);
    const source_tests = b.addRunArtifact(source_test);
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
        const suite_test = b.addTest(.{ .root_module = module, .use_llvm = true });
        if (test_exec.len > 0) suite_test.setExecCmd(test_exec_cmd);
        const run = b.addRunArtifact(suite_test);
        if (suite.runs_binary) {
            module.link_libc = true;
            run.step.dependOn(b.getInstallStep());
        }
        test_step.dependOn(&run.step);
        const name = std.fs.path.stem(suite.path)["test_".len..];
        b.step(b.fmt("test-{s}", .{name}), b.fmt("Run {s}", .{suite.path})).dependOn(&run.step);
    }
}

fn x86_64Entry(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    strip: bool,
    require_x86_64_v4: bool,
) *std.Build.Step.Compile {
    var query = target.query;
    query.cpu_model = .baseline;
    query.cpu_features_add = .empty;
    query.cpu_features_sub = .empty;
    const entry_options = b.addOptions();
    entry_options.addOption(bool, "simd_object", false);
    entry_options.addOption(bool, "avx2_object", false);
    entry_options.addOption(bool, "avx512_object", false);
    entry_options.addOption(bool, "require_x86_64_v4", require_x86_64_v4);
    const module = b.createModule(.{
        .root_source_file = b.path("src/x86_64_entry.zig"),
        .target = b.resolveTargetQuery(query),
        .optimize = optimize,
        .single_threaded = true,
        .strip = strip,
    });
    module.addOptions("build_options", entry_options);
    return b.addObject(.{ .name = "x86_64_entry", .root_module = module, .use_llvm = true });
}

fn kernelObject(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    strip: bool,
    name: []const u8,
    cpu_model: std.Target.Query.CpuModel,
    add: []const std.Target.x86.Feature,
    sub: []const std.Target.x86.Feature,
) *std.Build.Step.Compile {
    var query = target.query;
    query.cpu_model = cpu_model;
    query.cpu_features_add = .empty;
    query.cpu_features_sub = .empty;
    for (add) |feature| query.cpu_features_add.addFeature(@intFromEnum(feature));
    for (sub) |feature| query.cpu_features_sub.addFeature(@intFromEnum(feature));
    const object_options = b.addOptions();
    object_options.addOption(bool, "simd_object", true);
    object_options.addOption(bool, "avx2_object", false);
    object_options.addOption(bool, "avx512_object", false);
    const module = b.createModule(.{
        .root_source_file = b.path("src/simd.zig"),
        .target = b.resolveTargetQuery(query),
        .optimize = optimize,
        .single_threaded = true,
        .strip = strip,
    });
    module.addOptions("build_options", object_options);
    return b.addObject(.{ .name = name, .root_module = module, .use_llvm = true });
}
