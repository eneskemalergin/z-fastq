//! Build graph for z-fastq.

const std = @import("std");

pub fn build(b: *std.Build) void {
    const requested_target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    switch (optimize) {
        .Debug, .ReleaseFast => {},
        else => {
            std.debug.print(
                "error: z-fastq builds only in Debug or ReleaseFast\n",
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
    const static_release = optimize == .ReleaseFast;
    if (b.option(bool, "static", "Confirm the static ReleaseFast build")) |requested_static| {
        if (requested_static != static_release) {
            std.debug.print(
                "error: Debug is the development build and ReleaseFast is always static\n",
                .{},
            );
            std.process.exit(1);
        }
    }
    const target = if (static_release)
        b.resolveTargetQuery(.{
            .cpu_arch = .x86_64,
            .os_tag = .linux,
            .abi = .musl,
        })
    else
        requested_target;
    const strip = optimize == .ReleaseFast;
    const isa_l_supported = target.result.cpu.arch == .x86_64 and
        target.result.os.tag == .linux and
        switch (target.result.abi) {
            .gnu, .musl => true,
            else => false,
        };
    const use_isa_l = b.option(
        bool,
        "isa-l",
        "Use the vendored ISA-L gzip engine",
    ) orelse isa_l_supported;
    if (use_isa_l and !isa_l_supported) {
        std.debug.print("error: ISA-L requires a supported Linux x86-64 target\n", .{});
        std.process.exit(1);
    }
    const build_options = b.addOptions();
    build_options.addOption(bool, "use_isa_l", use_isa_l);

    const lib_module = b.createModule(.{
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

    if (use_isa_l) {
        const isal = addIsaL(b, target, optimize);
        lib_module.addIncludePath(b.path("vendor/ISA-L/include"));
        lib_module.linkLibrary(isal);
        exe.root_module.addIncludePath(b.path("vendor/ISA-L/include"));
        exe.root_module.linkLibrary(isal);
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

    const fastq_test_module = b.createModule(.{
        .root_source_file = b.path("src/fastq.zig"),
        .target = target,
        .optimize = optimize,
    });
    const fastq_test_options = b.addOptions();
    fastq_test_options.addOption(bool, "use_isa_l", false);
    fastq_test_module.addOptions("build_options", fastq_test_options);

    const reader_test_module = b.createModule(.{
        .root_source_file = b.path("tests/test_reader.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &import_lib,
    });

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

fn addIsaL(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
) *std.Build.Step.Compile {
    const module = b.createModule(.{
        .target = target,
        .optimize = optimize,
    });
    module.addIncludePath(b.path("vendor/ISA-L/include"));
    module.addIncludePath(b.path("vendor/ISA-L/igzip"));
    module.addIncludePath(b.path("vendor/ISA-L/crc"));
    module.link_libc = true;
    module.addCSourceFiles(.{
        .root = b.path("vendor/ISA-L"),
        .files = &.{
            "crc/crc_base.c",
            "crc/crc64_base.c",
            "igzip/adler32_base.c",
            "igzip/hufftables_c.c",
            "igzip/igzip_inflate.c",
            "igzip/inflate_helpers.c",
        },
        .flags = &.{ "-O2", "-DNDEBUG", "-Wall", "-ffunction-sections" },
    });
    module.addObjectFile(addIsaLAssembly(b, "isa-l-rfc1951-lookup.o", "igzip/rfc1951_lookup.asm"));
    module.addObjectFile(addIsaLAssembly(b, "isa-l-inflate-dispatch.o", "igzip/igzip_inflate_multibinary.asm"));
    module.addObjectFile(addIsaLAssembly(b, "isa-l-inflate-01.o", "igzip/igzip_decode_block_stateless_01.asm"));
    module.addObjectFile(addIsaLAssembly(b, "isa-l-inflate-04.o", "igzip/igzip_decode_block_stateless_04.asm"));
    const crc_assembly = [_][]const u8{
        "crc16_t10dif_01.asm",
        "crc16_t10dif_avx2.asm",
        "crc16_t10dif_by16_10.asm",
        "crc16_t10dif_copy_by4.asm",
        "crc16_t10dif_copy_by4_02.asm",
        "crc32_gzip_refl_avx2.asm",
        "crc32_gzip_refl_by16_10.asm",
        "crc32_gzip_refl_by8.asm",
        "crc32_ieee_01.asm",
        "crc32_ieee_avx2.asm",
        "crc32_ieee_by16_10.asm",
        "crc32_iscsi_01.asm",
        "crc32_iscsi_avx2.asm",
        "crc32_iscsi_by16_10.asm",
        "crc32_iscsi_by8_02.asm",
        "crc64_ecma_norm_avx2.asm",
        "crc64_ecma_norm_by16_10.asm",
        "crc64_ecma_norm_by8.asm",
        "crc64_ecma_refl_avx2.asm",
        "crc64_ecma_refl_by16_10.asm",
        "crc64_ecma_refl_by8.asm",
        "crc64_iso_norm_avx2.asm",
        "crc64_iso_norm_by16_10.asm",
        "crc64_iso_norm_by8.asm",
        "crc64_iso_refl_avx2.asm",
        "crc64_iso_refl_by16_10.asm",
        "crc64_iso_refl_by8.asm",
        "crc64_jones_norm_avx2.asm",
        "crc64_jones_norm_by16_10.asm",
        "crc64_jones_norm_by8.asm",
        "crc64_jones_refl_avx2.asm",
        "crc64_jones_refl_by16_10.asm",
        "crc64_jones_refl_by8.asm",
        "crc64_rocksoft_norm_avx2.asm",
        "crc64_rocksoft_norm_by16_10.asm",
        "crc64_rocksoft_norm_by8.asm",
        "crc64_rocksoft_refl_avx2.asm",
        "crc64_rocksoft_refl_by16_10.asm",
        "crc64_rocksoft_refl_by8.asm",
        "crc_const.asm",
        "crc_multibinary.asm",
    };
    for (crc_assembly) |file| {
        module.addObjectFile(addIsaLAssembly(
            b,
            b.fmt("isa-l-{s}.o", .{file[0 .. file.len - ".asm".len]}),
            b.fmt("crc/{s}", .{file}),
        ));
    }

    return b.addLibrary(.{
        .name = "isa-l",
        .root_module = module,
    });
}

fn addIsaLAssembly(b: *std.Build, output_name: []const u8, source: []const u8) std.Build.LazyPath {
    const assemble = b.addSystemCommand(&.{
        "nasm",
        "-f",
        "elf64",
        "-DINTEL_CET_ENABLED",
        "-Ivendor/ISA-L/",
        "-Ivendor/ISA-L/crc/",
        "-Ivendor/ISA-L/igzip/",
        "-Ivendor/ISA-L/include/",
        "-o",
    });
    const output = assemble.addOutputFileArg(output_name);
    assemble.addFileArg(b.path(b.fmt("vendor/ISA-L/{s}", .{source})));
    const includes = [_][]const u8{
        "vendor/ISA-L/include/multibinary.asm",
        "vendor/ISA-L/include/reg_sizes.asm",
        "vendor/ISA-L/crc/crc_const_extern.asm",
        "vendor/ISA-L/include/crc.inc",
        "vendor/ISA-L/include/memcpy.asm",
        "vendor/ISA-L/igzip/igzip_decode_block_stateless.asm",
        "vendor/ISA-L/igzip/inflate_data_structs.asm",
        "vendor/ISA-L/igzip/stdmac.asm",
    };
    for (includes) |path| assemble.addFileInput(b.path(path));
    return output;
}
