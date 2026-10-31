#!/usr/bin/env bash
# Release checks and packages for z-fastq.
#
#   release.sh version                     print the version build.zig.zon declares
#   release.sh notes [TAG]                 print the CHANGELOG.md entry for the version (TAG must be vVERSION)
#   release.sh package source OUTPUT_DIR   write z-fastq-VERSION-source.tar.gz and test it from an outside package
#   release.sh package SYSTEM OUTPUT_DIR   write z-fastq-VERSION-SYSTEM.tar.gz and smoke-test the unpacked binary
#
# SYSTEM is linux-x86_64 or linux-aarch64, and must be the machine's own: the binary it packages is the one it
# smoke-tests. Work directories are made under RUNNER_TEMP (or /var/tmp) and left in place; nothing is
# deleted.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

version=$(sed -n 's/^    \.version = "\([^"]*\)",$/\1/p' build.zig.zon)
[[ "$version" =~ ^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$ ]] ||
    die "build.zig.zon must declare one major.minor.patch version"

workdir() {
    mktemp -d "${RUNNER_TEMP:-/var/tmp}/z-fastq-release.XXXXXX"
}

package_files() {
    local paths=()
    mapfile -t paths < <(awk '/^    \.paths = \.\{$/ { inside = 1; next } inside && /^    \},$/ { exit } inside' build.zig.zon |
        sed -n 's/^ *"\([^"]*\)",$/\1/p')
    [[ ${#paths[@]} -gt 0 ]] || die "build.zig.zon has no .paths"
    git ls-files -- "${paths[@]}" | LC_ALL=C sort
}

package_source() {
    local output="$1" work name="z-fastq-$version-source"
    work=$(workdir)
    package_files > "$work/files"
    tar -czf "$work/$name.tar.gz" --transform "s,^,$name/," -T "$work/files"
    tar -tzf "$work/$name.tar.gz" | sed "s,^$name/,," | LC_ALL=C sort | cmp -s - "$work/files" ||
        die "the source archive does not hold exactly the tracked .paths files"

    mkdir "$work/consumer"
    (
        cd "$work/consumer"
        cat > build.zig.zon <<'ZON'
.{
    .name = .z_fastq_package_consumer,
    .version = "0.0.0",
    .fingerprint = 0x92feba6be5e30c69,
    .minimum_zig_version = "0.16.0",
    .dependencies = .{},
    .paths = .{ "build.zig", "build.zig.zon", "main.zig" },
}
ZON
        cat > build.zig <<'ZIG'
const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const z_fastq = b.dependency("z_fastq", .{ .target = target, .optimize = optimize });
    const options = b.addOptions();
    options.addOption([]const u8, "version", b.option([]const u8, "expect-version", "Expected version").?);
    const module = b.createModule(.{
        .root_source_file = b.path("main.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "z-fastq", .module = z_fastq.module("z-fastq") }},
    });
    module.addOptions("expected", options);
    const tests = b.addTest(.{ .root_module = module, .use_llvm = true });
    b.step("test", "Test the z-fastq module from an outside package").dependOn(&b.addRunArtifact(tests).step);
}
ZIG
        cat > main.zig <<'ZIG'
const std = @import("std");
const zfastq = @import("z-fastq");

test "every export is analyzable and reports the packaged version" {
    std.testing.refAllDecls(zfastq);
    try std.testing.expectEqualStrings(@import("expected").version, zfastq.VERSION);
}

test "a gzip member decodes through the packaged Zipir" {
    const member = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x03, 0xf3, 0x48, 0xcd, 0xc9, 0xc9, 0x57,
        0x28, 0xcf, 0x2f, 0xca, 0x49, 0xe1, 0x02, 0x00, 0xd5, 0xe0, 0x39, 0xb7, 0x0c, 0x00, 0x00, 0x00,
    };
    var input: std.Io.Reader = .fixed(&member);
    var gzip = zfastq.io.gzip.ReaderSource.init(&input);
    var output: [32]u8 = undefined;
    var len: usize = 0;
    while (true) {
        const n = try gzip.byteSource().read(output[len..]);
        if (n == 0) break;
        len += n;
    }
    try std.testing.expectEqualStrings("Hello world\n", output[0..len]);
}
ZIG
        zig fetch --save=z_fastq "$work/$name.tar.gz"
        zig build test -Dexpect-version="$version" --summary all
    )
    mkdir -p "$output"
    cp "$work/$name.tar.gz" "$output/"
    printf 'tested %s\n' "$output/$name.tar.gz"
}

package_binary() {
    local system="$1" output="$2" work stage name binary machine
    case "$system" in
        linux-x86_64) machine=x86-64 ;;
        linux-aarch64) machine='ARM aarch64' ;;
        *) die "unknown system: $system" ;;
    esac
    [[ "$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)" == "$system" ]] ||
        die "$system must be packaged on $system, not on $(uname -s)-$(uname -m)"
    work=$(workdir)
    name="z-fastq-$version-$system"
    zig build -Doptimize=ReleaseFast --prefix "$work/install"
    stage="$work/$name"
    mkdir "$stage"
    cp "$work/install/bin/z-fastq" LICENSE "$stage/"
    file -b "$stage/z-fastq" | grep -q "$machine.*statically linked.*stripped" ||
        die "release binary is not a static, stripped $machine executable: $(file -b "$stage/z-fastq")"
    tar -czf "$work/$name.tar.gz" -C "$work" "$name"

    mkdir "$work/smoke"
    tar -xzf "$work/$name.tar.gz" -C "$work/smoke"
    binary="$work/smoke/$name/z-fastq"
    cmp LICENSE "$work/smoke/$name/LICENSE"
    (
        cd "$work/smoke"
        export PATH=/usr/bin:/bin
        ! command -v zig > /dev/null || die "zig is still on PATH"
        [[ "$("$binary" --version)" == "z-fastq $version" ]] || die "--version is not 'z-fastq $version'"
        "$binary" --help > /dev/null
        awk 'BEGIN {
            for (i = 0; i < 2000; i++) {
                printf "@r%d/1\nACGTACGTNNACGTAGATCGGAAGAGCACACG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n", i > "r1.fastq"
                printf "@r%d/2\nTTGCAACGTTGCAAGATCGGAAGAGCGTCGTG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n", i > "r2.fastq"
            }
        }'
        [[ "$("$binary" count r1.fastq r2.fastq)" == $'2000\n2000' ]] || die "count"
        "$binary" stats --json r1.fastq | grep -q '"reads":2000,"bases":64000,' || die "stats --json"
        "$binary" check --paired r1.fastq r2.fastq
        [[ "$("$binary" sample --count 10 --seed 1 r1.fastq | wc -l)" -eq 40 ]] || die "sample --count"
        "$binary" interleave r1.fastq r2.fastq > pairs.fastq
        "$binary" deinterleave --out1 d1.fastq --out2 d2.fastq pairs.fastq
        cmp r1.fastq d1.fastq
        cmp r2.fastq d2.fastq
        "$binary" process --gzip --min-length 5 --report report.json \
            --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --out1 o1.fastq.gz --out2 o2.fastq.gz r1.fastq r2.fastq
        grep -q '"schema":"z-fastq/process-v1"' report.json || die "process report schema"
        "$binary" stats --json o1.fastq.gz | grep -q '"reads":2000,"bases":28000,"min_length":14,"max_length":14,' ||
            die "process read 1 trimming"
        "$binary" stats --json o2.fastq.gz | grep -q '"reads":2000,"bases":26000,"min_length":13,"max_length":13,' ||
            die "process read 2 trimming"
        head -n 7 r1.fastq > truncated.fastq
        status=0
        "$binary" check truncated.fastq 2> truncated.err || status=$?
        if [[ "$status" -ne 1 ]] || ! grep -q 'S004' truncated.err; then
            die "a truncated record must fail check with exit 1"
        fi
    )
    mkdir -p "$output"
    cp "$work/$name.tar.gz" "$output/"
    printf 'tested %s\n' "$output/$name.tar.gz"
}

case "${1:-}" in
    version)
        printf '%s\n' "$version"
        ;;
    notes)
        tag="${2:-v$version}"
        [[ "$tag" == "v$version" ]] || die "tag $tag does not match version v$version"
        awk -v version="$version" '
            /^## / {
                if (found) exit
                found = ($2 == "[" version "]")
                next
            }
            found {
                if ($0 ~ /[^[:space:]]/) nonempty = 1
                if (nonempty) print
            }
            END { if (!found || !nonempty) exit 1 }
        ' CHANGELOG.md || die "CHANGELOG.md needs a nonempty [$version] entry"
        ;;
    package)
        [[ $# -eq 3 ]] || die "usage: release.sh package source|SYSTEM OUTPUT_DIR"
        if [[ "$2" == source ]]; then package_source "$3"; else package_binary "$2" "$3"; fi
        ;;
    *)
        printf 'usage: release.sh version | notes [TAG] | package source|SYSTEM OUTPUT_DIR\n' >&2
        exit 2
        ;;
esac
