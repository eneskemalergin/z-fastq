//! Tests for the `z-fastq` library module exports.

const std = @import("std");
const zfastq = @import("z-fastq");

test {
    std.testing.refAllDecls(zfastq);
}

test "root: version string is non-empty" {
    try std.testing.expect(zfastq.VERSION.len > 0);
}
