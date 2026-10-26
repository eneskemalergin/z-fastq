//! CLI-private paired-read name parsing and matching.

const std = @import("std");
const simd = @import("simd.zig");

pub const NamePolicy = enum {
    illumina,
    exact,
};

pub const Name = struct {
    first_token: []const u8,
    normalized_id: []const u8,
    mate_markers: u2,
    first_mate_marker: ?u2,
};

pub fn parseName(header: []const u8, policy: NamePolicy) Name {
    const first_token = header[0..firstTokenEnd(header)];
    if (policy == .exact) {
        return .{
            .first_token = first_token,
            .normalized_id = first_token,
            .mate_markers = 0,
            .first_mate_marker = null,
        };
    }

    var normalized_id = first_token;
    var mate_markers: u2 = 0;
    const first_mate_marker = terminalMateMarker(first_token);
    if (first_mate_marker) |mate| {
        mate_markers |= mateMask(mate);
        normalized_id = first_token[0 .. first_token.len - 2];
    }

    const second_token = nextToken(header, first_token.len);
    if (leadingMateMarker(second_token)) |mate| {
        mate_markers |= mateMask(mate);
    }
    if (terminalMateMarker(second_token)) |mate| {
        mate_markers |= mateMask(mate);
    }
    return .{
        .first_token = first_token,
        .normalized_id = normalized_id,
        .mate_markers = mate_markers,
        .first_mate_marker = first_mate_marker,
    };
}

pub fn headersMatch(header1: []const u8, header2: []const u8, policy: NamePolicy) bool {
    return switch (policy) {
        .exact => exactHeadersMatch(header1, header2),
        .illumina => illuminaHeadersMatch(header1, header2),
    };
}

pub fn namesMatch(name1: Name, name2: Name) bool {
    if (name1.normalized_id.len == 0 or name2.normalized_id.len == 0) return false;
    if (!std.mem.eql(u8, name1.normalized_id, name2.normalized_id)) return false;
    if (name1.mate_markers == 0b11 or name2.mate_markers == 0b11) return false;
    if ((name1.mate_markers == 0) != (name2.mate_markers == 0)) return false;
    if (name1.mate_markers == 0) return true;
    return name1.mate_markers == 0b01 and name2.mate_markers == 0b10;
}

// Read names are short, so the token kernels have only the scalar and
// baseline tiers: no name would repay a call into the AVX2 object.
const TOKEN_LANES = simd.tierWidth(.baseline);
const TokenVector = @Vector(TOKEN_LANES, u8);
const TokenMask = @Int(.unsigned, TOKEN_LANES);

fn exactHeadersMatch(header1: []const u8, header2: []const u8) bool {
    return matchingTokenEnd(header1, header2) != null;
}

fn matchingTokenEnd(header1: []const u8, header2: []const u8) ?usize {
    const common_len = @min(header1.len, header2.len);
    const spaces: TokenVector = @splat(' ');
    const tabs: TokenVector = @splat('\t');

    var offset: usize = 0;
    while (common_len - offset >= TOKEN_LANES) : (offset += TOKEN_LANES) {
        const bytes1: TokenVector = header1[offset..][0..TOKEN_LANES].*;
        const bytes2: TokenVector = header2[offset..][0..TOKEN_LANES].*;
        const stop1: TokenMask = @bitCast((bytes1 == spaces) | (bytes1 == tabs));
        const stop2: TokenMask = @bitCast((bytes2 == spaces) | (bytes2 == tabs));
        const mismatch: TokenMask = @bitCast(bytes1 != bytes2);
        const event = stop1 | stop2 | mismatch;
        if (event != 0) {
            const lane = @as(usize, @intCast(@ctz(event)));
            const bit = @as(TokenMask, 1) << @intCast(lane);
            const token_end1 = stop1 & bit != 0;
            const token_end2 = stop2 & bit != 0;
            return if (offset + lane != 0 and token_end1 and token_end2) offset + lane else null;
        }
    }
    while (offset < common_len) : (offset += 1) {
        const token_end1 = header1[offset] == ' ' or header1[offset] == '\t';
        const token_end2 = header2[offset] == ' ' or header2[offset] == '\t';
        if (token_end1 or token_end2) {
            return if (offset != 0 and token_end1 and token_end2) offset else null;
        }
        if (header1[offset] != header2[offset]) return null;
    }

    if (common_len == 0) return null;
    if (header1.len == header2.len) return common_len;
    const longer = if (header1.len > header2.len) header1 else header2;
    return if (longer[common_len] == ' ' or longer[common_len] == '\t') common_len else null;
}

fn firstTokenEnd(header: []const u8) usize {
    return tokenEnd(header, 0);
}

fn tokenEnd(header: []const u8, start: usize) usize {
    const spaces: TokenVector = @splat(' ');
    const tabs: TokenVector = @splat('\t');

    var end = start;
    while (header.len - end >= TOKEN_LANES) : (end += TOKEN_LANES) {
        const bytes: TokenVector = header[end..][0..TOKEN_LANES].*;
        const mask: TokenMask = @bitCast((bytes == spaces) | (bytes == tabs));
        if (mask != 0) return end + @as(usize, @intCast(@ctz(mask)));
    }
    while (end < header.len and header[end] != ' ' and header[end] != '\t') {
        end += 1;
    }
    return end;
}

fn nextToken(header: []const u8, previous_end: usize) []const u8 {
    var start = previous_end;
    while (start < header.len and (header[start] == ' ' or header[start] == '\t')) {
        start += 1;
    }
    const end = tokenEnd(header, start);
    return header[start..end];
}

fn illuminaHeadersMatch(header1: []const u8, header2: []const u8) bool {
    return terminalPairHeadersMatch(header1, header2) or
        casavaHeadersMatch(header1, header2) or
        namesMatch(parseName(header1, .illumina), parseName(header2, .illumina));
}

fn casavaHeadersMatch(header1: []const u8, header2: []const u8) bool {
    const end = matchingTokenEnd(header1, header2) orelse return false;
    if (terminalMateMarker(header1[0..end]) != null) return false;
    const token1 = nextToken(header1, end);
    const token2 = nextToken(header2, end);
    return leadingMateMarker(token1) == 1 and leadingMateMarker(token2) == 2 and
        terminalMateMarker(token1) != 2 and terminalMateMarker(token2) != 1;
}

fn terminalPairHeadersMatch(header1: []const u8, header2: []const u8) bool {
    if (header1.len < 3 or header1.len != header2.len) return false;
    if (header1[header1.len - 2] != '/' or header2[header2.len - 2] != '/') return false;
    if (header1[header1.len - 1] != '1' or header2[header2.len - 1] != '2') return false;

    const prefix_len = header1.len - 1;
    var separator_plus_one: usize = 0;
    var offset: usize = 0;
    while (prefix_len - offset >= TOKEN_LANES) : (offset += TOKEN_LANES) {
        const bytes1: TokenVector = header1[offset..][0..TOKEN_LANES].*;
        const bytes2: TokenVector = header2[offset..][0..TOKEN_LANES].*;
        const spaces: TokenVector = @splat(' ');
        const tabs: TokenVector = @splat('\t');
        if (@reduce(.Or, (bytes1 != bytes2) | (bytes1 == tabs))) return false;

        const space_mask: TokenMask = @bitCast(bytes1 == spaces);
        if (space_mask != 0) {
            if (separator_plus_one != 0 or space_mask & (space_mask - 1) != 0) return false;
            separator_plus_one = offset + @as(usize, @intCast(@ctz(space_mask))) + 1;
        }
    }

    inline for (.{ 8, 4 }) |lanes| {
        if (prefix_len - offset >= lanes) {
            const Vector = @Vector(lanes, u8);
            const Mask = @Int(.unsigned, lanes);
            const bytes1: Vector = header1[offset..][0..lanes].*;
            const bytes2: Vector = header2[offset..][0..lanes].*;
            const spaces: Vector = @splat(' ');
            const tabs: Vector = @splat('\t');
            if (@reduce(.Or, (bytes1 != bytes2) | (bytes1 == tabs))) return false;

            const space_mask: Mask = @bitCast(bytes1 == spaces);
            if (space_mask != 0) {
                if (separator_plus_one != 0 or space_mask & (space_mask - 1) != 0) {
                    return false;
                }
                separator_plus_one = offset + @as(usize, @intCast(@ctz(space_mask))) + 1;
            }
            offset += lanes;
        }
    }
    while (offset < prefix_len) : (offset += 1) {
        if (header1[offset] != header2[offset] or header1[offset] == '\t') return false;
        if (header1[offset] == ' ') {
            if (separator_plus_one != 0) return false;
            separator_plus_one = offset + 1;
        }
    }

    if (separator_plus_one == 0) return true;
    const separator = separator_plus_one - 1;
    if (separator == 0) return false;

    const first_token = header1[0..separator];
    if (terminalMateMarker(first_token) != null) return false;

    return leadingMateMarker(header1[separator_plus_one..]) == null;
}

fn leadingMateMarker(token: []const u8) ?u2 {
    if (token.len < 2 or token[1] != ':') return null;
    return switch (token[0]) {
        '1' => 1,
        '2' => 2,
        else => null,
    };
}

fn terminalMateMarker(token: []const u8) ?u2 {
    if (token.len < 2 or token[token.len - 2] != '/') return null;
    return switch (token[token.len - 1]) {
        '1' => 1,
        '2' => 2,
        else => null,
    };
}

fn mateMask(mate: u2) u2 {
    return if (mate == 1) 0b01 else 0b10;
}

test "[property] - [paired names]: illumina and exact policies accept the listed header pairs" {
    const Row = struct { header1: []const u8, header2: []const u8, policy: NamePolicy, matches: bool };
    const rows = [_]Row{
        .{ .header1 = "cluster", .header2 = "cluster", .policy = .illumina, .matches = true },
        .{ .header1 = "cluster/1", .header2 = "cluster/2", .policy = .illumina, .matches = true },
        .{ .header1 = "cluster/2", .header2 = "cluster/1", .policy = .illumina, .matches = false },
        .{ .header1 = "cluster", .header2 = "cluster/1", .policy = .illumina, .matches = false },
        .{ .header1 = "cluster/1", .header2 = "other/2", .policy = .illumina, .matches = false },
        .{ .header1 = "cluster description/1", .header2 = "cluster description/2", .policy = .illumina, .matches = true },
        .{ .header1 = "cluster 1:N:0:index-a", .header2 = "cluster 2:N:0:index-b", .policy = .illumina, .matches = true },
        .{ .header1 = "cluster 2:N:0:index-b", .header2 = "cluster 1:N:0:index-a", .policy = .illumina, .matches = false },
        .{ .header1 = "cluster 1:N:0:index-a", .header2 = "cluster 1:N:0:index-a", .policy = .illumina, .matches = false },
        .{ .header1 = "cluster/1 1:N:0:x", .header2 = "cluster/2 2:N:0:x", .policy = .illumina, .matches = true },
        .{ .header1 = "cluster", .header2 = "cluster", .policy = .exact, .matches = true },
        .{ .header1 = "cluster", .header2 = "cluster/1", .policy = .exact, .matches = false },
        .{ .header1 = "cluster/1", .header2 = "cluster/2", .policy = .exact, .matches = false },
        .{ .header1 = "cluster description/1", .header2 = "cluster description/2", .policy = .exact, .matches = true },
        .{ .header1 = "cluster\tdescription/1", .header2 = "cluster description/2", .policy = .exact, .matches = true },
        .{ .header1 = "cluster 1:N:0", .header2 = "cluster 2:N:0", .policy = .exact, .matches = true },
        .{ .header1 = "", .header2 = "", .policy = .exact, .matches = false },
    };
    for (rows) |row| {
        try std.testing.expectEqual(row.matches, headersMatch(row.header1, row.header2, row.policy));
    }
}

test "[property] - [paired names]: illumina matching keeps mate markers and opaque suffixes" {
    const Row = struct { suffix1: []const u8, suffix2: []const u8, matches: bool };
    // Each row matches once a non-empty prefix precedes the suffixes; an empty
    // prefix leaves the first token empty, which never matches.
    const rows = [_]Row{
        .{ .suffix1 = " 1:N:0:index-a", .suffix2 = "\t2:N:0:index-b", .matches = true },
        .{ .suffix1 = " 2:N:0:index-b", .suffix2 = " 1:N:0:index-a", .matches = false },
        .{ .suffix1 = " 1:N:0:index-a", .suffix2 = " 1:N:0:index-a", .matches = false },
        .{ .suffix1 = " 1:/1", .suffix2 = " 2:/2", .matches = true },
        .{ .suffix1 = " 1:/2", .suffix2 = " 2:/1", .matches = false },
        .{ .suffix1 = " 1:", .suffix2 = " 2:", .matches = true },
        .{ .suffix1 = " 1", .suffix2 = " 2", .matches = true },
        .{ .suffix1 = "", .suffix2 = "", .matches = true },
        .{ .suffix1 = " 1:\x00\xff/1 extra/2", .suffix2 = " 2:\x00/2 extra/1", .matches = true },
        .{ .suffix1 = " description/1", .suffix2 = " description/2", .matches = true },
    };
    var storage1: [128]u8 = undefined;
    var storage2: [128]u8 = undefined;
    for ([_]usize{ 0, 1, 15, 16, 17, 33 }) |length| {
        @memset(storage1[0..length], 'x');
        @memset(storage2[0..length], 'x');
        for (rows) |row| {
            @memcpy(storage1[length..][0..row.suffix1.len], row.suffix1);
            @memcpy(storage2[length..][0..row.suffix2.len], row.suffix2);
            const header1 = storage1[0 .. length + row.suffix1.len];
            const header2 = storage2[0 .. length + row.suffix2.len];
            try std.testing.expectEqual(row.matches and length != 0, headersMatch(header1, header2, .illumina));
        }
    }

    const left = " 1:N:0:left";
    const right = "\t2:N:0:right/2 ignored/1";
    for ([_]usize{ 0, 1, 16, 17 }) |length| {
        @memset(storage1[0..length], 'x');
        @memset(storage2[0..length], 'x');
        @memcpy(storage1[length..][0..left.len], left);
        @memcpy(storage2[length..][0..right.len], right);
        const header1 = storage1[0 .. length + left.len];
        const header2 = storage2[0 .. length + right.len];
        try std.testing.expectEqual(length != 0, casavaHeadersMatch(header1, header2));
        for (0..length) |mismatch| {
            storage2[mismatch] = 'y';
            try std.testing.expect(!headersMatch(header1, header2, .illumina));
            storage2[mismatch] = 'x';
        }
    }
}

test "[edge] - [paired names]: token finder stops at delimiters across vector boundaries" {
    const lengths = [_]usize{ 0, 1, 15, 16, 17, 33 };
    const starts = [_]usize{ 0, 1, 16 };
    const stops = [_]usize{ 0, 1, 15, 16, 17, 32 };
    var bytes: [33]u8 = undefined;

    for (lengths) |length| {
        for (starts) |start| {
            if (start > length) continue;
            @memset(bytes[0..length], 'x');
            try std.testing.expectEqual(length, tokenEnd(bytes[0..length], start));

            for (stops) |stop| {
                if (stop < start or stop >= length) continue;
                for (" \t") |delimiter| {
                    @memset(bytes[0..length], 'x');
                    bytes[stop] = delimiter;
                    try std.testing.expectEqual(stop, tokenEnd(bytes[0..length], start));
                }
            }
        }
    }
}

test "[edge] - [paired names]: exact matching compares the first token across vector boundaries" {
    const Row = struct { header1: []const u8, header2: []const u8, matches: bool };
    const rows = [_]Row{
        .{ .header1 = "cluster", .header2 = "cluster", .matches = true },
        .{ .header1 = "cluster", .header2 = "clusters", .matches = false },
        .{ .header1 = "cluster", .header2 = "cluster/1", .matches = false },
        .{ .header1 = "cluster/1", .header2 = "cluster/2", .matches = false },
        .{ .header1 = "cluster 1", .header2 = "cluster\t2", .matches = true },
        .{ .header1 = "cluster\t", .header2 = "cluster ", .matches = true },
        .{ .header1 = " value", .header2 = " value", .matches = false },
        .{ .header1 = "0123456789abcdef a", .header2 = "0123456789abcdef\tb", .matches = true },
        .{ .header1 = "0123456789abcdef a", .header2 = "0123456789abcdeg a", .matches = false },
        .{ .header1 = "0123456789abcdef x", .header2 = "0123456789abcdef", .matches = true },
        .{ .header1 = "0123456789abcdef", .header2 = "0123456789abcdef1", .matches = false },
        .{ .header1 = "", .header2 = "", .matches = false },
        .{ .header1 = "", .header2 = " value", .matches = false },
    };
    for (rows) |row| {
        try std.testing.expectEqual(row.matches, exactHeadersMatch(row.header1, row.header2));
        try std.testing.expectEqual(row.matches, headersMatch(row.header1, row.header2, .exact));
    }
}

test "[property] - [paired names]: combined terminal proof covers vector boundaries" {
    const lengths = [_]usize{ 3, 4, 15, 16, 17, 33 };
    const separators = [_]usize{ 1, 15, 16, 17 };
    var header1 = [_]u8{'x'} ** 33;
    var header2 = [_]u8{'x'} ** 33;

    for (lengths) |length| {
        @memset(header1[0..length], 'x');
        @memset(header2[0..length], 'x');
        header1[length - 2] = '/';
        header2[length - 2] = '/';
        header1[length - 1] = '1';
        header2[length - 1] = '2';
        try std.testing.expect(terminalPairHeadersMatch(
            header1[0..length],
            header2[0..length],
        ));

        for (separators) |separator| {
            if (separator < 1 or separator >= length - 2) continue;
            header1[separator] = ' ';
            header2[separator] = ' ';
            try std.testing.expect(terminalPairHeadersMatch(
                header1[0..length],
                header2[0..length],
            ));
            header1[separator] = 'x';
            header2[separator] = 'x';
        }
    }

    const fallback_pairs = [_][2][]const u8{
        .{ "cluster  lane/1", "cluster  lane/2" },
        .{ "cluster\tlane/1", "cluster\tlane/2" },
        .{ "cluster left/1", "cluster right/2" },
        .{ "cluster/1 lane/1", "cluster/2 lane/2" },
    };
    for (fallback_pairs) |headers| {
        try std.testing.expect(!terminalPairHeadersMatch(headers[0], headers[1]));
        try std.testing.expect(headersMatch(headers[0], headers[1], .illumina));
    }

    try std.testing.expect(!terminalPairHeadersMatch(" 1:N:0/1", " 1:N:0/2"));
    try std.testing.expect(!terminalPairHeadersMatch("cluster 1:N:0/1", "cluster 1:N:0/2"));
}
