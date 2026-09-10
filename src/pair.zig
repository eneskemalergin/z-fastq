//! CLI-private paired-read name parsing and matching.

const std = @import("std");

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

const token_lanes = 16;
const TokenVector = @Vector(token_lanes, u8);

fn exactHeadersMatch(header1: []const u8, header2: []const u8) bool {
    const common_len = @min(header1.len, header2.len);
    const spaces: TokenVector = @splat(' ');
    const tabs: TokenVector = @splat('\t');

    var offset: usize = 0;
    while (common_len - offset >= token_lanes) : (offset += token_lanes) {
        const bytes1: TokenVector = header1[offset..][0..token_lanes].*;
        const bytes2: TokenVector = header2[offset..][0..token_lanes].*;
        const stop1: u16 = @bitCast((bytes1 == spaces) | (bytes1 == tabs));
        const stop2: u16 = @bitCast((bytes2 == spaces) | (bytes2 == tabs));
        const mismatch: u16 = @bitCast(bytes1 != bytes2);
        const event = stop1 | stop2 | mismatch;
        if (event != 0) {
            const lane = @as(usize, @intCast(@ctz(event)));
            const bit = @as(u16, 1) << @intCast(lane);
            const token_end1 = stop1 & bit != 0;
            const token_end2 = stop2 & bit != 0;
            return offset + lane != 0 and token_end1 and token_end2;
        }
    }
    while (offset < common_len) : (offset += 1) {
        const token_end1 = header1[offset] == ' ' or header1[offset] == '\t';
        const token_end2 = header2[offset] == ' ' or header2[offset] == '\t';
        if (token_end1 or token_end2) {
            return offset != 0 and token_end1 and token_end2;
        }
        if (header1[offset] != header2[offset]) return false;
    }

    if (common_len == 0) return false;
    if (header1.len == header2.len) return true;
    const longer = if (header1.len > header2.len) header1 else header2;
    return longer[common_len] == ' ' or longer[common_len] == '\t';
}

fn firstTokenEnd(header: []const u8) usize {
    return tokenEnd(header, 0);
}

fn tokenEnd(header: []const u8, start: usize) usize {
    const spaces: TokenVector = @splat(' ');
    const tabs: TokenVector = @splat('\t');

    var end = start;
    while (header.len - end >= token_lanes) : (end += token_lanes) {
        const bytes: TokenVector = header[end..][0..token_lanes].*;
        const mask: u16 = @bitCast((bytes == spaces) | (bytes == tabs));
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

// Keep the full matcher shared across paired command paths to avoid code growth.
noinline fn illuminaHeadersMatch(header1: []const u8, header2: []const u8) bool {
    return terminalPairHeadersMatch(header1, header2) or
        namesMatch(parseName(header1, .illumina), parseName(header2, .illumina));
}

fn terminalPairHeadersMatch(header1: []const u8, header2: []const u8) bool {
    if (header1.len < 3 or header1.len != header2.len) return false;
    if (header1[header1.len - 2] != '/' or header2[header2.len - 2] != '/') return false;
    if (header1[header1.len - 1] != '1' or header2[header2.len - 1] != '2') return false;

    const prefix_len = header1.len - 1;
    var separator_plus_one: usize = 0;
    var offset: usize = 0;
    while (prefix_len - offset >= token_lanes) : (offset += token_lanes) {
        const bytes1: TokenVector = header1[offset..][0..token_lanes].*;
        const bytes2: TokenVector = header2[offset..][0..token_lanes].*;
        const spaces: TokenVector = @splat(' ');
        const tabs: TokenVector = @splat('\t');
        if (@reduce(.Or, (bytes1 != bytes2) | (bytes1 == tabs))) return false;

        const space_mask: u16 = @bitCast(bytes1 == spaces);
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

test "[property] - [paired names]: success-first matching preserves established policies" {
    const headers = [_][]const u8{
        "cluster",
        "cluster/1",
        "cluster/2",
        "cluster description/1",
        "cluster description/2",
        "cluster\tdescription/1",
        "cluster\tdescription/2",
        "cluster  1:N:0:1",
        "cluster  2:N:0:1",
        "cluster 1:/1",
        "cluster 1:/2",
        "cluster/1 description/1",
        "cluster/1 description/2",
        "cluster description note/1",
        "cluster description note/2",
        "other/1",
        "other/2",
        "/1",
        "/2",
        "",
    };
    const policies = [_]NamePolicy{ .illumina, .exact };

    for (policies) |policy| {
        for (headers) |header1| {
            for (headers) |header2| {
                const expected = namesMatch(
                    parseName(header1, policy),
                    parseName(header2, policy),
                );
                try std.testing.expectEqual(
                    expected,
                    headersMatch(header1, header2, policy),
                );
            }
        }
    }
}

test "[property] - [paired names]: token finder matches scalar boundaries" {
    const lengths = [_]usize{ 0, 1, 15, 16, 17, 31, 32, 33, 65 };
    var bytes: [65]u8 = undefined;

    for (lengths) |length| {
        for (0..length + 1) |start| {
            @memset(bytes[0..length], 'x');
            try std.testing.expectEqual(length, tokenEnd(bytes[0..length], start));

            for (start..length) |stop| {
                for (" \t") |delimiter| {
                    @memset(bytes[0..length], 'x');
                    bytes[stop] = delimiter;
                    try std.testing.expectEqual(stop, tokenEnd(bytes[0..length], start));
                }
            }
        }
    }
}

test "[property] - [paired names]: exact comparison preserves parsed token equality" {
    const lengths = [_]usize{ 1, 15, 16, 17, 31, 32, 33, 63, 64, 65 };
    var header1 = [_]u8{'x'} ** 68;
    var header2 = [_]u8{'x'} ** 68;

    for (lengths) |length| {
        @memset(&header1, 'x');
        @memset(&header2, 'x');
        try std.testing.expectEqual(
            namesMatch(parseName(header1[0..length], .exact), parseName(header2[0..length], .exact)),
            exactHeadersMatch(header1[0..length], header2[0..length]),
        );

        header1[length] = ' ';
        header2[length] = '\t';
        header1[length + 1] = 'a';
        header2[length + 1] = 'b';
        try std.testing.expectEqual(
            namesMatch(parseName(header1[0 .. length + 2], .exact), parseName(header2[0 .. length + 2], .exact)),
            exactHeadersMatch(header1[0 .. length + 2], header2[0 .. length + 2]),
        );
        try std.testing.expectEqual(
            namesMatch(parseName(header1[0..length], .exact), parseName(header2[0 .. length + 2], .exact)),
            exactHeadersMatch(header1[0..length], header2[0 .. length + 2]),
        );

        header1[length] = 'x';
        header2[length] = 'x';
        try std.testing.expectEqual(
            namesMatch(parseName(header1[0..length], .exact), parseName(header2[0 .. length + 1], .exact)),
            exactHeadersMatch(header1[0..length], header2[0 .. length + 1]),
        );

        for (0..length) |mismatch| {
            header2[mismatch] = 'y';
            try std.testing.expectEqual(
                namesMatch(parseName(header1[0..length], .exact), parseName(header2[0..length], .exact)),
                exactHeadersMatch(header1[0..length], header2[0..length]),
            );
            header2[mismatch] = 'x';
        }
    }

    try std.testing.expect(!exactHeadersMatch("", ""));
    try std.testing.expect(!exactHeadersMatch("", " value"));
}

test "[property] - [paired names]: combined terminal proof covers vector boundaries" {
    const lengths = [_]usize{ 3, 4, 15, 16, 17, 31, 32, 33, 63, 64, 65 };
    var header1 = [_]u8{'x'} ** 65;
    var header2 = [_]u8{'x'} ** 65;

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

        for (1..length - 2) |separator| {
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
