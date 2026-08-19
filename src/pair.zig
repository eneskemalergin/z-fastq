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

    const second_token = nextToken(header, first_token.len).bytes;
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
        .exact => namesMatch(parseName(header1, policy), parseName(header2, policy)),
        .illumina => illuminaHeadersMatch(header1, header2),
    };
}

pub fn namesMatch(name1: Name, name2: Name) bool {
    if (!std.mem.eql(u8, name1.normalized_id, name2.normalized_id)) return false;
    if (name1.mate_markers == 0b11 or name2.mate_markers == 0b11) return false;
    if ((name1.mate_markers == 0) != (name2.mate_markers == 0)) return false;
    if (name1.mate_markers == 0) return true;
    return name1.mate_markers == 0b01 and name2.mate_markers == 0b10;
}

const HeaderToken = struct {
    bytes: []const u8,
    ends_header: bool,
};

fn firstTokenEnd(header: []const u8) usize {
    var end: usize = 0;
    while (end < header.len and header[end] != ' ' and header[end] != '\t') {
        end += 1;
    }
    return end;
}

fn nextToken(header: []const u8, previous_end: usize) HeaderToken {
    var start = previous_end;
    while (start < header.len and (header[start] == ' ' or header[start] == '\t')) {
        start += 1;
    }
    var end = start;
    while (end < header.len and header[end] != ' ' and header[end] != '\t') {
        end += 1;
    }
    return .{ .bytes = header[start..end], .ends_header = end == header.len };
}

// Keep the full matcher shared across paired command paths to avoid code growth.
noinline fn illuminaHeadersMatch(header1: []const u8, header2: []const u8) bool {
    return terminalPairHeadersMatch(header1, header2) or
        namesMatch(parseName(header1, .illumina), parseName(header2, .illumina));
}

fn terminalPairHeadersMatch(header1: []const u8, header2: []const u8) bool {
    if (header1.len < 2 or header1.len != header2.len) return false;
    if (header1[header1.len - 2] != '/' or header2[header2.len - 2] != '/') return false;
    if (header1[header1.len - 1] != '1' or header2[header2.len - 1] != '2') return false;
    if (!std.mem.eql(u8, header1[0 .. header1.len - 1], header2[0 .. header2.len - 1])) {
        return false;
    }

    const first_token = header1[0..firstTokenEnd(header1)];
    if (first_token.len == header1.len) return true;
    if (terminalMateMarker(first_token) != null) return false;

    const second_token = nextToken(header1, first_token.len);
    return second_token.ends_header and leadingMateMarker(second_token.bytes) == null;
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
