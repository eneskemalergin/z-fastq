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
    var first_end: usize = 0;
    while (first_end < header.len and header[first_end] != ' ' and header[first_end] != '\t') {
        first_end += 1;
    }
    const first_token = header[0..first_end];
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

    var second_start = first_end;
    while (second_start < header.len and
        (header[second_start] == ' ' or header[second_start] == '\t'))
    {
        second_start += 1;
    }
    var second_end = second_start;
    while (second_end < header.len and header[second_end] != ' ' and
        header[second_end] != '\t')
    {
        second_end += 1;
    }
    const second_token = header[second_start..second_end];
    if (second_token.len >= 2 and second_token[1] == ':' and
        (second_token[0] == '1' or second_token[0] == '2'))
    {
        mate_markers |= mateMask(@intCast(second_token[0] - '0'));
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

pub fn namesMatch(name1: Name, name2: Name) bool {
    if (!std.mem.eql(u8, name1.normalized_id, name2.normalized_id)) return false;
    if (name1.mate_markers == 0b11 or name2.mate_markers == 0b11) return false;
    if ((name1.mate_markers == 0) != (name2.mate_markers == 0)) return false;
    if (name1.mate_markers == 0) return true;
    return name1.mate_markers == 0b01 and name2.mate_markers == 0b10;
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
