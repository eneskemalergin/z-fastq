//! Stable format error codes and parse context.

pub const LintCode = enum {
    s001_invalid_plus_line,
    s002_invalid_sequence_alphabet,
    s003_invalid_header,
    s004_truncated_record,
    s005_length_mismatch,
    s006_invalid_quality_range,
    s007_duplicate_read_name,
    p001_paired_name_mismatch,
    p002_paired_count_mismatch,
    p003_paired_order_desync,
};

pub const ParseError = struct {
    code: LintCode,
    message: []const u8,
    record_index: u64,
    byte_offset: u64,
    line_in_record: u3,
};

pub fn codeTag(code: LintCode) []const u8 {
    return switch (code) {
        .s001_invalid_plus_line => "S001",
        .s002_invalid_sequence_alphabet => "S002",
        .s003_invalid_header => "S003",
        .s004_truncated_record => "S004",
        .s005_length_mismatch => "S005",
        .s006_invalid_quality_range => "S006",
        .s007_duplicate_read_name => "S007",
        .p001_paired_name_mismatch => "P001",
        .p002_paired_count_mismatch => "P002",
        .p003_paired_order_desync => "P003",
    };
}

pub const ReaderError = error{
    S003InvalidHeader,
    S004TruncatedRecord,
    S005LengthMismatch,
    LineTooLong,
    OutOfMemory,
    Io,
};
