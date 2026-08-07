//! Public library surface for z-fastq.

const std = @import("std");

pub const Record = @import("fastq/Record.zig").Record;
pub const OwnedRecord = @import("fastq/Record.zig").OwnedRecord;
pub const Reader = @import("fastq/Reader.zig").Reader;
pub const Writer = @import("fastq/Writer.zig").Writer;
pub const ParseError = @import("fastq/ParseError.zig").ParseError;
pub const LintCode = @import("fastq/ParseError.zig").LintCode;
pub const ReaderError = @import("fastq/ParseError.zig").ReaderError;
pub const codeTag = @import("fastq/ParseError.zig").codeTag;

pub const io = @import("io/mod.zig");
pub const limits = @import("io/limits.zig");
pub const count_scan = @import("fastq/count_scan.zig");

pub const VERSION = "0.0.1";
