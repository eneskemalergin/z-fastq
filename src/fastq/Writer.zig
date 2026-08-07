//! Streaming FASTQ writer with LF line endings.

const std = @import("std");
const Record = @import("Record.zig").Record;
const WriteError = @import("../io/ByteSink.zig").WriteError;
const ByteSink = @import("../io/ByteSink.zig").ByteSink;

pub const Writer = struct {
    sink: *const ByteSink,

    pub fn init(sink: *const ByteSink) Writer {
        return .{ .sink = sink };
    }

    pub fn writeRecord(self: *Writer, record: Record) WriteError!void {
        try self.sink.write("@");
        try self.sink.write(record.header);
        try self.sink.write("\n");
        try self.sink.write(record.sequence);
        try self.sink.write("\n+");
        try self.sink.write(record.plus);
        try self.sink.write("\n");
        try self.sink.write(record.quality);
        try self.sink.write("\n");
    }

    pub fn flush(self: *Writer) WriteError!void {
        try self.sink.flush();
    }
};
