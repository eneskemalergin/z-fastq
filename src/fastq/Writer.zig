//! Streaming FASTQ writer with LF line endings.

const std = @import("std");
const Record = @import("Record.zig").Record;
const WriteError = @import("../io/ByteSink.zig").WriteError;
const ByteSink = @import("../io/ByteSink.zig").ByteSink;

pub const WriterError = WriteError || error{InvalidRecord};

pub const Writer = struct {
    sink: ByteSink,

    /// The sink wrapper is copied; its referenced adapter must outlive the writer.
    pub fn init(sink: ByteSink) Writer {
        return .{ .sink = sink };
    }

    pub fn writeRecord(self: *Writer, record: Record) WriterError!void {
        if (record.sequence.len != record.quality.len or
            !isWritableField(record.header) or
            !isWritableField(record.sequence) or
            !isWritableField(record.plus) or
            !isWritableField(record.quality))
        {
            return error.InvalidRecord;
        }

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
        return self.sink.flush();
    }
};

fn isWritableField(bytes: []const u8) bool {
    return std.mem.indexOfScalar(u8, bytes, '\n') == null and
        (bytes.len == 0 or bytes[bytes.len - 1] != '\r');
}
