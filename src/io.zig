//! Byte interfaces, limits, and plain adapters for streaming FASTQ I/O.
//!
//! Adapters and borrowed backing state must stay alive and at stable addresses while wrapped.

const std = @import("std");

const ReadError = error{ReadFailed};
pub const WriteError = error{WriteFailed};

pub const DEFAULT_MAX_LINE_BYTES: usize = 16 * 1024 * 1024;
pub const DEFAULT_READER_BUFFER_BYTES: usize = 256 * 1024;
pub const COUNT_READ_BUFFER_BYTES: usize = DEFAULT_READER_BUFFER_BYTES;

/// Copied pull interface whose adapter must remain at a stable address and outlive it.
/// A read initializes the returned prefix of the destination and returns zero only at EOF.
pub const ByteSource = struct {
    vtable: *const VTable,
    ctx: *anyopaque,

    pub const VTable = struct {
        read: *const fn (ctx: *anyopaque, dest: []u8) ReadError!usize,
    };

    pub fn read(self: *const ByteSource, dest: []u8) ReadError!usize {
        return self.vtable.read(self.ctx, dest);
    }
};

/// Copied push interface whose adapter must remain at a stable address and outlive it.
/// Writes consume the complete slice; a missing flush callback is a successful no-op.
pub const ByteSink = struct {
    vtable: *const VTable,
    ctx: *anyopaque,

    pub const VTable = struct {
        write: *const fn (ctx: *anyopaque, data: []const u8) WriteError!void,
        flush: ?*const fn (ctx: *anyopaque) WriteError!void = null,
    };

    pub fn write(self: *const ByteSink, data: []const u8) WriteError!void {
        return self.vtable.write(self.ctx, data);
    }

    pub fn flush(self: *const ByteSink) WriteError!void {
        if (self.vtable.flush) |flush_fn| return flush_fn(self.ctx);
    }
};

pub const SliceSource = struct {
    data: []const u8,
    pos: usize = 0,

    pub fn init(data: []const u8) SliceSource {
        return .{ .data = data };
    }

    pub fn byteSource(self: *SliceSource) ByteSource {
        return .{
            .vtable = &SLICE_VTABLE,
            .ctx = self,
        };
    }
};

const SLICE_VTABLE = ByteSource.VTable{
    .read = sliceRead,
};

fn sliceRead(ctx: *anyopaque, dest: []u8) ReadError!usize {
    const self: *SliceSource = @ptrCast(@alignCast(ctx));
    const remaining = self.data[self.pos..];
    if (remaining.len == 0) return 0;
    const copy_len = @min(dest.len, remaining.len);
    @memcpy(dest[0..copy_len], remaining[0..copy_len]);
    self.pos += copy_len;
    return copy_len;
}

pub const ReaderSource = struct {
    reader: *std.Io.Reader,

    pub fn init(reader: *std.Io.Reader) ReaderSource {
        return .{ .reader = reader };
    }

    pub fn byteSource(self: *ReaderSource) ByteSource {
        return .{
            .vtable = &READER_VTABLE,
            .ctx = self,
        };
    }
};

const READER_VTABLE = ByteSource.VTable{
    .read = readerRead,
};

fn readerRead(ctx: *anyopaque, dest: []u8) ReadError!usize {
    const self: *ReaderSource = @ptrCast(@alignCast(ctx));
    return self.reader.readSliceShort(dest) catch error.ReadFailed;
}

pub const FileSource = struct {
    file_reader: std.Io.File.Reader,

    pub fn init(
        io: std.Io,
        file: std.Io.File,
        read_buf: []u8,
    ) FileSource {
        return .{ .file_reader = file.reader(io, read_buf) };
    }

    pub fn byteSource(self: *FileSource) ByteSource {
        return .{
            .vtable = &FILE_SOURCE_VTABLE,
            .ctx = self,
        };
    }
};

const FILE_SOURCE_VTABLE = ByteSource.VTable{
    .read = fileSourceRead,
};

fn fileSourceRead(ctx: *anyopaque, dest: []u8) ReadError!usize {
    const self: *FileSource = @ptrCast(@alignCast(ctx));
    return self.file_reader.interface.readSliceShort(dest) catch error.ReadFailed;
}

pub const SliceSink = struct {
    buffer: []u8,
    pos: usize = 0,

    pub fn init(buffer: []u8) SliceSink {
        return .{ .buffer = buffer };
    }

    pub fn written(self: *const SliceSink) []const u8 {
        return self.buffer[0..self.pos];
    }

    pub fn byteSink(self: *SliceSink) ByteSink {
        return .{
            .vtable = &SLICE_SINK_VTABLE,
            .ctx = self,
        };
    }
};

const SLICE_SINK_VTABLE = ByteSink.VTable{
    .write = sliceWrite,
    .flush = sliceFlush,
};

fn sliceWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *SliceSink = @ptrCast(@alignCast(ctx));
    const end = std.math.add(usize, self.pos, data.len) catch return error.WriteFailed;
    if (end > self.buffer.len) return error.WriteFailed;
    @memcpy(self.buffer[self.pos..end], data);
    self.pos = end;
}

fn sliceFlush(ctx: *anyopaque) WriteError!void {
    _ = ctx;
}

pub const WriterSink = struct {
    writer: *std.Io.Writer,

    pub fn init(writer: *std.Io.Writer) WriterSink {
        return .{ .writer = writer };
    }

    pub fn byteSink(self: *WriterSink) ByteSink {
        return .{
            .vtable = &WRITER_SINK_VTABLE,
            .ctx = self,
        };
    }
};

const WRITER_SINK_VTABLE = ByteSink.VTable{
    .write = writerWrite,
    .flush = writerFlush,
};

fn writerWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *WriterSink = @ptrCast(@alignCast(ctx));
    self.writer.writeAll(data) catch return error.WriteFailed;
}

fn writerFlush(ctx: *anyopaque) WriteError!void {
    const self: *WriterSink = @ptrCast(@alignCast(ctx));
    self.writer.flush() catch return error.WriteFailed;
}

pub const FileSink = struct {
    file_writer: std.Io.File.Writer,

    pub fn init(
        io: std.Io,
        file: std.Io.File,
        write_buf: []u8,
    ) FileSink {
        return .{ .file_writer = file.writer(io, write_buf) };
    }

    pub fn byteSink(self: *FileSink) ByteSink {
        return .{
            .vtable = &FILE_SINK_VTABLE,
            .ctx = self,
        };
    }
};

const FILE_SINK_VTABLE = ByteSink.VTable{
    .write = fileSinkWrite,
    .flush = fileSinkFlush,
};

fn fileSinkWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *FileSink = @ptrCast(@alignCast(ctx));
    self.file_writer.interface.writeAll(data) catch return error.WriteFailed;
}

fn fileSinkFlush(ctx: *anyopaque) WriteError!void {
    const self: *FileSink = @ptrCast(@alignCast(ctx));
    self.file_writer.interface.flush() catch return error.WriteFailed;
}
