//! Plain file and in-memory adapters for `ByteSource` and `ByteSink`.

const std = @import("std");
const ByteSource = @import("ByteSource.zig").ByteSource;
const ReadError = @import("ByteSource.zig").ReadError;
const ByteSink = @import("ByteSink.zig").ByteSink;
const WriteError = @import("ByteSink.zig").WriteError;

pub const SliceSource = struct {
    data: []const u8,
    pos: usize = 0,

    pub fn init(data: []const u8) SliceSource {
        return .{ .data = data };
    }

    pub fn byteSource(self: *SliceSource) ByteSource {
        return .{
            .vtable = &slice_vtable,
            .ctx = self,
        };
    }
};

const slice_vtable = ByteSource.VTable{
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
            .vtable = &reader_vtable,
            .ctx = self,
        };
    }
};

const reader_vtable = ByteSource.VTable{
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
        var self: FileSource = undefined;
        self.file_reader = file.reader(io, read_buf);
        return self;
    }

    pub fn byteSource(self: *FileSource) ByteSource {
        return .{
            .vtable = &file_source_vtable,
            .ctx = self,
        };
    }
};

const file_source_vtable = ByteSource.VTable{
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
            .vtable = &slice_sink_vtable,
            .ctx = self,
        };
    }
};

const slice_sink_vtable = ByteSink.VTable{
    .write = sliceWrite,
    .flush = sliceFlush,
};

fn sliceWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *SliceSink = @ptrCast(@alignCast(ctx));
    if (self.pos + data.len > self.buffer.len) return error.WriteFailed;
    @memcpy(self.buffer[self.pos..][0..data.len], data);
    self.pos += data.len;
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
            .vtable = &writer_sink_vtable,
            .ctx = self,
        };
    }
};

const writer_sink_vtable = ByteSink.VTable{
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
    adapter: WriterSink,

    pub fn init(
        io: std.Io,
        file: std.Io.File,
        write_buf: []u8,
    ) FileSink {
        var self: FileSink = undefined;
        self.file_writer = file.writer(io, write_buf);
        self.adapter = WriterSink.init(&self.file_writer.interface);
        return self;
    }

    pub fn byteSink(self: *FileSink) ByteSink {
        return self.adapter.byteSink();
    }

    pub fn flush(self: *FileSink, io: std.Io) WriteError!void {
        _ = io;
        return self.adapter.byteSink().flush();
    }
};
