//! Byte interfaces, limits, and plain or gzip adapters for streaming FASTQ I/O.
//!
//! Adapters and borrowed backing state must stay alive and at stable addresses while wrapped.

const std = @import("std");
const flate = std.compress.flate;

const ReadError = error{ReadFailed};
pub const WriteError = error{WriteFailed};

pub const DEFAULT_MAX_LINE_BYTES: usize = 16 * 1024 * 1024;
pub const DEFAULT_READER_BUFFER_BYTES: usize = 256 * 1024;
pub const COUNT_READ_BUFFER_BYTES: usize = DEFAULT_READER_BUFFER_BYTES;
const GZIP_OPTIONAL_HEADER_BYTES_MAX: usize = 64 * 1024;

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

// --- gzip input ---

/// Streams and validates complete RFC 1952 member sequences from a borrowed reader.
/// The reader needs at least ten buffer bytes and must share this adapter's stable lifetime.
pub const GzipSource = struct {
    input: *std.Io.Reader,
    decompressor: flate.Decompress = undefined,
    decompressor_buffer: [flate.max_window_len]u8 = undefined,
    crc: std.hash.Crc32 = .init(),
    size: u32 = 0,
    state: State = .between_members,
    member_seen: bool = false,

    const State = enum {
        between_members,
        payload,
        eof,
    };

    pub fn init(input: *std.Io.Reader) GzipSource {
        return .{ .input = input };
    }

    pub fn byteSource(self: *GzipSource) ByteSource {
        return .{
            .vtable = &GZIP_VTABLE,
            .ctx = self,
        };
    }

    fn read(self: *GzipSource, dest: []u8) ReadError!usize {
        var written: usize = 0;
        while (written < dest.len) {
            switch (self.state) {
                .eof => return written,
                .between_members => {
                    self.beginMember() catch |err| switch (err) {
                        error.EndOfStream => {
                            if (!self.member_seen) return error.ReadFailed;
                            self.state = .eof;
                            return written;
                        },
                        error.ReadFailed => return error.ReadFailed,
                    };
                },
                .payload => {
                    const n = self.decompressor.reader.readSliceShort(dest[written..]) catch
                        return error.ReadFailed;
                    const decoded = dest[written..][0..n];
                    self.crc.update(decoded);
                    self.size +%= @truncate(n);
                    written += n;
                    if (written == dest.len) return written;

                    self.finishMember() catch return error.ReadFailed;
                    self.state = .between_members;
                },
            }
        }
        return written;
    }

    fn beginMember(self: *GzipSource) std.Io.Reader.Error!void {
        _ = self.input.peekByte() catch |err| return err;
        self.parseHeader() catch return error.ReadFailed;

        self.decompressor = .init(self.input, .raw, &self.decompressor_buffer);
        self.crc = .init();
        self.size = 0;
        self.state = .payload;
        self.member_seen = true;
    }

    fn parseHeader(self: *GzipSource) std.Io.Reader.Error!void {
        const fixed = try self.input.takeArray(10);
        var header_crc: std.hash.Crc32 = .init();
        header_crc.update(fixed);
        if (fixed[0] != 0x1f or fixed[1] != 0x8b or fixed[2] != 8) {
            return error.ReadFailed;
        }

        const flags = fixed[3];
        if (flags & 0xe0 != 0) return error.ReadFailed;
        var optional_bytes: usize = 0;

        if (flags & 0x04 != 0) {
            try reserveOptionalBytes(&optional_bytes, 2);
            const length_bytes = try self.input.takeArray(2);
            header_crc.update(length_bytes);
            const length = std.mem.readInt(u16, length_bytes, .little);
            try self.consumeHeaderBytes(length, &optional_bytes, &header_crc);
        }
        if (flags & 0x08 != 0) {
            try self.consumeHeaderString(&optional_bytes, &header_crc);
        }
        if (flags & 0x10 != 0) {
            try self.consumeHeaderString(&optional_bytes, &header_crc);
        }
        if (flags & 0x02 != 0) {
            try reserveOptionalBytes(&optional_bytes, 2);
            const expected = try self.input.takeInt(u16, .little);
            if (expected != @as(u16, @truncate(header_crc.final()))) {
                return error.ReadFailed;
            }
        }
    }

    fn consumeHeaderBytes(
        self: *GzipSource,
        count: usize,
        optional_bytes: *usize,
        crc: *std.hash.Crc32,
    ) std.Io.Reader.Error!void {
        try reserveOptionalBytes(optional_bytes, count);
        for (0..count) |_| {
            const byte = try self.input.takeByte();
            crc.update(&.{byte});
        }
    }

    fn consumeHeaderString(
        self: *GzipSource,
        optional_bytes: *usize,
        crc: *std.hash.Crc32,
    ) std.Io.Reader.Error!void {
        while (true) {
            try reserveOptionalBytes(optional_bytes, 1);
            const byte = try self.input.takeByte();
            crc.update(&.{byte});
            if (byte == 0) return;
        }
    }

    fn finishMember(self: *GzipSource) std.Io.Reader.Error!void {
        const expected_crc = try self.input.takeInt(u32, .little);
        const expected_size = try self.input.takeInt(u32, .little);
        if (expected_crc != self.crc.final() or expected_size != self.size) {
            return error.ReadFailed;
        }
    }
};

const GZIP_VTABLE = ByteSource.VTable{
    .read = gzipRead,
};

fn gzipRead(ctx: *anyopaque, dest: []u8) ReadError!usize {
    const self: *GzipSource = @ptrCast(@alignCast(ctx));
    return self.read(dest);
}

fn reserveOptionalBytes(count: *usize, amount: usize) std.Io.Reader.Error!void {
    if (amount > GZIP_OPTIONAL_HEADER_BYTES_MAX - count.*) return error.ReadFailed;
    count.* += amount;
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
