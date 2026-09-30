//! Byte interfaces, limits, and plain or gzip adapters for streaming FASTQ I/O.
//!
//! Adapters, the gzip decoder workspace, and borrowed backing state must stay alive and at stable
//! addresses while wrapped.

const std = @import("std");
const zipir = @import("zipir");

const ReadError = error{ReadFailed};
pub const WriteError = error{WriteFailed};

pub const DEFAULT_MAX_LINE_BYTES: usize = 16 * 1024 * 1024;
pub const DEFAULT_READER_BUFFER_BYTES: usize = 256 * 1024;
pub const COUNT_READ_BUFFER_BYTES: usize = DEFAULT_READER_BUFFER_BYTES;
const GZIP_OPTIONAL_HEADER_BYTES_MAX: usize = 64 * 1024;

/// Copied pull interface whose adapter must remain at a stable address and outlive it.
/// A read initializes the returned prefix and rejects a count beyond the destination.
/// For a nonempty destination, zero reports EOF.
pub const ByteSource = struct {
    vtable: *const VTable,
    ctx: *anyopaque,

    pub const VTable = struct {
        read: *const fn (ctx: *anyopaque, dest: []u8) ReadError!usize,
    };

    pub fn read(self: *const ByteSource, dest: []u8) ReadError!usize {
        const count = try self.vtable.read(self.ctx, dest);
        if (count > dest.len) return error.ReadFailed;
        return count;
    }
};

/// Copied push interface whose adapter must remain at a stable address and outlive it.
/// Writes consume the complete slice; a missing flush callback is a successful no-op.
pub const ByteSink = struct {
    vtable: *const VTable,
    ctx: *anyopaque,

    pub const VTable = struct {
        write: *const fn (ctx: *anyopaque, data: []const u8) WriteError!void,
        writeVec: ?*const fn (ctx: *anyopaque, data: []const []const u8) WriteError!void = null,
        flush: ?*const fn (ctx: *anyopaque) WriteError!void = null,
    };

    pub fn write(self: *const ByteSink, data: []const u8) WriteError!void {
        return self.vtable.write(self.ctx, data);
    }

    pub fn writeVec(self: *const ByteSink, data: []const []const u8) WriteError!void {
        if (self.vtable.writeVec) |write_vec| return write_vec(self.ctx, data);
        for (data) |bytes| try self.write(bytes);
    }

    pub fn flush(self: *const ByteSink) WriteError!void {
        if (self.vtable.flush) |flush_fn| return flush_fn(self.ctx);
    }
};

/// Pull adapter over borrowed in-memory bytes.
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

/// Pull adapter over a borrowed standard reader.
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
    return readShort(self.reader, dest);
}

fn readShort(reader: *std.Io.Reader, dest: []u8) ReadError!usize {
    return reader.readSliceShort(dest) catch error.ReadFailed;
}

pub fn readChunk(reader: *std.Io.Reader) ReadError!?[]u8 {
    const chunk = reader.peekGreedy(1) catch |err| switch (err) {
        error.EndOfStream => return null,
        error.ReadFailed => return error.ReadFailed,
    };
    reader.toss(chunk.len);
    return chunk;
}

pub const ByteSourceReader = struct {
    source: ByteSource,
    interface: std.Io.Reader,

    pub fn init(source: ByteSource, buffer: []u8) ByteSourceReader {
        return .{
            .source = source,
            .interface = .{
                .vtable = &.{ .stream = byteSourceStream, .readVec = byteSourceReadVec },
                .buffer = buffer,
                .seek = 0,
                .end = 0,
            },
        };
    }
};

fn byteSourceStream(
    r: *std.Io.Reader,
    w: *std.Io.Writer,
    limit: std.Io.Limit,
) std.Io.Reader.StreamError!usize {
    const self: *ByteSourceReader = @alignCast(@fieldParentPtr("interface", r));
    const dest = limit.slice(try w.writableSliceGreedy(1));
    const n = self.source.read(dest) catch return error.ReadFailed;
    if (n == 0) return error.EndOfStream;
    w.advance(n);
    return n;
}

fn byteSourceReadVec(r: *std.Io.Reader, data: [][]u8) std.Io.Reader.Error!usize {
    _ = data;
    const self: *ByteSourceReader = @alignCast(@fieldParentPtr("interface", r));
    if (r.seek == r.end) {
        r.seek = 0;
        r.end = 0;
    }
    std.debug.assert(r.end < r.buffer.len);
    const n = self.source.read(r.buffer[r.end..]) catch return error.ReadFailed;
    if (n == 0) return error.EndOfStream;
    r.end += n;
    return 0;
}

/// Buffered pull adapter over a borrowed file handle and caller-owned buffer.
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
    return readShort(&self.file_reader.interface, dest);
}

// --- gzip input ---

/// Streams and validates complete RFC 1952 member sequences from a borrowed reader.
/// The reader accepts any buffer size and must share this adapter's stable lifetime.
pub const GzipSource = struct {
    input: *std.Io.Reader,
    decoder: zipir.gzip.Decompressor = undefined,
    started: bool = false,

    pub fn init(input: *std.Io.Reader) GzipSource {
        return .{ .input = input };
    }

    pub fn byteSource(self: *GzipSource) ByteSource {
        return .{
            .vtable = &GZIP_VTABLE,
            .ctx = self,
        };
    }

    fn decoded(self: *GzipSource) *std.Io.Reader {
        if (!self.started) {
            self.decoder.init(self.input, .{ .max_header_bytes = GZIP_OPTIONAL_HEADER_BYTES_MAX });
            self.started = true;
        }
        return &self.decoder.reader;
    }
};

pub fn gzipReader(self: *GzipSource) *std.Io.Reader {
    return self.decoded();
}

const GZIP_VTABLE = ByteSource.VTable{
    .read = gzipRead,
};

fn gzipRead(ctx: *anyopaque, dest: []u8) ReadError!usize {
    const self: *GzipSource = @ptrCast(@alignCast(ctx));
    return readShort(self.decoded(), dest);
}

/// Push adapter into a borrowed fixed-capacity byte slice.
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
};

fn sliceWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *SliceSink = @ptrCast(@alignCast(ctx));
    const end = std.math.add(usize, self.pos, data.len) catch return error.WriteFailed;
    if (end > self.buffer.len) return error.WriteFailed;
    @memcpy(self.buffer[self.pos..end], data);
    self.pos = end;
}

/// Push adapter over a borrowed standard writer.
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
    .writeVec = writerWriteVec,
    .flush = writerFlush,
};

fn writerWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *WriterSink = @ptrCast(@alignCast(ctx));
    self.writer.writeAll(data) catch return error.WriteFailed;
}

fn writerWriteVec(ctx: *anyopaque, data: []const []const u8) WriteError!void {
    const self: *WriterSink = @ptrCast(@alignCast(ctx));
    writeVecToWriter(self.writer, data) catch return error.WriteFailed;
}

fn writeVecToWriter(writer: *std.Io.Writer, data: []const []const u8) std.Io.Writer.Error!void {
    if (data.len == 0) return;

    var total: usize = 0;
    for (data) |bytes| {
        total = std.math.add(usize, total, bytes.len) catch {
            for (data) |fallback_bytes| try writer.writeAll(fallback_bytes);
            return;
        };
    }

    if (total <= writer.unusedCapacityLen()) {
        _ = try writer.writeVec(data);
        return;
    }

    for (data) |bytes| try writer.writeAll(bytes);
}

fn writerFlush(ctx: *anyopaque) WriteError!void {
    const self: *WriterSink = @ptrCast(@alignCast(ctx));
    self.writer.flush() catch return error.WriteFailed;
}

/// Buffered push adapter over a borrowed file handle and caller-owned buffer.
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
    .writeVec = fileSinkWriteVec,
    .flush = fileSinkFlush,
};

fn fileSinkWrite(ctx: *anyopaque, data: []const u8) WriteError!void {
    const self: *FileSink = @ptrCast(@alignCast(ctx));
    self.file_writer.interface.writeAll(data) catch return error.WriteFailed;
}

fn fileSinkWriteVec(ctx: *anyopaque, data: []const []const u8) WriteError!void {
    const self: *FileSink = @ptrCast(@alignCast(ctx));
    writeVecToWriter(&self.file_writer.interface, data) catch return error.WriteFailed;
}

fn fileSinkFlush(ctx: *anyopaque) WriteError!void {
    const self: *FileSink = @ptrCast(@alignCast(ctx));
    self.file_writer.interface.flush() catch return error.WriteFailed;
}

test "[property] - [gzip direct delivery]: preserves bytes across members" {
    const gzip_x = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x01, 0x01, 0x00, 0xfe, 0xff, 'x',
        0x83, 0x16, 0xdc, 0x8c, 0x01, 0x00, 0x00, 0x00,
    };
    const compressed = gzip_x ++ gzip_x;
    var input = std.Io.Reader.fixed(&compressed);
    var source = GzipSource.init(&input);
    var output: [2]u8 = undefined;
    var output_len: usize = 0;

    while (try readChunk(gzipReader(&source))) |decoded| {
        @memcpy(output[output_len..][0..decoded.len], decoded);
        output_len += decoded.len;
    }

    try std.testing.expectEqual(output.len, output_len);
    try std.testing.expectEqualStrings("xx", &output);
}

test "[property] - [gzip input]: compressed output may span many reads" {
    const compressed = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0xff,
        0x73, 0x74, 0x1c, 0x05, 0xa3, 0x60, 0x14, 0x8c, 0x54, 0x00,
        0x00, 0x1a, 0xfb, 0x37, 0xb7, 0x00, 0x04, 0x00, 0x00,
    };
    var input = std.Io.Reader.fixed(&compressed);
    var source = GzipSource.init(&input);
    const bytes = source.byteSource();
    var output: [1024]u8 = undefined;
    var written: usize = 0;
    while (written < output.len) {
        const end = @min(written + 17, output.len);
        const n = try bytes.read(output[written..end]);
        try std.testing.expect(n > 0);
        written += n;
    }

    try std.testing.expect(std.mem.allEqual(u8, &output, 'A'));
    try std.testing.expectEqual(@as(usize, 0), try bytes.read(output[0..17]));
}
