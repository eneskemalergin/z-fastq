//! Byte interfaces, limits, and plain or gzip adapters for streaming FASTQ I/O.
//!
//! Adapters and borrowed backing state must stay alive and at stable addresses while wrapped.

const std = @import("std");
const build_options = @import("build_options");
const flate = std.compress.flate;
const crc32 = @import("crc32.zig");
const Inflate = @import("inflate.zig");
const use_isa_l = build_options.use_isa_l;
const isal = if (use_isa_l) @cImport({
    @cInclude("igzip_lib.h");
}) else struct {};
const PayloadCrc32 = if (use_isa_l) void else crc32.Crc32;

const ReadError = error{ReadFailed};
pub const WriteError = error{WriteFailed};

pub const DEFAULT_MAX_LINE_BYTES: usize = 16 * 1024 * 1024;
pub const DEFAULT_READER_BUFFER_BYTES: usize = 256 * 1024;
pub const COUNT_READ_BUFFER_BYTES: usize = DEFAULT_READER_BUFFER_BYTES;
pub const COUNT_DECOMPRESS_BUFFER_BYTES: usize = if (use_isa_l)
    COUNT_READ_BUFFER_BYTES
else
    flate.history_len + COUNT_READ_BUFFER_BYTES;
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
        flush: ?*const fn (ctx: *anyopaque) WriteError!void = null,
    };

    pub fn write(self: *const ByteSink, data: []const u8) WriteError!void {
        return self.vtable.write(self.ctx, data);
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
    return self.reader.readSliceShort(dest) catch error.ReadFailed;
}

pub const PlainFileSource = struct {
    file_reader: *std.Io.File.Reader,

    pub fn init(file_reader: *std.Io.File.Reader) PlainFileSource {
        return .{ .file_reader = file_reader };
    }

    pub fn byteSource(self: *PlainFileSource) ByteSource {
        return .{
            .vtable = &PLAIN_FILE_VTABLE,
            .ctx = self,
        };
    }
};

const PLAIN_FILE_VTABLE = ByteSource.VTable{
    .read = plainFileRead,
};

fn plainFileRead(ctx: *anyopaque, dest: []u8) ReadError!usize {
    const self: *PlainFileSource = @ptrCast(@alignCast(ctx));
    const reader = &self.file_reader.interface;
    var written: usize = 0;

    const buffered = reader.buffered();
    if (buffered.len != 0) {
        const copy_len = @min(buffered.len, dest.len);
        @memcpy(dest[0..copy_len], buffered[0..copy_len]);
        reader.toss(copy_len);
        written = copy_len;
    }
    if (reader.bufferedLen() == 0) {
        reader.buffer = &.{};
        reader.seek = 0;
        reader.end = 0;
    }
    if (written == dest.len) return written;

    const direct = reader.readSliceShort(dest[written..]) catch return error.ReadFailed;
    return written + direct;
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
    return self.file_reader.interface.readSliceShort(dest) catch error.ReadFailed;
}

// --- gzip input ---

const IsalInflate = struct {
    state: isal.struct_inflate_state = undefined,
    ended: bool = false,

    fn init(self: *IsalInflate) void {
        self.* = .{};
        isal.isal_inflate_init(&self.state);
        self.state.crc_flag = isal.ISAL_GZIP_NO_HDR_VER;
    }

    fn read(
        self: *IsalInflate,
        input: *std.Io.Reader,
        output: []u8,
    ) ReadError!?[]const u8 {
        if (self.ended) return null;
        if (output.len == 0) return error.ReadFailed;

        while (true) {
            const compressed = input.peekGreedy(1) catch return error.ReadFailed;
            const input_len: u32 = @intCast(@min(compressed.len, std.math.maxInt(u32)));
            const output_len: u32 = @intCast(@min(output.len, std.math.maxInt(u32)));
            self.state.next_in = @ptrCast(@constCast(compressed.ptr));
            self.state.avail_in = input_len;
            self.state.next_out = @ptrCast(output.ptr);
            self.state.avail_out = output_len;

            if (isal.isal_inflate(&self.state) != isal.ISAL_DECOMP_OK) {
                return error.ReadFailed;
            }
            const consumed = input_len - self.state.avail_in;
            const produced = output_len - self.state.avail_out;
            input.toss(consumed);

            if (self.state.block_state == isal.ISAL_BLOCK_FINISH) {
                self.ended = true;
                return if (produced == 0) null else output[0..produced];
            }
            if (consumed == 0 and produced == 0) return error.ReadFailed;
            if (produced != 0) return output[0..produced];
        }
    }
};

const GzipInflate = if (use_isa_l) IsalInflate else Inflate;

/// Streams and validates complete RFC 1952 member sequences from a borrowed reader.
/// The reader needs at least ten buffer bytes and must share this adapter's stable lifetime.
pub const GzipSource = struct {
    input: *std.Io.Reader,
    decompressor: GzipInflate = if (use_isa_l) .{} else undefined,
    decompressor_buffer: [if (use_isa_l) 0 else flate.max_window_len]u8 = undefined,
    payload_crc: PayloadCrc32,
    size: u32 = 0,
    state: State = .between_members,
    member_seen: bool = false,

    const State = enum {
        between_members,
        payload,
        eof,
    };

    pub fn init(input: *std.Io.Reader) GzipSource {
        return .{
            .input = input,
            .payload_crc = if (use_isa_l) {} else crc32.Crc32.init(),
        };
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
                    if (!use_isa_l) {
                        const n = self.decompressor.reader.readSliceShort(dest[written..]) catch
                            return error.ReadFailed;
                        const decoded = dest[written..][0..n];
                        self.payload_crc.update(decoded);
                        self.size +%= @truncate(n);
                        written += n;
                        if (written == dest.len) return written;

                        self.finishMember() catch return error.ReadFailed;
                        self.state = .between_members;
                        continue;
                    }

                    const decoded = (self.decompressor.read(
                        self.input,
                        dest[written..],
                    ) catch return error.ReadFailed) orelse {
                        self.finishMember() catch return error.ReadFailed;
                        self.state = .between_members;
                        continue;
                    };
                    self.size +%= @truncate(decoded.len);
                    written += decoded.len;
                    if (written == dest.len) return written;
                },
            }
        }
        return written;
    }

    fn beginMember(self: *GzipSource) std.Io.Reader.Error!void {
        return self.beginMemberWithBuffer(&self.decompressor_buffer);
    }

    fn beginMemberWithBuffer(
        self: *GzipSource,
        decompressor_buffer: []u8,
    ) std.Io.Reader.Error!void {
        _ = self.input.peekByte() catch |err| return err;
        self.parseHeader() catch return error.ReadFailed;

        if (use_isa_l) {
            self.decompressor.init();
        } else {
            self.decompressor = .init(self.input, .raw, decompressor_buffer);
        }
        if (!use_isa_l) self.payload_crc.reset();
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
        if (use_isa_l) return;
        const expected_crc = try self.input.takeInt(u32, .little);
        const expected_size = try self.input.takeInt(u32, .little);
        if (expected_crc != self.payload_crc.final() or expected_size != self.size) {
            return error.ReadFailed;
        }
    }
};

pub fn readGzipChunk(
    self: *GzipSource,
    decompressor_buffer: []u8,
) ReadError!?[]const u8 {
    while (true) {
        switch (self.state) {
            .eof => return null,
            .between_members => {
                self.beginMemberWithBuffer(decompressor_buffer) catch |err| switch (err) {
                    error.EndOfStream => {
                        if (!self.member_seen) return error.ReadFailed;
                        self.state = .eof;
                        return null;
                    },
                    error.ReadFailed => return error.ReadFailed,
                };
            },
            .payload => {
                if (!use_isa_l) {
                    const reader = &self.decompressor.reader;
                    const decoded = reader.peekGreedy(1) catch |err| switch (err) {
                        error.EndOfStream => {
                            self.finishMember() catch return error.ReadFailed;
                            self.state = .between_members;
                            continue;
                        },
                        error.ReadFailed => return error.ReadFailed,
                    };
                    self.payload_crc.update(decoded);
                    self.size +%= @truncate(decoded.len);
                    reader.toss(decoded.len);
                    return decoded;
                }

                const decoded = (self.decompressor.read(self.input, decompressor_buffer) catch
                    return error.ReadFailed) orelse {
                    self.finishMember() catch return error.ReadFailed;
                    self.state = .between_members;
                    continue;
                };
                self.size +%= @truncate(decoded.len);
                return decoded;
            },
        }
    }
}

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

test "[property] - [gzip direct delivery]: preserves bytes across members" {
    const gzip_x = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x01, 0x01, 0x00, 0xfe, 0xff, 'x',
        0x83, 0x16, 0xdc, 0x8c, 0x01, 0x00, 0x00, 0x00,
    };
    const compressed = gzip_x ++ gzip_x;
    var input = std.Io.Reader.fixed(&compressed);
    var source = GzipSource.init(&input);
    var decompressor_buffer: [flate.max_window_len]u8 = undefined;
    var output: [2]u8 = undefined;
    var output_len: usize = 0;

    while (try readGzipChunk(&source, &decompressor_buffer)) |decoded| {
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
    var output: [1024]u8 = undefined;
    var written: usize = 0;
    while (written < output.len) {
        const end = @min(written + 17, output.len);
        const n = try source.read(output[written..end]);
        try std.testing.expect(n > 0);
        written += n;
    }

    try std.testing.expect(std.mem.allEqual(u8, &output, 'A'));
    try std.testing.expectEqual(@as(usize, 0), try source.read(output[0..17]));
}
