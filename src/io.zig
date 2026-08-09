//! Byte interfaces, limits, and plain or gzip adapters for streaming FASTQ I/O.
//!
//! Adapters and borrowed backing state must stay alive and at stable addresses while wrapped.

const std = @import("std");
const builtin = @import("builtin");
const flate = std.compress.flate;
const Inflate = @import("inflate.zig");

const ReadError = error{ReadFailed};
pub const WriteError = error{WriteFailed};

pub const DEFAULT_MAX_LINE_BYTES: usize = 16 * 1024 * 1024;
pub const DEFAULT_READER_BUFFER_BYTES: usize = 256 * 1024;
pub const COUNT_READ_BUFFER_BYTES: usize = DEFAULT_READER_BUFFER_BYTES;
const GZIP_OPTIONAL_HEADER_BYTES_MAX: usize = 64 * 1024;
const CRC32_POLYNOMIAL: u32 = 0xedb8_8320;
const CRC32_FOLD_MIN_BYTES = 128;

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

const GzipCrc32 = struct {
    value: u32 = 0,
    bulk_update: ?UpdateFn,

    const SLICES = 8;
    const TABLE = makeCrc32Table(SLICES);
    const UpdateFn = *const fn (start: u32, bytes: []const u8) u32;

    fn init() GzipCrc32 {
        return initWithPclmul(runtimeHasPclmul());
    }

    fn initWithPclmul(has_pclmul: bool) GzipCrc32 {
        if (comptime builtin.cpu.arch == .x86_64) {
            return .{ .bulk_update = if (has_pclmul) X86Crc32.update else null };
        }
        return .{ .bulk_update = null };
    }

    fn fresh(self: GzipCrc32) GzipCrc32 {
        return .{ .bulk_update = self.bulk_update };
    }

    fn reset(self: *GzipCrc32) void {
        self.value = 0;
    }

    fn update(self: *GzipCrc32, bytes: []const u8) void {
        if (bytes.len >= CRC32_FOLD_MIN_BYTES) {
            if (self.bulk_update) |bulk_update| {
                self.value = bulk_update(self.value, bytes);
                return;
            }
        }
        self.value = updatePortable(self.value, bytes);
    }

    fn updatePortable(start: u32, bytes: []const u8) u32 {
        var crc = ~start;
        var offset: usize = 0;
        while (bytes.len - offset >= SLICES) : (offset += SLICES) {
            const first = std.mem.readInt(u32, bytes[offset..][0..4], .little) ^ crc;
            crc = TABLE[7][@as(u8, @truncate(first))] ^
                TABLE[6][@as(u8, @truncate(first >> 8))] ^
                TABLE[5][@as(u8, @truncate(first >> 16))] ^
                TABLE[4][@as(u8, @truncate(first >> 24))] ^
                TABLE[3][bytes[offset + 4]] ^
                TABLE[2][bytes[offset + 5]] ^
                TABLE[1][bytes[offset + 6]] ^
                TABLE[0][bytes[offset + 7]];
        }

        while (offset < bytes.len) : (offset += 1) {
            crc = TABLE[0][@as(u8, @truncate(crc ^ bytes[offset]))] ^ (crc >> 8);
        }
        return ~crc;
    }

    fn final(self: GzipCrc32) u32 {
        return self.value;
    }
};

fn runtimeHasPclmul() bool {
    if (comptime builtin.cpu.arch != .x86_64) return false;
    return X86Crc32.isSupported();
}

// The folding schedule follows crc32fast 1.5.0; see THIRD_PARTY_NOTICES.md.
const X86Crc32 = struct {
    const Block = @Vector(2, u64);
    const Selector = enum {
        low_low,
        high_high,
        low_high,
    };

    fn isSupported() bool {
        const ecx = asm volatile (
            \\movl $1, %%eax
            \\xorl %%ecx, %%ecx
            \\cpuid
            : [ecx] "={ecx}" (-> u32),
            :
            : .{ .rax = true, .rbx = true, .rdx = true });
        return ecx & (1 << 1) != 0;
    }

    fn update(start: u32, bytes: []const u8) u32 {
        if (bytes.len < CRC32_FOLD_MIN_BYTES) {
            return GzipCrc32.updatePortable(start, bytes);
        }

        var offset: usize = 0;
        var x3 = readBlock(bytes, &offset);
        var x2 = readBlock(bytes, &offset);
        var x1 = readBlock(bytes, &offset);
        var x0 = readBlock(bytes, &offset);
        x3 ^= Block{ @as(u64, ~start), 0 };

        const k1k2 = Block{ 0x1_5444_2bd4, 0x1_c6e4_1596 };
        while (bytes.len - offset >= 64) {
            x3 = reduce128(x3, readBlock(bytes, &offset), k1k2);
            x2 = reduce128(x2, readBlock(bytes, &offset), k1k2);
            x1 = reduce128(x1, readBlock(bytes, &offset), k1k2);
            x0 = reduce128(x0, readBlock(bytes, &offset), k1k2);
        }

        const k3k4 = Block{ 0x1_7519_97d0, 0x0_ccaa_009e };
        var x = reduce128(x3, x2, k3k4);
        x = reduce128(x, x1, k3k4);
        x = reduce128(x, x0, k3k4);
        while (bytes.len - offset >= 16) {
            x = reduce128(x, readBlock(bytes, &offset), k3k4);
        }

        x = clmul(x, k3k4, .low_high) ^ shiftRightBytes(x, 8);
        x = clmul(
            x & Block{ 0xffff_ffff, 0 },
            Block{ 0x1_63cd_6124, 0 },
            .low_low,
        ) ^ shiftRightBytes(x, 4);

        const reduction = Block{ 0x1_db71_0641, 0x1_f701_1641 };
        const t1 = clmul(x & Block{ 0xffff_ffff, 0 }, reduction, .low_high);
        const t2 = clmul(t1 & Block{ 0xffff_ffff, 0 }, reduction, .low_low);
        const words: @Vector(4, u32) = @bitCast(x ^ t2);
        return GzipCrc32.updatePortable(~words[1], bytes[offset..]);
    }

    fn readBlock(bytes: []const u8, offset: *usize) Block {
        const block: Block = @bitCast(bytes[offset.*..][0..16].*);
        offset.* += 16;
        return block;
    }

    fn reduce128(value: Block, next: Block, keys: Block) Block {
        return next ^ clmul(value, keys, .low_low) ^ clmul(value, keys, .high_high);
    }

    fn clmul(value: Block, key: Block, comptime selector: Selector) Block {
        return switch (selector) {
            .low_low => asm (
                \\pclmulqdq $0x00, %[key], %[value]
                : [value] "=x" (-> Block),
                : [input] "0" (value),
                  [key] "x" (key),
            ),
            .high_high => asm (
                \\pclmulqdq $0x11, %[key], %[value]
                : [value] "=x" (-> Block),
                : [input] "0" (value),
                  [key] "x" (key),
            ),
            .low_high => asm (
                \\pclmulqdq $0x10, %[key], %[value]
                : [value] "=x" (-> Block),
                : [input] "0" (value),
                  [key] "x" (key),
            ),
        };
    }

    fn shiftRightBytes(value: Block, comptime count: u5) Block {
        return switch (count) {
            4 => @bitCast(@shuffle(
                u32,
                @as(@Vector(4, u32), @bitCast(value)),
                @as(@Vector(4, u32), @splat(0)),
                @as(@Vector(4, i32), .{ 1, 2, 3, -1 }),
            )),
            8 => @shuffle(
                u64,
                value,
                @as(Block, @splat(0)),
                @as(@Vector(2, i32), .{ 1, -1 }),
            ),
            else => @compileError("unsupported byte shift"),
        };
    }
};

fn makeCrc32Table(comptime slices: usize) [slices][256]u32 {
    @setEvalBranchQuota(100_000);

    var table: [slices][256]u32 = undefined;
    for (0..256) |index| {
        var value: u32 = @intCast(index);
        for (0..8) |_| {
            value = (value >> 1) ^ (CRC32_POLYNOMIAL & (0 -% (value & 1)));
        }
        table[0][index] = value;
    }
    for (1..slices) |slice| {
        for (0..256) |index| {
            const previous = table[slice - 1][index];
            table[slice][index] = table[0][@as(u8, @truncate(previous))] ^ (previous >> 8);
        }
    }
    return table;
}

/// Streams and validates complete RFC 1952 member sequences from a borrowed reader.
/// The reader needs at least ten buffer bytes and must share this adapter's stable lifetime.
pub const GzipSource = struct {
    input: *std.Io.Reader,
    decompressor: Inflate = undefined,
    decompressor_buffer: [flate.max_window_len]u8 = undefined,
    crc: GzipCrc32,
    size: u32 = 0,
    state: State = .between_members,
    member_seen: bool = false,

    const State = enum {
        between_members,
        payload,
        eof,
    };

    pub fn init(input: *std.Io.Reader) GzipSource {
        return .{ .input = input, .crc = .init() };
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
        return self.beginMemberWithBuffer(&self.decompressor_buffer);
    }

    fn beginMemberWithBuffer(
        self: *GzipSource,
        decompressor_buffer: []u8,
    ) std.Io.Reader.Error!void {
        _ = self.input.peekByte() catch |err| return err;
        self.parseHeader() catch return error.ReadFailed;

        self.decompressor = .init(self.input, .raw, decompressor_buffer);
        self.crc.reset();
        self.size = 0;
        self.state = .payload;
        self.member_seen = true;
    }

    fn parseHeader(self: *GzipSource) std.Io.Reader.Error!void {
        const fixed = try self.input.takeArray(10);
        var header_crc = self.crc.fresh();
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
        crc: *GzipCrc32,
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
        crc: *GzipCrc32,
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
                const reader = &self.decompressor.reader;
                const decoded = reader.peekGreedy(1) catch |err| switch (err) {
                    error.EndOfStream => {
                        self.finishMember() catch return error.ReadFailed;
                        self.state = .between_members;
                        continue;
                    },
                    error.ReadFailed => return error.ReadFailed,
                };
                self.crc.update(decoded);
                self.size +%= @truncate(decoded.len);
                reader.toss(decoded.len);
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

test "[unit] - [gzip CRC-32]: matches standard vectors" {
    var crc: GzipCrc32 = .init();
    try std.testing.expectEqual(@as(u32, 0), crc.final());

    crc.update("123456789");
    try std.testing.expectEqual(@as(u32, 0xcbf4_3926), crc.final());
}

test "[property] - [gzip CRC-32]: matches standard CRC across alignments and splits" {
    var storage: [521]u8 = undefined;
    for (&storage, 0..) |*byte, index| byte.* = @truncate(index *% 73 +% 19);

    for (0..16) |alignment| {
        const input = storage[alignment..];
        var expected: std.hash.Crc32 = .init();
        expected.update(input);

        var whole: GzipCrc32 = .init();
        whole.update(input);
        try std.testing.expectEqual(expected.final(), whole.final());

        var bytewise: GzipCrc32 = .init();
        for (input) |byte| bytewise.update(&.{byte});
        try std.testing.expectEqual(expected.final(), bytewise.final());

        for (0..input.len + 1) |split| {
            var split_expected: std.hash.Crc32 = .init();
            split_expected.update(input[0..split]);
            split_expected.update(input[split..]);

            var actual: GzipCrc32 = .init();
            actual.update(input[0..split]);
            actual.update(input[split..]);
            try std.testing.expectEqual(split_expected.final(), actual.final());
        }
    }
}

test "[unit] - [gzip CRC-32 dispatch]: matches CPUID and preserves forced fallback" {
    var bytes: [257]u8 = undefined;
    for (&bytes, 0..) |*byte, index| byte.* = @truncate(index *% 41 +% 7);

    var expected: std.hash.Crc32 = .init();
    expected.update(&bytes);
    var fallback = GzipCrc32.initWithPclmul(false);
    fallback.update(&bytes);
    try std.testing.expect(fallback.bulk_update == null);
    try std.testing.expectEqual(expected.final(), fallback.final());

    const selected = GzipCrc32.init();
    if (comptime builtin.cpu.arch == .x86_64) {
        try std.testing.expectEqual(X86Crc32.isSupported(), selected.bulk_update != null);
    } else {
        try std.testing.expect(selected.bulk_update == null);
    }
}

test "[property] - [gzip CRC-32]: PCLMUL matches portable folding boundaries" {
    if (comptime builtin.cpu.arch != .x86_64) return error.SkipZigTest;
    if (!X86Crc32.isSupported()) return error.SkipZigTest;

    var storage: [545]u8 = undefined;
    for (&storage, 0..) |*byte, index| byte.* = @truncate(index *% 73 +% 19);

    for (0..16) |alignment| {
        for (0..storage.len - alignment + 1) |length| {
            const input = storage[alignment..][0..length];
            const expected = GzipCrc32.updatePortable(0x1234_5678, input);
            try std.testing.expectEqual(expected, X86Crc32.update(0x1234_5678, input));
        }

        const input = storage[alignment..];
        for (0..input.len + 1) |split| {
            const expected = GzipCrc32.updatePortable(0x89ab_cdef, input);
            var actual = X86Crc32.update(0x89ab_cdef, input[0..split]);
            actual = X86Crc32.update(actual, input[split..]);
            try std.testing.expectEqual(expected, actual);
        }
    }
}
