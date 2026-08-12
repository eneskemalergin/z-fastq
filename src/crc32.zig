//! Private CRC-32/ISO-HDLC with portable and runtime-selected x86-64 paths.

const std = @import("std");
const builtin = @import("builtin");

const POLYNOMIAL: u32 = 0xedb8_8320;
const FOLD_MIN_BYTES = 128;
// Zig's self-hosted x86 backend only encodes PCLMUL when the target enables it.
const CAN_COMPILE_PCLMUL = builtin.cpu.arch == .x86_64 and
    (builtin.zig_backend != .stage2_x86_64 or builtin.cpu.has(.x86, .pclmul));

pub const Crc32 = struct {
    value: u32 = 0,
    bulk_update: ?UpdateFn = null,

    const SLICES = 8;
    const TABLE = makeTable(SLICES);
    const UpdateFn = *const fn (start: u32, bytes: []const u8) u32;

    pub fn init() Crc32 {
        return initWithPclmul(runtimeHasPclmul());
    }

    fn initPortable() Crc32 {
        return .{};
    }

    pub fn update(self: *Crc32, bytes: []const u8) void {
        if (bytes.len >= FOLD_MIN_BYTES) {
            if (self.bulk_update) |bulk_update| {
                self.value = bulk_update(self.value, bytes);
                return;
            }
        }
        self.value = updatePortable(self.value, bytes);
    }

    pub fn reset(self: *Crc32) void {
        self.value = 0;
    }

    pub fn final(self: Crc32) u32 {
        return self.value;
    }

    fn usesPclmul(self: Crc32) bool {
        return self.bulk_update != null;
    }

    fn initWithPclmul(has_pclmul: bool) Crc32 {
        if (comptime CAN_COMPILE_PCLMUL) {
            return .{ .bulk_update = if (has_pclmul) X86Crc32.update else null };
        }
        return .{};
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
};

fn runtimeHasPclmul() bool {
    if (comptime !CAN_COMPILE_PCLMUL) return false;
    return X86Crc32.isSupported();
}

// --- x86-64 PCLMUL ---

// The folding schedule and constants are adapted from crc32fast 1.5.0.
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
        if (bytes.len < FOLD_MIN_BYTES) return Crc32.updatePortable(start, bytes);

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
        return Crc32.updatePortable(~words[1], bytes[offset..]);
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

fn makeTable(comptime slices: usize) [slices][256]u32 {
    @setEvalBranchQuota(100_000);

    var table: [slices][256]u32 = undefined;
    for (0..256) |index| {
        var value: u32 = @intCast(index);
        for (0..8) |_| {
            value = (value >> 1) ^ (POLYNOMIAL & (0 -% (value & 1)));
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

test "[unit] - [CRC-32]: matches standard vectors" {
    var crc = Crc32.init();
    try std.testing.expectEqual(@as(u32, 0), crc.final());

    crc.update("123456789");
    try std.testing.expectEqual(@as(u32, 0xcbf4_3926), crc.final());
}

test "[property] - [CRC-32]: matches std across alignments and splits" {
    var storage: [521]u8 = undefined;
    for (&storage, 0..) |*byte, index| byte.* = @truncate(index *% 73 +% 19);

    for (0..16) |alignment| {
        const input = storage[alignment..];
        var expected: std.hash.Crc32 = .init();
        expected.update(input);

        var whole = Crc32.init();
        whole.update(input);
        try std.testing.expectEqual(expected.final(), whole.final());

        var bytewise = Crc32.init();
        for (input) |byte| bytewise.update(&.{byte});
        try std.testing.expectEqual(expected.final(), bytewise.final());

        for (0..input.len + 1) |split| {
            var actual = Crc32.init();
            actual.update(input[0..split]);
            actual.update(input[split..]);
            try std.testing.expectEqual(expected.final(), actual.final());
        }
    }
}

test "[unit] - [CRC-32 dispatch]: selects safely and permits portable use" {
    var bytes: [257]u8 = undefined;
    for (&bytes, 0..) |*byte, index| byte.* = @truncate(index *% 41 +% 7);

    var expected: std.hash.Crc32 = .init();
    expected.update(&bytes);
    var portable = Crc32.initPortable();
    portable.update(&bytes);
    try std.testing.expect(!portable.usesPclmul());
    try std.testing.expectEqual(expected.final(), portable.final());

    const selected = Crc32.init();
    if (comptime CAN_COMPILE_PCLMUL) {
        try std.testing.expectEqual(X86Crc32.isSupported(), selected.usesPclmul());
    } else {
        try std.testing.expect(!selected.usesPclmul());
    }
}

test "[unit] - [CRC-32 reset]: preserves the selected implementation" {
    var crc = Crc32.init();
    const used_pclmul = crc.usesPclmul();
    crc.update("123456789");
    crc.reset();

    try std.testing.expectEqual(@as(u32, 0), crc.final());
    try std.testing.expectEqual(used_pclmul, crc.usesPclmul());
}

test "[property] - [CRC-32 PCLMUL]: matches portable folding boundaries" {
    if (comptime !CAN_COMPILE_PCLMUL) return error.SkipZigTest;
    if (!X86Crc32.isSupported()) return error.SkipZigTest;

    var storage: [545]u8 = undefined;
    for (&storage, 0..) |*byte, index| byte.* = @truncate(index *% 73 +% 19);

    for (0..16) |alignment| {
        for (0..storage.len - alignment + 1) |length| {
            const input = storage[alignment..][0..length];
            const expected = Crc32.updatePortable(0x1234_5678, input);
            try std.testing.expectEqual(expected, X86Crc32.update(0x1234_5678, input));
        }

        const input = storage[alignment..];
        for (0..input.len + 1) |split| {
            const expected = Crc32.updatePortable(0x89ab_cdef, input);
            var actual = X86Crc32.update(0x89ab_cdef, input[0..split]);
            actual = X86Crc32.update(actual, input[split..]);
            try std.testing.expectEqual(expected, actual);
        }
    }
}
