//! Deterministic fraction parsing and record selection for the private sample command.

const std = @import("std");

pub const FractionError = error{InvalidFraction};
pub const SeedError = error{ InvalidSeed, Overflow };

pub const Fraction = union(enum) {
    none,
    all,
    probability: f64,

    pub fn parse(text: []const u8) FractionError!Fraction {
        if (std.mem.eql(u8, text, "0")) return .none;
        if (std.mem.eql(u8, text, "1")) return .all;
        if (text.len < 3 or text[1] != '.') return error.InvalidFraction;

        if (text[0] == '1') {
            for (text[2..]) |byte| {
                if (byte != '0') return error.InvalidFraction;
            }
            return .all;
        }
        if (text[0] != '0') return error.InvalidFraction;
        var nonzero = false;
        for (text[2..]) |byte| {
            if (byte < '0' or byte > '9') return error.InvalidFraction;
            nonzero = nonzero or byte != '0';
        }
        if (!nonzero) return .none;
        const probability = std.fmt.parseFloat(f64, text) catch
            return error.InvalidFraction;
        if (probability == 0.0) return .none;
        if (probability == 1.0) return .all;
        return .{ .probability = probability };
    }
};

pub fn parseSeed(text: []const u8) SeedError!u64 {
    if (text.len == 0) return error.InvalidSeed;
    for (text) |byte| {
        if (byte < '0' or byte > '9') return error.InvalidSeed;
    }
    return std.fmt.parseInt(u64, text, 10) catch |err| switch (err) {
        error.Overflow => error.Overflow,
        error.InvalidCharacter => error.InvalidSeed,
    };
}

pub const Selector = union(enum) {
    none,
    all,
    probability: Probability,

    const Probability = struct {
        fraction: f64,
        generator: Mt19937_64,
    };

    pub fn init(fraction: Fraction, seed: u64) Selector {
        return switch (fraction) {
            .none => .none,
            .all => .all,
            .probability => |probability| .{ .probability = .{
                .fraction = probability,
                .generator = Mt19937_64.init(seed),
            } },
        };
    }

    pub fn selectRecord(self: *Selector) bool {
        return switch (self.*) {
            .none => false,
            .all => true,
            .probability => |*state| state.generator.nextUnitFloat() < state.fraction,
        };
    }
};

pub const Mt19937_64 = struct {
    state: [STATE_WORDS]u64,
    index: usize,

    const STATE_WORDS = 312;
    const MIDDLE_WORD = 156;
    const MATRIX_A: u64 = 0xb5026f5aa96619e9;
    const UPPER_MASK: u64 = 0xffffffff80000000;
    const LOWER_MASK: u64 = 0x000000007fffffff;

    pub fn init(seed: u64) Mt19937_64 {
        var generator: Mt19937_64 = undefined;
        generator.state[0] = seed;
        for (1..STATE_WORDS) |index| {
            const previous = generator.state[index - 1];
            generator.state[index] = 6364136223846793005 *%
                (previous ^ (previous >> 62)) +% index;
        }
        generator.index = STATE_WORDS;
        return generator;
    }

    pub fn nextU64(self: *Mt19937_64) u64 {
        if (self.index == STATE_WORDS) self.twist();

        var value = self.state[self.index];
        self.index += 1;
        value ^= (value >> 29) & 0x5555555555555555;
        value ^= (value << 17) & 0x71d67fffeda60000;
        value ^= (value << 37) & 0xfff7eee000000000;
        value ^= value >> 43;
        return value;
    }

    pub fn nextUnitFloat(self: *Mt19937_64) f64 {
        const top_53 = self.nextU64() >> 11;
        return @as(f64, @floatFromInt(top_53)) * (1.0 / 9007199254740992.0);
    }

    fn twist(self: *Mt19937_64) void {
        for (0..STATE_WORDS - MIDDLE_WORD) |index| {
            const joined = (self.state[index] & UPPER_MASK) |
                (self.state[index + 1] & LOWER_MASK);
            self.state[index] = self.state[index + MIDDLE_WORD] ^
                (joined >> 1) ^ matrixTerm(joined);
        }
        for (STATE_WORDS - MIDDLE_WORD..STATE_WORDS - 1) |index| {
            const joined = (self.state[index] & UPPER_MASK) |
                (self.state[index + 1] & LOWER_MASK);
            self.state[index] = self.state[index + MIDDLE_WORD - STATE_WORDS] ^
                (joined >> 1) ^ matrixTerm(joined);
        }
        const joined = (self.state[STATE_WORDS - 1] & UPPER_MASK) |
            (self.state[0] & LOWER_MASK);
        self.state[STATE_WORDS - 1] = self.state[MIDDLE_WORD - 1] ^
            (joined >> 1) ^ matrixTerm(joined);
        self.index = 0;
    }

    fn matrixTerm(value: u64) u64 {
        return if (value & 1 == 0) 0 else MATRIX_A;
    }
};
