//! I/O layer re-exports (`ByteSource`, `ByteSink`, plain adapters, limits).

pub const ByteSource = @import("ByteSource.zig").ByteSource;
pub const ByteSink = @import("ByteSink.zig").ByteSink;
pub const limits = @import("limits.zig");
pub const plain = @import("plain.zig");
