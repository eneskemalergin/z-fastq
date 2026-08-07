//! Push interface for byte-oriented FASTQ output.

pub const WriteError = error{
    WriteFailed,
};

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
        if (self.vtable.flush) |flush_fn| {
            return flush_fn(self.ctx);
        }
    }
};
