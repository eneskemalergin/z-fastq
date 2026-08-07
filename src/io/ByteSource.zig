//! Pull interface for byte-oriented FASTQ input.

pub const ReadError = error{
    ReadFailed,
};

pub const ByteSource = struct {
    vtable: *const VTable,
    ctx: *anyopaque,

    pub const VTable = struct {
        read: *const fn (ctx: *anyopaque, dest: []u8) ReadError!usize,
    };

    /// Return bytes written to `dest`, or zero at end of input.
    pub fn read(self: *const ByteSource, dest: []u8) ReadError!usize {
        return self.vtable.read(self.ctx, dest);
    }
};
