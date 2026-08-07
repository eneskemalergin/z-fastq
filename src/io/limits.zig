//! Hard caps for streaming FASTQ I/O.

pub const default_max_line_bytes: usize = 16 * 1024 * 1024;
pub const default_reader_buffer_bytes: usize = 256 * 1024;
/// `count` uses one stack buffer for streaming reads and the scan window.
pub const count_read_buffer_bytes: usize = default_reader_buffer_bytes;
