//! Hard caps for streaming FASTQ I/O.

pub const default_max_line_bytes: usize = 16 * 1024 * 1024;
pub const default_reader_buffer_bytes: usize = 256 * 1024;
/// `count` uses one stack buffer (`count_read_buffer_bytes`) split between I/O
/// syscall batching and the scan window.
pub const count_read_buffer_bytes: usize = default_reader_buffer_bytes;
pub const file_io_buffer_bytes: usize = 8 * 1024;
