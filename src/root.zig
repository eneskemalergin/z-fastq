//! Public library surface for z-fastq.

const fastq = @import("fastq.zig");
const io_layer = @import("io.zig");
const statistics = @import("stats.zig");

pub const Record = fastq.Record;
pub const OwnedRecord = fastq.OwnedRecord;
pub const toOwned = fastq.toOwned;
pub const Reader = fastq.Reader;
pub const RecordOffsets = fastq.RecordOffsets;
pub const Writer = fastq.Writer;
pub const WriterError = fastq.WriterError;
pub const ParseError = fastq.ParseError;
pub const LintCode = fastq.LintCode;
pub const ReaderError = fastq.ReaderError;
pub const codeTag = fastq.codeTag;
pub const Alphabet = fastq.Alphabet;
pub const ValidationOptions = fastq.ValidationOptions;
pub const SemanticField = fastq.SemanticField;
pub const SemanticError = fastq.SemanticError;
pub const validateRecord = fastq.validateRecord;

pub const Stats = statistics.Stats;
pub const StatsResult = statistics.StatsResult;
pub const StatsError = statistics.StatsError;
pub const QualityError = statistics.QualityError;
pub const decodePhred33 = fastq.decodePhred33;

pub const io = struct {
    pub const ByteSource = io_layer.ByteSource;
    pub const ByteSink = io_layer.ByteSink;
    pub const limits = Limits;
    pub const plain = struct {
        pub const SliceSource = io_layer.SliceSource;
        pub const ReaderSource = io_layer.ReaderSource;
        pub const FileSource = io_layer.FileSource;
        pub const SliceSink = io_layer.SliceSink;
        pub const WriterSink = io_layer.WriterSink;
        pub const FileSink = io_layer.FileSink;
    };
    pub const gzip = struct {
        pub const ReaderSource = io_layer.GzipSource;
    };
};

const Limits = struct {
    pub const DEFAULT_MAX_LINE_BYTES = io_layer.DEFAULT_MAX_LINE_BYTES;
    pub const DEFAULT_READER_BUFFER_BYTES = io_layer.DEFAULT_READER_BUFFER_BYTES;
    pub const COUNT_READ_BUFFER_BYTES = io_layer.COUNT_READ_BUFFER_BYTES;
};

pub const limits = Limits;
pub const count_scan = @import("count_scan.zig");

pub const VERSION = "0.0.10";
