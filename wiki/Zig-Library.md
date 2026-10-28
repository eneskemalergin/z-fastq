# Zig library

The command line is the main product, but the reader, writer, validation, and statistics behind it are also a small Zig module, for programs that want to handle FASTQ directly.

## Import it

With the repository as a Zig dependency, the module is named `z-fastq`:

```zig
const zfastq = @import("z-fastq");
```

It is registered from `src/root.zig` and targets Zig 0.16.0 on Linux x86-64.

## What it exports

| Name                                                                  | For                                  |
| --------------------------------------------------------------------- | ------------------------------------ |
| `Record`, `OwnedRecord`, `toOwned`                                    | records, borrowed and owned          |
| `Reader`, `RecordOffsets`, `ParseError`, `ReaderError`                | reading records from a byte source   |
| `Writer` and its errors                                               | writing four-line records            |
| `LintCode`, `codeTag`, `Alphabet`, `ValidationOptions`                | validation rules and codes           |
| `SemanticField`, `SemanticError`, `validateRecord`                    | checking one record                  |
| `Stats`, `StatsResult`, `StatsError`, `QualityError`, `decodePhred33` | statistics                           |
| `io.ByteSource`, `io.ByteSink`, and adapters                          | plain, file, writer, and gzip I/O    |
| `count_scan`                                                          | counting records without allocating  |
| `limits`, `VERSION`                                                   | default sizes and the version string |

Option parsing, pair-name matching, sampling, JSON output, and `process` are part of the command line, not of the library.

## Borrowed and owned records

`Reader.next` returns a borrowed `Record`: its fields point into the reader's buffer and stay valid until the next call or `deinit`. That is what keeps streaming cheap. When a record must live longer, copy it:

```zig
const record = (try reader.next()) orelse return;
var owned = try zfastq.toOwned(allocator, record);
defer owned.deinit();
```

A record's `id` is the first space- or tab-separated word of its `header`, in the same bytes.

## Writing

`Writer` writes four-line records with LF endings. It refuses a record before writing any of it when the sequence and quality lengths differ, the header is empty or starts with a space or tab, or a field holds a line break or ends in a carriage return. The output is complete only after `flush`, and the sink it writes to must outlive it.

## Validation and statistics

`validateRecord` checks the sequence letters and Phred+33 qualities of one record without allocating; the record structure is the reader's job. `decodePhred33` refuses bytes outside ASCII 33 to 126.

`Stats` allocates nothing. `addRecord` validates a record and leaves every counter unchanged if it refuses it; `takeLastQualityError` gives the details. In `StatsResult`, the length fields are `null` when no read was added, and a ratio is `null` when its denominator is zero.

## Sources, sinks, and errors

The caller opens and closes files; the adapters only borrow them. `ByteSource` and `ByteSink` are small values that point at their adapter, so the adapter must stay at one address and outlive every copy.

- A read fills part of the destination and returns its length; zero for a nonempty destination means the end.
- A write consumes the whole slice; a sink without a flush callback treats `flush` as a no-op.
- `FileSink` writes whole pages and keeps a partial last page until the next write or `flush`.

After an error, read the diagnostic with `takeLastError` and stop: do not call `next`, `feed`, or `finishEof` again, because the position may reflect a half-read record.

Pin the package version if your program depends on exact names or error shapes, and read the changelog before upgrading.

## Related pages

- [Validation](Validation): the record rules the reader and `validateRecord` apply.
- [Counting and statistics](Counting-and-Statistics#stats): the meaning of each `Stats` field.
- [For contributors](For-Contributors): where each part lives in the source.
