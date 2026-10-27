# Limits and supported formats

## Supported today

- Linux x86-64 builds.
- Plain four-line FASTQ input.
- gzip-compressed FASTQ input, including concatenated gzip members.
- gzip output from `process` (`--gzip`).
- Explicit standard input with `-`.
- Phred+33 quality bytes from ASCII 33 through 126.
- IUPAC sequence validation, or the narrower ACGTN policy.

## Not supported today

- FASTA input.
- Wrapped sequence or quality lines.
- gzip output from commands other than `process`.
- BGZF, ZIP, tar, bzip2, xz, and zstd containers.
- Remote object stores or network paths.
- Memory-mapped or worker-thread CLI execution.
- Automatic repair or record recovery after malformed input.
- fastp's average-quality, maximum-length, and low-complexity filters, and writing the surviving mate of a dropped pair.
- fastp's automatic adapter detection and automatic poly-G mode.

These are current boundaries, not promises about the final project. I would rather state them clearly than make a command appear broader than its tested behavior.

## Logical line limit

The default maximum logical line length is 16 MiB. Override it per command:

```bash
z-fastq check --max-line-bytes 33554432 reads.fastq
```

The limit applies to the line content, not the LF terminator. Raising it allows longer single-line fields; it does not enable wrapped FASTQ.

## Memory and the line limit

Memory follows what a command keeps, and the line limit bounds that:

- `count`, `check`, `stats`, and sampling that selects nothing keep no record bytes: they stream lines of any length, so memory stays at a few MiB.
- `check --paired` and `check --interleaved` keep each mate's header line while they compare names.
- Commands that write records (`sample`, `interleave`, `deinterleave`, `process`) keep a record that does not fit the read buffer until it is written: up to four lines per input, plus mate 1 for interleaved input.
- `process` also keeps its quality tables, up to 800 KiB per mate for the reads as read and as many again in a trimmed or filtered run, plus 48 KiB of adapter counters per mate and population. Trimming itself uses under 16 KiB, the same at any read length.

With 15 MiB header, sequence, and quality lines, `stats` peaks near 5 MiB, `sample` near 50 MiB, and the pair commands near 100 MiB. Raise `--max-line-bytes` only as far as the data needs, because the worst case grows with it.

## gzip limits

- An optional gzip header field (file name, comment, or extra field) may be at most 64 KiB; a longer one fails with an I/O error. `gzip -dc` accepts it.
- Bytes after the last gzip member fail, including zero padding that `gzip -dc` accepts.

## Location semantics

Validation locations use:

- A zero-based record index.
- A zero-based byte offset in the decompressed FASTQ stream.
- A one-based line number within the record.

For gzip input, the byte offset is not a compressed-file offset. It points into the decoded FASTQ stream that the parser sees.

## Record shape

The current parser expects one header line, one sequence line, one plus line, and one quality line per record. The sequence and quality fields must have equal byte lengths. Empty fields are valid when both are empty.

