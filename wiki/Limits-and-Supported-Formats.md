# Limits and supported formats

## Supported today

- Linux x86-64 builds.
- Plain four-line FASTQ input.
- gzip-compressed FASTQ input, including concatenated gzip members.
- Explicit standard input with `-`.
- Phred+33 quality bytes from ASCII 33 through 126.
- IUPAC sequence validation, or the narrower ACGTN policy.

## Not supported today

- FASTA input.
- Wrapped sequence or quality lines.
- gzip output from the CLI.
- BGZF, ZIP, tar, bzip2, xz, and zstd containers.
- Remote object stores or network paths.
- Memory-mapped or worker-thread CLI execution.
- Automatic repair or record recovery after malformed input.

These are current boundaries, not promises about the final project. I would rather state them clearly than make a command appear broader than its tested behavior.

## Logical line limit

The default maximum logical line length is 16 MiB. Override it per command:

```bash
z-fastq check --max-line-bytes 33554432 reads.fastq
```

The limit applies to the line content, not the LF terminator. Raising it allows longer single-line fields; it does not enable wrapped FASTQ.

## Location semantics

Validation locations use:

- A zero-based record index.
- A zero-based byte offset in the decompressed FASTQ stream.
- A one-based line number within the record.

For gzip input, the byte offset is not a compressed-file offset. It points into the decoded FASTQ stream that the parser sees.

## Record shape

The current parser expects one header line, one sequence line, one plus line, and one quality line per record. The sequence and quality fields must have equal byte lengths. Empty fields are valid when both are empty.

