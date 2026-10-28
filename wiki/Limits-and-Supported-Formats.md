# Limits and supported formats

I would rather state the edges clearly than let a command look broader than what has been tested. These are today's boundaries, not promises about the final project.

## Supported

- Linux on x86-64.
- Four-line FASTQ, plain or gzip, including concatenated gzip members; standard input as `-`.
- Phred+33 qualities, ASCII 33 to 126.
- IUPAC or ACGTN sequence letters.
- gzip output from `process --gzip`.

## Not supported

- FASTA, and wrapped sequence or quality lines.
- BGZF, ZIP, tar, bzip2, xz, and zstd.
- gzip output from commands other than `process`.
- Remote paths and object stores.
- Several threads per command.
- Repairing or skipping malformed records.
- From fastp: automatic adapter detection, automatic poly-G mode, and the average-quality, maximum-length, and low-complexity filters; writing the surviving mate of a dropped pair.

## Line length

The longest line accepted is 16 MiB by default. Raise it per command when a file really has longer lines, as long-read files can:

```bash
z-fastq check --max-line-bytes 33554432 reads.fastq
```

The limit counts the line's content, not its line ending. Raising it does not allow wrapped records.

## Memory

What a command keeps decides its memory, and the line limit bounds it.

| Command                                 | Keeps                                                       |
| --------------------------------------- | ----------------------------------------------------------- |
| `count`, `check`, `stats`               | no record bytes: a few MiB at any line length               |
| `check --paired`, `check --interleaved` | each mate's header, to compare names                        |
| `sample`, `interleave`, `deinterleave`  | a record too large for the read buffer, until it is written |
| `sample --count K`                      | about 4 bytes per chosen record                             |
| `process`                               | the same, plus the QC tables below when QC is requested     |

With 15 MiB lines, `count`, `check`, and `stats` peak near 5 MiB, `check --paired` near 37 MiB, `sample` near 50 MiB, and `interleave` near 100 MiB, so raise `--max-line-bytes` only as far as the data needs.

When `process` measures QC, per mate and per population (raw and retained) it keeps up to 800 KiB of quality tables, 48 KiB of adapter counters, up to 5 MB for the 100,000 sequence keys behind duplication, and 8 bytes per tile and cycle for per-tile quality: about 0.7 MB for 600 tiles of 151-base reads, at most 21 MB. `--pair-duplication` adds about 10 bytes per pair, at most 1 GiB. Trimming itself uses under 16 KiB at any read length.

## gzip

- A gzip header field (file name, comment, or extra field) may be at most 64 KiB.
- Bytes after the last gzip member are an error, even zero padding that `gzip -dc` accepts.

## Related pages

- [Input and output](Input-and-Output): how gzip and standard input are read.
- [Choosing a FASTQ tool](Choosing-a-FASTQ-Tool): tools that cover what z-fastq does not.
- [Troubleshooting](Troubleshooting#a-line-is-too-long): raising the line limit.
