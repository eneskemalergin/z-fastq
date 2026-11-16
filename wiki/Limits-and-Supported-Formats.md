# Limits and supported formats

I would rather state the edges clearly than let a command look broader than what has been tested. These are today's boundaries, not promises about the final project.

## Supported

- Linux on x86-64 or aarch64, and macOS 13 or later on Apple Silicon or Intel.
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
| `count`, `check`, `stats`               | no record bytes                                             |
| `check --paired`, `check --interleaved` | each mate's header, to compare names                        |
| `sample`, `interleave`, `deinterleave`  | one record too large for the read buffer                    |
| `sample --count K`                      | about 4 bytes per chosen record                             |
| `process`                               | as above, plus the QC tables below                          |

With 15 MiB lines, `count`, `check`, and `stats` stay under 1 MiB. `check --paired` peaks near 31 MiB when the headers are that long, and near 1 MiB when only sequence and quality are. When header, sequence, and quality are all 15 MiB, `sample` peaks near 46 MiB, and `interleave` and `deinterleave` near 91 MiB (61 MiB for `interleave` with short headers). Raise `--max-line-bytes` only as far as the data needs.

When `process` measures QC, it keeps, per mate and per population (raw and retained):

- Up to 800 KiB of quality tables and 48 KiB of adapter counters.
- Up to 4.3 MiB for the 100,000 sequence keys behind duplication, or 3.2 MiB when the keys hold only A, C, G, and T, which take two bits per base. When trimming keeps the start of each read, the retained population reuses the raw keys and needs about 1.4 MiB.
- 8 bytes per tile and cycle for per-tile quality: about 0.7 MiB for 600 tiles of 151-base reads, at most about 20 MiB.

`--pair-duplication` adds about 4 to 6 bytes per pair, at most 1 GiB. Each stored sequence key counts up to 4,294,967,295 reads; a run where one key passes that stops with status `4`, a statistics arithmetic limit. Trimming itself uses under 16 KiB at any read length.

## gzip

- The optional fields of a gzip header (file name, comment, and extra field) may hold at most 64 KiB together, counting the extra field's two length bytes and the zero byte that ends a name or comment.
- Bytes after the last gzip member are an error, even zero padding that `gzip -dc` accepts.

## Related pages

- [Input and output](Input-and-Output): how gzip and standard input are read.
- [Choosing a FASTQ tool](Choosing-a-FASTQ-Tool): tools that cover what z-fastq does not.
- [Troubleshooting](Troubleshooting#a-line-is-too-long): raising the line limit.
