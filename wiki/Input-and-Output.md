# Input and output

## Plain and gzip input

z-fastq accepts plain FASTQ and gzip-compressed FASTQ. It chooses the decoder from the first two input bytes, not from the filename suffix.

That means both of these are valid situations:

- `reads.fastq` contains gzip bytes.
- `reads.gz` contains plain FASTQ.

If the first bytes identify gzip, z-fastq commits to gzip decoding. A damaged gzip stream is an input error. The tool does not retry it as plain text.

The gzip reader checks the gzip header, DEFLATE stream, CRC32, and uncompressed size. Concatenated gzip members are supported.

CRC32 and uncompressed size are stored at the end of each gzip member, so they are checked after that member has been decoded. A command can process and write records from a member before it finds that the member's check fails. It still exits with status `3`, and the FASTQ it wrote is the beginning of the output a valid input would give. See [Automation](Automation) for handling partial output.

## Standard input

Use `-` explicitly:

```bash
gzip -dc reads.fastq.gz | z-fastq stats -
```

An omitted path does not mean stdin. A command accepts stdin at most once. When stdin is a pipe or FIFO, that includes other names for it, such as `/dev/stdin`, because a stream can be read only once; a regular file on stdin may also be named by its path. Exact-count sampling needs a regular file because it reads the input twice.

## Standard output

`count`, `stats`, and successful JSON output go to standard output. Diagnostics go to standard error. Sampling and interleaving commands write FASTQ to standard output.

FASTQ output uses four logical lines per record and LF line endings. z-fastq does not currently write gzip output. Use an external compressor:

```bash
z-fastq sample --fraction 0.1 reads.fastq | gzip > sample.fastq.gz
```

`deinterleave` is the exception: it writes two requested output files.

Standard output may not be one of the inputs: `z-fastq count a.fastq b.fastq >> b.fastq` exits with status `2` and `input file is output file` before reading anything, so the input is never changed. Every command that writes to standard output applies this check, including when stdin is redirected from the output file.

If standard output fails, including a pipe whose reader has exited, the command reports `error: standard output: I/O error` and exits with status `3`.

## Multiple inputs

`count` and `stats` accept multiple independent inputs. They process them in argument order and keep their results separate. Paired commands have fixed shapes:

- `check --paired` and `sample --paired` take two inputs.
- `check --interleaved` and `sample --interleaved` take one input.
- `interleave` takes two inputs.
- `deinterleave` takes one input and two output paths.

## Streaming boundary

The parser reads through bounded byte sources. It does not map or load a complete FASTQ file into memory for normal count, stats, validation, sampling, or pair operations. Exact-count sampling keeps selection state for the requested count and makes a second pass over the file.

See [Limits and supported formats](Limits-and-Supported-Formats) for the line and format boundaries.

