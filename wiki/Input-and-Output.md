# Input and output

## Plain and gzip input

z-fastq accepts plain FASTQ and gzip-compressed FASTQ. It chooses the decoder from the first two input bytes, not from the filename suffix.

That means both of these are valid situations:

- `reads.fastq` contains gzip bytes.
- `reads.gz` contains plain FASTQ.

If the first bytes identify gzip, z-fastq commits to gzip decoding. A damaged gzip stream is an input error. The tool does not retry it as plain text.

The gzip reader checks the gzip header, DEFLATE stream, CRC32, and uncompressed size. Concatenated gzip members are supported.

## Standard input

Use `-` explicitly:

```bash
gzip -dc reads.fastq.gz | z-fastq stats -
```

An omitted path does not mean stdin. A command accepts stdin at most once. Exact-count sampling needs a regular file because it reads the input twice.

## Standard output

`count`, `stats`, and successful JSON output go to standard output. Diagnostics go to standard error. Sampling and interleaving commands write FASTQ to standard output.

FASTQ output uses four logical lines per record and LF line endings. z-fastq does not currently write gzip output. Use an external compressor:

```bash
z-fastq sample --fraction 0.1 reads.fastq | gzip > sample.fastq.gz
```

`deinterleave` is the exception: it writes two requested output files.

## Multiple inputs

`count` and `stats` accept multiple independent inputs. They process them in argument order and keep their results separate. Paired commands have fixed shapes:

- `check --paired` and `sample --paired` take two inputs.
- `check --interleaved` and `sample --interleaved` take one input.
- `interleave` takes two inputs.
- `deinterleave` takes one input and two output paths.

## Streaming boundary

The parser reads through bounded byte sources. It does not map or load a complete FASTQ file into memory for normal count, stats, validation, sampling, or pair operations. Exact-count sampling keeps selection state for the requested count and makes a second pass over the file.

See [Limits and supported formats](Limits-and-Supported-Formats) for the line and format boundaries.

