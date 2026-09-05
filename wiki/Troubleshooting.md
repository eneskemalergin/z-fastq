# Troubleshooting

## `z-fastq` is not found

Run the binary from the build output:

```bash
./zig-out/bin/z-fastq --version
```

If that works, place a copy or link in a directory on your `PATH`.

## The build fails while assembling ISA-L

The default Linux x86-64 release build needs NASM `2.14.01` or newer. Install NASM or use the Zig fallback:

```bash
zig build -Dstatic=true -Doptimize=ReleaseFast -Disa-l=false
```

## A gzip input fails even though the name is unusual

z-fastq checks the first two bytes, not `.gz`. If those bytes identify gzip, the complete stream must have valid framing, DEFLATE data, CRC32, and uncompressed size.

Check it independently:

```bash
gzip -t reads.fastq.gz
```

## A line-limit error appears

The default logical line limit is 16 MiB. Increase it only when the input really uses longer single-line fields:

```bash
z-fastq check --max-line-bytes 33554432 reads.fastq
```

This does not enable wrapped FASTQ.

## Exact sampling rejects `-`

That mode needs two passes. Give it a regular path:

```bash
z-fastq sample --count 1000 reads.fastq
```

Use fraction sampling when the source is a pipe:

```bash
gzip -dc reads.fastq.gz | z-fastq sample --fraction 0.1 -
```

## Paired names do not match

Start with the default `illumina` policy. If your headers use a different convention, inspect the first tokens and try the strict policy:

```bash
z-fastq check --paired --pair-names exact R1.fastq R2.fastq
```

Changing policy does not make genuinely different read identifiers equal. See [Paired reads](Paired-Reads).

## `deinterleave` refuses an output path

Output files are created exclusively. The command does not truncate or overwrite an existing file. Choose a new path after checking that the old output is safe to keep or remove.

## A FASTQ-producing command failed after writing data

Sampling and interleaving stream output. A late parser, pair, input, or output failure can leave a partial stream. Write to a temporary path and move it only after a zero exit status.

## JSON is missing

The current pre-alpha build supports `--json` only for `stats` and `check`. `count` prints decimal counts. Sampling and pair-layout commands write FASTQ. The JSON document shape is provisional; see [Automation](Automation) before depending on it.

## The command returns status `4`

Status `4` means z-fastq reached a configured or representable limit: a line limit, sample-selection or staging limit, or arithmetic limit. The error message identifies the relevant boundary when the command can report it.
