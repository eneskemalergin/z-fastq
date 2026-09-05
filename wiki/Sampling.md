# Sampling

`sample` writes selected FASTQ records to standard output. It supports a reproducible fraction mode and an exact-count mode. The two modes have different input requirements.

## Select a fraction

```bash
z-fastq sample --fraction 0.1 --seed 42 reads.fastq.gz > sample.fastq
```

The selector makes one deterministic decision per record. The output count is therefore approximate, not an exact tenth of the input.

Accepted fraction forms are `0`, `1`, decimals beginning with `0.`, and decimals beginning with `1.` whose remaining digits are zero. Percentages and scientific notation are not accepted.

The default seed is `11`. Set it explicitly when you need the current selection to repeat:

```bash
z-fastq sample --fraction 0.01 --seed 2025 reads.fastq > sample-a.fastq
z-fastq sample --fraction 0.01 --seed 2025 reads.fastq > sample-b.fastq
cmp sample-a.fastq sample-b.fastq
```

The same seed and input record stream repeat the current selection. The project may change its generator before the `0.1.0` interface settles.

## Select an exact count

```bash
z-fastq sample --count 1000 --seed 42 reads.fastq > sample.fastq
```

The command selects `min(1000, record_count)` records. It uses two passes over a regular file, so exact-count mode rejects stdin. The selection state grows with the requested count, not with the size of the FASTQ file.

## Sample pairs

For two mate files, the output is an interleaved stream:

```bash
z-fastq sample --paired --fraction 0.1 --seed 42 R1.fastq R2.fastq > sample.interleaved.fastq
```

For an already interleaved input:

```bash
z-fastq sample --interleaved --count 1000 --seed 42 reads.interleaved.fastq > sample.interleaved.fastq
```

The pair is the sampling unit. z-fastq never selects R1 and R2 independently. Paired exact-count sampling requires regular file paths for both mates.

## Validation and partial output

z-fastq validates records before writing selected records. A later malformed record can still leave earlier output on stdout. Do not replace a trusted file until the command succeeds:

```bash
if z-fastq sample --fraction 0.1 --seed 42 reads.fastq > sample.tmp.fastq; then
    mv sample.tmp.fastq sample.fastq
else
    rm -f sample.tmp.fastq
    exit 1
fi
```

This is a stream-oriented command. It does not write gzip output. Pipe the completed stream through `gzip` when needed.

