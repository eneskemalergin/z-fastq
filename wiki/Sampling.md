# Sampling

`sample` writes selected FASTQ records to standard output. It supports a reproducible fraction mode and an exact-count mode. The two modes have different input requirements.

## Select a fraction

```bash
z-fastq sample --fraction 0.1 --seed 42 reads.fastq.gz > sample.fastq
```

The selector makes one deterministic decision per record. The output count is therefore approximate, not an exact tenth of the input.

Accepted fraction forms are `0`, `1`, decimals beginning with `0.`, and decimals beginning with `1.` whose remaining digits are zero. Percentages and scientific notation are not accepted.

`--fraction 0` and `--count 0` only validate the input: they write nothing, so a record that could not be written (a field ending in CR) does not fail them, and a closed standard output still exits `0`, as with `check`.

The default seed is `11`. Set it explicitly when you need the current selection to repeat:

```bash
z-fastq sample --fraction 0.01 --seed 2025 reads.fastq > sample-a.fastq
z-fastq sample --fraction 0.01 --seed 2025 reads.fastq > sample-b.fastq
cmp sample-a.fastq sample-b.fastq
```

The same seed, input record stream, and sampling options repeat the selection in `0.1.0`. Pin the package version when exact selected records must remain reproducible across upgrades.

## Select an exact count

```bash
z-fastq sample --count 1000 --seed 42 reads.fastq > sample.fastq
```

The command selects `min(1000, record_count)` records: when the input has no more records than the count, every record is written and no note is printed. It uses two passes over a regular file, so exact-count mode rejects stdin. The file must not change between the passes; z-fastq compares its inode, size, and modification time before and after each pass and fails with `input changed during exact sampling` when they differ, but a rewrite that keeps all three is not detected.

The selection state grows with the requested count, not with the size of the FASTQ file: once the input has more records than the count, z-fastq reserves about 4 bytes per requested record (8 bytes above 2^32 records), so `--count 10000000` holds about 40 MB.

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

