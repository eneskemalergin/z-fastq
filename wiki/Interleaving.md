# Interleaving and deinterleaving

These commands change pair layout. They do not trim, reorder, or repair records.

## Combine R1 and R2

```bash
z-fastq interleave R1.fastq.gz R2.fastq.gz > reads.interleaved.fastq
```

For each pair, z-fastq:

1. Reads one record from R1 and one from R2.
2. Validates both records.
3. Checks their names using the selected pair policy.
4. Writes R1, then R2.

The output is plain FASTQ with LF line endings. Use `--pair-names exact` when the complete first tokens must match:

```bash
z-fastq interleave --pair-names exact R1.fastq R2.fastq > reads.interleaved.fastq
```

To compress the output, use an external compressor:

```bash
z-fastq interleave R1.fastq R2.fastq | gzip > reads.interleaved.fastq.gz
```

## Split an interleaved stream

```bash
z-fastq deinterleave \
  --out1 R1.fastq \
  --out2 R2.fastq \
  reads.interleaved.fastq.gz
```

The first record in each pair goes to `--out1`; the second goes to `--out2`.

Output paths are created exclusively. z-fastq does not truncate or overwrite an existing file. It opens the input before creating the output paths, but a later validation or write failure can still leave empty or partial output files.

## Input and failure behavior

Both input commands accept plain FASTQ, gzip FASTQ, and explicit stdin. At most one input may be `-`. A record error, pair-name error, or unequal pair count stops the operation at the first failing pair. See [Paired reads](Paired-Reads) for name matching and [Automation](Automation) for exit status classes.
