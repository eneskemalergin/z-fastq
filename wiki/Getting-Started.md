# Getting started

This is the short path I recommend when you first try z-fastq. It starts with a tiny file so you can see what each command owns.

## 1. Create a FASTQ file

Save this as `reads.fastq`:

```fastq
@read-1
ACGTN
+
IIIII
@read-2
GGC
+
III
```

The file has two records. The first sequence has five bases and the second has three.

## 2. Check the file

```bash
z-fastq check reads.fastq
```

There is no success message. A zero exit status means the input passed structural, sequence, and quality validation.

Try the narrower alphabet policy too:

```bash
z-fastq check --alphabet acgtn reads.fastq
```

This file passes because it contains only A, C, G, T, and N.

## 3. Count the records

```bash
z-fastq count reads.fastq
```

Output:

```text
2
```

`count` is the command to use when you do not need to materialize or summarize the records.

## 4. Read the summary

```bash
z-fastq stats reads.fastq
```

The report includes read count, total bases, length range, nucleotide counts, GC fraction, and Phred quality summaries. [Counting and statistics](Counting-and-Statistics) defines each field.

## 5. Read compressed input through the same command

```bash
gzip -c reads.fastq > reads.fastq.gz
z-fastq stats reads.fastq.gz
```

z-fastq detects gzip from the input header. The `.gz` suffix is not what selects the decoder.

## 6. Read from a pipe

```bash
gzip -dc reads.fastq.gz | z-fastq count -
```

The `-` is explicit. z-fastq does not treat a missing path as standard input.

## 7. Make a repeatable sample

```bash
z-fastq sample --fraction 0.5 --seed 42 reads.fastq > sample.fastq
```

Use the same seed and input when you need the current selection to repeat. Read [Sampling](Sampling) before using exact-count mode or paired data.

## What to read next

- [Command cheat sheet](Command-Cheat-Sheet) for the complete invocation shapes.
- [Validation](Validation) for error codes and alphabet policies.
- [Paired reads](Paired-Reads) before processing R1 and R2 files.
- [Input and output](Input-and-Output) for gzip, stdin, and output boundaries.

