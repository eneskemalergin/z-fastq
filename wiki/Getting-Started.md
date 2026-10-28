# Getting started

This is the path I recommend for a first look: a tiny file, and one command at a time, so you can see exactly what each one does. It takes a few minutes, and assumes `z-fastq` is on your `PATH` ([Installation](Installation)).

## 1. Make a file

Save this as `reads.fastq`. It holds two records, one of five bases and one of three:

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

## 2. Check it

```bash
z-fastq check reads.fastq
```

Nothing is printed. For `check`, silence and exit status `0` mean the file is valid, which is exactly what a script wants. [Validation](Validation) explains what is checked.

## 3. Count and summarize it

```bash
z-fastq count reads.fastq
z-fastq stats reads.fastq
```

`count` prints `2`. `stats` adds bases, lengths, base composition, GC, and quality; [Counting and statistics](Counting-and-Statistics) explains each number.

## 4. Read gzip and pipes

```bash
gzip -c reads.fastq > reads.fastq.gz
z-fastq stats reads.fastq.gz
gzip -dc reads.fastq.gz | z-fastq count -
```

z-fastq recognizes gzip by its first bytes rather than by the `.gz` name, and reads standard input only when you write `-` ([Input and output](Input-and-Output)).

## 5. Take a sample

```bash
z-fastq sample --fraction 0.5 --seed 42 reads.fastq > sample.fastq
```

Run it twice: the same seed and input give the same records. [Sampling](Sampling) also covers exact counts and pairs.

## Where next

- **Paired reads?** Read [Paired reads](Paired-Reads), then [Processing paired reads](Processing).
- **A script or workflow?** Read [Automation](Automation) for exit statuses and partial output.
- **An option?** The [Command cheat sheet](Command-Cheat-Sheet) has them all.
