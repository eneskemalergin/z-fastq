# Command cheat sheet

Use this page to choose a command. Open the linked page when the command's edge cases matter.

## Pick the command

`count` prints one record count for each independent input:

```bash
z-fastq count reads.fastq.gz
```

`stats` prints aggregate read, base, composition, and quality fields:

```bash
z-fastq stats reads.fastq.gz
```

`check` validates one input:

```bash
z-fastq check reads.fastq.gz
```

`sample` writes a selected FASTQ stream to standard output:

```bash
z-fastq sample --fraction 0.1 --seed 11 reads.fastq.gz > sample.fastq
```

`interleave` validates two mates and writes R1, then R2, to standard output:

```bash
z-fastq interleave R1.fastq.gz R2.fastq.gz > interleaved.fastq
```

`deinterleave` validates one interleaved input and creates two output files:

```bash
z-fastq deinterleave --out1 R1.fastq --out2 R2.fastq interleaved.fastq
```

`process` validates two mates into two new files, plain or with `--gzip`, and can write a JSON report:

```bash
z-fastq process --gzip --out1 clean_R1.fastq.gz --out2 clean_R2.fastq.gz --report sample.json R1.fastq.gz R2.fastq.gz
```

With adapters, poly-G, and a 3' quality window, `process` trims each pair by fastp 1.1.0's rules:

```bash
z-fastq process --adapter1 CTGTCTCTTATACACATCT --adapter2 CTGTCTCTTATACACATCT --trim-poly-g 10 --trim-tail-window 4 --trim-tail-quality 20 --out1 trimmed_R1.fastq --out2 trimmed_R2.fastq R1.fastq.gz R2.fastq.gz
```

## Input shapes

```text
count         <path|-> [<path|-> ...]
stats         <path|-> [<path|-> ...]
check         <path|-> [<path|-> ...]
check --paired R1 R2
check --interleaved interleaved.fastq
sample        <path|->
sample --paired R1 R2
sample --interleaved interleaved.fastq
interleave    R1 R2
deinterleave --out1 R1 --out2 R2 interleaved.fastq
process --out1 R1 --out2 R2 R1-input R2-input
```

`sample` requires exactly one ordinary input, two inputs for `--paired`, and one input for `--interleaved`. Exact-count forms require regular file paths because they use two passes.

## Shared options

```text
-h, --help
--max-line-bytes N
```

`-h` and `--help` work before or after a command. `-V` and `--version` work only before a command:

```bash
z-fastq --version
```

The default maximum logical line length is 16 MiB. An option that takes a value may be given once; repeating it fails with status `2` instead of letting the last value win. Use `--` before a path that begins with a hyphen:

```bash
z-fastq count -- ./-named.fastq
```

## Validation and pair options

```text
--alphabet iupac
--alphabet acgtn
--paired
--interleaved
--pair-names illumina
--pair-names exact
```

The pair flags apply to `check` and `sample`. The standalone pair commands accept `--pair-names` directly.

## Sample options

```text
--fraction P
--count K
--seed S
```

`--fraction` and `--count` are mutually exclusive. The default seed is `11`.

## Machine output

```bash
z-fastq stats --json reads.fastq.gz
z-fastq check --json reads.fastq.gz
```

The current build supports `--json` for `stats` and `check`. Version `0.1.0` freezes the `z-fastq/stats-v1` and `z-fastq/check-v1` document shapes. See [Automation](Automation).
