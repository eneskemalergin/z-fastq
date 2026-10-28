# Command cheat sheet

The whole command line on one page: what each command is for, how to call it, and every option. Each section links to the page that explains it in depth.

`<path|->` means a file, or `-` for standard input. Every command reads plain FASTQ and gzip alike.

## The commands

### count

Counts records, as fast as the file can be read. See [Counting and statistics](Counting-and-Statistics).

```bash
z-fastq count <path|-> [<path|-> ...]
z-fastq count lane1.fastq.gz lane2.fastq.gz
```

### stats

Summarizes reads, bases, length, GC, and quality. See [Counting and statistics](Counting-and-Statistics#stats).

```bash
z-fastq stats [--json] <path|-> [<path|-> ...]
z-fastq stats --json reads.fastq.gz > stats.json
```

### check

Validates every record, or every pair with `--paired` or `--interleaved`. Prints nothing on success. See [Validation](Validation) and [Paired reads](Paired-Reads).

```bash
z-fastq check [--json] <path|-> [<path|-> ...]
z-fastq check --paired R1.fastq.gz R2.fastq.gz
z-fastq check --interleaved interleaved.fastq.gz
```

### sample

Writes a reproducible subset to standard output, by fraction or by exact count. Exact counts need regular files. See [Sampling](Sampling).

```bash
z-fastq sample (--fraction P | --count K) [--seed S] <path|->
z-fastq sample --paired --count 100000 --seed 7 R1.fastq.gz R2.fastq.gz > sample.fastq
```

### interleave and deinterleave

Join R1 and R2 into one interleaved stream, or split one into two new files. See [Interleaving](Interleaving).

```bash
z-fastq interleave R1.fastq.gz R2.fastq.gz > interleaved.fastq
z-fastq deinterleave --out1 R1.fastq --out2 R2.fastq interleaved.fastq.gz
```

### process

Validates, trims, filters, and measures pairs into two new files, with an optional report and MultiQC files. See [Processing paired reads](Processing) and [Quality control](Quality-Control).

```bash
z-fastq process --out1 OUT1 --out2 OUT2 [options] R1 R2
z-fastq process --gzip --adapter1 CTGTCTCTTATACACATCT --adapter2 CTGTCTCTTATACACATCT \
  --trim-poly-g 10 --min-length 30 --report sample.json \
  --out1 trimmed_R1.fastq.gz --out2 trimmed_R2.fastq.gz R1.fastq.gz R2.fastq.gz
```

## Options

### Everywhere

- `-h`, `--help`: help, before or after a command. `-V`, `--version`: the version, before a command.
- `--max-line-bytes N`: the longest line accepted, 16 MiB by default ([Limits](Limits-and-Supported-Formats#line-length)).
- `--`: ends the options, for a path that starts with `-`, as in `z-fastq count -- ./-odd.fastq`.

An option that takes a value may be given only once; a repeat is an error rather than the last one winning.

### Records and pairs

- `--alphabet iupac|acgtn`: which sequence letters pass, `iupac` by default. Every command but `count` and `stats` takes it.
- `--paired`, `--interleaved`: for `check` and `sample`, read two mate files or one interleaved file.
- `--pair-names illumina|exact`: how mate names are matched, `illumina` by default ([Paired reads](Paired-Reads#how-names-are-matched)). For every pair command and `process`.
- `--json`: machine-readable output from `stats` and `check` ([Automation](Automation#json)).

### Sampling

- `--fraction P`: keep each record with chance P, in one pass.
- `--count K`: keep exactly K records, in two passes.
- `--seed S`: the seed, 11 by default.

### Processing

Outputs:

- `--out1 PATH`, `--out2 PATH`: the two new FASTQ files.
- `--gzip`, `--gzip-preset fast|even|dense`: compress them, with the `fast` preset by default.
- `--report PATH`: the JSON report. `--qc-tables DIR` and `--qc-summary PATH`: files for MultiQC.
- `--max-pairs N`: stop after N pairs.

Trims, each off unless given:

- `--adapter1 A1 --adapter2 A2`: adapters and mate overlaps.
- `--trim-poly-g N`: poly-G tails, scanning at least N bases.
- `--trim-tail-window W --trim-tail-quality Q`: low-quality 3' ends.

Filters, each off unless given:

- `--low-quality Q --max-low-quality-percent P`: too many bases below quality Q.
- `--max-n N`: more than N `N` bases.
- `--min-length N`: shorter than N bases.

Measurements:

- `--contaminants PATH`: name the likely source of overrepresented sequences.
- `--pair-duplication`: the pair duplication rate.

[Processing](Processing) gives each option's range, and [Quality control](Quality-Control) says which options need each other.

## Related pages

- [Getting started](Getting-Started): the main commands, one at a time.
- [Automation](Automation): exit statuses and JSON output.
- [Recipes](Recipes): complete commands for everyday jobs.
