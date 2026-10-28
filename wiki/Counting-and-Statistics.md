# Counting and statistics

Two commands answer the first question asked of any FASTQ file: what is in it? `count` gives the number of records and nothing else, as fast as the file can be read. `stats` reads every base and quality score and summarizes the whole file.

## count

```bash
z-fastq count lane1.fastq.gz lane2.fastq.gz
```

It prints one number per input, in the order given. `count` checks that every record is well formed: four lines, the `@` and `+`, and a quality line as long as the sequence. It does not check the letters or the quality bytes; [`check`](Validation) does.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/eneskemalergin/z-fastq/main/bench/count/results/figures/file-overview-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/eneskemalergin/z-fastq/main/bench/count/results/figures/file-overview-light.png">
    <img alt="Wall time and peak memory of five other record counters relative to z-fastq, for plain and gzip input." src="https://raw.githubusercontent.com/eneskemalergin/z-fastq/main/bench/count/results/figures/file-overview-light.png" width="760">
  </picture>
</p>
<p align="center"><sub><b>Figure 1.</b> Time and peak memory of other record counters relative to z-fastq, on 1 GiB inputs and one CPU; lower left is better. The <a href="https://github.com/eneskemalergin/z-fastq/blob/main/bench/count/README.md">count benchmark</a> has the inputs, the tools, and every number.</sub></p>

## stats

```bash
z-fastq stats reads.fastq.gz
```

Each result starts with an `input` line naming the file, followed by these fields.

**Reads and bases.** `reads` counts records and `bases` counts sequence letters. `min_length`, `max_length`, and `mean_length` describe read lengths. An empty record still counts as a read, so the two can tell different stories.

**Composition.** `a`, `c`, `g`, `t`, and `n` count each base, upper or lower case, and `other_bases` counts any other letter. GC is measured only among the four definite bases, because z-fastq does not guess what an ambiguous base would have been:

```math
\text{gc\_fraction} = \frac{g + c}{a + c + g + t}
```

**Quality.** Scores are read as Phred+33. `quality_sum` adds them up, `mean_quality` divides by `bases`, and `q20_bases` and `q30_bases` count bases scoring at least 20 and 30, with `q20_fraction` and `q30_fraction` as their shares of `bases`. A quality byte outside ASCII 33 to 126 is an error.

A value that cannot be computed, such as a mean over no reads, prints as `-`, or `null` in JSON.

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/eneskemalergin/z-fastq/main/bench/stats/results/figures/file-overview-dark.png">
    <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/eneskemalergin/z-fastq/main/bench/stats/results/figures/file-overview-light.png">
    <img alt="Wall time and peak memory of other statistics tools relative to z-fastq, for plain and gzip input." src="https://raw.githubusercontent.com/eneskemalergin/z-fastq/main/bench/stats/results/figures/file-overview-light.png" width="760">
  </picture>
</p>
<p align="center"><sub><b>Figure 2.</b> The same comparison for <code>stats</code>. Details in the <a href="https://github.com/eneskemalergin/z-fastq/blob/main/bench/stats/README.md">stats benchmark</a>.</sub></p>

## Several inputs at once

Both commands take any number of inputs and keep each result separate. When one input fails, the rest are still read, and the exit status reports the worst failure ([Automation](Automation#exit-status)).

`stats` summarizes a file as a whole. Quality by cycle, adapter content, duplication, and the rest of a QC report come from `process` ([Quality control](Quality-Control)).

## Related pages

- [Validation](Validation): the full record checks, and what their error codes mean.
- [Automation](Automation#json): the `stats --json` document.
- [Quality control](Quality-Control): per-cycle and per-read quality, adapters, and duplication.
