# Counting and statistics

## Count records

Use `count` when the only fact you need is how many FASTQ records an input contains:

```bash
z-fastq count reads.fastq.gz
```

It prints one decimal count for each successful input, in argument order:

```bash
z-fastq count lane1.fastq.gz lane2.fastq.gz
```

This command uses a specialized count path. It does not need to retain the header, sequence, or quality fields after it recognizes the record boundaries.

## Report aggregate statistics

```bash
z-fastq stats reads.fastq.gz
```

The human-readable report contains these fields:

| Field | Meaning |
| --- | --- |
| `reads` | Number of FASTQ records. |
| `bases` | Total sequence bases. |
| `min_length`, `max_length` | Shortest and longest sequence lengths. |
| `mean_length` | `bases / reads` when reads exist. |
| `a`, `c`, `g`, `t`, `n` | Case-insensitive counts for those symbols. |
| `other_bases` | Symbols outside A, C, G, T, and N. |
| `gc_fraction` | `(g + c) / (a + c + g + t)`. |
| `quality_sum` | Sum of decoded Phred scores. |
| `mean_quality` | `quality_sum / bases` when bases exist. |
| `q20_bases`, `q20_fraction` | Bases with Phred score at least 20. |
| `q30_bases`, `q30_fraction` | Bases with Phred score at least 30. |

The quality encoding is Phred+33. Quality bytes must be ASCII 33 through 126.

## Interpret empty and ambiguous data

Undefined values print as `-` in human-readable output. This includes ratios with an empty denominator. A zero-length record still counts as one read, so `reads` and `bases` answer different questions. The current machine-readable representation is provisional; see [Automation](Automation) before building an integration around it.

`gc_fraction` excludes N and other symbols from its denominator. z-fastq does not guess whether an ambiguous symbol should count as A, C, G, or T.

## Multiple inputs

`count` and `stats` accept independent inputs. They keep each result separate and process inputs in argument order. If one input fails, the command can still attempt later inputs. The final exit status reports the highest failure class. See [Automation](Automation).

## What this report does not provide

The current report is an aggregate summary. It does not include quality histograms, per-cycle summaries, N50, adapter content, duplicate estimates, or platform classification. Those belong to a broader QC workflow rather than this command's current contract.
