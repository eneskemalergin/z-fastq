# Automation

This page is for scripts and workflow engines. The one rule to remember: **trust the exit status, not the output.** A run can write records and then fail, and only the status tells the two apart.

## Exit status

| Status       | Meaning                                           |
| ------------ | ------------------------------------------------- |
| `0`          | success                                           |
| `1`          | a record or pair failed validation                |
| `2`          | the command line is invalid                       |
| `3`          | I/O, gzip, output, or memory failure              |
| `4`          | a line, sampling, or arithmetic limit was reached |
| `130`, `143` | interrupted by SIGINT or SIGTERM                  |

`count`, `stats`, and `check` on several independent inputs carry on past a failed input and exit with the highest status seen. Pair commands and `process` stop at the first failure.

## Data and messages

Data goes to standard output and messages to standard error, so capture them separately:

```bash
if ! z-fastq stats reads.fastq > reads.stats.txt 2> reads.stats.err; then
    cat reads.stats.err >&2
    exit 1
fi
```

An invalid command line (status `2`) is reported on standard error before anything else is written, including the case of an input that is also the output.

## Partial output

Commands that write FASTQ stream it: records go out as they are read. A failure later in the input, such as a malformed record or a damaged gzip member, leaves the earlier records written. An interrupted run leaves its output in place too, and `deinterleave` keeps both files.

Write to a temporary path and move it into place only after a zero status:

```bash
if z-fastq sample --count 10000 --seed 2025 reads.fastq > sample.fastq.tmp; then
    mv sample.fastq.tmp sample.fastq
else
    rm -f sample.fastq.tmp
    exit 1
fi
```

## JSON

| Command                | Schema                  | Fields                                                   |
| ---------------------- | ----------------------- | -------------------------------------------------------- |
| `stats --json`         | `z-fastq/stats-v1`      | [Counting and statistics](Counting-and-Statistics#stats) |
| `check --json`         | `z-fastq/check-v1`      | [Validation](Validation)                                 |
| `process --report`     | `z-fastq/process-v1`    | [Report reference](Report-Reference#report)              |
| `process --qc-summary` | `z-fastq/qc-summary-v1` | [Report reference](Report-Reference#qc-summary)          |

Check the `schema` field before reading a document; `tool.name` and `tool.version` say which build wrote it. Version 0.1.0 froze the `stats-v1` and `check-v1` shapes. Results follow the order of the inputs, and an input that failed has its error inside its own entry. Still check the exit status: a document can hold failed inputs, and an output failure can cut it short.

## Related pages

- [Input and output](Input-and-Output): where data and messages go.
- [Report reference](Report-Reference): the fields of the `process` documents.
- [Recipes](Recipes): short patterns built on these rules.
