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

z-fastq does not catch SIGINT or SIGTERM: the signal ends it, and the shell reports `130` or `143`. A program that waits for z-fastq directly sees the signal, not an exit status.

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

Commands that write FASTQ stream it: records go out as they are read. A failure later in the input, such as a malformed record or a damaged gzip member, leaves the earlier records written; the two files of `deinterleave` then end with the same complete pairs. A failed write, such as a full disk, can leave an output ending inside a record, and the two files of `deinterleave` with different numbers of records. An interrupted run leaves its output in place too, and `deinterleave` keeps both files.

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

### The `stats` and `check` documents

Each is one JSON object on one line. After `schema`, `tool`, and `byte_strings` ([Report reference](Report-Reference#shared-conventions)), `results` holds one entry per input:

```json
{"schema":"z-fastq/check-v1","tool":{"name":"z-fastq","version":"0.2.1"},"byte_strings":"escaped-bytes-v1","results":[{"input":"ok.fq","status":"ok"},{"input":"bad.fq","status":"error","error":{"code":"S002","message":"sequence byte is outside the selected alphabet","record_index":0,"byte_offset":5,"line_in_record":2}}]}
```

- **`input`**: the path as given, or **`inputs`**, the two paths, for `check --paired`.
- **`status`**: `ok` or `error`. A `stats-v1` entry with `ok` also holds every field of [`stats`](Counting-and-Statistics#stats), with `null` where the text output prints `-`.
- **`error`**: the failure's `code`, `message`, and position. A record failure (`S001` to `S006`) gives `record_index`, `byte_offset`, and `line_in_record`; a failure with no position in the records, such as `io_error` for a missing file, gives the same fields as `null`.
- **`failed_input`**: in a failed `check --paired` entry, the input whose record failed.
- A pair failure in `check --paired` or `--interleaved` gives `pair_index` and `record_indexes`, one per mate: in `--paired`, the last record read on each side; in `--interleaved`, `null` for a mate that does not exist. `P001` (names differ) adds `byte_offsets`, and `first_tokens` and `normalized_ids`, each mate's as `prefix`, `length`, and `truncated`, and `mate_markers`; `P002` (a missing mate) adds `remaining_side`, `R1` or `R2`.

## Related pages

- [Input and output](Input-and-Output): where data and messages go.
- [Report reference](Report-Reference): the fields of the `process` documents.
- [Recipes](Recipes): short patterns built on these rules.
