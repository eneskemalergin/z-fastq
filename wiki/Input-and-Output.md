# Input and output

## Reading

Every command reads plain FASTQ and gzip through the same path. z-fastq looks at the first two bytes to decide: a file called `reads.fastq` that holds gzip is read as gzip, and `reads.gz` that holds plain text is read as text.

Once the first bytes say gzip, the whole stream must be valid gzip; a damaged stream is an error and is never retried as text. Concatenated gzip members are read in order, and each one's CRC-32 and length are checked. Those checks sit at the end of a member, so a command may already have written records from a member before it finds the member damaged ([Automation](Automation#partial-output)).

## Standard input

Standard input is read only when you write `-`; leaving out the path is an error, not a request for stdin.

```bash
gzip -dc reads.fastq.gz | z-fastq stats -
```

A command reads standard input at most once. A pipe can be read only once, so other names for it, such as `/dev/stdin`, count as the same input. Exact-count sampling needs a regular file, because it reads twice.

## Writing

Results, JSON, and FASTQ go to standard output; messages go to standard error. FASTQ is always written as four lines per record with LF endings.

| Command                                       | Reads              | Writes                       |
| --------------------------------------------- | ------------------ | ---------------------------- |
| `count`, `stats`, `check`                     | one or more inputs | standard output              |
| `check --paired`, `sample --paired`           | two inputs         | standard output              |
| `check --interleaved`, `sample --interleaved` | one input          | standard output              |
| `interleave`                                  | two inputs         | standard output              |
| `deinterleave`                                | one input          | two new files                |
| `process`                                     | two inputs         | two new files, plain or gzip |

Only `process --gzip` writes gzip; for the other commands, pipe the output to `gzip`. New output files are never overwritten: an existing path fails before anything is written.

Standard output may not be one of the inputs. `z-fastq count a.fastq b.fastq >> b.fastq` fails with `input file is output file` (status `2`) before reading anything, so the input is never touched. If standard output fails, for example because the reader of a pipe has exited, the command reports `error: standard output: I/O error` and exits with `3`.

## Memory

z-fastq reads through fixed buffers; no command maps or loads a whole file. The [limits page](Limits-and-Supported-Formats#memory) lists what each command keeps.

## Related pages

- [Automation](Automation): exit statuses, and what to do with partial output.
- [Limits and supported formats](Limits-and-Supported-Formats): the formats read and the memory each command keeps.
- [Troubleshooting](Troubleshooting#a-gzip-input-fails): when a gzip input fails.
