# Processing paired reads

`process` reads two mate files in lockstep, validates every pair as `check --paired` does, and writes each accepted pair to two new files, plain or gzip. With `--report` it also writes a JSON report with per-mate statistics. It changes no read.

```bash
z-fastq process --gzip --out1 clean_R1.fastq.gz --out2 clean_R2.fastq.gz --report sample.json R1.fastq.gz R2.fastq.gz
```

## Outputs

- `--out1` and `--out2` are created as new files. An existing path, including one of the inputs, fails with `output path already exists` (status `3`) before anything is written, and outputs created by the same run are removed.
- Records are written with LF line endings, in input order, mate 1 to `--out1` and mate 2 to `--out2`.
- `--gzip` compresses both outputs. A name ending in `.gz` without `--gzip` is rejected (status `2`), so a plain file never carries a gzip name. The gzip header holds no time stamp, so the same inputs and options give the same bytes.
- Standard input may be one of the two inputs; outputs cannot be standard output.

## Stop after N pairs

```bash
z-fastq process --max-pairs 1000000 --out1 sub_R1.fastq --out2 sub_R2.fastq R1.fastq.gz R2.fastq.gz
```

`--max-pairs N` writes the first N complete pairs in input order and stops reading. Nothing after pair N is checked, so a missing mate or a damaged gzip trailer later in the files goes unseen; the report says `"input_complete": false`. Use `check --paired` to validate whole files.

## The report

`--report PATH` writes a `z-fastq/process-v1` document after both outputs are finished; `--report -` writes it to standard output. It carries:

- `options`: the effective options.
- `inputs`: each input's path, and its size and modification time in nanoseconds when it is a regular file (`null` for a pipe).
- `outputs`: the output paths.
- `status`: `ok` or `error`; `input_complete`: whether both inputs were read to the end.
- `pairs`: pairs `read` and `written`.
- `mates`: one statistics object per mate, with the same fields as `stats --json`.
- `failed_input`, `failed_output`, and `error`: the failure, in the form `check --paired --json` uses, or `null`.

The report holds no time of the run or host name, so a rerun gives the same bytes.

## Failure and completion

A failure stops the run. Pairs already written stay in the outputs; if mate 2 cannot be written after mate 1, `--out1` holds one record more than `--out2`. The report, when requested, is still written with `"status": "error"`.

A run is complete only when the status is `0`, and, when a report was requested, the report says `"status": "ok"`. Treat any other result as a failed sample. Exit statuses are the usual ones: `1` validation, `2` usage, `3` I/O or allocation, `4` limits.
