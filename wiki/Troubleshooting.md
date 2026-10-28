# Troubleshooting

Find the symptom, then the fix.

## `z-fastq` is not found

The binary stays in the build tree: run `./zig-out/bin/z-fastq`, or link it into a directory on your `PATH` ([Installation](Installation)).

## The build cannot download zipir

The first build fetches the zipir release named in `build.zig.zon` and checks its hash. Build once with network access; later builds reuse the download.

## A gzip input fails

z-fastq decided it was gzip from its first bytes, so the whole stream must be valid. Test it on its own:

```bash
gzip -t reads.fastq.gz
```

Trailing bytes after the last member, even zero padding, also fail ([Limits](Limits-and-Supported-Formats#gzip)).

## A line is too long

The default limit is 16 MiB per line. Raise it only as far as the data needs, for example with `--max-line-bytes 33554432` for long reads. It does not allow wrapped records.

## Exact sampling refuses `-`

Exact counts read the file twice, so they need a regular file. On a pipe, use `--fraction`.

## Pair names do not match

Look at the first words of the two headers. The default `illumina` policy understands `/1` and `/2` and `1:N:...` markers; `exact` compares the first words byte for byte ([Paired reads](Paired-Reads)). No policy makes different read names equal.

## An output path is refused

`deinterleave` and `process` never overwrite a file. Remove the old output yourself once you are sure it is safe, or choose a new path.

## Output was written, but the command failed

Streaming commands write records as they go, and a later failure, such as a broken record or a gzip member whose checksum fails, leaves the earlier records written. Trust the exit status ([Automation](Automation#partial-output)).

## There is no JSON output

`--json` exists for `stats` and `check`; `process` writes its JSON with `--report` and `--qc-summary`. `count` prints plain numbers ([Automation](Automation#json)).

## The command exits with status 4

A limit was reached: the line limit, a sampling limit, or an arithmetic limit. The message names which.

## Related pages

- [FAQ](FAQ): short answers to common questions.
- [Automation](Automation): exit statuses in full.
- [Validation](Validation#when-a-record-fails): what each error code means.
